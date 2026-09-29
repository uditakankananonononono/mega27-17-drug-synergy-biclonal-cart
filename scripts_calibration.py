"""Amendment-queue Tier-3 #16: calibration of Lane A probability models.
Prespec: results/calibration_prespec.json (committed 1655bb1 BEFORE this run).
5-fold CV out-of-fold probabilities from the #6 missingness logistic models;
Brier (+ prevalence null), ECE (10 bins), reliability tables; consequence rule
applied to the #6 IPW weights verbatim from the prespec.
Run: python3 scripts_calibration.py -> results/calibration.json"""
import json, csv, gzip, io, zipfile
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
from sklearn.metrics import roc_auc_score, brier_score_loss

SEED = 20260929

# ---- same covariate rebuild + outcomes as scripts_missingness.py ----
with gzip.open("data/propensity/gtex_v8_median_tpm.gct.gz", "rt") as f:
    f.readline(); f.readline()
    header = f.readline().rstrip("\n").split("\t")
    expr, ensg_of_sym = {}, {}
    for line in f:
        p = line.rstrip("\n").split("\t")
        ensg = p[0].split(".")[0]
        expr[ensg] = np.array([float(x) for x in p[2:]], dtype=np.float32)
        ensg_of_sym.setdefault(p[1], ensg)
hpa_class, hpa_prog = {}, {}
with zipfile.ZipFile("data/propensity/proteinatlas.tsv.zip") as z:
    with z.open("proteinatlas.tsv") as fb:
        rd = csv.reader(io.TextIOWrapper(fb), delimiter="\t", quoting=csv.QUOTE_ALL)
        cols = next(rd)
        i_class = cols.index("Protein class")
        prog_idx = [i for i, c in enumerate(cols) if c.startswith("Cancer prognostics -")]
        for row in rd:
            hpa_class[row[2]] = row[i_class]
            vals = [row[i] for i in prog_idx]
            hpa_prog[row[2]] = (sum(1 for v in vals if "favorable" in v and "unfavorable" not in v),
                                sum(1 for v in vals if "unfavorable" in v))
lengths = json.load(open("data/propensity/mygene_lengths.json"))
lengths_ensg = {ensg_of_sym[s]: L for s, L in lengths.items() if s in ensg_of_sym}

def feats_ensg(g):
    v = expr[g]
    pc = hpa_class[g].lower()
    fav, unfav = hpa_prog.get(g, (0, 0))
    return [np.log1p(v.mean()), (v >= 1).mean(), np.log10(lengths_ensg[g]),
            float("membrane" in pc), float("secreted" in pc), float(fav), float(unfav)]

prows = list(csv.DictReader(open("results/pharos_per_gene.csv")))
gated_ensg = {ensg_of_sym[r["gene"]] for r in prows if r["group"] == "gated" and r["gene"] in ensg_of_sym}
excl = gated_ensg | {ensg_of_sym[r["gene"]] for r in prows if r["group"] in ("references", "positive_controls") and r["gene"] in ensg_of_sym}
universe = sorted(set(expr) & set(hpa_class) & set(lengths_ensg))
pool = [g for g in universe if g not in excl]

subcell = set()
with open("data/subcellular_location.tsv") as f:
    rd = csv.reader(f, delimiter="\t"); next(rd)
    for row in rd:
        subcell.add(row[0].split(".")[0])
uniprot_mem = set()
with open("data/uniprot_cellmembrane.tsv") as f:
    rd = csv.reader(f, delimiter="\t")
    cols = next(rd); gi = cols.index("Gene Names")
    for row in rd:
        toks = row[gi].split()
        if toks and toks[0] in ensg_of_sym:
            uniprot_mem.add(ensg_of_sym[toks[0]])

Xp = np.array([feats_ensg(g) for g in pool])
mu, sd = Xp.mean(0), Xp.std(0); sd[sd == 0] = 1.0
Zp = (Xp - mu) / sd
YS = {"has_subcellular": np.array([float(g in subcell) for g in pool]),
      "has_uniprot_membrane": np.array([float(g in uniprot_mem) for g in pool])}

kf = KFold(n_splits=5, shuffle=True, random_state=SEED)
out = {"seed": SEED, "prespec": "results/calibration_prespec.json", "n_pool": len(pool), "per_model": {}}
for oc, y in YS.items():
    oof = np.zeros_like(y)
    for tr, te in kf.split(Zp):
        m = LogisticRegression(max_iter=2000).fit(Zp[tr], y[tr])
        oof[te] = m.predict_proba(Zp[te])[:, 1]
    brier = brier_score_loss(y, oof)
    prev = y.mean()
    prev_brier = float(np.mean((y - prev) ** 2))
    bins = np.linspace(0, 1, 11)
    rel, ece = [], 0.0
    for i in range(10):
        m_ = (oof >= bins[i]) & (oof < bins[i + 1] if i < 9 else oof <= bins[i + 1])
        if m_.sum() == 0:
            continue
        mp, mf = float(oof[m_].mean()), float(y[m_].mean())
        rel.append({"bin": [float(bins[i]), float(bins[i + 1])], "n": int(m_.sum()),
                    "mean_pred": mp, "empirical": mf, "gap": mp - mf})
        ece += (m_.sum() / len(y)) * abs(mp - mf)
    auc = roc_auc_score(y, oof)
    calibrated = bool(brier < prev_brier and ece < 0.05)
    over = [r for r in rel if r["gap"] > 0.05]
    under = [r for r in rel if r["gap"] < -0.05]
    direction = ("over-confident in " + str(len(over)) + " bin(s)" if over else "") + ("; " if over and under else "") + ("under-confident in " + str(len(under)) + " bin(s)" if under else "") or "no bin deviates >5pp"
    out["per_model"][oc] = {"cv_auc_oof": float(auc), "brier": float(brier),
                            "prevalence_brier": prev_brier, "ece_10bin": float(ece),
                            "usefully_calibrated_prespec": calibrated,
                            "miscalibration_direction": direction,
                            "reliability": rel}
    print(f"{oc}: AUC {auc:.3f} | Brier {brier:.4f} vs prev {prev_brier:.4f} | ECE {ece:.4f} | calibrated={calibrated} | {direction}")

worst = max(m["ece_10bin"] for m in out["per_model"].values())
out["consequence_for_6_ipw"] = ("all models ECE < 0.05 - #6 IPW weights stand as computed" if worst < 0.05
                                 else f"max ECE {worst:.4f} >= 0.05 - #6 P(observation) miscalibrated; IPW survival numbers are approximate per prespec consequence rule")
json.dump(out, open("results/calibration.json", "w"), indent=1)
print("wrote results/calibration.json |", out["consequence_for_6_ipw"])

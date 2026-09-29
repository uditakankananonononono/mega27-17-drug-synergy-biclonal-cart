"""Amendment-queue Tier-3 #7: MS protein-detectability model.
Prespec: results/detectability_prespec.json (committed 00951b4 BEFORE this run).
Predicts MS normal-tissue detection breadth (n_normal_tissues, committed
proteomicsdb_per_gene.csv) from the 7 locked propensity features; tests whether
the H1 breadth contrast (gated 20 vs background 5 normal tissues) is a
detectability-expectation artifact; lists under-detected antigens (safety-relevant:
their normal breadth is likely undercounted).
Run: python3 scripts_detectability.py -> results/detectability_model.json"""
import json, csv, gzip, io, zipfile
import numpy as np
from sklearn.linear_model import PoissonRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import mean_poisson_deviance

SEED = 20260929
FN = ["expr_log1p_mean", "expr_breadth_tpm1", "gene_len_log10", "is_membrane",
      "is_secreted", "prog_favourable_n", "prog_unfavourable_n"]

# ---- covariate table (identical rebuild to #3/#6) ----
with gzip.open("data/propensity/gtex_v8_median_tpm.gct.gz", "rt") as f:
    f.readline(); f.readline()
    header = f.readline().rstrip("\n").split("\t")
    expr, ensg_of_sym, sym_of_ensg = {}, {}, {}
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

usable = set(expr) & set(hpa_class) & set(lengths_ensg)

# ---- study genes with breadth outcome ----
rows = [r for r in csv.DictReader(open("results/proteomicsdb_per_gene.csv")) if r["mapped"] == "True"]
X, Y, G, SYM, dropped = [], [], [], [], 0
for r in rows:
    e = ensg_of_sym.get(r["gene"])
    if e is None or e not in usable:
        dropped += 1
        continue
    X.append(feats_ensg(e)); Y.append(float(r["n_normal_tissues"])); G.append(r["group"]); SYM.append(r["gene"])
X, Y = np.array(X), np.array(Y)
print(f"study genes usable {len(SYM)}/{len(rows)} (dropped {dropped})")

kf = KFold(n_splits=5, shuffle=True, random_state=SEED)
devs, pr2 = [], []
for tr, te in kf.split(X):
    m = PoissonRegressor(alpha=1.0).fit(X[tr], Y[tr])
    d = mean_poisson_deviance(Y[te], m.predict(X[te]))
    nd = mean_poisson_deviance(Y[te], np.full_like(Y[te], Y[tr].mean()))
    devs.append(d); pr2.append(1.0 - d / nd if nd > 0 else float("nan"))
mfull = PoissonRegressor(alpha=1.0).fit(X, Y)
EXP = mfull.predict(X)
print(f"CV pseudo-R2 {np.nanmean(pr2):.3f} deviance {np.mean(devs):.3f}")

def grp_med(mask, arr):
    return float(np.median(arr[mask])) if mask.any() else float("nan")
G = np.array(G)
gated_m, bg_m = G == "gated", G == "background"
obs_gap = grp_med(gated_m, Y) - grp_med(bg_m, Y)
exp_gap = grp_med(gated_m, EXP) - grp_med(bg_m, EXP)
explained = exp_gap / obs_gap if obs_gap != 0 else float("nan")
h1_reading = ("observed breadth gap substantially explained by detectability expectation (>=50%)"
              if abs(explained) >= 0.5 and np.sign(explained) == np.sign(obs_gap)
              else "observed breadth gap exceeds detectability expectation - residual real difference")
print(f"gated obs med {grp_med(gated_m, Y):.1f} exp med {grp_med(gated_m, EXP):.1f} | bg obs med {grp_med(bg_m, Y):.1f} exp med {grp_med(bg_m, EXP):.1f}")
print(f"observed gap {obs_gap:+.1f} vs expected gap {exp_gap:+.1f} (explained {explained:.2f}) -> {h1_reading}")

resid = Y - EXP
order = np.argsort(resid)
under = [{"gene": SYM[i], "group": G[i], "observed": float(Y[i]), "expected": float(EXP[i]), "residual": float(resid[i])} for i in order[:12]]
out = {"seed": SEED, "prespec": "results/detectability_prespec.json", "features": FN,
       "n_genes": len(SYM), "dropped_no_covariates": dropped,
       "cv_poisson_deviance_mean": float(np.mean(devs)), "cv_pseudo_r2_mean": float(np.nanmean(pr2)),
       "group_medians": {"gated_observed": grp_med(gated_m, Y), "gated_expected": grp_med(gated_m, EXP),
                          "background_observed": grp_med(bg_m, Y), "background_expected": grp_med(bg_m, EXP)},
       "observed_gap": obs_gap, "expected_gap": exp_gap, "gap_fraction_explained": float(explained),
       "h1_reading": h1_reading,
       "most_underdetected": under,
       "per_gene_residuals": [{"gene": SYM[i], "group": G[i], "observed": float(Y[i]),
                                "expected": float(EXP[i]), "residual": float(resid[i])} for i in range(len(SYM))],
       "thin_n_caveat": "n~196 with 7 features; point estimates with that caveat, no tuning. Sequence-derived features unavailable (gene length proxy). Breadth from committed CSV, raw wiped."}
json.dump(out, open("results/detectability_model.json", "w"), indent=1)
print("wrote results/detectability_model.json")

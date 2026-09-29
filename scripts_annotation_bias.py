"""Amendment-queue Tier-3 #3: annotation-bias correction model.
Prespec: results/annotation_bias_prespec.json (committed 52fecf8 BEFORE this run).
Idea (judge verdict): candidates are over-annotated vs matched nulls on literature
metrics (Tier-1 #2: publications 118 vs 42.5, p<0.001). Is that biology or residual
annotation bias? Fit an EXPECTED-ANNOTATION model: Poisson regression (log link) of
each annotation outcome on the 7 locked propensity features, trained on the 831
unique matched controls (the covariate-equivalent reference class). Corrected
enrichment = observed - expected per candidate. Pre-declared endpoints and
interpretation rule in the prespec; results reported as measured either direction.
Run: python3 scripts_annotation_bias.py -> results/annotation_bias_correction.json"""
import json, csv, gzip, io, zipfile
import numpy as np
from sklearn.linear_model import PoissonRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import mean_poisson_deviance

SEED = 20260929
OUTCOMES = ["publication_count", "generif_count", "gwas_total", "ppi_total", "n_drugs"]
FN = ["expr_log1p_mean", "expr_breadth_tpm1", "gene_len_log10", "is_membrane",
      "is_secreted", "prog_favourable_n", "prog_unfavourable_n"]

# ---- rebuild the covariate table (identical logic to scripts_propensity.py) ----
with gzip.open("data/propensity/gtex_v8_median_tpm.gct.gz", "rt") as f:
    f.readline(); f.readline()
    header = f.readline().rstrip("\n").split("\t")
    expr, ensg_of_sym, sym_of_ensg = {}, {}, {}
    for line in f:
        p = line.rstrip("\n").split("\t")
        ensg = p[0].split(".")[0]
        expr[ensg] = np.array([float(x) for x in p[2:]], dtype=np.float32)
        sym_of_ensg[ensg] = p[1]
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

# ---- controls (831 symbols) + candidates (gated) with outcomes ----
ctrl_stats = json.load(open("data/propensity/matched_pharos_stats.json"))
prows = list(csv.DictReader(open("results/pharos_per_gene.csv")))
cand_stats = {r["gene"]: {k: float(r[k]) for k in OUTCOMES} for r in prows if r["group"] == "gated" and r["status"] == "ok"}

def assemble(stats_by_sym):
    X, Y, names, dropped = [], [], [], 0
    for sym, oc in stats_by_sym.items():
        e = ensg_of_sym.get(sym)
        if e is None or e not in usable:
            dropped += 1
            continue
        X.append(feats_ensg(e)); Y.append([oc[k] for k in OUTCOMES]); names.append(sym)
    return np.array(X), np.array(Y), names, dropped

Xc, Yc, ctrl_names, ctrl_drop = assemble(ctrl_stats)
Xg, Yg, cand_names, cand_drop = assemble(cand_stats)
print(f"controls usable {len(ctrl_names)}/{len(ctrl_stats)} (dropped {ctrl_drop}); candidates usable {len(cand_names)}/{len(cand_stats)} (dropped {cand_drop})")

# ---- per-outcome Poisson expectation model ----
kf = KFold(n_splits=5, shuffle=True, random_state=SEED)
out = {"seed": SEED, "prespec": "results/annotation_bias_prespec.json",
       "features": FN, "n_controls_train": len(ctrl_names), "n_candidates": len(cand_names),
       "controls_dropped_no_covariates": ctrl_drop, "candidates_dropped_no_covariates": cand_drop,
       "per_outcome": {}, "per_candidate": []}
for j, oc in enumerate(OUTCOMES):
    y = Yc[:, j]
    devs, pr2 = [], []
    for tr, te in kf.split(Xc):
        m = PoissonRegressor(alpha=1.0).fit(Xc[tr], y[tr])
        pred = m.predict(Xc[te])
        devs.append(mean_poisson_deviance(y[te], pred))
        null_dev = mean_poisson_deviance(y[te], np.full_like(y[te], y[tr].mean()))
        pr2.append(1.0 - devs[-1] / null_dev if null_dev > 0 else float("nan"))
    mfull = PoissonRegressor(alpha=1.0).fit(Xc, y)
    exp_g = mfull.predict(Xg)
    obs_g = Yg[:, j]
    corr = obs_g - exp_g
    frac = float((obs_g > exp_g).mean())
    out["per_outcome"][oc] = {
        "cv_poisson_deviance_mean": float(np.mean(devs)),
        "cv_pseudo_r2_mean": float(np.nanmean(pr2)),
        "candidate_expected_median": float(np.median(exp_g)),
        "candidate_observed_median": float(np.median(obs_g)),
        "median_corrected_enrichment": float(np.median(corr)),
        "fraction_observed_gt_expected": frac,
        "control_observed_median": float(np.median(y))}
    print(f"{oc}: pseudoR2 {np.nanmean(pr2):.3f} | obs med {np.median(obs_g):.1f} vs exp med {np.median(exp_g):.1f} | corrected med {np.median(corr):+.1f} | O>E {frac:.2f}")

# ---- per-candidate detail (all 5 outcomes) ----
models = {oc: PoissonRegressor(alpha=1.0).fit(Xc, Yc[:, j]) for j, oc in enumerate(OUTCOMES)}
for i, sym in enumerate(cand_names):
    rec = {"gene": sym}
    for j, oc in enumerate(OUTCOMES):
        e = float(models[oc].predict(Xg[i:i+1])[0])
        rec[oc] = {"observed": float(Yg[i, j]), "expected": e, "corrected": float(Yg[i, j] - e)}
    out["per_candidate"].append(rec)

# ---- pre-declared endpoint verdicts (rule fixed in prespec) ----
verdicts = {}
for oc in ["publication_count", "generif_count"]:
    po = out["per_outcome"][oc]
    residual = po["median_corrected_enrichment"] > 0 and po["fraction_observed_gt_expected"] > 0.5
    verdicts[oc] = ("residual over-annotation survives covariate correction - annotation bias on axes outside the 7 matched features (literature-covariate gap from #2), reported as measured"
                    if residual else
                    "raw gap confounded by the matched features - corrected enrichment ~ 0")
out["endpoint_verdicts"] = verdicts
json.dump(out, open("results/annotation_bias_correction.json", "w"), indent=1)
print("wrote results/annotation_bias_correction.json")

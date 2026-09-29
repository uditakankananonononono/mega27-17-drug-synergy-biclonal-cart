"""Amendment-queue Tier-3 #8: RNA->protein translation model (closes PARTIAL gap).
Prespec: results/translation_prespec.json (committed 123b332 BEFORE this run).
Per-cancer Huber fit protein_median ~ log1p(HPA tumor RNA mean) across the study
gene set; per-gene residuals; cross-cancer same-sign consistency (>=3/7 cancers,
top-decile |median residual|) = unusual translation behavior, listed with direction.
Reuses committed pairing helpers from scripts_cptac_protein.py; all local.
Run: python3 scripts_translation.py -> results/translation_model.json"""
import json
import numpy as np
from sklearn.linear_model import HuberRegressor
from scripts_cptac_protein import (ROWS, META, SETS, TP, CANCER_MAP, load_rows,
                                   gene_study_metrics, symbol_to_ensg)

SEED = 20260929
df = load_rows()
meta = json.load(open(META))
sets = json.load(open(SETS))
tp = json.load(open(TP))["data"]
sc = meta["sample_counts"]
gsm = gene_study_metrics(df, sc)
sym2ensg = symbol_to_ensg()

per_cancer, gene_res = {}, {}
for study, cancer in CANCER_MAP.items():
    xs, ys, names = [], [], []
    for sym, met in gsm.items():
        ensg = sym2ensg.get(sym)
        if ensg is None or ensg not in tp or study not in met:
            continue
        rec = tp[ensg].get(cancer)
        if rec is None:
            continue
        xs.append(float(rec[3])); ys.append(met[study]["median"]); names.append(sym)
    if len(xs) < 20:
        continue
    X = np.log1p(np.array(xs)).reshape(-1, 1)
    y = np.array(ys)
    m = HuberRegressor().fit(X, y)
    pred = m.predict(X)
    r2 = m.score(X, y)
    per_cancer[cancer] = {"n": len(xs), "r2": float(r2),
                          "slope": float(m.coef_[0]), "intercept": float(m.intercept_)}
    for i, sym in enumerate(names):
        gene_res.setdefault(sym, {})[cancer] = float(y[i] - pred[i])
    print(f"{cancer}: n={len(xs)} R2={r2:.3f} slope={m.coef_[0]:.3f}")

# cross-cancer consistency
med_res = {s: float(np.median(list(r.values()))) for s, r in gene_res.items() if len(r) >= 3}
vals = np.array(list(med_res.values()))
p90 = np.percentile(np.abs(vals), 90) if len(vals) else float("nan")
outliers = []
gated_set, bg_set = set(sets["gated"]), set(sets["background"])
for sym, mr in med_res.items():
    r = gene_res[sym]
    signs = {np.sign(v) for v in r.values()}
    if len(r) >= 3 and len(signs) == 1 and abs(mr) >= p90:
        outliers.append({"gene": sym, "group": "gated" if sym in gated_set else ("background" if sym in bg_set else "other"),
                         "n_cancers": len(r), "median_residual": mr,
                         "direction": "protein higher than RNA predicts" if mr > 0 else "protein lower than RNA predicts"})
outliers.sort(key=lambda d: -abs(d["median_residual"]))
gated_out = [o for o in outliers if o["group"] == "gated"]
print(f"outliers: {len(outliers)} total, {len(gated_out)} gated")
out = {"seed": SEED, "prespec": "results/translation_prespec.json",
       "per_cancer_model": per_cancer,
       "n_genes_with_3plus_cancers": len(med_res),
       "abs_residual_p90": float(p90),
       "outliers": outliers,
       "gated_outliers": gated_out,
       "residuals_all": {s: {"n_cancers": len(gene_res[s]), "median_residual": med_res[s],
                              "per_cancer": gene_res[s]} for s in med_res},
       "caveats": ["cross-cohort pairing (HPA RNA tumor mean vs CPTAC protein median), declared in prespec",
                    "tumor tissue only - no normal-tissue claim"]}
json.dump(out, open("results/translation_model.json", "w"), indent=1)
print("wrote results/translation_model.json")

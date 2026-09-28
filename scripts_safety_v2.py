"""Tier-1 #13 gap: fold GWAS/Mendelian/essentiality terms into the safety axis.
ADDITIVE analysis - rankscore v1 is NOT modified. safety_v2 = equal-weight mean of
4 min-max-normalized terms (higher = safer): (1) v1 normal_safety carried over,
(2) essentiality safety = -DepMap median Chronos (24Q4, 1,178 models),
(3) Mendelian safety = -MONARCH causal_total, (4) GWAS safety = -n_sig (p<=5e-8),
each normalized across the 7-antigen panel. Equal weights are arbitrary and
documented; reported as measured. Consistency check vs #15 GEO external findings.
Writes results/safety_v2.json."""
import json, csv
import numpy as np
from scipy import stats

PANEL = ["CA9","MSLN","SLC39A6","PSCA","CLDN18","CLDN6","CA12"]
rk = json.load(open("results/rankscore_v1.json"))
ns_v1 = {g: rk["rows"][g]["components"]["normal_safety"] for g in PANEL}

gwas = {r["gene"]: int(r["n_sig"]) for r in csv.DictReader(open("results/gwas_per_gene.csv"))}
mon = {r["gene"]: (int(r["causal_total"]) if r["causal_total"] not in ("", None) else 0) for r in csv.DictReader(open("results/monarch_per_gene.csv"))}

rows = list(csv.DictReader(open("results/depmap_subset.csv")))
dep = {}
for g in PANEL:
    vals = [float(r[g]) for r in rows if r.get(g) not in (None, "")]
    dep[g] = float(np.median(vals)) if vals else None

def norm(d, keys, invert=False):
    v = np.array([d[k] for k in keys], dtype=float)
    if invert: v = -v
    lo, hi = v.min(), v.max()
    return {k: (0.5 if hi == lo else (x - lo) / (hi - lo)) for k, x in zip(keys, v)}

t_ns   = ns_v1
t_dep  = norm(dep, PANEL, invert=True)   # more negative Chronos = more essential = less safe
t_mon  = norm(mon, PANEL, invert=True)
t_gwas = norm(gwas, PANEL, invert=True)

out = {"design": "equal-weight mean of 4 min-max-normalized safety terms; additive; v1 untouched",
       "terms": {}, "safety_v2": {}, "raw_inputs": {"depmap_median_chronos": dep,
       "monarch_causal_total": {g: mon[g] for g in PANEL}, "gwas_n_sig": {g: gwas[g] for g in PANEL}}}
for g in PANEL:
    out["terms"][g] = {"normal_safety_v1": t_ns[g], "essentiality": t_dep[g],
                       "mendelian": t_mon[g], "gwas": t_gwas[g]}
    out["safety_v2"][g] = float(np.mean([t_ns[g], t_dep[g], t_mon[g], t_gwas[g]]))
v2_order = sorted(PANEL, key=lambda g: -out["safety_v2"][g])
v1_order = sorted(PANEL, key=lambda g: -t_ns[g])
out["safety_ranking_v2"] = v2_order
out["safety_ranking_v1"] = v1_order
out["spearman_v1_vs_v2"] = float(stats.spearmanr(
    [v2_order.index(g) for g in PANEL], [v1_order.index(g) for g in PANEL]).statistic)
geo = json.load(open("results/external_bulk_validation.json"))
out["consistency_check_geo15"] = {g: {"safety_v2_rank": v2_order.index(g) + 1,
    "geo_enriched": geo["per_antigen"][g].get("enriched")} for g in PANEL}
json.dump(out, open("results/safety_v2.json", "w"), indent=1)
print("v2 order (safest first):", v2_order)
print("v1 order (safest first):", v1_order)
print("spearman v1 vs v2:", round(out["spearman_v1_vs_v2"], 3))
for g in PANEL:
    print(f"{g:9s} v2={out['safety_v2'][g]:.3f}  terms ns={t_ns[g]:.2f} dep={t_dep[g]:.2f} mon={t_mon[g]:.2f} gwas={t_gwas[g]:.2f}  raw dep={dep[g]:.3f} mon={mon[g]} gwas={gwas[g]}")

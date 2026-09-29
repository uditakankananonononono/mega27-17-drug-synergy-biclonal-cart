"""Tier-4 #10 v2: pooled section-level occupancy test.
Prespec: results/spatial_occupancy_v2_prespec.json (committed 2b9a75c BEFORE this run).
Reuses v1 per-section numbers VERBATIM from results/spatial_occupancy.json (87eb0cb);
no null recomputation. Endpoints: paired Wilcoxon on (obs - null median), Fisher
combined p of v1 per-section co p-values, bootstrap 95% CI on median effect.
Run: python3 scripts_spatial_occupancy_v2.py -> results/spatial_occupancy_v2.json"""
import json
import numpy as np
from scipy import stats

v1 = json.load(open("results/spatial_occupancy.json"))
rng = np.random.default_rng(20260929)
out = {"prespec": "results/spatial_occupancy_v2_prespec.json",
       "source": "results/spatial_occupancy.json (v1, verbatim per-section reuse)",
       "pairs": {}}
control_pass = {}
for pair, rec in v1["pairs"].items():
    secs = rec.get("sections", [])
    n = len(secs)
    cls = rec["class"]
    if n == 0:
        out["pairs"][pair] = {"class": cls, "n_sections": 0, "verdict": "NO QUALIFYING SECTIONS"}
        continue
    eff = np.array([s["observed_jaccard"] - s["null_median"] for s in secs])
    p_co = np.array([s["p_co"] for s in secs])
    if n < 8:
        out["pairs"][pair] = {"class": cls, "n_sections": n,
                              "median_effect": float(np.median(eff)),
                              "verdict": "UNDERPOWERED (n<8 per v2 prespec) - descriptive only"}
        continue
    w = stats.wilcoxon(eff, alternative="two-sided")
    fisher = stats.combine_pvalues(p_co, method="fisher")
    boots = [float(np.median(eff[rng.integers(0, n, n)])) for _ in range(10000)]
    lo, hi = float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))
    verdict = ("CO-OCCUPANCY" if (w.pvalue < 0.05 and fisher.pvalue < 0.05 and np.median(eff) > 0)
               else "NOT CO-OCCUPANCY (pooled)")
    if cls == "positive_controls":
        control_pass[pair] = verdict == "CO-OCCUPANCY"
    out["pairs"][pair] = {"class": cls, "n_sections": n,
                          "median_effect": float(np.median(eff)),
                          "effect_ci95": [lo, hi],
                          "wilcoxon_p": float(w.pvalue), "fisher_p": float(fisher.pvalue),
                          "verdict_pre_gate": verdict}
    print(f"{pair:16s} [{cls:18s}] n={n:3d} eff={np.median(eff):+.4f} [{lo:+.4f},{hi:+.4f}] W_p={w.pvalue:.2e} F_p={fisher.pvalue:.2e} -> {verdict}")
gate_ok = all(control_pass.values()) and len(control_pass) == 3
out["method_check"] = {"controls": control_pass,
                       "status": "METHOD OK - verdicts stand" if gate_ok else "METHOD FAILED v2 gate - pair verdicts WITHHELD"}
if not gate_ok:
    for pair, rec in out["pairs"].items():
        if rec.get("verdict_pre_gate"):
            rec["verdict"] = "WITHHELD (v2 control gate failed)"
else:
    for pair, rec in out["pairs"].items():
        if "verdict_pre_gate" in rec:
            rec["verdict"] = rec.pop("verdict_pre_gate")
json.dump(out, open("results/spatial_occupancy_v2.json", "w"), indent=1)
print("v2 gate:", out["method_check"]["status"])

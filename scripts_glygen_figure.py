#!/usr/bin/env python3
"""Glycan-shielding figure: reported glycosylation sites for AND-gate antigens + approved references (reported_with_glycan structures annotated), plus gated-vs-background density strip within the surface stratum. Reads results/glygen_audit.json."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = json.load(open("results/glygen_audit.json"))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.4), gridspec_kw={"width_ratios": [3, 2]})
ag = sorted(A["and_gate"], key=lambda r: -r["n_shield_sites"])
rf = sorted(A["references"], key=lambda r: -r["n_shield_sites"])
names = [r["gene"] for r in ag + rf]
vals = [r["n_shield_sites"] for r in ag + rf]
cols = ["#2c7fb8"] * len(ag) + ["#999999"] * len(rf)
ax1.bar(range(len(names)), vals, color=cols)
ax1.set_xticks(range(len(names)))
ax1.set_xticklabels(names, rotation=45, ha="right", fontsize=7)
ax1.set_ylabel("unique reported glycosylation sites")
ax1.set_title("AND-gate antigens vs approved references")
for i, r in enumerate(ag + rf):
    if r["n_shield_structures"]:
        ax1.text(i, vals[i] + 0.3, str(r["n_shield_structures"]), ha="center", fontsize=6, color="#444444")
s = A["surface_stratum"]
ax2.bar([0, 1], [100 * s["gated_with_shield"] / s["n_gated"], 100 * s["background_with_shield"] / s["n_background"]], color=["#2c7fb8", "#f03b20"], width=0.55)
ax2.set_xticks([0, 1])
ax2.set_xticklabels(["gated\n(n=%d)" % s["n_gated"], "background\n(n=%d)" % s["n_background"]], fontsize=8)
ax2.set_ylabel("% surface proteins with reported sites")
ax2.set_title("surface stratum, Fisher p=%.3f" % s["fisher_p"], fontsize=8)
fig.tight_layout()
fig.savefig("paper/glygen_fig.pdf")
print("figure written")

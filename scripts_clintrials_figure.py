"""Figure for the ClinicalTrials.gov landscape audit (paper/clintrials_fig.pdf).
Panel a: alias-tier trial counts for the 7 AND-gate antigens + 4 references (log scale;
filled = total, hatched = active/recruiting subset). Panel b: strict-tier counts for the
99 gated vs 98 background surfaceome genes (log1p strip + median bars, MW p annotated).
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

a = json.load(open("results/clintrials_audit.json"))
tab = a["alias_table"]
genes = [t["gene"] for t in tab][::-1]
tot = [t["total"] for t in tab][::-1]
act = [t["active"] for t in tab][::-1]
cls = [t["class"] for t in tab][::-1]
colors = ["#b2182b" if c == "and_gate" else "#2166ac" for c in cls]

import csv
rows = list(csv.DictReader(open("results/clintrials_per_gene.csv")))
gated = [int(r["verified"]) for r in rows if r["class"] == "gated"]
bg = [int(r["verified"]) for r in rows if r["class"] == "background"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw={"width_ratios": [1.25, 1]})
y = np.arange(len(genes))
ax1.barh(y, [t + 0.5 for t in tot], color=colors, alpha=0.85, height=0.7)
ax1.barh(y, [t + 0.5 for t in act], color=colors, alpha=1.0, height=0.7,
         hatch="//", edgecolor="white", linewidth=0)
ax1.set_yticks(y); ax1.set_yticklabels(genes, fontsize=7)
ax1.set_xscale("log"); ax1.set_xlim(0.5, 3000)
ax1.set_xlabel("ClinicalTrials.gov targeted-therapy trials (alias tier)", fontsize=7)
ax1.tick_params(labelsize=7)
caps = [t["capped"] for t in tab][::-1]
for yi, t, cp in zip(y, tot, caps):
    ax1.text(t + 0.7 if t > 0 else 0.6, yi, ("\u2265" if cp else "") + str(t), va="center", fontsize=6.5)

rng = np.random.default_rng(7)
for i, (vals, lab) in enumerate([(gated, f"gated (n={len(gated)})"), (bg, f"background (n={len(bg)})")]):
    ax2.scatter(rng.normal(i, 0.06, len(vals)), np.log1p(vals), s=6, alpha=0.45,
                color="#b2182b" if i == 0 else "#2166ac")
    med = np.median(np.log1p(vals))
    ax2.hlines(med, i - 0.25, i + 0.25, color="black", lw=1.5)
ax2.set_xticks([0, 1]); ax2.set_xticklabels(["gated", "background"], fontsize=7)
ax2.set_ylabel("log1p curated strict-tier trials", fontsize=7)
p = a["stats"]["verified_mw"]["p"]
ax2.set_title(f"curated strict tier, MW p = {p:.3g}", fontsize=7)
ax2.tick_params(labelsize=7)
fig.tight_layout()
fig.savefig("paper/clintrials_fig.pdf")
print("wrote paper/clintrials_fig.pdf")

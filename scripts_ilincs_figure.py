#!/usr/bin/env python3
"""Figure for the iLINCS visibility/connectivity audit -> paper/figs/fig_ilincs.pdf."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

S = json.load(open("results/ilincs_summary.json"))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.6))

lm, cg = S["landmark"], S["cgs"]
cats = ["L1000 landmark\npanel (978 genes)", "LINCS CGS knockdown\nsignatures"]
gv = [lm["gated"][0]/lm["gated"][1], cg["gated"][0]/cg["gated"][1]]
bv = [lm["background"][0]/lm["background"][1], cg["background"][0]/cg["background"][1]]
x = range(2); w = 0.36
ax1.bar([i-w/2 for i in x], gv, w, label="gated (n=99)", color="#c44e52")
ax1.bar([i+w/2 for i in x], bv, w, label="background (n=98)", color="#4c72b0")
for i, p in enumerate([lm["fisher_p"], cg["fisher_p"]]):
    ax1.text(i, max(gv[i], bv[i])+0.012, f"p={p:.3f}" if p >= 0.001 else f"p={p:.1e}", ha="center", fontsize=8)
ax1.set_xticks(list(x)); ax1.set_xticklabels(cats, fontsize=8)
ax1.set_ylabel("fraction of genes covered"); ax1.set_ylim(0, 0.42)
ax1.set_title("LINCS L1000 visibility of the gated surfaceome\n(AND-gate antigens: 0/7 landmark, 2/7 CGS)")
ax1.legend(fontsize=8, frameon=False)

genes = ["CA12", "SLC39A6", "ERBB2"]
cd = S["connectivity"]["per_gene"]
import csv
rows = list(csv.DictReader(open("results/ilincs_connectivity_rows.csv")))
data = [[float(r[g]) for r in rows] for g in genes]
parts = ax2.violinplot(data, positions=range(3), showmedians=True, widths=0.7)
for b in parts["bodies"]: b.set_facecolor("#4c72b0"); b.set_alpha(0.6)
for k in ("cmedians", "cbars", "cmins", "cmaxes"):
    if k in parts: parts[k].set_color("black")
ax2.axhline(0, color="grey", lw=0.6, ls="--")
ax2.set_xticks(range(3)); ax2.set_xticklabels([f"{g}-KD\nconsensus" for g in genes], fontsize=8)
ax2.set_ylabel("Spearman connectivity to\n198 GDSC compounds")
ax2.set_title("Connectivity-map pilot: positive control NOT recovered\n(EGFR/HER2 inhibitors 1/20 top ERBB2 mimics, p=1.0)")
fig.tight_layout()
fig.savefig("paper/figs/fig_ilincs.pdf")
print("figure written")

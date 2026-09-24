#!/usr/bin/env python3
"""Figure for the tumor-vs-adjacent window audit (committed JSON only)."""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

D = json.load(open("results/firebrowse_window.json"))
G = D["per_gene"]
ORDER = ["CLDN18", "MSLN", "CA9", "CA12", "CLDN6", "PSCA", "SLC39A6", "ERBB2", "MS4A1", "TNFRSF17"]
floor = D["controls"]["hk_adj_max_abs_delta"]
med = [G[g]["median_delta"] for g in ORDER]
lo = [min(v["delta_adj"] for v in G[g]["cells"].values()) for g in ORDER]
hi = [max(v["delta_adj"] for v in G[g]["cells"].values()) for g in ORDER]
cols = ["#c0392b" if g in ORDER[:7] else "#7f8c8d" for g in ORDER]
x = np.arange(len(ORDER))
fig, ax = plt.subplots(figsize=(8.6, 3.4))
ax.axhspan(-floor, floor, color="0.9", zorder=0)
ax.axhline(0, c="k", lw=0.8)
ax.bar(x, med, color=cols, zorder=3)
ax.errorbar(x, med, yerr=[np.array(med) - np.array(lo), np.array(hi) - np.array(med)],
            fmt="none", ecolor="k", elinewidth=0.8, capsize=2, zorder=4)
ax.text(len(ORDER) - 0.4, floor + 0.12, f"normalization noise floor (+/-{floor:.1f} log2, housekeeping residual)", fontsize=7, ha="right")
ax.set_xticks(x); ax.set_xticklabels(ORDER, rotation=45, ha="right", fontsize=8)
ax.set_ylabel("tumor - adjacent normal (log2 RSEM, housekeeping-adjusted)")
ax.set_title("Within-organ therapeutic window: only CA9 clears the noise floor", fontsize=9)
fig.tight_layout()
fig.savefig("paper/figs/fig_firebrowse_window.pdf")
print("wrote fig")

#!/usr/bin/env python3
"""Two-panel HPA figure: membrane annotation and IF-data coverage by gene group."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = json.load(open("results/hpa_audit.json"))
GROUPS = ["gated", "background", "references", "positive_controls"]
LBL = ["gated\n(n=%d)" % 99, "background\n(n=%d)" % 98, "approved\nreferences", "intracellular\ncontrols"]

def frac(entry):
    k, n = entry
    return (k / n) if n else 0.0

pm = [frac(A["predicted_membrane"][g]) if g in ("gated", "background") else
      frac(A["calibration_gates"]["G1_refs_predicted_membrane"]) if g == "references" else
      frac(A["calibration_gates"]["G2_controls_predicted_membrane"]) for g in GROUPS]
ifc = [frac(A["if_coverage"][g]) for g in GROUPS]
ifp = [frac(A["if_plasma_membrane"][g]) for g in GROUPS]

fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.4))
c_main, c_alt = "#2166ac", "#b2182b"
b = axes[0].bar(range(4), pm, color=[c_main, c_main, "#67a9cf", "#969696"], edgecolor="k", lw=0.5)
axes[0].set_ylabel("fraction predicted membrane")
axes[0].set_title("A  HPA protein-class membrane prediction")
for i, v in enumerate(pm):
    axes[0].text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=8)
axes[0].set_ylim(0, 1.08)
w = 0.38
axes[1].bar([i - w / 2 for i in range(4)], ifc, width=w, color="#67a9cf", edgecolor="k", lw=0.5, label="has HPA IF data")
axes[1].bar([i + w / 2 for i in range(4)], ifp, width=w, color=c_alt, edgecolor="k", lw=0.5, label="IF plasma membrane")
axes[1].set_ylabel("fraction of genes")
axes[1].set_title("B  IF subcellular coverage vs. plasma membrane")
axes[1].legend(fontsize=7, frameon=False)
for ax in axes:
    ax.set_xticks(range(4)); ax.set_xticklabels(LBL, fontsize=7.5)
    ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig("paper/hpa_fig.pdf")
print("wrote paper/hpa_fig.pdf")

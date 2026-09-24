#!/usr/bin/env python3
"""Figure for the epitope structural audit (committed results only)."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

a = json.load(open("results/epitope_audit.json"))
df = pd.read_csv("results/epitope_per_gene.csv")
ann = df[(df.found == 1) & (df.topology != "unannotated")]
fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.2))
ax = axs[0]
for i, grp in enumerate(("gated", "background", "reference")):
    s = ann[ann.group == grp].exp_cov_ecto.dropna()
    ax.scatter([i] * len(s), s, s=14, alpha=0.55, color=["#c0392b", "#7f8c8d", "#2980b9"][i])
    ax.plot([i - 0.15, i + 0.15], [s.median(), s.median()], c="k", lw=1.5)
ax.set_xticks(range(3)); ax.set_xticklabels(["gated", "background", "CAR-T ref"], fontsize=8)
ax.set_ylabel("experimental ectodomain coverage")
ax.set_title("A. Structure coverage", fontsize=9)
ax = axs[1]
ag = pd.DataFrame(a["and_gate"])
ax.barh(ag.gene, ag.n_ab_structures, color="#c0392b")
ax.set_xlabel("antibody co-structures (PDB)")
ax.set_title("B. AND-gate antigens", fontsize=9)
ax = axs[2]
g = ann[ann.group == "gated"]; b = ann[ann.group == "background"]
vals = [g.has_structure.mean(), b.has_structure.mean(), g.has_ab_structure.mean(), b.has_ab_structure.mean()]
ax.bar([0, 1, 2.5, 3.5], vals, color=["#c0392b", "#7f8c8d", "#c0392b", "#7f8c8d"])
ax.set_xticks([0.5, 3]); ax.set_xticklabels(["any ectodomain\nstructure", "antibody\nco-structure"], fontsize=8)
ax.set_ylabel("fraction of annotated genes"); ax.set_ylim(0, 1)
ax.set_title("C. Gated vs background (null)", fontsize=9)
plt.tight_layout(); plt.savefig("paper/figs/fig_epitope.pdf"); print("ok")

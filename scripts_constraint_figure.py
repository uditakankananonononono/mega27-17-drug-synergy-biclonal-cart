#!/usr/bin/env python3
"""Figure: germline LoF tolerance (gnomAD v4 LOEUF) vs DepMap CRISPR dependency."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

AND_GATE = ["CLDN18", "MSLN", "CA9", "CA12", "CLDN6", "PSCA", "SLC39A6"]
df = pd.read_csv("results/constraint_per_gene.csv").dropna(subset=["loeuf", "depmap_frac_dep"])
fig, ax = plt.subplots(figsize=(6.2, 4.4))
sty = {"background": ("#bbbbbb", "o", "background surfaceome"), "gated": ("#1f77b4", "o", "gated antigens"),
       "references": ("#d62728", "s", "CAR-T references"), "positive_controls": ("#2ca02c", "^", "essential controls")}
for g, (c, m, lab) in sty.items():
    d = df[df.group == g]
    ax.scatter(d.loeuf, d.depmap_frac_dep, c=c, marker=m, s=22 if g in ("gated", "background") else 46, label=lab, alpha=0.8)
for _, r in df.iterrows():
    if r.gene in AND_GATE or r.group in ("references", "positive_controls"):
        ax.annotate(r.gene, (r.loeuf, r.depmap_frac_dep), fontsize=6.5, xytext=(3, 2), textcoords="offset points")
ax.axvline(0.6, ls="--", c="k", lw=0.8)
ax.axhline(0.1, ls=":", c="k", lw=0.8)
ax.set_xlabel("LOEUF (gnomAD v4; right = LoF-tolerant)")
ax.set_ylabel("fraction of DepMap lines dependent")
ax.set_title("Antigen-loss cost: lower-right quadrant = free loss", fontsize=10)
ax.legend(fontsize=7, loc="upper right")
fig.tight_layout()
fig.savefig("paper/figs/constraint_escape.pdf")
print("wrote paper/figs/constraint_escape.pdf")

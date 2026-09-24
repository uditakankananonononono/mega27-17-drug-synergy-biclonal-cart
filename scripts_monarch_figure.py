#!/usr/bin/env python3
"""Monarch audit figure: (a) fraction of genes with >=1 causal gene-disease
association, gated vs background; (b) HPO phenotype-breadth (symlog), antigens ringed."""
import csv, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ANTS = {"CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"}


def main():
    rows = [r for r in csv.DictReader(open("results/monarch_per_gene.csv"))
            if r["resolved"] == "True"]
    a = json.load(open("results/monarch_audit.json"))
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.4))
    for xi, (grp, col) in enumerate([("gated", "#c0392b"), ("background", "#2c3e50")]):
        sub = [r for r in rows if r["group"] == grp]
        fr = sum(1 for r in sub if int(r["causal_total"]) >= 1) / len(sub)
        ax[0].bar(xi, fr, 0.5, color=col)
    ax[0].set_xticks([0, 1]); ax[0].set_xticklabels(["gated", "background"])
    ax[0].set_ylabel("fraction with >=1 causal disease")
    ax[0].set_title(f"(a) Mendelian causal-disease carriers (Fisher p={a['H1_causal_carriage']['fisher_p']:.3g})",
                    fontsize=9)
    rng = np.random.default_rng(7)
    for xi, (grp, col) in enumerate([("gated", "#c0392b"), ("background", "#2c3e50")]):
        vals = [int(r["pheno_total"]) for r in rows if r["group"] == grp]
        ax[1].boxplot([vals], positions=[xi], widths=0.45, showfliers=False,
                      medianprops=dict(color=col, lw=2), patch_artist=True,
                      boxprops=dict(facecolor="white", edgecolor=col))
        ax[1].scatter(np.full(len(vals), xi) + rng.uniform(-0.13, 0.13, len(vals)), vals,
                      s=4, alpha=0.3, color=col)
    for r in rows:
        if r["gene"] in ANTS:
            ax[1].scatter([0], [int(r["pheno_total"])], s=26, facecolor="none",
                          edgecolor="#c0392b", lw=1.4, zorder=5)
    ax[1].set_yscale("symlog")
    ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["gated", "background"])
    ax[1].set_ylabel("HPO phenotype annotations (symlog)")
    ax[1].set_title(f"(b) Phenotypic breadth (MW p={a['H2_pheno_breadth']['mw_p']:.3g}; antigens ringed)",
                    fontsize=9)
    fig.tight_layout()
    fig.savefig("paper/monarch_fig.pdf")
    print("wrote paper/monarch_fig.pdf")


if __name__ == "__main__":
    main()

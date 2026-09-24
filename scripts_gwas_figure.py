#!/usr/bin/env python3
"""GWAS audit figure: (a) fraction of genes with >=1 significant association,
any-trait vs cancer-trait, gated vs background; (b) significant-association
burden (symlog) with AND-gate antigens ringed."""
import csv, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ANTS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]

def main():
    rows = list(csv.DictReader(open("results/gwas_per_gene.csv")))
    a = json.load(open("results/gwas_audit.json"))
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.4))
    x = np.arange(2); w = 0.36
    for k, (grp, col) in enumerate([("gated", "#c0392b"), ("background", "#2c3e50")]):
        sub = [r for r in rows if r["group"] == grp]
        any_fr = sum(1 for r in sub if int(r["n_sig"]) >= 1) / len(sub)
        can_fr = sum(1 for r in sub if int(r["n_sig_cancer"]) >= 1) / len(sub)
        ax[0].bar(x + (k - 0.5) * w, [any_fr, can_fr], w, color=col, label=grp)
    ax[0].set_xticks(x); ax[0].set_xticklabels(["any trait", "cancer trait"])
    ax[0].set_ylabel("fraction with >=1 sig association")
    ax[0].set_title(f"(a) Sig-association carriers (p={a['H1_any_sig']['fisher_p']:.3g} / {a['H2_cancer_sig']['fisher_p']:.3g})")
    ax[0].legend(fontsize=8, frameon=False)
    rng = np.random.default_rng(5)
    for xi, (grp, col) in enumerate([("gated", "#c0392b"), ("background", "#2c3e50")]):
        vals = [int(r["n_sig"]) for r in rows if r["group"] == grp]
        ax[1].boxplot([vals], positions=[xi], widths=0.45, showfliers=False,
                      medianprops=dict(color=col, lw=2), patch_artist=True,
                      boxprops=dict(facecolor="white", edgecolor=col))
        ax[1].scatter(np.full(len(vals), xi) + rng.uniform(-0.13, 0.13, len(vals)), vals,
                      s=4, alpha=0.3, color=col)
    for r in rows:
        if r["gene"] in ANTS:
            ax[1].scatter([0], [int(r["n_sig"])], s=26, facecolor="none",
                          edgecolor="#c0392b", lw=1.4, zorder=5)
    ax[1].set_yscale("symlog")
    ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["gated", "background"])
    ax[1].set_ylabel("significant associations per gene (symlog)")
    ax[1].set_title(f"(b) Burden (MW p={a['H3_burden']['mw_p']:.3g}; antigens ringed)")
    fig.tight_layout()
    fig.savefig("paper/gwas_fig.pdf")
    print("wrote paper/gwas_fig.pdf")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""BioPlex audit figure: (a) detection (presence) rate gated vs background per
network; (b) AP-MS degree distributions among detected genes, per network."""
import csv, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

def main():
    rows = list(csv.DictReader(open("results/bioplex_per_gene.csv")))
    a = json.load(open("results/bioplex_audit.json"))
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.4))
    nets = ["293T", "HCT116"]
    x = np.arange(2)
    w = 0.36
    for k, (grp, col) in enumerate([("gated", "#c0392b"), ("background", "#2c3e50")]):
        fr = []
        for n in nets:
            sub = [r for r in rows if r["group"] == grp]
            fr.append(sum(1 for r in sub if r[n + "_present"] == "True") / len(sub))
        ax[0].bar(x + (k - 0.5) * w, fr, w, color=col, label=grp)
    for i, n in enumerate(nets):
        p = a["H1_per_network"][n]["presence_fisher_p"]
        ax[0].text(i, 1.0, f"p={p:.2g}", ha="center", fontsize=8)
    ax[0].set_xticks(x); ax[0].set_xticklabels(["HEK293T 10K", "HCT116 5.5K"])
    ax[0].set_ylabel("fraction of genes detected")
    ax[0].set_ylim(0, 1.12)
    ax[0].set_title("(a) AP-MS detection rate")
    ax[0].legend(fontsize=8, frameon=False)
    rng = np.random.default_rng(3)
    pos = [0, 1]
    for i, n in enumerate(nets):
        for k, (grp, col) in enumerate([("gated", "#c0392b"), ("background", "#2c3e50")]):
            vals = [int(r[n + "_degree"]) for r in rows
                    if r["group"] == grp and r[n + "_degree"] not in (None, "", "None")]
            xi = i + (k - 0.5) * 0.42
            ax[1].boxplot([vals], positions=[xi], widths=0.3, showfliers=False,
                          medianprops=dict(color=col, lw=2), patch_artist=True,
                          boxprops=dict(facecolor="white", edgecolor=col))
            ax[1].scatter(np.full(len(vals), xi) + rng.uniform(-0.08, 0.08, len(vals)), vals,
                          s=3, alpha=0.3, color=col)
    ax[1].set_yscale("symlog")
    ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["HEK293T 10K", "HCT116 5.5K"])
    ax[1].set_ylabel("AP-MS degree (symlog)")
    p1 = a["H1_per_network"]["293T"]["degree_mw_p"]
    p2 = a["H1_per_network"]["HCT116"]["degree_mw_p"]
    ax[1].set_title(f"(b) Degree among detected (MW p={p1:.3g} / {p2:.3g})")
    fig.tight_layout()
    fig.savefig("paper/bioplex_fig.pdf")
    print("wrote paper/bioplex_fig.pdf")

if __name__ == "__main__":
    main()

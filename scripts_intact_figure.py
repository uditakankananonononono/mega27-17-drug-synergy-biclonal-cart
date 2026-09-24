#!/usr/bin/env python3
"""IntAct audit figure: (a) curated interaction degree gated vs background;
(b) 7x7 AND-gate antigen pairwise physical-interaction matrix."""
import csv, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ANTIGENS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]

def main():
    rows = list(csv.DictReader(open("results/intact_per_gene.csv")))
    audit = json.load(open("results/intact_audit.json"))
    gd = [int(r["count"] or 0) for r in rows if r["group"] == "gated" and r["count"] not in (None, "None", "")]
    bd = [int(r["count"] or 0) for r in rows if r["group"] == "background" and r["count"] not in (None, "None", "")]
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.6))
    rng = np.random.default_rng(7)
    for xi, (vals, lab, col) in enumerate([(gd, "gated", "#c0392b"), (bd, "background", "#2c3e50")]):
        bp = ax[0].boxplot([vals], positions=[xi], widths=0.45, showfliers=False,
                           medianprops=dict(color=col, lw=2), patch_artist=True,
                           boxprops=dict(facecolor="white", edgecolor=col))
        ax[0].scatter(np.full(len(vals), xi) + rng.uniform(-0.13, 0.13, len(vals)), vals,
                      s=4, alpha=0.35, color=col)
    ax[0].set_yscale("symlog")
    ax[0].set_xticks([0, 1]); ax[0].set_xticklabels(["gated", "background"])
    ax[0].set_ylabel("IntAct curated interactions (count)")
    ax[0].set_title(f"(a) Interaction degree (MW p={audit['H2_degree']['mw_p']})")
    mat = np.zeros((7, 7))
    hits = {tuple(h.split("|")) for h in audit["H1_pairwise_antigen_complex"]["edges_detail"]}
    for i, a in enumerate(ANTIGENS):
        for j, b in enumerate(ANTIGENS):
            if tuple(sorted((a, b))) in hits:
                mat[i, j] = 1
            if i == j:
                mat[i, j] = np.nan
    masked = np.ma.masked_invalid(mat)
    cmap = matplotlib.colors.ListedColormap(["#f5f6fa", "#c0392b"])
    cmap.set_bad("#dddddd")
    ax[1].imshow(masked, cmap=cmap, vmin=0, vmax=1)
    ax[1].set_xticks(range(7)); ax[1].set_xticklabels(ANTIGENS, rotation=45, ha="right", fontsize=8)
    ax[1].set_yticks(range(7)); ax[1].set_yticklabels(ANTIGENS, fontsize=8)
    ax[1].set_title(f"(b) Antigen-pair curated physical edges ({len(hits)}/21)")
    for i in range(7):
        for j in range(7):
            if i != j:
                ax[1].text(j, i, "Y" if mat[i, j] == 1 else "-", ha="center", va="center",
                           fontsize=8, color="white" if mat[i, j] == 1 else "#999999")
    fig.tight_layout()
    fig.savefig("paper/intact_fig.pdf")
    print("wrote paper/intact_fig.pdf")

if __name__ == "__main__":
    main()

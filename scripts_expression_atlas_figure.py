#!/usr/bin/env python3
"""Expression Atlas figure: (a) tau by group, antigens ringed; (b) Body Map
breadth vs ProteomicsDB normal-tissue breadth."""
import csv, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ANTS = {"CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"}
COLS = [("gated", "#c0392b"), ("background", "#2c3e50")]


def main():
    rows = [r for r in csv.DictReader(open("results/expression_atlas_per_gene.csv"))
            if r["matched"] == "True" and r["tau"] not in ("", "None")]
    a = json.load(open("results/expression_atlas_audit.json"))
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.4))
    rng = np.random.default_rng(7)
    for xi, (grp, col) in enumerate(COLS):
        v = [float(r["tau"]) for r in rows if r["group"] == grp]
        ax[0].boxplot([v], positions=[xi], widths=0.45, showfliers=False,
                      medianprops=dict(color=col, lw=2), patch_artist=True,
                      boxprops=dict(facecolor="white", edgecolor=col))
        ax[0].scatter(np.full(len(v), xi) + rng.uniform(-0.13, 0.13, len(v)), v,
                      s=4, alpha=0.35, color=col)
    for r in rows:
        if r["gene"] in ANTS:
            ax[0].scatter([0], [float(r["tau"])], s=26, facecolor="none",
                          edgecolor="#c0392b", lw=1.4, zorder=5)
    ax[0].set_xticks([0, 1]); ax[0].set_xticklabels(["gated", "background"])
    ax[0].set_ylabel("tissue specificity tau (Body Map)")
    ax[0].set_title("(a) RNA specificity (one-sided p=%.2g; antigens ringed)"
                    % a["H1_tau"]["p"], fontsize=9)
    pdb = {r["gene"]: r for r in csv.DictReader(open("results/proteomicsdb_per_gene.csv"))
           if r["mapped"] == "True"}
    for grp, col in COLS:
        pts = [(int(r["n_tissues_tpm1"]), int(pdb[r["gene"]]["n_normal_tissues"]))
               for r in rows if r["group"] == grp and r["gene"] in pdb]
        ax[1].scatter([x + rng.uniform(-0.2, 0.2) for x, _ in pts], [y for _, y in pts],
                      s=7, alpha=0.5, color=col, label=grp)
    ax[1].set_xlabel("Body Map tissues with TPM >= 1 (of 16)")
    ax[1].set_ylabel("ProteomicsDB normal tissues detected")
    ax[1].set_title("(b) RNA vs protein breadth (Spearman rho=%.2f)"
                    % a["H2_cross_platform"]["spearman_rho"], fontsize=9)
    ax[1].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig("paper/expression_atlas_fig.pdf")
    print("wrote paper/expression_atlas_fig.pdf")


if __name__ == "__main__":
    main()

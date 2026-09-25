#!/usr/bin/env python3
"""ProteomicsDB figure: (a) fraction of genes detected by MS in any sample;
(b) distinct normal human tissues with detection, antigens ringed."""
import csv, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ANTS = {"CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"}
COLS = [("gated", "#c0392b"), ("background", "#2c3e50")]


def main():
    rows = [r for r in csv.DictReader(open("results/proteomicsdb_per_gene.csv"))
            if r["mapped"] == "True"]
    a = json.load(open("results/proteomicsdb_audit.json"))
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.4))
    for xi, (grp, col) in enumerate(COLS):
        sub = [r for r in rows if r["group"] == grp]
        frac = sum(r["detected_any"] == "True" for r in sub) / len(sub)
        ax[0].bar(xi, frac, 0.5, color=col)
    ax[0].set_xticks([0, 1]); ax[0].set_xticklabels(["gated", "background"])
    ax[0].set_ylabel("fraction detected in any sample")
    ax[0].set_title("(a) MS detection (Fisher p=%.2g)"
                    % a["H2_detected_any"]["fisher_p"], fontsize=9)
    rng = np.random.default_rng(7)
    for xi, (grp, col) in enumerate(COLS):
        v = [int(r["n_normal_tissues"]) for r in rows if r["group"] == grp]
        ax[1].boxplot([v], positions=[xi], widths=0.45, showfliers=False,
                      medianprops=dict(color=col, lw=2), patch_artist=True,
                      boxprops=dict(facecolor="white", edgecolor=col))
        ax[1].scatter(np.full(len(v), xi) + rng.uniform(-0.13, 0.13, len(v)), v,
                      s=4, alpha=0.35, color=col)
    for r in rows:
        if r["gene"] in ANTS:
            ax[1].scatter([0], [int(r["n_normal_tissues"])], s=26, facecolor="none",
                          edgecolor="#c0392b", lw=1.4, zorder=5)
    ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["gated", "background"])
    ax[1].set_ylabel("normal human tissues with detection")
    ax[1].set_title("(b) Normal-tissue breadth (one-sided 'less' p=%.2g; antigens ringed)"
                    % a["H1_normal_breadth"]["p"], fontsize=9)
    fig.tight_layout()
    fig.savefig("paper/proteomicsdb_fig.pdf")
    print("wrote paper/proteomicsdb_fig.pdf")


if __name__ == "__main__":
    main()

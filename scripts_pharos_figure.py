#!/usr/bin/env python3
"""Pharos audit figure: (a) TDL mix gated vs background (stacked fractions);
(b) TIN-X novelty (log10) gated vs background with AND-gate antigens marked."""
import csv, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ANTIGENS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]
TDLS = ["Tclin", "Tchem", "Tbio", "Tdark"]
COLS = {"Tclin": "#c0392b", "Tchem": "#e67e22", "Tbio": "#2980b9", "Tdark": "#7f8c8d"}

def main():
    rows = [r for r in csv.DictReader(open("results/pharos_per_gene.csv")) if r["status"] == "ok"]
    audit = json.load(open("results/pharos_audit.json"))
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.4))
    for xi, grp in enumerate(["gated", "background"]):
        sub = [r for r in rows if r["group"] == grp]
        n = len(sub)
        bottom = 0.0
        for t in TDLS:
            frac = sum(1 for r in sub if r["tdl"] == t) / n
            ax[0].bar(xi, frac, bottom=bottom, color=COLS[t], width=0.55,
                      label=t if xi == 0 else None)
            bottom += frac
    ax[0].set_xticks([0, 1]); ax[0].set_xticklabels(["gated", "background"])
    ax[0].set_ylabel("fraction of genes")
    h1 = audit["H1_tclin"]
    ax[0].set_title(f"(a) Target Development Level (Tclin Fisher p={h1['fisher_p']:.3g})")
    ax[0].legend(fontsize=7, frameon=False)
    rng = np.random.default_rng(11)
    for xi, (grp, col) in enumerate([("gated", "#c0392b"), ("background", "#2c3e50")]):
        vals = [np.log10(float(r["novelty"])) for r in rows
                if r["group"] == grp and r["novelty"] not in (None, "", "None") and float(r["novelty"]) > 0]
        ax[1].boxplot([vals], positions=[xi], widths=0.45, showfliers=False,
                      medianprops=dict(color=col, lw=2), patch_artist=True,
                      boxprops=dict(facecolor="white", edgecolor=col))
        ax[1].scatter(np.full(len(vals), xi) + rng.uniform(-0.13, 0.13, len(vals)), vals,
                      s=4, alpha=0.3, color=col)
    for r in rows:
        if r["gene"] in ANTIGENS and r["novelty"] not in (None, "", "None") and float(r["novelty"]) > 0:
            ax[1].scatter([0], [np.log10(float(r["novelty"]))], s=26, facecolor="none",
                          edgecolor="#c0392b", lw=1.4, zorder=5)
    ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["gated", "background"])
    ax[1].set_ylabel("log10 TIN-X novelty (higher = less studied)")
    h2 = audit["H2_novelty"]
    ax[1].set_title(f"(b) Novelty (MW p={h2['mw_p']:.3g}; antigens ringed)")
    fig.tight_layout()
    fig.savefig("paper/pharos_fig.pdf")
    print("wrote paper/pharos_fig.pdf")

if __name__ == "__main__":
    main()

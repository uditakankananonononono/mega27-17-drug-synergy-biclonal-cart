#!/usr/bin/env python3
"""Harmonizome audit figure: per-gene association counts (symlog) for gated vs
background in curated, systematic-expression and systematic-non-expression
datasets; AND-gate antigens ringed in the curated panel."""
import csv, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scripts_harmonizome_audit as H

ANTS = set(H.ANTS)


def per_gene_split():
    groups = H.gene_set()
    out = {}
    for g, grp in groups.items():
        if grp not in ("gated", "background"):
            continue
        p = os.path.join(H.GDIR, g + ".json")
        if not os.path.exists(p):
            continue
        rec = json.load(open(p))
        if "_error" in rec or not rec.get("n_associations"):
            continue
        c = {"curated": 0, "expr": 0, "nonexpr": 0}
        for ds, (up, dn) in rec["per_dataset"].items():
            k = H.classify(ds)
            if k == "curated":
                c["curated"] += up + dn
            elif k == "systematic":
                c["expr" if H.is_expression(ds) else "nonexpr"] += up + dn
        out[g] = (grp, c)
    return out


def main():
    a = json.load(open("results/harmonizome_audit.json"))
    d = per_gene_split()
    panels = [("curated", "(a) curated/literature", a["tests"]["curated"]["p"]),
              ("expr", "(b) systematic: expression", a["H2_control"]["expression"]["p"]),
              ("nonexpr", "(c) systematic: non-expression", a["H2_control"]["nonexpression"]["p"])]
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.4))
    rng = np.random.default_rng(7)
    for i, (key, title, p) in enumerate(panels):
        for xi, (grp, col) in enumerate([("gated", "#c0392b"), ("background", "#2c3e50")]):
            vals = [c[key] for g, (gg, c) in d.items() if gg == grp]
            ax[i].boxplot([vals], positions=[xi], widths=0.45, showfliers=False,
                          medianprops=dict(color=col, lw=2), patch_artist=True,
                          boxprops=dict(facecolor="white", edgecolor=col))
            ax[i].scatter(np.full(len(vals), xi) + rng.uniform(-0.13, 0.13, len(vals)),
                          vals, s=4, alpha=0.3, color=col)
        if key == "curated":
            for g, (gg, c) in d.items():
                if g in ANTS:
                    ax[i].scatter([0], [c[key]], s=26, facecolor="none",
                                  edgecolor="#c0392b", lw=1.4, zorder=5)
        ax[i].set_yscale("symlog")
        ax[i].set_xticks([0, 1]); ax[i].set_xticklabels(["gated", "background"])
        ax[i].set_title(f"{title}\nMW p={p:.2g}", fontsize=9)
    ax[0].set_ylabel("Harmonizome associations per gene")
    fig.tight_layout()
    fig.savefig("paper/harmonizome_fig.pdf")
    print("wrote paper/harmonizome_fig.pdf")


if __name__ == "__main__":
    main()

"""DepMap CRISPR dependency audit of gated CAR-T antigens (antigen-loss escape risk).

A CAR-T antigen the tumor does not need can be lost under immune pressure.
Question: are the 99 gated antigens more essential than a random surfaceome
background? Reads committed results/depmap_subset.csv (DepMap 24Q4 Chronos
gene effect; 0 = no effect, -1 = median common-essential) and
results/depmap_gene_sets.json. Metrics per gene: median effect across models and
fraction of models with effect < -0.5. Mann-Whitney gated vs background.
"""
import json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

THR = -0.5


def gene_metrics(df, genes):
    out = {}
    for g in genes:
        if g in df.columns:
            v = df[g].dropna().to_numpy(float)
            out[g] = {"median": float(np.median(v)), "frac_dep": float((v < THR).mean()),
                      "n": int(len(v))}
    return out


def main():
    df = pd.read_csv("results/depmap_subset.csv")
    sets = json.load(open("results/depmap_gene_sets.json"))
    res = {k: gene_metrics(df, sets[k]) for k in ["gated", "background", "references", "positive_controls"]}
    gm = np.array([m["median"] for m in res["gated"].values()])
    bm = np.array([m["median"] for m in res["background"].values()])
    gf = np.array([m["frac_dep"] for m in res["gated"].values()])
    bf = np.array([m["frac_dep"] for m in res["background"].values()])
    out = {
        "design": __doc__.strip().splitlines()[0],
        "source": "DepMap Public 24Q4 CRISPRGeneEffect.csv (figshare 27993248, file 51064667) + Model.csv (51065297)",
        "threshold": THR, "n_models": int(len(df)),
        "n_gated": len(gm), "n_background": len(bm),
        "gated_median_of_medians": float(np.median(gm)),
        "background_median_of_medians": float(np.median(bm)),
        "mw_p_median_effect_gated_lt_background": float(mannwhitneyu(gm, bm, alternative="less").pvalue),
        "mw_p_median_effect_two_sided": float(mannwhitneyu(gm, bm).pvalue),
        "gated_n_any_dependency_ge_10pct": int((gf >= 0.10).sum()),
        "background_n_any_dependency_ge_10pct": int((bf >= 0.10).sum()),
        "gated_top_dependencies": sorted(({"gene": g, **m} for g, m in res["gated"].items()),
                                         key=lambda r: r["median"])[:8],
        "references": res["references"], "positive_controls": res["positive_controls"],
        "missing_in_depmap": sets["missing_in_depmap"],
        "per_gene": res,
    }
    json.dump(out, open("results/depmap_dependency.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k not in ("per_gene",)}, indent=1)[:3500])


if __name__ == "__main__":
    main()

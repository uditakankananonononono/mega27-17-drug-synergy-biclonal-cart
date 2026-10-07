#!/usr/bin/env python3
"""P4-CC1 (prespec results/catalog_annotation_adjusted_prespec.json): annotation-depth-stratified
permutation test of the catalog->genuine enrichment. Inputs are committed CSV/JSON only."""
import json, numpy as np, pandas as pd
SEED, NPERM = 20261007, 20000
rng = np.random.default_rng(SEED)
c = pd.read_csv("results/clintrials_per_gene.csv")
p = pd.read_csv("results/pharos_per_gene.csv")[["gene", "publication_count", "generif_count", "status"]]
cat = set(json.load(open("results/catalog_comparator.json"))["catalog_gene_verdicts"])
m = c.merge(p, on="gene", how="left")
bad = m[(m.status != "ok") | m.publication_count.isna()]
m = m.drop(bad.index).copy()
m["cat"] = m.gene.isin(cat); m["gen"] = m.verdict == "genuine"

def perm_test(df, col):
    df = df.copy(); df["s"] = pd.qcut(np.log1p(df[col]), 5, labels=False, duplicates="drop")
    obs = int((df.cat & df.gen).sum())
    groups = [(g.cat.values.copy(), g.gen.values) for _, g in df.groupby("s")]
    null = np.empty(NPERM, int)
    for i in range(NPERM):
        null[i] = sum(int(rng.permutation(cv)[gv].sum()) if len(cv) else 0 for cv, gv in groups)
    return obs, null, df

def mh_or(df):
    num = den = 0.0
    for _, g in df.groupby("s"):
        n = len(g); a = (g.cat & g.gen).sum(); b = (g.cat & ~g.gen).sum(); cc = (~g.cat & g.gen).sum(); d = (~g.cat & ~g.gen).sum()
        num += a * d / n; den += b * cc / n
    return float(num / den) if den > 0 else "undefined"

obs, null, df = perm_test(m, "publication_count")
pval = float((1 + (null >= obs).sum()) / (1 + NPERM))
obs2, null2, df2 = perm_test(m, "generif_count")
p2 = float((1 + (null2 >= obs2).sum()) / (1 + NPERM))
verdict = "SURVIVES" if pval < 0.05 else ("WEAKENED" if pval < 0.20 else "EXPLAINED_BY_ANNOTATION")
res = {"prespec": "results/catalog_annotation_adjusted_prespec.json", "seed": SEED, "n_perm": NPERM,
       "excluded_genes": bad.gene.tolist(), "n_retained": int(len(m)), "n_catalog_retained": int(m.cat.sum()),
       "observed_catalog_genuine": obs, "stratified_null_mean": float(null.mean()), "stratified_null_p95": float(np.percentile(null, 95)),
       "stratified_null_max": int(null.max()), "p_publication_strata": pval, "mh_or_publication_strata": mh_or(df),
       "secondary_p_generif_strata": p2, "secondary_observed_generif": obs2, "secondary_null_mean_generif": float(null2.mean()),
       "pubcount_median_catalog": float(m[m.cat].publication_count.median()), "pubcount_median_noncatalog": float(m[~m.cat].publication_count.median()),
       "unstratified_fisher_reference": "see results/catalog_comparator.json", "verdict": verdict}
json.dump(res, open("results/catalog_annotation_adjusted.json", "w"), indent=1)
print(json.dumps(res, indent=1))

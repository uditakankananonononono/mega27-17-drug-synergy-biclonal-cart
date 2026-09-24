"""CAR-T antigen ranking from Human Protein Atlas real data.

Score = tumor expression (max nTPM across patients of a cancer) penalized by
expression in vital normal tissues (heart, liver, lung, kidney, brain).
Therapeutic window: tumor_nTPM / (1 + vital_normal_nTPM). Known benchmarks:
CD19, BCMA, HER2/ERBB2, MSLN, GPC3, GD2 (GD2 is a glycolipid, excluded),
which lets us validate the ranking against clinical reality.
"""
import pandas as pd

VITAL = {"heart muscle", "liver", "lung", "kidney", "cerebral cortex"}
KNOWN_TARGETS = ["CD19", "MS4A1", "BCMA", "TNFRSF17", "ERBB2", "MSLN", "GPC3", "CD22", "EGFR", "FOLR1"]


def load(data_dir="data"):
    norm = pd.read_csv(f"{data_dir}/rna_tissue_consensus.tsv", sep="\t")
    path = pd.read_csv(f"{data_dir}/pathology.tsv", sep="\t")
    return norm, path


def vital_normal_expression(norm):
    v = norm[norm["Tissue"].isin(VITAL)]
    return v.groupby("Gene name")["nTPM"].max()


def tumor_expression(path):
    """Weighted detection score per gene x cancer from HPA patient counts:
    (3*High + 2*Medium + 1*Low) / n_patients."""
    for col in ("High", "Medium", "Low", "Not detected"):
        path[col] = pd.to_numeric(path[col], errors="coerce").fillna(0)
    n = path[["High", "Medium", "Low", "Not detected"]].sum(axis=1)
    path = path.assign(_score=(3 * path["High"] + 2 * path["Medium"] + path["Low"]) / n.clip(lower=1))
    return path.groupby(["Gene name", "Cancer"])["_score"].max().unstack(fill_value=0.0)


def rank(norm, path, cancer=None, top=30):
    vital = vital_normal_expression(norm)
    tum = tumor_expression(path)
    if cancer:
        tum = tum[[cancer]] if cancer in tum.columns else tum.iloc[:, :1]
    score = tum.div(1 + vital, axis=0).dropna(how="all")
    ranked = score.max(axis=1).sort_values(ascending=False)
    return ranked.head(top), vital


def surface_genes(loc):
    """Genes with plasma-membrane / extracellular annotation (CAR-T addressable)."""
    mask = loc["Main location"].fillna("").str.contains("Plasma Membrane", case=False)          | loc["Additional location"].fillna("").str.contains("Plasma Membrane", case=False)
    return set(loc.loc[mask, "Gene name"])

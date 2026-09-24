"""Per-cancer CAR-T antigen windows with patient coverage, surface-filtered."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "src"))
import pandas as pd
from cart import antigen_rank as ar

def per_cancer_table(norm, path, loc):
    surf = ar.surface_genes(loc)
    vital = ar.vital_normal_expression(norm)
    for col in ("High", "Medium", "Low", "Not detected"):
        path[col] = pd.to_numeric(path[col], errors="coerce").fillna(0)
    n = path[["High", "Medium", "Low", "Not detected"]].sum(axis=1)
    path = path.assign(score=(3*path["High"] + 2*path["Medium"] + path["Low"]) / n.clip(lower=1),
                       frac_detected=(path["High"] + path["Medium"] + path["Low"]) / n.clip(lower=1))
    g = path.groupby(["Gene name", "Cancer"]).agg(score=("score", "max"), coverage=("frac_detected", "max")).reset_index()
    g = g[g["Gene name"].isin(surf)]
    g["window"] = g["score"] / (1 + g["Gene name"].map(vital).fillna(0))
    return g

if __name__ == "__main__":
    norm, path = ar.load()
    loc = pd.read_csv("data/subcellular_location.tsv", sep="\t")
    t = per_cancer_table(norm, path, loc)
    for cancer in ["pancreatic cancer", "glioma", "ovarian cancer", "lung cancer"]:
        sub = t[t["Cancer"] == cancer].nlargest(8, "window")
        print(f"\n{cancer}:")
        for _, r in sub.iterrows():
            print(f"  {r['Gene name']:10s} window={r['window']:.2f} coverage={r['coverage']:.2f}")
    # validation: where do known targets land in their clinical indication?
    checks = [("MSLN", "pancreatic cancer"), ("ERBB2", "breast cancer"), ("GPC3", "liver cancer"), ("CD19", "lymphoma")]
    print("\nknown target in its clinical indication (rank within that cancer):")
    for gene, cancer in checks:
        sub = t[t["Cancer"] == cancer].sort_values("window", ascending=False).reset_index(drop=True)
        hit = sub[sub["Gene name"] == gene]
        if len(hit):
            r = hit.index[0]
            print(f"  {gene:6s} in {cancer:18s} rank #{r+1} of {len(sub)}  window={hit.iloc[0]['window']:.2f} coverage={hit.iloc[0]['coverage']:.2f}")

"""Extract the DepMap 24Q4 CRISPRGeneEffect columns needed for the dependency audit.

Source: DepMap Public 24Q4 (figshare article 27993248): CRISPRGeneEffect.csv
(file 51064667) and Model.csv (file 51065297). Writes a committed subset:
results/depmap_subset.csv (ModelID, OncotreeLineage, one column per gene symbol).
Genes: the 99 gated antigens + surfaceome background from results/ot_tractability.json,
the 4 CAR-T references, and 4 pan-essential positive controls.
"""
import json, sys
import pandas as pd

CRISPR = sys.argv[1] if len(sys.argv) > 1 else "/tmp/CRISPRGeneEffect.csv"
MODEL = sys.argv[2] if len(sys.argv) > 2 else "/tmp/Model.csv"
POS_CTRL = ["POLR2A", "RPS3", "PCNA", "PSMA1"]


def gene_sets(ot):
    gated = sorted(v["symbol"] for v in ot["gated"]["rows"].values())
    background = sorted(ot["background"]["rows"].keys())
    refs = sorted(ot["references"].keys())
    return gated, background, refs


def main():
    ot = json.load(open("results/ot_tractability.json"))
    gated, background, refs = gene_sets(ot)
    want = set(gated) | set(background) | set(refs) | set(POS_CTRL)
    header = pd.read_csv(CRISPR, nrows=0).columns
    sym = {c: c.split(" (")[0] for c in header[1:]}
    cols = [header[0]] + [c for c in header[1:] if sym[c] in want]
    df = pd.read_csv(CRISPR, usecols=cols).rename(columns={header[0]: "ModelID", **sym})
    model = pd.read_csv(MODEL, usecols=["ModelID", "OncotreeLineage"])
    df = model.merge(df, on="ModelID", how="inner")
    df.to_csv("results/depmap_subset.csv", index=False)
    missing = sorted(want - set(df.columns))
    json.dump({"gated": gated, "background": background, "references": refs,
               "positive_controls": POS_CTRL, "missing_in_depmap": missing,
               "n_models": int(len(df))}, open("results/depmap_gene_sets.json", "w"), indent=1)
    print(len(df), "models;", len(cols) - 1, "genes; missing", missing)


if __name__ == "__main__":
    main()

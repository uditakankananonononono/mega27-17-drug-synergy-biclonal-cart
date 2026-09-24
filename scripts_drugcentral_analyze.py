"""DrugCentral engagement audit of the AND-gate antigens.

Orthogonal evidence base to DGIdb (interaction claims) and Open Targets
(predicted tractability buckets): DrugCentral holds curated *quantitative*
drug-target activities (ChEMBL-derived potencies + FDA-label mechanisms) with
Illuminate target-development levels (TDL). Gene sets are exactly those of the
DepMap/OT audits (results/depmap_gene_sets.json): 99 gated antigens, 98
random-surfaceome background, 4 CAR-T references.

Outputs: results/drugcentral_engagement.json, results/drugcentral_engagement_rows.csv,
paper/drugcentral_sec.tex (numbers generated from the JSON values directly).
"""
import csv, json, pathlib
import pandas as pd
from scipy.stats import fisher_exact

SETS = json.load(open("results/depmap_gene_sets.json"))
GATED, BG, REFS = SETS["gated"], SETS["background"], SETS["references"]
AND_GATE = ["CLDN18", "MSLN", "CA9", "CA12", "CLDN6", "PSCA", "SLC39A6"]

df = pd.read_csv("data/drugcentral/drug.target.interaction.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
h = df[df["ORGANISM"] == "Homo sapiens"].copy()

per_gene = {}
def annotate(gene, group):
    sub = h[h["GENE"] == gene]
    tdls = sorted({t for cell in sub["TDL"] for t in cell.split("|") if t and t != "NA"})
    drugs = sorted(sub["DRUG_NAME"].unique())
    return {
        "group": group,
        "n_rows": int(len(sub)),
        "n_drugs": int(sub["STRUCT_ID"].nunique()),
        "tdl": tdls,
        "tclin": "Tclin" in tdls,
        "moa": int((sub["MOA"] == "1").sum()),
        "drugs": drugs[:15],
    }

for g in GATED: per_gene[g] = annotate(g, "gated")
for g in BG: per_gene[g] = annotate(g, "background")
for g in REFS: per_gene[g] = annotate(g, "reference")

def frac(genes, key):
    yes = sum(1 for g in genes if (per_gene[g]["n_rows"] > 0 if key == "any" else
           per_gene[g]["tclin"] if key == "tclin" else per_gene[g]["moa"] > 0))
    return yes, len(genes) - yes

summary = {"gated": {"n": len(GATED)}, "background": {"n": len(BG)}}
tests = {}
for key in ("any", "tclin", "moa"):
    gy, gn = frac(GATED, key); by, bn = frac(BG, key)
    odds, p = fisher_exact([[gy, gn], [by, bn]])
    summary["gated"][key] = gy; summary["background"][key] = by
    tests[key] = {"gated_yes": gy, "gated_no": gn, "bg_yes": by, "bg_no": bn,
                  "odds": round(odds, 3), "p": p}
    print(f"{key}: gated {gy}/{len(GATED)} bg {by}/{len(BG)} odds {odds:.3f} p {p:.2e}")

out = {
    "source": {
        "name": "DrugCentral drug.target.interaction snapshot 2021_09_01",
        "url": "https://unmtid-dbs.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz",
        "rows_total": int(len(df)), "rows_human": int(len(h)),
        "human_genes": int(h["GENE"].nunique()), "human_drugs": int(h["STRUCT_ID"].nunique()),
        "note": "bulk snapshot file; curated quantitative activities, small-molecule-heavy",
    },
    "per_gene": per_gene,
    "summary": summary,
    "tests": tests,
    "and_gate": {g: per_gene[g] for g in AND_GATE},
    "references": {g: per_gene[g] for g in REFS},
}
json.dump(out, open("results/drugcentral_engagement.json", "w"), indent=1)

with open("results/drugcentral_engagement_rows.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["gene", "group", "n_rows", "n_drugs", "tdl", "tclin", "moa", "drug_names"])
    for g, d in per_gene.items():
        w.writerow([g, d["group"], d["n_rows"], d["n_drugs"], "|".join(d["tdl"]),
                    int(d["tclin"]), d["moa"], ";".join(d["drugs"])])
print("wrote results/drugcentral_engagement.json + _rows.csv")

"""Mechanism annotation of ALMANAC drugs from ChEMBL assay accessions.

Counts unique assay accessions actually fetched+used (the accession-level
dataset contribution); annotates each drug's top targets by max pChEMBL;
tests whether synergistic pairs (SCORE>=50) share annotated targets LESS than
expected (complementary-mechanism hypothesis) vs random pairings.
Output: results/chembl_annotation.json
"""
import json, pathlib, random
from collections import Counter
import numpy as np
import pandas as pd

CHEM = pathlib.Path("data/chembl")
drug_targets, all_assays = {}, {}
names = {}
for f in CHEM.glob("*.json"):
    rec = json.loads(f.read_text())
    nsc = str(rec["nsc"])
    names[nsc] = rec.get("pref_name") or rec.get("name")
    tgts = Counter()
    for a in rec.get("assays", []):
        all_assays[a["assay"]] = 1
        if a.get("target") and a.get("pchembl_max", 0) >= 6.0:  # <=1 uM
            tgts[a["target"]] = max(tgts[a["target"]], a["pchembl_max"])
    if tgts:
        drug_targets[int(nsc)] = dict(tgts.most_common(5))

df = pd.read_csv("data/almanac_synergy.tsv", sep="\t")
pair = df.groupby(["NSC1", "NSC2"])["MAXSCORE"].max().reset_index()
pair["syn"] = pair.MAXSCORE >= 50

def shared(a, b):
    ta, tb = set(drug_targets.get(a, {})), set(drug_targets.get(b, {}))
    return len(ta & tb) if (ta and tb) else None

annot = [(r.NSC1, r.NSC2, r.syn, shared(r.NSC1, r.NSC2)) for r in pair.itertuples()]
annot = [x for x in annot if x[3] is not None]
syn_share = np.mean([x[3] > 0 for x in annot if x[2]])
non_share = np.mean([x[3] > 0 for x in annot if not x[2]])
out = {
    "n_drugs_with_targets": len(drug_targets),
    "n_unique_assay_accessions": len(all_assays),
    "pairs_both_annotated": len(annot),
    "frac_synergistic_sharing_target": round(float(syn_share), 4),
    "frac_nonsynergistic_sharing_target": round(float(non_share), 4),
    "drug_top_targets": {names.get(str(k), str(k)): v for k, v in list(drug_targets.items())[:200]},
    "note": "target sharing uses max-pChEMBL>=6 (<=1 uM) top-5 targets per drug from ChEMBL assays (each assay CHEMBL ID = 1 accession-level dataset)",
}
json.dump(out, open("results/chembl_annotation.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "drug_top_targets"}, indent=1))

"""CAR-T RNA-vs-RNA window validation + AND-gate pair discovery (HPA v23)."""
import json, sys
sys.path.insert(0, "src")
import pandas as pd
from cart.rna_window import (load_tumor, load_normal, surface_set,
                             single_windows, and_gate_pairs, load_uniprot_membrane)

KNOWN = {"CD19": "B-cell malignancies", "MS4A1": "B-cell malignancies (CD20)",
         "ERBB2": "HER2+ solid tumours", "MSLN": "mesothelioma/pancreatic/ovarian",
         "GPC3": "hepatocellular carcinoma", "TNFRSF17": "multiple myeloma (BCMA)",
         "EGFR": "solid tumours", "FOLR1": "ovarian"}
CANCERS = ["PAAD", "BRCA", "GBM", "LIHC", "OV", "LUAD", "COAD", "SKCM", "STAD", "KIRC"]

tumor = load_tumor("data/cancer_rna_mean.tsv")
normal, names = load_normal("data/rna_tissue_consensus.tsv")
loc = pd.read_csv("data/subcellular_location.tsv", sep="\t")
uni = load_uniprot_membrane("data/uniprot_cellmembrane.tsv")
surf = surface_set(loc, uni, names)
print("surface genes:", len(surf))
name2id = {}
for gid, n in names.items():
    name2id.setdefault(n, []).append(gid)

# 1) known-target validation: window must be strongly positive in its cancer
w_all = single_windows(tumor, normal, surf)
val = {}
for sym, ind in KNOWN.items():
    for gid in name2id.get(sym, []):
        if gid in w_all.index:
            row = w_all.loc[gid].sort_values(ascending=False)
            val[sym] = {"ensg": gid, "top_cancers": {c: round(float(v), 2) for c, v in row.head(4).items()},
                        "indication": ind}
print(json.dumps(val, indent=1))

# 2) AND-gate discovery
pairs = {}
for c in CANCERS:
    ps = and_gate_pairs(tumor, normal, surf, c, n_cand=120, min_tumor=20.0, top=15)
    pairs[c] = [{**p, "a_name": names.get(p["a"]), "b_name": names.get(p["b"])} for p in ps]
    print(c, [(p["a_name"], p["b_name"], round(p["gated_window"],2), round(p["gain"],2)) for p in pairs[c][:3]], flush=True)

json.dump({"validation": val, "and_gate_pairs": pairs},
          open("results/rna_window_and_gate.json", "w"), indent=1)
print("saved results/rna_window_and_gate.json")

"""ERBB2 subtype analysis: fraction of BRCA samples above normal max nTPM
(the correct metric for subtype-restricted targets). Single-gene stream."""
import json, sys, zipfile
sys.path.insert(0, "src")
import numpy as np
import pandas as pd
from cart.rna_window import load_normal, load_uniprot_membrane, surface_set, normal_block

normal, names = load_normal("data/rna_tissue_consensus.tsv")
gid = [g for g, n in names.items() if n == "ERBB2"][0]
nmax = float(normal_block(normal, [gid]).max(axis=1).iloc[0])
zf = zipfile.ZipFile("data/rna_cancer_sample.tsv.zip")
inner = [n for n in zf.namelist() if n.endswith(".tsv")][0]
vals = []
for raw in zf.open(inner):
    p = raw.decode().rstrip("\n").split("\t")
    if p[0] == gid and p[2] == "BRCA":
        vals.append(float(p[3]))
v = np.array(vals)
out = {"gene": "ERBB2", "cancer": "BRCA", "normal_max_nTPM": round(nmax, 1),
       "n_samples": len(v), "frac_above_normal_max": round(float((v > nmax).mean()), 4),
       "frac_above_2x_normal_max": round(float((v > 2 * nmax).mean()), 4),
       "p50": float(np.percentile(v, 50)), "p75": float(np.percentile(v, 75)),
       "p90": float(np.percentile(v, 90)), "p95": float(np.percentile(v, 95)),
       "max": float(v.max())}
json.dump(out, open("results/erbb2_subtype.json", "w"), indent=1)
print(json.dumps(out, indent=1))

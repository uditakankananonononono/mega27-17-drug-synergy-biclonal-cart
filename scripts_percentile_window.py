"""Percentile-based tumor metric (p50/p75/p90) to fix mean-expression dilution
of subtype-restricted targets (documented ERBB2 -0.19 artifact). Streams the
7GB per-sample TCGA file (sorted by Gene) with O(samples-per-gene) memory.
Output: results/rna_window_percentile.json"""
import json, sys, zipfile
sys.path.insert(0, "src")
import numpy as np
import pandas as pd
from cart.rna_window import load_normal, surface_set, normal_block, load_uniprot_membrane

CANCERS = ["PAAD", "BRCA", "GBM", "LIHC", "OV", "LUAD", "COAD", "SKCM", "STAD", "KIRC"]
PCTS = [50, 75, 90]

normal, names = load_normal("data/rna_tissue_consensus.tsv")
loc = pd.read_csv("data/subcellular_location.tsv", sep="\t")
uni = load_uniprot_membrane("data/uniprot_cellmembrane.tsv")
surf = surface_set(loc, uni, names)
cset = set(CANCERS)
print("surface genes:", len(surf), flush=True)

zf = zipfile.ZipFile("data/rna_cancer_sample.tsv.zip")
inner = [n for n in zf.namelist() if n.endswith(".tsv")][0]

result = {}
cur, vals, nlines = None, {}, 0
def finalize():
    if cur in surf and vals:
        result[cur] = {c: [float(np.percentile(v, p)) for p in PCTS] + [float(np.mean(v)), len(v)]
                       for c, v in vals.items()}

stream = zf.open(inner)
header = stream.readline()
for raw in stream:
    nlines += 1
    g, s, c, f = raw.decode().rstrip("\n").split("\t")
    if g != cur:
        finalize(); cur, vals = g, {}
        if nlines % 20_000_000 < 20:
            print("lines", nlines, "genes done", len(result), flush=True)
    if g in surf and c in cset:
        vals.setdefault(c, []).append(float(f))
finalize()
print("lines", nlines, "surface genes with data", len(result), flush=True)
json.dump({"percentiles": PCTS + ["mean", "n"], "data": result},
          open("results/tumor_percentiles.json", "w"))
print("saved results/tumor_percentiles.json", flush=True)

# windows + validation + AND-gate under p75
genes = sorted(result)
N = normal_block(normal, genes).to_numpy()
nmax = N.max(axis=1)
def win(metric_idx):
    W = {}
    for g in genes:
        W[g] = {c: float(np.log2(result[g][c][metric_idx] + 1) - np.log2(nmax[genes.index(g)] + 1))
                for c in result[g]}
    return W
Wmean, W75 = win(3), win(1)
name2id = {}
for gid, nm in names.items():
    name2id.setdefault(nm, []).append(gid)
KNOWN = {"ERBB2": "BRCA", "MSLN": "OV", "GPC3": "LIHC", "EGFR": "GBM", "FOLR1": "OV"}
valid = {}
for sym, can in KNOWN.items():
    for gid in name2id.get(sym, []):
        if gid in result and can in result[gid]:
            valid[sym] = {"mean_window": round(Wmean[gid][can], 2), "p75_window": round(W75[gid][can], 2),
                          "p50": round(result[gid][can][0], 1), "p75": round(result[gid][can][1], 1),
                          "p90": round(result[gid][can][2], 1), "mean": round(result[gid][can][3], 1),
                          "n_samples": result[gid][can][4]}
print(json.dumps(valid, indent=1), flush=True)

pairs = {}
for c in CANCERS:
    have = [(g, result[g][c][1]) for g in genes if c in result[g]]
    have = [(g, v) for g, v in have if v >= 20.0]
    have.sort(key=lambda x: -x[1]); have = have[:120]
    ids = [g for g, _ in have]
    Tv = np.array([v for _, v in have])
    Nv = np.array([normal_block(normal, [g]).to_numpy()[0] for g in ids])
    single = np.log2(Tv + 1) - np.log2(Nv.max(axis=1) + 1)
    out = []
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            tum = min(Tv[a], Tv[b]); nrm = np.minimum(Nv[a], Nv[b]).max()
            w = np.log2(tum + 1) - np.log2(nrm + 1)
            out.append({"a": ids[a], "b": ids[b], "gated_window": float(w),
                        "gain": float(w - max(single[a], single[b]))})
    out.sort(key=lambda d: -d["gated_window"])
    pairs[c] = [{**p, "a_name": names.get(p["a"]), "b_name": names.get(p["b"])} for p in out[:15]]
    print(c, [(p["a_name"], p["b_name"], round(p["gated_window"], 2), round(p["gain"], 2)) for p in pairs[c][:3]], flush=True)

json.dump({"metric": "p75 of per-sample TCGA FPKM (subtype-robust) vs per-cancer mean; normal max nTPM unchanged",
           "known_target_windows": valid, "and_gate_pairs_p75": pairs},
          open("results/rna_window_percentile.json", "w"), indent=1)
print("saved results/rna_window_percentile.json")

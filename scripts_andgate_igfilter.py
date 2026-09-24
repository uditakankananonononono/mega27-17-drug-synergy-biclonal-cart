"""p75 AND-gate rerun with plasma-cell infiltrate filter.

Paper flagged: p75 pair rankings dominated by immunoglobulin genes (IGHA2/
IGHG2/IGHV*) = plasma-cell infiltrate, not malignant-cell surface expression.
Filter: exclude Ig-locus genes (IGH*/IGK*/IGL* name prefixes + JCHAIN) from
the candidate pool, rerun p75 AND-gate from committed tumor_percentiles.json
(no re-stream of the 7GB file). Reports before/after top-15 deltas.
Output: results/andgate_p75_igfiltered.json
"""
import json, re, sys
sys.path.insert(0, "src")
import numpy as np
import pandas as pd
from cart.rna_window import load_normal, surface_set, normal_block, load_uniprot_membrane

CANCERS = ["PAAD", "BRCA", "GBM", "LIHC", "OV", "LUAD", "COAD", "SKCM", "STAD", "KIRC"]
IG_RE = re.compile(r"^(IGH|IGK|IGL|JCHAIN)")

d = json.load(open("results/tumor_percentiles.json"))
result = d["data"]
normal, names = load_normal("data/rna_tissue_consensus.tsv")
loc = pd.read_csv("data/subcellular_location.tsv", sep="\t")
uni = load_uniprot_membrane("data/uniprot_cellmembrane.tsv")
surf = surface_set(loc, uni, names)
genes = sorted(result)
N = normal_block(normal, genes).to_numpy()
nmax = {g: float(N[i].max()) for i, g in enumerate(genes)}

ig_genes = [g for g in genes if IG_RE.match(names.get(g) or "")]
print(f"surface genes with p75 data: {len(genes)}; Ig-locus filtered: {len(ig_genes)}", flush=True)

def andgate(pool):
    pairs = {}
    for c in CANCERS:
        have = [(g, result[g][c][1]) for g in pool if c in result[g] and result[g][c][1] >= 20.0]
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
        out.sort(key=lambda x: -x["gated_window"])
        pairs[c] = [{**p, "a_name": names.get(p["a"]), "b_name": names.get(p["b"])} for p in out[:15]]
    return pairs

before = andgate(genes)
after = andgate([g for g in genes if g not in set(ig_genes)])

def top3(pairs, c):
    return [(p["a_name"], p["b_name"], round(p["gated_window"], 2)) for p in pairs[c][:3]]

delta = {}
for c in CANCERS:
    b_ig = sum(1 for p in before[c][:15] if IG_RE.match(p["a_name"] or "") or IG_RE.match(p["b_name"] or ""))
    delta[c] = {"ig_pairs_in_top15_before": b_ig,
                "top3_before": top3(before, c), "top3_after": top3(after, c)}
    print(c, f"Ig pairs in top15 before: {b_ig}", "top3 after:", top3(after, c), flush=True)

json.dump({"filter": "Ig-locus name filter (^(IGH|IGK|IGL|JCHAIN)) - plasma-cell infiltrate exclusion",
           "n_filtered": len(ig_genes), "filtered_gene_names": sorted({names[g] for g in ig_genes}),
           "metric": "p75 per-sample tumor FPKM vs normal max nTPM",
           "delta": delta, "and_gate_pairs_p75_igfiltered": after},
          open("results/andgate_p75_igfiltered.json", "w"), indent=1)
print("saved results/andgate_p75_igfiltered.json", flush=True)

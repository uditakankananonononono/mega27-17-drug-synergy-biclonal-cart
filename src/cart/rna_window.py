"""RNA-vs-RNA therapeutic window for CAR-T antigens and AND-gate antigen pairs
(HPA v23). T[g,c] = per-cancer mean FPKM (TCGA); N[g,t] = consensus nTPM.

Single antigen:  W(g,c)   = log2(T[g,c]+1) - log2(max_t N[g,t]+1)
AND gate (a,b):  W(a,b,c) = log2(min(T[a,c],T[b,c])+1)
                          - log2(max_t min(N[a,t],N[b,t])+1)
Gain(a,b,c) = W(a,b,c) - max(W(a,c), W(b,c)).
"""
from __future__ import annotations
import numpy as np
import pandas as pd

VITAL = {"cerebral cortex", "heart muscle", "lung", "liver", "kidney",
         "pancreas", "stomach", "colon", "small intestine", "esophagus",
         "bone marrow", "spleen", "skeletal muscle", "smooth muscle",
         "adrenal gland", "thyroid gland", "cerebellum", "basal ganglia",
         "hypothalamus", "spinal cord", "retina", "pituitary gland"}
SURFACE = ("Plasma membrane", "Cell Junctions")


def load_tumor(path):
    df = pd.read_csv(path, sep="\t")
    df = df[df["Cancer"] != "Cancer"]
    df["mean_FPKM"] = pd.to_numeric(df["mean_FPKM"], errors="coerce")
    return df.pivot_table(index="Gene", columns="Cancer", values="mean_FPKM")


def load_normal(path):
    df = pd.read_csv(path, sep="\t")
    names = dict(zip(df["Gene"], df["Gene name"]))
    return df.pivot_table(index="Gene", columns="Tissue", values="nTPM"), names


def surface_set(loc, uniprot_syms=None, names=None):
    """HPA subcellular plasma-membrane/cell-junction genes, optionally unioned
    with an external reviewed cell-membrane list (UniProt KW-1003 gene symbols)
    mapped to ENSG via the HPA names dict. HPA subcellular annotation misses
    canonical CAR-T antigens (CD19, BCMA/TNFRSF17, FOLR1 have no HPA row)."""
    txt = loc["Main location"].fillna("") + ";" + loc["Additional location"].fillna("")
    s = set(loc.loc[txt.str.contains("|".join(SURFACE)), "Gene"])
    if uniprot_syms and names:
        name2id = {}
        for gid, nm in names.items():
            name2id.setdefault(nm, []).append(gid)
        for sym in uniprot_syms:
            for gid in name2id.get(sym, []):
                s.add(gid)
    return s


def load_uniprot_membrane(path):
    """Primary + alias gene symbols from a UniProt TSV (Entry, Gene Names)."""
    import csv
    out = set()
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            for tok in (row.get("Gene Names") or "").split():
                out.add(tok)
    return out


def normal_block(normal, genes, vital_only=True):
    cols = [c for c in normal.columns if (c in VITAL or not vital_only)]
    return normal.loc[genes, cols].fillna(0.0)


def single_windows(tumor, normal, genes, vital_only=True):
    g = tumor.index.intersection(normal.index).intersection(pd.Index(sorted(genes)))
    nmax = normal_block(normal, g, vital_only).max(axis=1)
    return np.log2(tumor.loc[g].fillna(0.0) + 1).sub(np.log2(nmax + 1), axis=0)


def and_gate_pairs(tumor, normal, genes, cancer, n_cand=150, min_tumor=10.0, top=25):
    """Exhaustive AND-gate pair search over the top tumour-expressed surface
    genes in one cancer. Returns list of dicts sorted by gated window."""
    g = tumor.index.intersection(normal.index).intersection(pd.Index(sorted(genes)))
    t = tumor.loc[g, cancer].fillna(0.0)
    t = t[t >= min_tumor].sort_values(ascending=False).head(n_cand)
    ids = list(t.index)
    N = normal_block(normal, ids).to_numpy()
    T = t.to_numpy()
    single = np.log2(T + 1) - np.log2(N.max(axis=1) + 1)
    out = []
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            tum = min(T[i], T[j])
            nrm = np.minimum(N[i], N[j]).max()
            w = np.log2(tum + 1) - np.log2(nrm + 1)
            out.append({"a": ids[i], "b": ids[j], "cancer": cancer,
                        "gated_window": float(w),
                        "gain": float(w - max(single[i], single[j])),
                        "tumor_min_fpk": float(tum),
                        "normal_max_nTPM": float(nrm)})
    out.sort(key=lambda d: -d["gated_window"])
    return out[:top]

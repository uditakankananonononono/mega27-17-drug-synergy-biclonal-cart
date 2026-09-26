"""AMENDMENT-1 (Judge round 1): SINGLE-CELL AND-GATE FEASIBILITY ATLAS.

Question (from ChatGPT judge r1 critique): bulk co-expression does not prove the
two antigens co-occur on the SAME tumor cell, and bulk normal-tissue signal does
not prove the AND-gate spares normal cells. This arm measures, per public
scRNA-seq dataset (TISCH2), on malignant vs non-malignant (and normal-source) cells:

  coverage_AND  = fraction of malignant cells positive for BOTH antigens
  leakage_AND   = fraction of non-malignant (or normal-source) cells double+
  vs the same for each single antigen, at three positivity thresholds:
  (a) detected (>0), (b) >=2 UMI-ish units, (c) top-25% among expressors
Controls: per real pair, 200 expression-matched random gene pairs (detection
rate matched within +/-5 pct points per gene) -> null distribution of the AND
separation score; report percentile of the real pair.
Holdout: for cancers with 2 datasets, rank pairs on dataset A, test on B;
report Spearman of pair scores across datasets.

Datasets live OUTSIDE the repo (/home/sandbox/singlecell_atlas/<CANCER>/)
because of size; this script is the durable IP. Outputs:
  results/scatlas_per_dataset.json  (all per-dataset pair metrics)
  results/scatlas_holdout.json      (cross-dataset rank tests)
  results/scatlas_summary.md        (paper-ready table)
"""
import json, os, sys, csv, math, random
import numpy as np

ATLAS_DIR = "/home/sandbox/singlecell_atlas"
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(RESULTS, exist_ok=True)

GENES = ["CA9","CA12","CLDN18","CLDN3","CLDN6","MSLN","PSCA","SPP1","KLK7",
         "EGFR","GAP43","SLC34A2","SLC17A3","CLTRN","DPEP1","NOX1","APOC3",
         "OR2I1P","PRAME","DCT","SLC39A6","HLA-DRA","TACSTD2","GPC3","FOLR1",
         "EPCAM","CEACAM5","ERBB2"]

PAIRS = [  # judge-specified
 ("CA9","CLDN18"),("CA9","MSLN"),("CLDN6","MSLN"),("PSCA","CA9"),
 # bulk Table-1 winners (andgate_p75_igfiltered.json)
 ("CA9","SLC17A3"),("CA9","CA12"),("CLDN18","CLDN3"),("SPP1","PSCA"),
 ("KLK7","CLDN6"),("SLC34A2","SPP1"),("GAP43","EGFR"),("DPEP1","NOX1"),
 ("APOC3","OR2I1P"),("HLA-DRA","SLC39A6"),("PRAME","DCT"),
]

def load_meta(path):
    """Header-driven TISCH2 meta parse -> dict cell -> (is_malignant, is_normal_source, lineage)."""
    out = {}
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            cell = row.get("Cell")
            mal = (row.get("Celltype (malignancy)") or "").strip()
            src = (row.get("Source") or row.get("Tissue") or "").strip()
            lin = (row.get("Celltype (major-lineage)") or "").strip()
            out[cell] = (mal == "Malignant cells", src.lower() == "normal", lin)
    return out

def load_matrix(h5path):
    """Return (genes list, cells list, scipy.sparse.csr matrix cells x genes).
    TISCH2 h5 layout discovered by inspection; fill in after first file opens."""
    import h5py
    from scipy import sparse
    f = h5py.File(h5path, "r")
    keys = list(f.keys())
    # common layouts: 10x-style 'matrix' group or flat datasets
    def dec(a):
        return [x.decode() if isinstance(x, bytes) else str(x) for x in a]
    if "matrix" in keys:
        g = f["matrix"]
        data, indices, indptr = g["data"][:], g["indices"][:], g["indptr"][:]
        shape = tuple(g["shape"][:])
        barcodes = dec(g["barcodes"][:])
        genes = dec(g["features"]["name"][:] if "features" in g else g["genes"][:])
        M = sparse.csr_matrix((data, indices, indptr), shape=shape)
    else:
        # flat: X + obs/var style or simple datasets
        ks = {k.lower(): k for k in keys}
        if "x" in ks:
            X = f[ks["x"]][:]
            genes = dec(f[ks.get("var", ks.get("genes"))][:] if (ks.get("var") or ks.get("genes")) else [])
            barcodes = dec(f[ks.get("obs", ks.get("barcodes"))][:] if (ks.get("obs") or ks.get("barcodes")) else [])
            M = sparse.csr_matrix(X)
        else:
            raise RuntimeError(f"unknown h5 layout: {keys}")
    f.close()
    # ensure cells x genes orientation: barcodes length == M.shape[0]
    if M.shape[0] != len(barcodes) and M.shape[1] == len(barcodes):
        M = M.T.tocsr()
    return genes, barcodes, M

def thresholds(vec):
    """vec: 1D np array of expression. Return dict of boolean masks."""
    pos = vec > 0
    out = {"detected": pos}
    out["ge2"] = vec >= 2
    if pos.sum() >= 4:
        q75 = np.quantile(vec[pos], 0.75)
        out["top25"] = vec >= q75 if q75 > 0 else (vec > 0)
    else:
        out["top25"] = vec >= vec.max() if vec.max() > 0 else pos
    return out

def pair_metrics(mA, mB, mal, norm):
    """masks: boolean arrays over cells; mal: malignant mask; norm: normal-source mask.
    Returns coverage/leakage for singles + AND."""
    res = {}
    for name, m in [("A", mA), ("B", mB), ("AND", mA & mB)]:
        res[f"cov_{name}"] = float((m & mal).sum() / max(mal.sum(), 1))
        nonmal = ~mal
        res[f"leak_nonmal_{name}"] = float((m & nonmal).sum() / max(nonmal.sum(), 1))
        if norm.sum() > 0:
            res[f"leak_normsrc_{name}"] = float((m & norm).sum() / max(norm.sum(), 1))
    # separation: coverage / (leakage + pseudocount)
    for leakkey in ["leak_nonmal", "leak_normsrc"]:
        if f"{leakkey}_AND" in res:
            for name in ["A", "B", "AND"]:
                res[f"sep_{leakkey}_{name}"] = res[f"cov_{name}"] / (res[f"{leakkey}_{name}"] + 0.001)
    return res

def analyze_dataset(cancer, dsname, rng):
    d = os.path.join(ATLAS_DIR, cancer)
    h5 = os.path.join(d, f"{dsname}_expression.h5")
    meta = os.path.join(d, f"{dsname}_CellMetainfo_table.tsv")
    genes, cells, M = load_matrix(h5)
    meta_d = load_meta(meta)
    cell_arr = np.array(cells)
    mal = np.array([meta_d.get(c, (None, None, None))[0] for c in cells], dtype=object)
    keep = np.array([m is not None for m in mal])
    M = M[keep]; cell_arr = cell_arr[keep]
    mal = np.array([bool(m) for m in mal[keep]])
    norm = np.array([bool(meta_d.get(c, (None, None, None))[1]) for c in cell_arr])
    gidx = {g: i for i, g in enumerate(genes)}
    present = [g for g in GENES if g in gidx]
    cols = {g: np.asarray(M[:, gidx[g]].todense()).ravel() for g in present}
    det_rate = {g: float((cols[g] > 0).mean()) for g in present}
    out = {"dataset": dsname, "cancer": cancer,
           "n_cells": int(keep.sum()), "n_malignant": int(mal.sum()),
           "n_normal_source": int(norm.sum()),
           "genes_present": present, "detection_rates": det_rate, "pairs": {}}
    all_genes = list(genes)
    for a, b in PAIRS:
        if a not in cols or b not in cols:
            out["pairs"][f"{a}-{b}"] = {"skipped": "gene absent"}; continue
        va, vb = cols[a], cols[b]
        rec = {}
        for thr in ["detected", "ge2", "top25"]:
            ma, mb = thresholds(va)[thr], thresholds(vb)[thr]
            rec[thr] = pair_metrics(ma, mb, mal, norm)
        # expression-matched random-pair null on 'detected' threshold
        seps = []
        da, db = det_rate[a], det_rate[b]
        cand_a = [g for g in all_genes if abs(float((np.asarray(M[:, gidx[g]].todense()).ravel() > 0).mean()) - da) <= 0.05 and g not in PAIRS_GENES]
        # cheaper: precompute detection for candidate pool once outside loop (see below)
        out["pairs"][f"{a}-{b}"] = rec
    return out, M, gidx, genes, mal, norm, cols

PAIRS_GENES = set(g for p in PAIRS for g in p)

def null_distribution(M, gidx, genes, mal, da, db, n=200, seed=1):
    """Random pairs matched to detection rates da, db; AND separation on detected."""
    rng = np.random.default_rng(seed)
    # precompute detection vector lazily for sampled genes
    def det(g):
        return np.asarray(M[:, gidx[g]].todense()).ravel() > 0
    pool = [g for g in genes if g not in PAIRS_GENES]
    rng.shuffle(pool)
    seps = []
    made = 0
    for i in range(0, len(pool) - 1, 2):
        if made >= n: break
        ga, gb = pool[i], pool[i + 1]
        ma, mb = det(ga), det(gb)
        ra, rb = ma.mean(), mb.mean()
        if abs(ra - da) > 0.05 or abs(rb - db) > 0.05: continue
        mand = ma & mb
        cov = (mand & mal).sum() / max(mal.sum(), 1)
        nonmal = ~mal
        leak = (mand & nonmal).sum() / max(nonmal.sum(), 1)
        seps.append(float(cov / (leak + 0.001)))
        made += 1
    return seps

if __name__ == "__main__":
    # phase 1: inspect only
    if len(sys.argv) > 1 and sys.argv[1] == "inspect":
        import h5py
        p = sys.argv[2]
        f = h5py.File(p, "r")
        def walk(n, o):
            print(n, type(o).__name__, getattr(o, "shape", ""), getattr(o, "dtype", ""))
        f.visititems(walk)
        sys.exit(0)
    print("run per-dataset after downloads complete")

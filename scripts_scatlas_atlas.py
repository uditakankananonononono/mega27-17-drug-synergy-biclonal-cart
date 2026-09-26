"""AMENDMENT-1 (Judge round 1): SINGLE-CELL AND-GATE FEASIBILITY ATLAS.

Question (judge r1): bulk co-expression does not prove two antigens co-occur on
the SAME tumor cell; bulk normal signal does not prove the AND-gate spares
normal cells. This arm measures, per public scRNA-seq dataset (TISCH2, uniform
malignant annotations), on malignant vs non-malignant (and normal-source) cells:

  coverage  = fraction of malignant cells positive (single or AND)
  leakage   = fraction of non-malignant (or normal-source) cells positive
  separation = coverage / (leakage + 0.001)

Thresholds (TISCH2 matrices are log-normalized, verified: continuous values,
min>0; UMI thresholds meaningless): (a) detected (>0), (b) above median among
expressors, (c) top quartile among expressors.

Controls: per real pair, up to 200 random gene pairs matched per-gene on
detection rate (+/-5 pct points) -> null distribution of AND separation;
report the real pair's percentile.

Holdout: for cancers with 2 datasets, rank all pairs by AND separation on
dataset A, report Spearman vs dataset B ranking (cross-dataset reproducibility).

Datasets (size) live OUTSIDE the repo at /home/sandbox/singlecell_atlas/<CANCER>/;
this script + results JSONs are the durable IP. Outputs:
  results/scatlas_per_dataset.json, results/scatlas_holdout.json,
  results/scatlas_summary.md
"""
import json, os, sys, csv
import numpy as np
from scipy import sparse
import h5py

ATLAS_DIR = "/home/sandbox/singlecell_atlas"
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

GENES = ["CA9","CA12","CLDN18","CLDN3","CLDN6","MSLN","PSCA","SPP1","KLK7",
         "EGFR","GAP43","SLC34A2","SLC17A3","CLTRN","DPEP1","NOX1","APOC3",
         "OR2I1P","PRAME","DCT","SLC39A6","HLA-DRA","TACSTD2","GPC3","FOLR1",
         "EPCAM","CEACAM5","ERBB2"]
PAIRS = [("CA9","CLDN18"),("CA9","MSLN"),("CLDN6","MSLN"),("PSCA","CA9"),
         ("CA9","SLC17A3"),("CA9","CA12"),("CLDN18","CLDN3"),("SPP1","PSCA"),
         ("KLK7","CLDN6"),("SLC34A2","SPP1"),("GAP43","EGFR"),("DPEP1","NOX1"),
         ("APOC3","OR2I1P"),("HLA-DRA","SLC39A6"),("PRAME","DCT")]
PAIR_GENES = set(g for p in PAIRS for g in p)

DATASETS = {  # cancer -> [dataset names]
 "PAAD": ["PAAD_GSE111672","PAAD_CRA001160"],
 "STAD": ["STAD_GSE134520","STAD_GSE167297"],
 "OV":   ["OV_GSE154763","OV_EMTAB8107"],
 "LIHC": ["LIHC_GSE140228_10X"],
 "KIRC": ["KIRC_GSE159115"],
}

def load_meta(path):
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            mal = (row.get("Celltype (malignancy)") or "").strip()
            src = (row.get("Source") or row.get("Tissue") or "").strip()
            out[row["Cell"]] = (mal == "Malignant cells", src.lower() == "normal")
    return out

def load_matrix(h5path):
    """10x-style h5 (genes x cells CSC) -> (genes, cells, csr cells x genes)."""
    f = h5py.File(h5path, "r")
    g = f["matrix"]
    dec = lambda a: [x.decode() if isinstance(x, bytes) else str(x) for x in a]
    M = sparse.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]),
                          shape=tuple(g["shape"][:]))  # 10x stores CSC (genes x cells)
    barcodes = dec(g["barcodes"][:]); genes = dec(g["features"]["name"][:])
    f.close()
    if M.shape[0] != len(barcodes):
        M = M.T.tocsr()
    assert M.shape[0] == len(barcodes) and M.shape[1] == len(genes)
    return genes, barcodes, M

def thr_masks(v):
    pos = v > 0
    out = {"detected": pos}
    if pos.sum() >= 8:
        med = np.median(v[pos]); q75 = np.quantile(v[pos], 0.75)
        out["above_med"] = v >= med
        out["top25"] = v >= q75 if q75 > 0 else pos
    else:
        out["above_med"] = pos; out["top25"] = pos
    return out

def sep(m, mal, nonmal, norm):
    cov = float((m & mal).sum() / max(mal.sum(), 1))
    leak_nm = float((m & nonmal).sum() / max(nonmal.sum(), 1))
    r = {"cov": cov, "leak_nonmal": leak_nm,
         "sep_nonmal": cov / (leak_nm + 0.001)}
    if norm.sum() > 0:
        leak_n = float((m & norm).sum() / max(norm.sum(), 1))
        r["leak_normsrc"] = leak_n
        r["sep_normsrc"] = cov / (leak_n + 0.001)
    return r

def analyze(cancer, ds, rng):
    d = os.path.join(ATLAS_DIR, cancer)
    genes, cells, M = load_matrix(os.path.join(d, f"{ds}_expression.h5"))
    meta = load_meta(os.path.join(d, f"{ds}_CellMetainfo_table.tsv"))
    flags = np.array([meta.get(c) for c in cells], dtype=object)
    keep = np.array([f is not None for f in flags])
    M = M[keep]
    mal = np.array([bool(f[0]) for f in flags[keep]])
    norm = np.array([bool(f[1]) for f in flags[keep]])
    nonmal = ~mal
    gidx = {g: i for i, g in enumerate(genes)}
    present = [g for g in GENES if g in gidx]
    col = {g: np.asarray(M[:, gidx[g]].todense()).ravel() for g in present}
    det = {g: float((col[g] > 0).mean()) for g in present}
    rec = {"dataset": ds, "cancer": cancer, "n_cells": int(keep.sum()),
           "n_malignant": int(mal.sum()), "n_normal_source": int(norm.sum()),
           "genes_present": present,
           "genes_absent": [g for g in GENES if g not in gidx],
           "detection_rates": det, "pairs": {}}
    # binarize once; pool = dense bool matrix of non-pair genes (sampled order)
    pool = [g for g in genes if g not in PAIR_GENES]
    rng.shuffle(pool)
    pool = pool[:4000]
    pidx = [gidx[g] for g in pool]
    PB = np.asarray(M[:, pidx].todense()) > 0           # cells x pool
    pdr = PB.mean(axis=0)
    for a, b in PAIRS:
        key = f"{a}-{b}"
        if a not in col or b not in col:
            rec["pairs"][key] = {"skipped": "gene absent"}; continue
        pa = {}
        for thr in ["detected", "above_med", "top25"]:
            ma, mb = thr_masks(col[a])[thr], thr_masks(col[b])[thr]
            pa[thr] = {"A": sep(ma, mal, nonmal, norm),
                       "B": sep(mb, mal, nonmal, norm),
                       "AND": sep(ma & mb, mal, nonmal, norm)}
        # expression-matched random-pair null on 'detected'
        da, db = det[a], det[b]
        nulls = []
        for i in range(0, len(pool) - 1, 2):
            if len(nulls) >= 200: break
            if abs(pdr[i] - da) > 0.05 or abs(pdr[i + 1] - db) > 0.05: continue
            s = sep(PB[:, i] & PB[:, i + 1], mal, nonmal, norm)
            nulls.append(s["sep_nonmal"])
        real = pa["detected"]["AND"]["sep_nonmal"]
        pa["null_sep_nonmal_detected"] = {
            "n": len(nulls),
            "median": float(np.median(nulls)) if nulls else None,
            "pctl_of_real": float(np.mean([x <= real for x in nulls])) if nulls else None}
        rec["pairs"][key] = pa
    return rec

def main():
    rng = np.random.default_rng(7)
    only = sys.argv[1] if len(sys.argv) > 1 else None
    out_path = os.path.join(RESULTS, "scatlas_per_dataset.json")
    allrec = {}
    if os.path.exists(out_path):
        allrec = json.load(open(out_path))
    for cancer, dss in DATASETS.items():
        for ds in dss:
            if only and ds != only: continue
            h5 = os.path.join(ATLAS_DIR, cancer, f"{ds}_expression.h5")
            meta = os.path.join(ATLAS_DIR, cancer, f"{ds}_CellMetainfo_table.tsv")
            if not (os.path.exists(h5) and os.path.exists(meta)):
                print(f"SKIP {ds}: files missing"); continue
            print(f"analyzing {ds} ...", flush=True)
            allrec[ds] = analyze(cancer, ds, rng)
            json.dump(allrec, open(out_path, "w"), indent=1)
            print(f"  saved {ds}", flush=True)
    print("done")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "inspect":
        f = h5py.File(sys.argv[2], "r")
        f.visititems(lambda n, o: print(n, type(o).__name__, getattr(o, "shape", ""), getattr(o, "dtype", "")))
    else:
        main()

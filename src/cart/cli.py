"""cart-target-rank - rank CAR-T target candidates by tumor-vs-normal window.

For each gene symbol: surface membership (HPA subcellular + UniProt KW-1003),
normal vital-tissue max nTPM, and per-sample tumor stats streamed from the HPA
per-sample TCGA file: FANM (fraction above normal max), p75, mean - the metrics
validated in results/erbb2_subtype.json and results/rna_window_percentile.json.

Usage:
  cart-target-rank GENE [GENE...] --cancer BRCA [--data-dir data] [--json out.json]
"""
import argparse, csv, json, pathlib, sys, zipfile

import numpy as np
import pandas as pd


def vital_normal_max(normal, names, symbols, vital_only=True):
    from .rna_window import normal_block
    sym2gid = {}
    for gid, nm in names.items():
        sym2gid.setdefault(nm, gid)
    gids = {s: sym2gid.get(s) for s in symbols}
    found = {s: g for s, g in gids.items() if g}
    nmax = normal_block(normal, list(found.values()), vital_only).max(axis=1) if found else pd.Series(dtype=float)
    return gids, {g: float(nmax[g]) for g in found.values()}


def stream_tumor_stats(zip_path, gid_set, cancer, normal_max):
    """Single pass over the per-sample tumor TSV (sorted by Gene) collecting
    FANM/percentiles for the requested genes in one cancer."""
    vals = {g: [] for g in gid_set}
    with zipfile.ZipFile(zip_path) as zf:
        inner = [n for n in zf.namelist() if n.endswith(".tsv")][0]
        for raw in zf.open(inner):
            p = raw.decode().rstrip("\n").split("\t")
            if p[0] in vals and p[2] == cancer:
                try:
                    vals[p[0]].append(float(p[3]))
                except (ValueError, IndexError):
                    continue
    out = {}
    for g, v in vals.items():
        v = np.array(v)
        if len(v) == 0:
            out[g] = {"n_samples": 0}
            continue
        nm = normal_max.get(g)
        row = {"n_samples": int(len(v)), "mean": round(float(v.mean()), 1),
               "p75": round(float(np.percentile(v, 75)), 1),
               "p95": round(float(np.percentile(v, 95)), 1), "max": round(float(v.max()), 1)}
        if nm is not None:
            row["normal_max_nTPM"] = round(nm, 1)
            row["fanm"] = round(float((v > nm).mean()), 4)
            row["fanm_2x"] = round(float((v > 2 * nm).mean()), 4)
        out[g] = row
    return out


def rank_genes(symbols, cancer, data_dir="data", vital_only=True):
    from .rna_window import load_normal, load_uniprot_membrane, surface_set
    data_dir = pathlib.Path(data_dir)
    normal, names = load_normal(data_dir / "rna_tissue_consensus.tsv")
    loc = pd.read_csv(data_dir / "subcellular_location.tsv", sep="\t")
    uni = load_uniprot_membrane(data_dir / "uniprot_cellmembrane.tsv")
    surf = surface_set(loc, uni, names)
    gids, nmax = vital_normal_max(normal, names, symbols, vital_only)
    stats = stream_tumor_stats(data_dir / "rna_cancer_sample.tsv.zip",
                               {g for g in gids.values() if g}, cancer, nmax)
    rows = {}
    for sym in symbols:
        gid = gids.get(sym)
        if not gid:
            rows[sym] = {"error": "symbol not in HPA consensus"}
            continue
        row = dict(stats.get(gid, {}))
        row["gene_id"] = gid
        row["surface"] = gid in surf
        rows[sym] = row
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(prog="cart-target-rank",
        description="Rank CAR-T targets: surface membership + FANM/percentile tumor-vs-normal stats.")
    p.add_argument("genes", nargs="+", help="gene symbols, e.g. ERBB2 MSLN GPC3")
    p.add_argument("--cancer", required=True, help="TCGA cancer code, e.g. BRCA")
    p.add_argument("--data-dir", default="data")
    p.add_argument("--all-tissues", action="store_true", help="use all normal tissues, not just vital")
    p.add_argument("--json", help="write results JSON here")
    a = p.parse_args(argv)
    rows = rank_genes(a.genes, a.cancer, a.data_dir, vital_only=not a.all_tissues)
    out = {"tool": "cart-target-rank", "cancer": a.cancer, "genes": rows}
    if a.json:
        json.dump(out, open(a.json, "w"), indent=1)
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

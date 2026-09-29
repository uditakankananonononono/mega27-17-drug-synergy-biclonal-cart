"""cart-target-rank - rank CAR-T target candidates by tumor-vs-normal window.

--report mode: joins the live ranking axes with the committed per-gene audit
corpus (results/*.csv: Pharos TDL/novelty/drugs, BioPlex interactome presence,
ClinicalTrials CAR-T status, gnomAD constraint, DrugCentral engagement) for the
188-gene study corpus. Genes outside the corpus get live ranking axes only plus
an explicit boundary note - the full battery requires the fetch scripts
(scripts_*.py, documented in docs/TOOLS.md). No network access in report mode.

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


REPORT_SOURCES = {
    "pharos": ("pharos_per_gene.csv", ["tdl", "novelty", "n_drugs", "n_ligands",
                                       "publication_count", "gwas_total", "ppi_total"]),
    "bioplex": ("bioplex_per_gene.csv", ["293T_present", "293T_degree", "HCT116_present", "HCT116_degree"]),
    "clintrials": ("clintrials_per_gene.csv", ["strict_total", "strict_onc_total", "strict_onc_active",
                                               "verdict", "n_nct_union"]),
    "constraint": ("constraint_per_gene.csv", ["loeuf", "pli", "mis_z", "depmap_median",
                                               "depmap_frac_dep", "escape_quadr"]),
    "drugcentral": ("drugcentral_engagement_rows.csv", ["n_drugs", "tdl", "moa", "drug_names"]),
}

BOUNDARY_NOTE = ("outside the 188-gene study corpus - live ranking axes only; "
                 "full battery needs scripts_*.py fetch runs (docs/TOOLS.md); "
                 "CPTAC protein, DepMap matrix and DailyMed label tables are "
                 "study-level artifacts not joined by this tool")


def load_report_corpus(results_dir):
    corpus = {}
    for src, (fn, cols) in REPORT_SOURCES.items():
        path = pathlib.Path(results_dir) / fn
        if not path.exists():
            continue
        with open(path) as fh:
            for row in csv.DictReader(fh):
                gene = row.get("gene")
                if not gene:
                    continue
                entry = corpus.setdefault(gene, {})
                entry[src] = {c: row[c] for c in cols if c in row and row[c] != ""}
    return corpus


def main(argv=None):
    p = argparse.ArgumentParser(prog="cart-target-rank",
        description="Rank CAR-T targets: surface membership + FANM/percentile tumor-vs-normal stats.")
    p.add_argument("genes", nargs="+", help="gene symbols, e.g. ERBB2 MSLN GPC3")
    p.add_argument("--cancer", required=True, help="TCGA cancer code, e.g. BRCA")
    p.add_argument("--data-dir", default="data")
    p.add_argument("--all-tissues", action="store_true", help="use all normal tissues, not just vital")
    p.add_argument("--json", help="write results JSON here")
    p.add_argument("--report", action="store_true",
                   help="join committed per-gene audit corpus (results/*.csv); no network")
    p.add_argument("--results-dir", default="results")
    a = p.parse_args(argv)
    rows = rank_genes(a.genes, a.cancer, a.data_dir, vital_only=not a.all_tissues)
    out = {"tool": "cart-target-rank", "cancer": a.cancer, "genes": rows}
    if a.report:
        corpus = load_report_corpus(a.results_dir)
        for sym in a.genes:
            if sym in corpus:
                rows[sym]["audit"] = corpus[sym]
                rows[sym]["corpus"] = "188-gene study corpus"
            else:
                rows[sym]["audit"] = {}
                rows[sym]["corpus"] = BOUNDARY_NOTE
        out["report_sources"] = {s: spec[0] for s, spec in REPORT_SOURCES.items()}
    if a.json:
        json.dump(out, open(a.json, "w"), indent=1)
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

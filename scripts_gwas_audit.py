#!/usr/bin/env python3
"""GWAS Catalog germline-association audit for the 205-gene study set.
Question: do gated tumor-antigen loci carry excess germline association
signal - overall and for cancer traits - versus the random surfaceome
background? GWAS captures locus-level population relevance, orthogonal to
the somatic/clinical axes already audited (constraint, trials, Pharos).
Source: GWAS Catalog full ontology-annotated associations TSV (latest
release, bulk file). Association assigned to a gene when the gene appears in
the MAPPED_GENE column (multi-gene LD blocks assign to each mapped gene -
caveat kept). Significance: PVALUE_MLOG >= 7.301 (p <= 5e-8).
Pre-registered gates:
  G1: PSCA has >=1 significant CANCER-trait association (gastric/bladder
      carcinoma locus rs2294008 documented).
  G2: parsed association rows >= 400,000 (catalog depth / file integrity).
  G3: EGFR sentinel has >=1 significant cancer-trait association (glioma
      locus documented).
Analyses:
  H1: fraction of genes with >=1 significant association, gated vs bg.
  H2: fraction with >=1 significant CANCER association, gated vs bg;
      per-antigen table.
  H3: significant-association burden per gene (MW)."""
import csv, glob, json, re
from collections import defaultdict
from scipy.stats import mannwhitneyu, fisher_exact

ANTS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]
SIG_MLOG = 7.301  # p <= 5e-8
CANCER = re.compile(r"carcinoma|cancer\b|neoplasm|melanoma|leukemia|leukaemia|"
                    r"lymphoma|sarcoma|myeloma|glioma|glioblastoma|mesothelioma|"
                    r"blastoma|adenocarcinoma|malignan", re.I)


def gene_set():
    gs = json.load(open("results/depmap_gene_sets.json"))
    genes = {}
    for k in ("gated", "background", "references", "positive_controls"):
        for g in gs[k]:
            genes[g] = k
    return genes


def split_mapped(s):
    toks = re.split(r",\s*|\s+-\s+", s.strip())
    return [t.strip() for t in toks if t.strip()]


def main():
    genes = gene_set()
    universe = set(genes)
    path = glob.glob("data/gwas/gwas-catalog-download-associations-*.tsv")
    assert path, "associations TSV missing - run scripts_gwas_fetch.py"
    per = {g: {"n_assoc": 0, "n_sig": 0, "n_sig_cancer": 0,
               "sig_snps": set(), "cancer_traits": set()} for g in universe}
    n_rows = 0
    with open(path[0], newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            n_rows += 1
            hits = [g for g in split_mapped(r["MAPPED_GENE"]) if g in universe]
            if not hits:
                continue
            try:
                mlog = float(r["PVALUE_MLOG"])
            except (ValueError, KeyError):
                continue
            sig = mlog >= SIG_MLOG
            cancer = bool(CANCER.search(r["MAPPED_TRAIT"] or ""))
            for g in hits:
                p = per[g]
                p["n_assoc"] += 1
                if sig:
                    p["n_sig"] += 1
                    p["sig_snps"].add(r["SNPS"])
                    if cancer:
                        p["n_sig_cancer"] += 1
                        p["cancer_traits"].add(r["MAPPED_TRAIT"])
    # EGFR sentinel
    egfr = {"n_assoc": 0, "n_sig": 0, "n_sig_cancer": 0}
    with open(path[0], newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            if "EGFR" not in split_mapped(r["MAPPED_GENE"]):
                continue
            try:
                mlog = float(r["PVALUE_MLOG"])
            except (ValueError, KeyError):
                continue
            egfr["n_assoc"] += 1
            if mlog >= SIG_MLOG:
                egfr["n_sig"] += 1
                if CANCER.search(r["MAPPED_TRAIT"] or ""):
                    egfr["n_sig_cancer"] += 1

    gates = {
        "G1_psca_sig_cancer": {"pass": per["PSCA"]["n_sig_cancer"] >= 1,
                               "detail": {"n_sig_cancer": per["PSCA"]["n_sig_cancer"],
                                          "traits": sorted(per["PSCA"]["cancer_traits"])}},
        "G2_catalog_rows_ge_400k": {"pass": n_rows >= 400000, "n_rows": n_rows},
        "G3_egfr_sig_cancer": {"pass": egfr["n_sig_cancer"] >= 1, "detail": egfr},
    }
    grp = lambda k: [g for g in universe if genes[g] == k]
    gated, bg = grp("gated"), grp("background")

    def frac_test(field):
        a = sum(1 for g in gated if per[g][field] >= 1)
        b = len(gated) - a
        c = sum(1 for g in bg if per[g][field] >= 1)
        d = len(bg) - c
        orr, p = fisher_exact([[a, b], [c, d]])
        return {"gated": a, "gated_n": len(gated), "bg": c, "bg_n": len(bg),
                "fisher_or": orr, "fisher_p": p}

    import statistics as st
    sg = [per[g]["n_sig"] for g in gated]
    sb = [per[g]["n_sig"] for g in bg]
    _, p_burden = mannwhitneyu(sg, sb)
    out = {
        "tool": "GWAS Catalog full ontology-annotated associations (ftp.ebi.ac.uk releases/latest)",
        "n_genes_total": len(genes), "n_assoc_rows": n_rows,
        "sig_threshold": "p<=5e-8 (PVALUE_MLOG>=7.301)",
        "gates": gates, "all_gates_pass": all(g["pass"] for g in gates.values()),
        "H1_any_sig": frac_test("n_sig"),
        "H2_cancer_sig": frac_test("n_sig_cancer"),
        "H3_burden": {"gated_median": st.median(sg), "bg_median": st.median(sb),
                      "mw_p": p_burden},
        "and_gate_antigens": {a: {"n_assoc": per[a]["n_assoc"], "n_sig": per[a]["n_sig"],
                                  "n_sig_cancer": per[a]["n_sig_cancer"],
                                  "n_sig_snps": len(per[a]["sig_snps"]),
                                  "cancer_traits": sorted(per[a]["cancer_traits"])}
                              for a in ANTS},
        "egfr_sentinel": egfr,
    }
    with open("results/gwas_audit.json", "w") as f:
        json.dump(out, f, indent=1)
    with open("results/gwas_per_gene.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "group", "n_assoc", "n_sig", "n_sig_snps", "n_sig_cancer", "cancer_traits"])
        for g in sorted(per):
            p = per[g]
            w.writerow([g, genes[g], p["n_assoc"], p["n_sig"], len(p["sig_snps"]),
                        p["n_sig_cancer"], "; ".join(sorted(p["cancer_traits"]))])
    print(json.dumps({"gates": {k: v["pass"] for k, v in gates.items()},
                      "n_rows": n_rows, "H1": out["H1_any_sig"],
                      "H2": out["H2_cancer_sig"], "H3": out["H3_burden"],
                      "antigens": {a: v["n_sig_cancer"] for a, v in out["and_gate_antigens"].items()}},
                     indent=1))


if __name__ == "__main__":
    main()

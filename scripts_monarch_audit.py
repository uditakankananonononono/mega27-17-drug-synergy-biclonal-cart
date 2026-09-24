#!/usr/bin/env python3
"""Monarch v3 Mendelian-disease / phenotype audit: do gated antigens differ from
background in causal gene-disease associations (OMIM/Orphanet), correlated
associations, or HPO phenotypic breadth? Rare/Mendelian complement to the GWAS
common-variant germline audit. Reads committed caches in data/monarch/ only.

Pre-registered calibration gates:
  G1: >= 180/205 genes resolved to an HGNC id with a clean Monarch record
  G2: CD19 sentinel carries >= 1 causal gene-disease association (OMIM CVID3)
  G3: CD19 sentinel carries >= 10 HPO phenotype annotations
Gate values calibrated from a recorded connectivity probe of the live API
(CD19 observed: causal 1, phenotypes 42) before the full fetch.
Outputs: results/monarch_audit.json, results/monarch_per_gene.csv
"""
import csv, json, os, re
from scipy.stats import mannwhitneyu, fisher_exact

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "monarch")
ANTS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]
CANCER = re.compile(r"carcinoma|cancer|neoplasm|tumou?r|sarcoma|lymphoma|"
                    r"leukemia|leukaemia|melanoma|blastoma|malignan|adenocarcinoma", re.I)


def gene_set():
    gs = json.load(open("results/depmap_gene_sets.json"))
    return {g: grp for grp in ("gated", "background", "references", "positive_controls")
            for g in gs[grp]}


def load_records():
    recs = {}
    for f in os.listdir(OUT):
        if f.endswith(".json") and f != "manifest.json":
            r = json.load(open(os.path.join(OUT, f)))
            recs[r["gene"]] = r
    return recs


def cancer_causal(r):
    return [it["label"] for it in r.get("CausalGeneToDiseaseAssociation_items", [])
            if it.get("label") and CANCER.search(it["label"])]


def main():
    groups = gene_set()
    recs = load_records()
    rows = []
    for g, grp in groups.items():
        r = recs.get(g)
        if not r or not r.get("hgnc_id") or "_error" in r:
            rows.append({"gene": g, "group": grp, "hgnc_id": (r or {}).get("hgnc_id"),
                         "resolved": False, "causal_total": None, "correlated_total": None,
                         "pheno_total": None, "n_cancer_causal": None, "causal_labels": ""})
            continue
        cc = cancer_causal(r)
        labels = [it["label"] for it in r.get("CausalGeneToDiseaseAssociation_items", [])]
        rows.append({"gene": g, "group": grp, "hgnc_id": r["hgnc_id"], "resolved": True,
                     "causal_total": r.get("CausalGeneToDiseaseAssociation_total", 0),
                     "correlated_total": r.get("CorrelatedGeneToDiseaseAssociation_total", 0),
                     "pheno_total": r.get("GeneToPhenotypicFeatureAssociation_total", 0),
                     "n_cancer_causal": len(cc),
                     "causal_labels": "; ".join(labels[:12])})

    res = [r for r in rows if r["resolved"]]
    gated = [r for r in res if r["group"] == "gated"]
    bg = [r for r in res if r["group"] == "background"]
    refs = {r["gene"]: r for r in res if r["group"] == "references"}

    cd19 = recs.get("CD19", {})
    gates = {
        "G1_coverage_ge_180": {"pass": len(res) >= 180, "resolved": len(res),
                               "note": "pre-registered"},
        "G2_cd19_causal_ge_1": {"pass": cd19.get("CausalGeneToDiseaseAssociation_total", 0) >= 1,
                                "cd19_causal": cd19.get("CausalGeneToDiseaseAssociation_total"),
                                "note": "pre-registered; OMIM CVID3 expected"},
        "G3_cd19_pheno_ge_10": {"pass": cd19.get("GeneToPhenotypicFeatureAssociation_total", 0) >= 10,
                                "cd19_pheno": cd19.get("GeneToPhenotypicFeatureAssociation_total"),
                                "note": "pre-registered; probe observed 42"},
    }
    all_pass = all(g["pass"] for g in gates.values())

    def carriage(rs, key):
        return sum(1 for r in rs if (r[key] or 0) >= 1)

    a, c = carriage(gated, "causal_total"), carriage(bg, "causal_total")
    or1, p1 = fisher_exact([[a, len(gated) - a], [c, len(bg) - c]])
    _, p2 = mannwhitneyu([r["pheno_total"] for r in gated],
                         [r["pheno_total"] for r in bg])
    a3, c3 = carriage(gated, "correlated_total"), carriage(bg, "correlated_total")
    or3, p3 = fisher_exact([[a3, len(gated) - a3], [c3, len(bg) - c3]])
    a4, c4 = carriage(gated, "n_cancer_causal"), carriage(bg, "n_cancer_causal")
    or4, p4 = fisher_exact([[a4, len(gated) - a4], [c4, len(bg) - c4]])

    import statistics as st
    audit = {
        "source": "Monarch Initiative v3 API (api.monarchinitiative.org/v3/api/association; "
                  "biolink:CausalGeneToDiseaseAssociation, CorrelatedGeneToDiseaseAssociation, "
                  "GeneToPhenotypicFeatureAssociation)",
        "n_genes_total": len(rows), "n_resolved": len(res),
        "gates": gates, "all_gates_pass": all_pass,
        "H1_causal_carriage": {"gated": a, "gated_n": len(gated), "bg": c, "bg_n": len(bg),
                               "fisher_or": or1, "fisher_p": p1},
        "H2_pheno_breadth": {"gated_median": st.median(r["pheno_total"] for r in gated),
                             "bg_median": st.median(r["pheno_total"] for r in bg), "mw_p": p2},
        "H3_correlated_carriage": {"gated": a3, "gated_n": len(gated), "bg": c3, "bg_n": len(bg),
                                   "fisher_or": or3, "fisher_p": p3},
        "H4_cancer_causal_carriage": {"gated": a4, "gated_n": len(gated), "bg": c4,
                                      "bg_n": len(bg), "fisher_or": or4, "fisher_p": p4},
        "and_gate_antigens": {},
        "reference_genes": {g: {"causal": refs[g]["causal_total"],
                                "pheno": refs[g]["pheno_total"]} for g in refs},
    }
    for s in ANTS:
        r = recs.get(s, {})
        audit["and_gate_antigens"][s] = {
            "causal_total": r.get("CausalGeneToDiseaseAssociation_total", 0),
            "causal_labels": [it["label"] for it in
                              r.get("CausalGeneToDiseaseAssociation_items", [])],
            "correlated_total": r.get("CorrelatedGeneToDiseaseAssociation_total", 0),
            "pheno_total": r.get("GeneToPhenotypicFeatureAssociation_total", 0),
            "cancer_causal_labels": cancer_causal(r)}

    with open("results/monarch_audit.json", "w") as f:
        json.dump(audit, f, indent=1)
    with open("results/monarch_per_gene.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "group", "hgnc_id", "resolved",
                                          "causal_total", "correlated_total", "pheno_total",
                                          "n_cancer_causal", "causal_labels"])
        w.writeheader()
        w.writerows(rows)
    print("H1 causal carriage: gated %d/%d vs bg %d/%d OR=%.2f p=%.3g" % (a, len(gated), c, len(bg), or1, p1))
    print("H2 pheno breadth MW p=%.3g" % p2)
    print("H3 correlated carriage: gated %d/%d vs bg %d/%d OR=%.2f p=%.3g" % (a3, len(gated), c3, len(bg), or3, p3))
    print("H4 cancer-causal carriage: gated %d/%d vs bg %d/%d p=%.3g" % (a4, len(gated), c4, len(bg), p4))
    print("gates:", {k: v["pass"] for k, v in gates.items()}, "ALL", all_pass)


if __name__ == "__main__":
    main()

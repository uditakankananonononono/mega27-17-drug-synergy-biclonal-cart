#!/usr/bin/env python3
"""Harmonizome curated-vs-systematic annotation-density audit.

Meta-finding under test (from IntAct / HPA-IF / BioPlex audits): the gated
panel looks 'study-biased upward' in curated resources, while systematic
measurement erases the difference. Harmonizome's ~174 datasets let us test this
on one independent integrator with a fixed, pre-registered name-based split.

Pre-registered hypothesis H1 (two parts, BOTH required to confirm):
  H1a gated > background in curated/literature association count
      (Mann-Whitney one-sided, p < 0.05)
  H1b no gated-vs-background difference in systematic association count
      (Mann-Whitney two-sided, p >= 0.05)
Falsified if H1a fails or H1b fails.
Calibration gates (from the recorded CD19 probe, 11,116 associations):
  G1 >= 180/205 genes resolved with a clean record
  G2 CD19 sentinel total associations >= 1,000
  G3 all 4 housekeeping positive controls resolved with >= 1,000 associations
Datasets matching neither keyword list are excluded and listed.
Outputs: results/harmonizome_audit.json, results/harmonizome_per_gene.csv
"""
import csv, json, os
from scipy.stats import mannwhitneyu

BASE = os.path.dirname(os.path.abspath(__file__))
GDIR = os.path.join(BASE, "data", "harmonizome", "genes")
ANTS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]
CURATED_KEYS = ["Curated", "Text-mining", "Textmining", "GeneRIF", "GO ", "KEGG",
                "Reactome", "WikiPathways", "Biocarta", "PID Pathways", "PANTHER",
                "HumanCyc", "OMIM", "ClinVar", "HPO ", "GAD ", "HuGE", "CTD ",
                "DisGeNET", "DrugBank", "DGIdb", "Guide to Pharmacology", "CORUM",
                "Pathway Commons", "GeneSigDB", "PFOCR", "PhosphoSitePlus", "MGI ",
                "MPO ", "dbGAP", "GWAS", "CellMarker", "HuBMAP ASCT", "SynGO", "NURSA",
                "Hub Proteins", "DEPOD", "KEA ", "HMDB", "MW Enzyme", "Virus MINT",
                "GlyGen", "DeepCoverMOA"]
SYSTEMATIC_KEYS = ["Expression Profiles", "CNV", "Mutation", "Achilles", "DepMap",
                   "Dependency", "ENCODE", "Roadmap", "Signatures", "Perturb", "L1000",
                   "CMAP", "Proteomics", "eQTL", "Experimental", "Tabula", "BioGPS",
                   "Predicted", "MotifMap", "TargetScan", "ChEA", "KinomeScan",
                   "Kinativ", "SILAC", "Kinase Library", "CM4AI", "Azimuth", "MoTrPAC",
                   "IMPC", "KnockTF", "Sci-Plex", "Tahoe", "Rummagene", "RummaGEO",
                   "Carcinogenome", "ESCAPE", "InterPro", "GTEx", "HPA ", "HPM ",
                   "Allen Brain"]


def classify(ds):
    # curated first: 'Curated'/'Text-mining' variants of TISSUES/COMPARTMENTS/DISEASES
    # must not be caught by a generic systematic key
    if any(k in ds for k in CURATED_KEYS):
        return "curated"
    if any(k in ds for k in SYSTEMATIC_KEYS):
        return "systematic"
    return "unclassified"


def gene_set():
    gs = json.load(open(os.path.join(BASE, "results", "depmap_gene_sets.json")))
    return {g: grp for grp in ("gated", "background", "references", "positive_controls")
            for g in gs[grp]}


def per_gene(rec):
    out = {"curated": 0, "systematic": 0, "unclassified": 0,
           "n_datasets_curated": 0, "n_datasets_systematic": 0}
    for ds, (up, dn) in rec.get("per_dataset", {}).items():
        c = classify(ds)
        out[c] += up + dn
        if c != "unclassified" and up + dn:
            out["n_datasets_" + c] += 1
    tot = out["curated"] + out["systematic"]
    out["curated_frac"] = out["curated"] / tot if tot else None
    return out


def mw(a, b, alt):
    a = [x for x in a if x is not None]; b = [x for x in b if x is not None]
    if not a or not b:
        return None
    return float(mannwhitneyu(a, b, alternative=alt).pvalue)


def median(v):
    v = sorted(x for x in v if x is not None)
    if not v:
        return None
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def main():
    groups = gene_set()
    rows, dsets = [], {}
    for g, grp in groups.items():
        p = os.path.join(GDIR, g + ".json")
        rec = json.load(open(p)) if os.path.exists(p) else {"_error": "missing"}
        ok = "_error" not in rec and rec.get("n_associations", 0) > 0
        row = {"gene": g, "group": grp, "resolved": ok,
               "n_associations": rec.get("n_associations") if ok else None}
        if ok:
            row.update(per_gene(rec))
            for ds in rec["per_dataset"]:
                dsets[ds] = classify(ds)
        rows.append(row)
    res = {r["gene"]: r for r in rows}
    ok = [r for r in rows if r["resolved"]]
    G = [r for r in ok if r["group"] == "gated"]
    B = [r for r in ok if r["group"] == "background"]
    cd19 = res.get("CD19", {})
    pcs = [res[g] for g in groups if groups[g] == "positive_controls"]
    gates = {
        "G1_resolved": {"obs": len(ok), "threshold": 180, "pass": len(ok) >= 180},
        "G2_CD19_total": {"obs": cd19.get("n_associations"), "threshold": 1000,
                          "pass": (cd19.get("n_associations") or 0) >= 1000},
        "G3_posctrl_resolved_ge1000": {
            "obs": {r["gene"]: r.get("n_associations") for r in pcs}, "threshold": 1000,
            "pass": all((r.get("n_associations") or 0) >= 1000 for r in pcs)}}
    tests = {}
    for key, alt in (("curated", "greater"), ("systematic", "two-sided"),
                     ("n_datasets_curated", "greater"), ("n_datasets_systematic", "two-sided"),
                     ("curated_frac", "two-sided"), ("n_associations", "two-sided")):
        a = [r.get(key) for r in G]; b = [r.get(key) for r in B]
        tests[key] = {"alternative": alt, "median_gated": median(a),
                      "median_background": median(b), "n_gated": len(G),
                      "n_background": len(B), "p": mw(a, b, alt)}
    h1a = tests["curated"]["p"] is not None and tests["curated"]["p"] < 0.05
    h1b = tests["systematic"]["p"] is not None and tests["systematic"]["p"] >= 0.05
    verdict = ("CONFIRMED" if h1a and h1b else
               "FALSIFIED (" + ", ".join(x for x, ok_ in (("H1a curated not elevated", h1a),
                                                          ("H1b systematic differs", h1b))
                                         if not ok_) + ")")
    out = {"source": "Harmonizome API 1.0 (maayanlab.cloud/Harmonizome)",
           "calibration_gates": gates, "all_gates_pass": all(v["pass"] for v in gates.values()),
           "hypothesis": "H1a gated>background curated (one-sided p<0.05) AND "
                         "H1b systematic no difference (two-sided p>=0.05)",
           "tests": tests, "H1a": h1a, "H1b": h1b, "verdict": verdict,
           "dataset_classes": {c: sum(1 for v in dsets.values() if v == c)
                               for c in ("curated", "systematic", "unclassified")},
           "unclassified_datasets": sorted(d for d, c in dsets.items() if c == "unclassified"),
           "antigens": {a: {k: res[a].get(k) for k in ("n_associations", "curated",
                                                        "systematic", "curated_frac")}
                        for a in ANTS if a in res},
           "unresolved": sorted(r["gene"] for r in rows if not r["resolved"])}
    json.dump(out, open(os.path.join(BASE, "results", "harmonizome_audit.json"), "w"), indent=1)
    cols = ["gene", "group", "resolved", "n_associations", "curated", "systematic",
            "unclassified", "n_datasets_curated", "n_datasets_systematic", "curated_frac"]
    with open(os.path.join(BASE, "results", "harmonizome_per_gene.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(json.dumps({k: out[k] for k in ("calibration_gates", "verdict", "dataset_classes")}))
    return out


if __name__ == "__main__":
    main()

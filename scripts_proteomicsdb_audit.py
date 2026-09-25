#!/usr/bin/env python3
"""ProteomicsDB protein-level normal-tissue detection audit.

Gated antigens were chosen for tumour-restricted RNA expression. Does mass
spectrometry agree at the protein level, i.e. are gated proteins detected in
fewer distinct normal human tissues than background proteins?

Normal-tissue sample rule (fixed before any group comparison): TAXCODE 9606,
empty DISEASE field, and TISSUE in the explicit NORMAL_TISSUES whitelist of
human organ/tissue names (cell lines, body fluids, stem cells and sorted
immune/neural cell types excluded).

Pre-registered hypothesis H1: gated < background in distinct normal tissues
with detection (Mann-Whitney one-sided 'less', p < 0.05).
Secondary H2 (two-sided): fraction of genes detected in any sample at all.
Calibration gates:
  G1 >= 180/205 genes mapped to a ProteomicsDB target protein
  G2 CD19 detected in >= 1 B-lineage sample (recorded probe: 1,119 rows,
     mostly lymphoblastoid / RAMOS / SU-DHL-4 / ARH-77)
  G3 each of the 4 housekeeping controls detected in >= 10 normal tissues
Outputs: results/proteomicsdb_audit.json, results/proteomicsdb_per_gene.csv
"""
import csv, json, os, re
from scipy.stats import mannwhitneyu, fisher_exact

BASE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(BASE, "data", "proteomicsdb")
ANTS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]
BLINEAGE = re.compile(r"lymphoblastoid|B-lymph|RAMOS|SU-DHL|ARH-77|lymph node|"
                      r"OCI-LY|JeKo|NU-DHL|Reh|SEM", re.I)
NORMAL_TISSUES = {
    "adipose tissue", "adrenal gland", "anus", "bone", "bone marrow", "brain",
    "brain stem", "breast", "cardia", "cerebellum", "cerebral cortex",
    "cervical mucosa", "colon", "colon ascendens", "colon descendens",
    "colon muscle", "corpus callosum", "corpus striatum", "duodenum",
    "epididymis", "esophagus", "eye", "frontal lobe", "gall bladder", "gut",
    "hair follicle", "heart", "hindbrain", "hippocampus", "kidney",
    "large intestine", "liver", "lung", "lymph node", "mammary gland",
    "motor cortex", "auditory cortex", "somatosensory cortex", "nasopharynx",
    "occipital lobe", "olfactory bulb", "optic nerve", "oral epithelium",
    "ovary", "oviduct", "pancreas", "pancreatic islet", "placenta",
    "prefrontal cortex", "prostate gland", "rectum", "retina",
    "salivary gland", "seminal vesicle", "skeletal muscle", "skin",
    "small intestine", "smooth muscle", "spinal cord", "spleen", "stomach",
    "temporal lobe", "testis", "thalamus", "thymus", "thyroid gland",
    "tongue", "tonsil", "trachea", "urinary bladder", "uterine cervix",
    "uterine endometrium", "uterus", "vagina", "vermiform appendix"}


def is_normal(s):
    return (s.get("TAXCODE") == 9606 and not (s.get("DISEASE") or "")
            and (s.get("TISSUE") or "") in NORMAL_TISSUES)


def gene_set():
    gs = json.load(open(os.path.join(BASE, "results", "depmap_gene_sets.json")))
    return {g: grp for grp in ("gated", "background", "references", "positive_controls")
            for g in gs[grp]}


def median(v):
    v = sorted(x for x in v if x is not None)
    if not v:
        return None
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def per_gene(rec, samples):
    """Detection = any ProteinExpression row with PEPTIDES >= 1 for the sample."""
    det = {sid for sid, _expr, _m, pep in rec.get("expression", []) if (pep or 0) >= 1}
    norm = {samples[s]["TISSUE"] for s in det if s in samples and is_normal(samples[s])}
    blin = sum(1 for s in det if s in samples and BLINEAGE.search(samples[s]["TISSUE"] or ""))
    return {"n_samples_detected": len(det), "n_normal_tissues": len(norm),
            "n_blineage_samples": blin, "detected_any": len(det) > 0,
            "normal_tissues": ";".join(sorted(norm))}


def main():
    groups = gene_set()
    samples = {s["SAMPLE_ID"]: s for s in json.load(open(os.path.join(D, "samples.json")))}
    rows = []
    for g, grp in groups.items():
        p = os.path.join(D, "genes", g + ".json")
        rec = json.load(open(p)) if os.path.exists(p) else {"_error": "missing"}
        ok = "_error" not in rec and bool(rec.get("protein_ids"))
        row = {"gene": g, "group": grp, "mapped": ok, "uniprot": rec.get("uniprot")}
        if ok:
            row.update(per_gene(rec, samples))
        rows.append(row)
    res = {r["gene"]: r for r in rows}
    G = [r for r in rows if r["mapped"] and r["group"] == "gated"]
    B = [r for r in rows if r["mapped"] and r["group"] == "background"]
    pcs = [res[g] for g, grp in groups.items() if grp == "positive_controls"]
    n_mapped = sum(r["mapped"] for r in rows)
    gates = {
        "G1_mapped": {"obs": n_mapped, "threshold": 180, "pass": n_mapped >= 180},
        "G2_CD19_blineage": {"obs": res["CD19"].get("n_blineage_samples"), "threshold": 1,
                             "pass": (res["CD19"].get("n_blineage_samples") or 0) >= 1},
        "G3_posctrl_normal_tissues": {"obs": {r["gene"]: r.get("n_normal_tissues") for r in pcs},
                                      "threshold": 10,
                                      "pass": all((r.get("n_normal_tissues") or 0) >= 10 for r in pcs)}}
    a = [r["n_normal_tissues"] for r in G]; b = [r["n_normal_tissues"] for r in B]
    p1 = float(mannwhitneyu(a, b, alternative="less").pvalue) if a and b else None
    tab = [[sum(r["detected_any"] for r in G), sum(not r["detected_any"] for r in G)],
           [sum(r["detected_any"] for r in B), sum(not r["detected_any"] for r in B)]]
    or2, p2 = fisher_exact(tab) if G and B else (None, None)
    out = {"source": "ProteomicsDB api_v2 OData (proteomicsdb.org)",
           "n_genes_total": len(rows), "calibration_gates": gates,
           "all_gates_pass": all(v["pass"] for v in gates.values()),
           "n_normal_tissue_whitelist": len(NORMAL_TISSUES),
           "H1_normal_breadth": {"alternative": "less", "median_gated": median(a),
                                 "median_background": median(b), "n_gated": len(a),
                                 "n_background": len(b), "p": p1,
                                 "verdict": "CONFIRMED" if p1 is not None and p1 < 0.05 else "FALSIFIED"},
           "H2_detected_any": {"table_gated_bg_detected_not": tab,
                               "fisher_or": float(or2) if or2 is not None else None,
                               "fisher_p": float(p2) if p2 is not None else None},
           "antigens": {s: {k: res[s].get(k) for k in ("n_samples_detected",
                                                        "n_normal_tissues", "normal_tissues")}
                        for s in ANTS},
           "unmapped": sorted(r["gene"] for r in rows if not r["mapped"])}
    json.dump(out, open(os.path.join(BASE, "results", "proteomicsdb_audit.json"), "w"), indent=1)
    cols = ["gene", "group", "mapped", "uniprot", "detected_any", "n_samples_detected",
            "n_normal_tissues", "n_blineage_samples", "normal_tissues"]
    with open(os.path.join(BASE, "results", "proteomicsdb_per_gene.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(json.dumps({"gates": gates, "H1": out["H1_normal_breadth"], "H2": out["H2_detected_any"]}))
    return out


if __name__ == "__main__":
    main()


# ---- Follow-up control H3 (declared after H1 came out reversed) ----
# Explanation under test: background contains many proteins MS rarely sees
# (IG/TCR segments, olfactory receptors, small secreted proteins), so a gated
# excess may be pure detectability. H3 compares normal-tissue breadth only
# among genes detected in >= 1 sample (two-sided MW); a remaining gated excess
# means the reversal is not just detectability.
def h3():
    rows = list(csv.DictReader(open(os.path.join(BASE, "results", "proteomicsdb_per_gene.csv"))))
    det = [r for r in rows if r["detected_any"] == "True"]
    a = [int(r["n_normal_tissues"]) for r in det if r["group"] == "gated"]
    b = [int(r["n_normal_tissues"]) for r in det if r["group"] == "background"]
    p = float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
    return {"hypothesis": "H3 among detected genes, gated vs background normal-tissue "
                          "breadth (two-sided MW)",
            "median_gated": median(a), "median_background": median(b),
            "n_gated": len(a), "n_background": len(b), "p": p,
            "gated_excess_survives": p < 0.05 and median(a) > median(b)}


if __name__ == "__main__":
    _p = os.path.join(BASE, "results", "proteomicsdb_audit.json")
    _o = json.load(open(_p))
    _o["H3_detected_only"] = h3()
    json.dump(_o, open(_p, "w"), indent=1)
    print(json.dumps(_o["H3_detected_only"]))

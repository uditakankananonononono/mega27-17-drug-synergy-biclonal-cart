#!/usr/bin/env python3
"""Expression Atlas (E-MTAB-513, Illumina Body Map 2.0) tissue-specificity audit.

Gated antigens were selected for tumour-restricted RNA using HPA/GTEx windows.
ProteomicsDB showed the proteins are MS-detected in MORE normal tissues than
background. Does an independent RNA-seq atlas reproduce RNA-level restriction?

Pre-registered H1: tissue specificity tau (log2 TPM+1, 16 tissues) is higher
in gated than background genes (Mann-Whitney one-sided, p < 0.05).
Pre-registered H2 (cross-platform): across mapped genes, Body Map breadth
(tissues with TPM >= 1) correlates with ProteomicsDB normal-tissue breadth
(Spearman rho > 0.3 and p < 0.05).
Calibration gates (G2 threshold set from a recorded probe row, CD19 lymph node 70 TPM):
  G1 >= 180/205 genes matched by Ensembl gene id (HGNC cache) or symbol
  G2 CD19 lymph-node TPM >= 10
  G3 all 4 housekeeping controls tau < 0.5 (broadly expressed)
Outputs: results/expression_atlas_audit.json, results/expression_atlas_per_gene.csv
"""
import csv, json, math, os, re
from scipy.stats import mannwhitneyu, spearmanr

BASE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(BASE, "data", "expression_atlas")
ANTS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]


def tissues():
    x = open(os.path.join(D, "E-MTAB-513-configuration.xml")).read()
    return dict(re.findall(r'<assay_group id="(g\d+)" label="([^"]+)"', x))


def load_tpm():
    """{ensg: (symbol, {tissue: tpm})}; replicate cells 'a,b,c' -> median."""
    lab = tissues()
    out = {}
    with open(os.path.join(D, "E-MTAB-513-tpms.tsv")) as f:
        rows = (l for l in f if not l.startswith("#"))
        hdr = next(rows).rstrip("\n").split("\t")
        for l in rows:
            c = l.rstrip("\n").split("\t")
            vals = {}
            for h, v in zip(hdr[2:], c[2:]):
                xs = sorted(float(t) for t in v.split(",") if t != "")
                vals[lab[h]] = xs[len(xs) // 2] if xs else 0.0
            out[c[0]] = (c[1], vals)
    return out


def tau(vals):
    """Yanai tissue-specificity index on log2(TPM+1); 0 broad, 1 specific."""
    x = [math.log2(v + 1) for v in vals]
    m = max(x)
    if m == 0:
        return None
    return sum(1 - xi / m for xi in x) / (len(x) - 1)


def ensg_for(sym):
    p = os.path.join(BASE, "data", "constraint", "hgnc_%s.json" % sym)
    try:
        return json.load(open(p))["response"]["docs"][0].get("ensembl_gene_id")
    except Exception:
        return None


def gene_set():
    gs = json.load(open(os.path.join(BASE, "results", "depmap_gene_sets.json")))
    return {g: grp for grp in ("gated", "background", "references", "positive_controls")
            for g in gs[grp]}


def median(v):
    v = sorted(x for x in v if x is not None)
    n = len(v)
    return None if not n else (v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2]))


def main():
    groups = gene_set()
    tpm = load_tpm()
    by_sym = {s: e for e, (s, _v) in tpm.items()}
    rows = []
    for g, grp in groups.items():
        e = ensg_for(g)
        e = e if e in tpm else by_sym.get(g)
        row = {"gene": g, "group": grp, "ensembl": e, "matched": e is not None}
        if e is not None:
            v = tpm[e][1]
            row["tau"] = tau(list(v.values()))
            row["n_tissues_tpm1"] = sum(1 for x in v.values() if x >= 1)
            row["max_tissue"] = max(v, key=v.get)
            row["max_tpm"] = max(v.values())
            row["lymph_node_tpm"] = v.get("lymph node")
        rows.append(row)
    res = {r["gene"]: r for r in rows}
    G = [r for r in rows if r["matched"] and r["group"] == "gated"]
    B = [r for r in rows if r["matched"] and r["group"] == "background"]
    pcs = [res[g] for g, grp in groups.items() if grp == "positive_controls"]
    nm = sum(r["matched"] for r in rows)
    gates = {"G1_matched": {"obs": nm, "threshold": 180, "pass": nm >= 180},
             "G2_CD19_lymph_node_tpm": {"obs": res["CD19"].get("lymph_node_tpm"),
                                        "threshold": 10,
                                        "pass": (res["CD19"].get("lymph_node_tpm") or 0) >= 10},
             "G3_posctrl_tau_lt_0.5": {"obs": {r["gene"]: r.get("tau") for r in pcs},
                                       "threshold": 0.5,
                                       "pass": all(r.get("tau") is not None and r["tau"] < 0.5
                                                   for r in pcs)}}
    a = [r["tau"] for r in G if r.get("tau") is not None]
    b = [r["tau"] for r in B if r.get("tau") is not None]
    p1 = float(mannwhitneyu(a, b, alternative="greater").pvalue)
    pdb = {r["gene"]: r for r in csv.DictReader(open(os.path.join(
        BASE, "results", "proteomicsdb_per_gene.csv"))) if r["mapped"] == "True"}
    pairs = [(r["n_tissues_tpm1"], int(pdb[r["gene"]]["n_normal_tissues"]))
             for r in rows if r["matched"] and r["gene"] in pdb]
    rho, p2 = spearmanr([x for x, _ in pairs], [y for _, y in pairs])
    rho, p2 = float(rho), float(p2)
    out = {"source": "EBI Expression Atlas E-MTAB-513 (Illumina Body Map 2.0, 16 tissues)",
           "n_genes_total": len(rows), "calibration_gates": gates,
           "all_gates_pass": all(v["pass"] for v in gates.values()),
           "H1_tau": {"alternative": "greater", "median_gated": median(a),
                      "median_background": median(b), "n_gated": len(a),
                      "n_background": len(b), "p": p1,
                      "verdict": "CONFIRMED" if p1 < 0.05 else "FALSIFIED"},
           "H2_cross_platform": {"n_pairs": len(pairs), "spearman_rho": rho, "p": p2,
                                 "verdict": "CONFIRMED" if rho > 0.3 and p2 < 0.05
                                 else "FALSIFIED"},
           "antigens": {s: {k: res[s].get(k) for k in ("tau", "n_tissues_tpm1",
                                                        "max_tissue", "max_tpm")}
                        for s in ANTS},
           "unmatched": sorted(r["gene"] for r in rows if not r["matched"])}
    json.dump(out, open(os.path.join(BASE, "results", "expression_atlas_audit.json"), "w"),
              indent=1)
    cols = ["gene", "group", "ensembl", "matched", "tau", "n_tissues_tpm1",
            "max_tissue", "max_tpm", "lymph_node_tpm"]
    with open(os.path.join(BASE, "results", "expression_atlas_per_gene.csv"), "w",
              newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(json.dumps({k: out[k] for k in ("calibration_gates", "H1_tau",
                                          "H2_cross_platform", "unmatched")}))
    return out


if __name__ == "__main__":
    main()

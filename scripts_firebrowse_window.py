#!/usr/bin/env python3
"""Tumor vs adjacent-normal therapeutic window of the AND-gate CAR-T antigens
(FireBrowse RSEM log2, 5 TCGA cohorts matching the CPTAC protein audit).

Question: within the organ itself, does the tumor overexpress each antigen
versus its own adjacent-normal tissue? GTEx/HPA answer "other organs"; this
answers "same organ, field-effect included" - the hardest background a local
CAR must beat.

Per cohort x gene: median log2 RSEM in primary tumor (TP) and solid-tissue
normal (NT), delta = median(TP) - median(NT), Mann-Whitney p. Kept when
n_nt >= 4 and n_tp >= 20 (PAAD and GBM normals are n=4-5; reported, not hidden).
Controls: GAPDH must be flat (|delta| < 0.75 everywhere); ERBB2 must be up in
BRCA; CA9 (hypoxia) broadly up. Outputs: results/firebrowse_window.json,
results/firebrowse_window_rows.csv (participant-level records)."""
import csv, json, statistics
from scipy.stats import mannwhitneyu

GENES = {"targets": ["CLDN18", "MSLN", "CA9", "CA12", "CLDN6", "PSCA", "SLC39A6"],
         "references": ["MS4A1", "TNFRSF17", "ERBB2"], "housekeeping": ["GAPDH", "ACTB", "RPLP0", "TBP"]}
ALL = GENES["targets"] + GENES["references"] + GENES["housekeeping"]
COHORTS = ["PAAD", "COAD", "LUAD", "GBM", "BRCA"]

vals = {}   # (cohort, gene, stype) -> [expr]
rows = []
for c in COHORTS:
    for r in json.load(open(f"data/firebrowse/{c}.json")) + json.load(open(f"data/firebrowse/hk_{c}.json")):
        if r.get("expression_log2") is None:
            continue
        key = (r["cohort"], r["gene"], r["sample_type"])
        vals.setdefault(key, []).append(r["expression_log2"])
        rows.append({"cohort": r["cohort"], "gene": r["gene"],
                     "sample_type": r["sample_type"],
                     "participant": r["tcga_participant_barcode"],
                     "expression_log2_rsem": round(r["expression_log2"], 4)})

cells = {}
for c in COHORTS:
    for g in ALL:
        tp, nt = vals.get((c, g, "TP"), []), vals.get((c, g, "NT"), [])
        if len(nt) < 4 or len(tp) < 20:
            cells[f"{c}:{g}"] = None
            continue
        mt, mn = statistics.median(tp), statistics.median(nt)
        p = float(mannwhitneyu(tp, nt, alternative="greater").pvalue)
        cells[f"{c}:{g}"] = {"n_tp": len(tp), "n_nt": len(nt),
                             "median_tp": round(mt, 3), "median_nt": round(mn, 3),
                             "delta": round(mt - mn, 3), "mw_p_greater": p}

offset = {}
for c in COHORTS:
    hk_ds = [cells[f"{c}:{g}"]["delta"] for g in GENES["housekeeping"] if cells.get(f"{c}:{g}")]
    offset[c] = round(statistics.median(hk_ds), 3) if hk_ds else 0.0
for k, v in cells.items():
    if v:
        c = k.split(":")[0]
        v["delta_adj"] = round(v["delta"] - offset[c], 3)

per_gene = {}
for g in ALL:
    ok = {c: cells[f"{c}:{g}"] for c in COHORTS if cells[f"{c}:{g}"]}
    ds = [v["delta_adj"] for v in ok.values()]
    per_gene[g] = {"n_cohorts": len(ok),
                   "cohorts_delta_ge1": sum(1 for d in ds if d >= 1.0),
                   "cohorts_delta_ge1_sig": sum(1 for c, v in ok.items() if v["delta_adj"] >= 1.0 and v["mw_p_greater"] < 0.05),
                   "median_delta": round(statistics.median(ds), 3) if ds else None,
                   "min_nt_median": round(min(v["median_nt"] for v in ok.values()), 3) if ok else None,
                   "cells": ok}

sample_records = {(r["cohort"], r["participant"], r["sample_type"]) for r in rows}
summary = {
    "source": "FireBrowse (Broad Firehose) REST API, /Samples/mRNASeq, RSEM log2, TP+NT",
    "genes": GENES, "cohorts": COHORTS,
    "n_value_records": len(rows),
    "n_sample_records": len(sample_records),
    "per_gene": per_gene,
    "housekeeping_offset_per_cohort": offset,
    "controls": {
        "hk_raw_delta_range": [round(min(v["delta"] for k, v in cells.items() if v and k.split(":")[1] in GENES["housekeeping"]), 3),
                               round(max(v["delta"] for k, v in cells.items() if v and k.split(":")[1] in GENES["housekeeping"]), 3)],
        "hk_adj_max_abs_delta": max(abs(v["delta_adj"]) for k, v in cells.items() if v and k.split(":")[1] in GENES["housekeeping"]),
        "erbb2_brca_delta_adj": cells["BRCA:ERBB2"]["delta_adj"],
        "ca9_median_delta_adj": per_gene["CA9"]["median_delta"],
    },
}
json.dump(summary, open("results/firebrowse_window.json", "w"), indent=1)
with open("results/firebrowse_window_rows.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

print("sample records:", len(sample_records), " value records:", len(rows))
for g in ALL:
    pg = per_gene[g]
    print(f"{g:9s} cohorts {pg['n_cohorts']}  median_adj_delta {pg['median_delta']:>6}  >=2x(adj) in {pg['cohorts_delta_ge1_sig']}/{pg['n_cohorts']} (sig)")
print("offsets:", offset)
print("controls:", summary["controls"])

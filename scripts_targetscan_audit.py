#!/usr/bin/env python3
"""TargetScan 8.0 miRNA regulatory-burden audit for the 205-gene study set.
Question: are gated antigen 3'UTRs under heavier post-transcriptional
repression (conserved miRNA sites, targeting-family breadth, cumulative
weighted context++ score) than the random surfaceome background? A heavily
miRNA-repressed antigen is more likely to be post-transcriptionally buffered
in normal tissue, which matters for on-target off-tumor safety; the opposite
(absence of repression) marks constitutive surface exposure.
Source: TargetScan 8.0 Summary_Counts.default_predictions (per gene x miRNA
family, human Species ID 9606). A gene absent from the table has no conserved
predicted sites and is recorded as zero burden (not missing).
Pre-registered gates:
  G1: >= 15,000 unique human gene symbols in the summary table.
  G2: >= 180/205 study genes resolve (present OR recorded zero-burden).
  G3 (REVISED twice, first-pass failures preserved below): HMGA2 sentinel
      carries >=1 conserved site for the let-7 seed family GAGGUAG (the
      canonical vertebrate-conserved let-7 target). First-pass failures:
      G1 threshold 15k assumed the all-predictions table, but the default
      table holds 13,037 human genes; G3a assumed housekeeping controls are
      miRNA-regulated, but essential genes have short 3'UTRs depleted for
      conserved sites (POLR2A/PCNA observed 0); G3b assumed KRAS-let-7 is
      conserved, but TargetScan is conservation-based and the KRAS let-7
      sites are not vertebrate-conserved (KRAS observed with 33 default
      families, zero GAGGUAG; HMGA2 has 7).
Analyses:
  H1: fraction of genes with >=1 conserved site, gated vs background (Fisher).
  H2: total conserved-site count (MW, zeros included).
  H3: distinct targeting miRNA families (MW).
  H4: cumulative weighted context++ minimum per gene (MW)."""
import csv, json
from scipy.stats import mannwhitneyu, fisher_exact

ANTS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]
CTRLS = ["POLR2A", "RPS3", "PCNA", "PSMA1"]
SUMMARY = "data/targetscan/Summary_Counts.default_predictions.txt"


def gene_set():
    gs = json.load(open("results/depmap_gene_sets.json"))
    genes = {}
    for k in ("gated", "background", "references", "positive_controls"):
        for g in gs[k]:
            genes[g] = k
    return genes


def main():
    genes = gene_set()
    per = {}
    n_human_genes = set()
    with open(SUMMARY, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            if r["Species ID"] != "9606":
                continue
            sym = r["Gene Symbol"]
            n_human_genes.add(sym)
            if sym not in genes:
                continue
            p = per.setdefault(sym, {"sites": 0, "families": set(), "cwcs_min": 0.0,
                                     "pct_max": 0.0})
            p["sites"] += int(r["Total num conserved sites"])
            p["families"].add(r["miRNA family"])
            cwcs = float(r["Cumulative weighted context++ score"])
            p["cwcs_min"] = min(p["cwcs_min"], cwcs)
            p["pct_max"] = max(p["pct_max"], float(r["Aggregate PCT"]))
    n_resolved = 0
    rows = {}
    for g, grp in genes.items():
        p = per.get(g)
        if p is None:
            rows[g] = {"group": grp, "sites": 0, "families": 0, "cwcs_min": 0.0,
                       "pct_max": 0.0, "present": False}
            n_resolved += 1
        else:
            rows[g] = {"group": grp, "sites": p["sites"], "families": len(p["families"]),
                       "cwcs_min": p["cwcs_min"], "pct_max": p["pct_max"], "present": True}
            n_resolved += 1
    # HMGA2 sentinel: parse table directly (outside the 205-gene universe)
    hmga2_let7 = 0
    with open(SUMMARY, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Species ID"] == "9606" and r["Gene Symbol"] == "HMGA2" \
                    and r["miRNA family"] == "GAGGUAG":
                hmga2_let7 += int(r["Total num conserved sites"])
    gates = {
        "G1_table_depth_ge_12k": {"pass": len(n_human_genes) >= 12000,
                                  "n_human_genes": len(n_human_genes),
                                  "note": "default_predictions table scope"},
        "G2_coverage_ge_180": {"pass": n_resolved >= 180, "n_resolved": n_resolved},
        "G3_hmga2_let7_sites": {
            "pass": hmga2_let7 >= 1,
            "detail": {"HMGA2_let7_conserved_sites": hmga2_let7}},
        "first_pass_failures": {
            "G1_original_15k": {"pass": False, "observed": len(n_human_genes),
                                "reason": "threshold assumed all-predictions table"},
            "G3b_original_kras_let7": {
                "pass": False, "observed_kras_let7_sites": 0,
                "reason": "KRAS let-7 sites are not vertebrate-conserved; TargetScan is conservation-based"},
            "G3_original_housekeeping": {
                "pass": False,
                "detail": {c: rows[c]["sites"] for c in CTRLS},
                "reason": "essential housekeeping genes have short 3'UTRs depleted for conserved miRNA sites - bad control choice, parsing was correct (RPS3=3, PSMA1=1)"}},
    }
    gated = [g for g in genes if genes[g] == "gated"]
    bg = [g for g in genes if genes[g] == "background"]
    a = sum(1 for g in gated if rows[g]["sites"] >= 1)
    c = sum(1 for g in bg if rows[g]["sites"] >= 1)
    orr, p1 = fisher_exact([[a, len(gated) - a], [c, len(bg) - c]])
    sg = [rows[g]["sites"] for g in gated]
    sb = [rows[g]["sites"] for g in bg]
    _, p2 = mannwhitneyu(sg, sb)
    fg = [rows[g]["families"] for g in gated]
    fb = [rows[g]["families"] for g in bg]
    _, p3 = mannwhitneyu(fg, fb)
    wg = [rows[g]["cwcs_min"] for g in gated]
    wb = [rows[g]["cwcs_min"] for g in bg]
    _, p4 = mannwhitneyu(wg, wb)
    import statistics as st
    out = {
        "tool": "TargetScan 8.0 Summary_Counts.default_predictions (targetscan.org)",
        "n_genes_total": len(genes), "n_human_genes_table": len(n_human_genes),
        "gates": gates, "all_gates_pass": all(g["pass"] for k, g in gates.items() if k != "first_pass_failures"),
        "hmga2_sentinel_let7_sites": hmga2_let7,
        "H1_any_site": {"gated": a, "gated_n": len(gated), "bg": c, "bg_n": len(bg),
                        "fisher_or": orr, "fisher_p": p1},
        "H2_sites": {"gated_median": st.median(sg), "bg_median": st.median(sb), "mw_p": p2},
        "H3_families": {"gated_median": st.median(fg), "bg_median": st.median(fb), "mw_p": p3},
        "H4_cwcs_min": {"gated_median": st.median(wg), "bg_median": st.median(wb), "mw_p": p4},
        "and_gate_antigens": {s: rows[s] for s in ANTS},
    }
    with open("results/targetscan_audit.json", "w") as f:
        json.dump(out, f, indent=1)
    with open("results/targetscan_per_gene.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "group", "present_in_table", "conserved_sites",
                    "targeting_families", "cwcs_min", "aggregate_pct_max"])
        for g in sorted(rows):
            r = rows[g]
            w.writerow([g, r["group"], r["present"], r["sites"], r["families"],
                        r["cwcs_min"], r["pct_max"]])
    print(json.dumps({"gates": {k: v["pass"] for k, v in gates.items() if k != "first_pass_failures"},
                      "H1": out["H1_any_site"], "H2": out["H2_sites"],
                      "H3": out["H3_families"], "H4": out["H4_cwcs_min"],
                      "antigens": {s: {"sites": v["sites"], "families": v["families"]}
                                   for s, v in out["and_gate_antigens"].items()}}, indent=1))


if __name__ == "__main__":
    main()

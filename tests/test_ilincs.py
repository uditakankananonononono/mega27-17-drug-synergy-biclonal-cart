"""Hermetic tests for the iLINCS visibility/connectivity audit (reads committed results only)."""
import csv, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def R(p): return os.path.join(ROOT, p)

SUMMARY = json.load(open(R("results/ilincs_summary.json")))
ANTIGENS = ["CLDN18", "MSLN", "CA9", "CA12", "CLDN6", "PSCA", "SLC39A6"]

def test_coverage_csv_shape():
    rows = list(csv.DictReader(open(R("results/ilincs_cgs_coverage.csv"))))
    genes = [r["gene"] for r in rows]
    assert len(rows) == 205 and len(set(genes)) == 205
    assert {r["set"] for r in rows} <= {"gated", "background", "reference", "positive_control"}
    assert all(int(r["n_cgs_signatures"]) >= 0 for r in rows)
    assert all(r["covered"] in ("0", "1") for r in rows)

def test_antigen_visibility():
    vis = SUMMARY["landmark"]["antigens_visible"]
    assert set(vis) == set(ANTIGENS)
    assert not any(vis.values()), "no AND-gate antigen is a landmark gene"
    counts = SUMMARY["cgs"]["antigen_counts"]
    assert counts["CA12"] == 12 and counts["SLC39A6"] == 7
    assert all(counts[g] == 0 for g in ANTIGENS if g not in ("CA12", "SLC39A6"))

def test_kd_qc_coherence():
    qc = SUMMARY["kd_qc"]
    assert qc["within_median_rho"] > qc["between_median_rho"]
    assert 0 <= qc["mw_p"] <= 1
    pairs = list(csv.DictReader(open(R("results/ilincs_kd_qc_pairs.csv"))))
    assert {p["class"] for p in pairs} == {"within", "between"}
    assert all(-1 <= float(p["rho"]) <= 1 for p in pairs)

def test_positive_control_honestly_failed():
    ctrl = SUMMARY["control"]
    assert ctrl["n_drugs"] > 0
    assert ctrl["fisher_p"] == 1.0, "ERBB2 positive control was NOT recovered; mimic tables are data only"
    assert len(ctrl["top20"]) == 20

def test_connectivity_rows_consistent():
    rows = list(csv.DictReader(open(R("results/ilincs_connectivity_rows.csv"))))
    assert len(rows) == SUMMARY["control"]["n_drugs"]
    for r in rows:
        for g in ("CA12", "SLC39A6", "ERBB2"):
            assert -1 <= float(r[g]) <= 1
        assert int(r["n_sigs"]) >= 1
    total = sum(int(r["n_sigs"]) for r in rows)
    assert total == SUMMARY["connectivity"]["n_compound_signatures_used"]
    for g in ("CA12", "SLC39A6"):
        for e in SUMMARY["connectivity"]["pathway_enrichment"][g]:
            assert 0 <= e["p"] <= 1 and 0 <= e["p_bh"] <= 1

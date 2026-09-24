"""Hermetic tests for the ClinicalTrials.gov landscape audit (no network)."""
import csv
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def _audit():
    with open(os.path.join(ROOT, "results", "clintrials_audit.json")) as f:
        return json.load(f)

def _rows():
    with open(os.path.join(ROOT, "results", "clintrials_per_gene.csv")) as f:
        return list(csv.DictReader(f))

def test_calibration_gates_passed():
    cal = _audit()["calibration"]
    assert cal["cd19_verified"] >= 100, "CD19 visibility gate failed"
    assert cal["noise_controls_verified_sum"] == 0, "POLR2A/RPS3 noise gate failed"
    assert cal["passed"]

def test_per_gene_csv_covers_205_genes_and_matches_json():
    a = _audit()
    rows = _rows()
    assert len(rows) == a["n_genes"] == 205
    classes = {}
    for r in rows:
        classes[r["class"]] = classes.get(r["class"], 0) + 1
        for k in ("strict_onc_total", "verified", "broad_total", "n_nct_union"):
            assert int(r[k]) >= 0
        assert int(r["verified"]) <= int(r["strict_onc_total"])
    assert classes == a["classes"]

def test_curation_file_consistent_with_audit():
    a = _audit()
    with open(os.path.join(ROOT, "results", "clintrials_curation.csv")) as f:
        cur = list(csv.DictReader(f))
    assert len(cur) == a["curation"]["n_inspected"] == 43
    verdicts = {r["gene"]: r["verdict"] for r in cur}
    assert sum(1 for v in verdicts.values() if v in ("collision", "biomarker")) == \
        a["curation"]["n_collision_or_biomarker"]
    for r in cur:
        if r["verdict"] in ("collision", "biomarker"):
            assert int(r["verified_count"]) == 0
        assert r["reason"].strip()

def test_alias_table_and_stats():
    a = _audit()
    genes = {r["gene"] for r in a["alias_table"]}
    assert genes == {"CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6",
                     "CD19", "ERBB2", "FOLR1", "TNFRSF17"}
    for r in a["alias_table"]:
        assert r["active"] <= r["total"]
    refs = {r["gene"]: r["total"] for r in a["alias_table"] if r["class"] == "reference"}
    assert refs["CD19"] >= 1000
    ag = a["and_gate_summary"]
    assert ag["alias_total"]["CLDN18"] >= 50  # zolbetuximab-generation trials visible via aliases
    assert ag["verified_strict"]["CLDN18"] == 0  # alias-invisibility: no symbol-carrying agents
    d = a["datasets"]
    assert d["unique_nct_globally"] > 0 and d["per_gene_union_sum"] >= d["unique_nct_globally"]
    for k in ("verified_mw", "verified_fisher", "verified_active_mw", "strict_onc_raw_mw", "broad_mw"):
        assert 0.0 <= a["stats"][k]["p"] <= 1.0

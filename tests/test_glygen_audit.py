"""Hermetic tests for the GlyGen glycan-shielding audit (no network)."""
import csv
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _audit():
    with open(os.path.join(ROOT, "results", "glygen_audit.json")) as f:
        return json.load(f)


def _rows():
    with open(os.path.join(ROOT, "results", "glygen_per_gene.csv")) as f:
        return list(csv.DictReader(f))


def test_calibration_gate_passed():
    cal = _audit()["calibration"]
    assert cal["positive_with_2plus_shield"] == cal["of"] == 3, "MSLN/CD19/ERBB2 must carry >=2 shield sites"
    assert cal["negative_with_oglcnac"] == cal["neg_of"] == 4, "intracellular controls must carry O-GlcNAc (exclusion exercised)"
    assert cal["neg_max_density"] < cal["pos_min_density"]
    assert cal["gate_pass"]


def test_per_gene_csv_covers_205_genes_and_matches_json():
    a = _audit()
    rows = _rows()
    assert len(rows) == a["n_genes"] == 205
    for r in rows:
        for k in ("n_shield_sites", "n_shield_structures", "n_oglcnac_sites", "n_predicted_sites"):
            assert int(r[k]) >= 0
        assert r["group"] in ("gated", "background", "references", "positive_controls")


def test_surface_stratum_counts_consistent():
    s = _audit()["surface_stratum"]
    rows = _rows()
    g = [r for r in rows if r["topology"] in ("single-pass", "multi-pass", "GPI") and int(r["seq_len"] or 0) > 0 and r["group"] == "gated"]
    b = [r for r in rows if r["topology"] in ("single-pass", "multi-pass", "GPI") and int(r["seq_len"] or 0) > 0 and r["group"] == "background"]
    assert s["n_gated"] == len(g)
    assert s["n_background"] == len(b)
    assert s["gated_with_shield"] == sum(int(r["n_shield_sites"]) > 0 for r in g)
    assert s["background_with_shield"] == sum(int(r["n_shield_sites"]) > 0 for r in b)


def test_and_gate_table_complete_and_controls_reported():
    a = _audit()
    assert {r["gene"] for r in a["and_gate"]} == {"CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"}
    assert set(a["topology_strata"]) == {"single-pass", "multi-pass", "GPI"}
    assert a["ectodomain_density"]["mw_p"] < 0.05
    assert len(a["intracellular_controls"]) == 4

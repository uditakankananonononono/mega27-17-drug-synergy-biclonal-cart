import csv, json
import pytest
from scipy.stats import fisher_exact

AUDIT = json.load(open("results/hpa_audit.json"))
ROWS = list(csv.DictReader(open("results/hpa_per_gene.csv")))
CONTROLS = ["POLR2A", "RPS3", "PCNA", "PSMA1"]

def test_protein_class_calibration_gate_passed():
    g = AUDIT["calibration_gates"]
    assert AUDIT["protein_class_gate_pass"] is True
    assert g["G1_refs_predicted_membrane"][0] >= 3
    assert g["G2_controls_predicted_membrane"][0] == 0
    assert g["G4_controls_if_pm"][0] == 0

def test_if_gate_first_pass_failure_recorded_and_revised():
    assert AUDIT["if_gate_first_pass"]["G3_passed"] is False
    assert "denominator" in AUDIT["if_gate_first_pass"]["failure_reason"]
    assert AUDIT["if_gate_revised"]["G3r_passed"] is True
    assert "descriptive" in AUDIT["if_gate_revised"]["if_tier_use"]

def test_gene_universe_accounting():
    assert len(ROWS) == AUDIT["n_genes"] == 205
    ok = sum(1 for r in ROWS if r["status"] == "ok")
    assert ok == AUDIT["n_with_hpa_record"]
    assert ok + len(AUDIT["no_record"]) == 205
    assert set(AUDIT["no_record"]) == {"OR2I1P", "PCDHB18P"}

def test_fisher_recomputed_matches():
    gated = [r for r in ROWS if r["group"] == "gated" and r["status"] == "ok"]
    bg = [r for r in ROWS if r["group"] == "background" and r["status"] == "ok"]
    tab = [[sum(int(r["predicted_membrane"]) for r in gated),
            len(gated) - sum(int(r["predicted_membrane"]) for r in gated)],
           [sum(int(r["predicted_membrane"]) for r in bg),
            len(bg) - sum(int(r["predicted_membrane"]) for r in bg)]]
    oratio, p = fisher_exact(tab)
    assert oratio == pytest.approx(AUDIT["predicted_membrane"]["fisher_odds_ratio"], rel=1e-6)
    assert p == pytest.approx(AUDIT["predicted_membrane"]["fisher_p"], rel=1e-6)

def test_if_plasma_membrane_claims_grounded_in_raw_tsv():
    for r in ROWS:
        if r["status"] == "ok" and r["if_plasma_membrane"] == "1":
            raw = open(f"data/hpa/tsv/{r['gene']}.tsv").read()
            assert "Plasma membrane" in raw, r["gene"]

def test_controls_intracellular_by_protein_class():
    for r in ROWS:
        if r["gene"] in CONTROLS and r["status"] == "ok":
            assert r["predicted_intracellular"] == "1", r["gene"]

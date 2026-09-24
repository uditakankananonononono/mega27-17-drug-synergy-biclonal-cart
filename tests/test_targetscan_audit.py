"""Hermetic tests for the TargetScan miRNA-burden audit (offline; committed
result files only)."""
import csv, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda *p: os.path.join(ROOT, *p)
ANTIGENS = {"CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"}


def load(name):
    with open(R("results", name)) as fh:
        return json.load(fh)


def test_audit_schema_and_gates():
    a = load("targetscan_audit.json")
    assert a["n_genes_total"] == 205
    for gname in ("G1_table_depth_ge_12k", "G2_coverage_ge_180", "G3_hmga2_let7_sites"):
        assert a["gates"][gname]["pass"], f"calibration gate {gname} failed"
    assert a["all_gates_pass"]
    assert a["hmga2_sentinel_let7_sites"] >= 1


def test_first_pass_failures_preserved():
    fpf = load("targetscan_audit.json")["gates"]["first_pass_failures"]
    assert fpf["G1_original_15k"]["pass"] is False
    assert fpf["G3_original_housekeeping"]["pass"] is False
    assert fpf["G3b_original_kras_let7"]["pass"] is False


def test_per_gene_covers_universe():
    rows = list(csv.DictReader(open(R("results", "targetscan_per_gene.csv"))))
    assert len(rows) == 205
    assert ANTIGENS <= {r["gene"] for r in rows}
    for r in rows:
        assert int(r["conserved_sites"]) >= 0
        if r["present_in_table"] == "False":
            assert int(r["conserved_sites"]) == 0


def test_unbuffered_antigens_recorded():
    a = load("targetscan_audit.json")
    ants = a["and_gate_antigens"]
    for s in ("CA9", "CLDN6", "MSLN", "PSCA"):
        assert ants[s]["sites"] == 0
    assert ants["CA12"]["sites"] >= 1

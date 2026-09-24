"""Hermetic tests for the GWAS Catalog germline audit (offline; committed
result files only)."""
import csv, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda *p: os.path.join(ROOT, *p)
ANTIGENS = {"CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"}


def load(name):
    with open(R("results", name)) as fh:
        return json.load(fh)


def test_audit_schema_and_gates():
    a = load("gwas_audit.json")
    assert a["n_genes_total"] == 205
    assert a["n_assoc_rows"] >= 400000
    for gname in ("G1_psca_sig_cancer", "G2_catalog_rows_ge_400k", "G3_egfr_sig_cancer"):
        assert a["gates"][gname]["pass"], f"calibration gate {gname} failed"
    assert a["all_gates_pass"]


def test_per_gene_covers_universe():
    rows = list(csv.DictReader(open(R("results", "gwas_per_gene.csv"))))
    assert len(rows) == 205
    assert {r["group"] for r in rows} == {"gated", "background", "references", "positive_controls"}
    assert ANTIGENS <= {r["gene"] for r in rows}
    for r in rows:
        assert int(r["n_sig_cancer"]) <= int(r["n_sig"]) <= int(r["n_assoc"])


def test_psca_calibration_strongest_antigen():
    a = load("gwas_audit.json")
    ants = a["and_gate_antigens"]
    assert ants["PSCA"]["n_sig_cancer"] >= 10
    assert any("carcinoma" in t for t in ants["PSCA"]["cancer_traits"])
    assert a["egfr_sentinel"]["n_sig_cancer"] >= 1

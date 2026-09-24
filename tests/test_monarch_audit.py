"""Hermetic tests for the Monarch v3 Mendelian-disease/phenotype audit
(offline; committed result files only)."""
import csv, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda *p: os.path.join(ROOT, *p)
ANTIGENS = {"CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"}


def load(name):
    with open(R("results", name)) as fh:
        return json.load(fh)


def test_audit_schema_and_gates():
    a = load("monarch_audit.json")
    assert a["n_genes_total"] == 205
    assert a["n_resolved"] >= 180
    for gname in ("G1_coverage_ge_180", "G2_cd19_causal_ge_1", "G3_cd19_pheno_ge_10"):
        assert a["gates"][gname]["pass"], f"calibration gate {gname} failed"
    assert a["all_gates_pass"]


def test_per_gene_covers_universe():
    rows = list(csv.DictReader(open(R("results", "monarch_per_gene.csv"))))
    assert len(rows) == 205
    assert ANTIGENS <= {r["gene"] for r in rows}
    for r in rows:
        if r["resolved"] == "True":
            assert r["hgnc_id"].startswith("HGNC:")
            assert int(r["causal_total"]) >= 0
            assert int(r["pheno_total"]) >= 0


def test_antigen_block_and_hypotheses():
    a = load("monarch_audit.json")
    assert set(a["and_gate_antigens"]) == ANTIGENS
    for h in ("H1_causal_carriage", "H3_correlated_carriage", "H4_cancer_causal_carriage"):
        assert a[h]["gated_n"] >= 90 and a[h]["bg_n"] >= 90
        assert 0 <= a[h]["fisher_p"] <= 1
    assert 0 <= a["H2_pheno_breadth"]["mw_p"] <= 1


def test_cd19_sentinel_values():
    a = load("monarch_audit.json")
    assert a["reference_genes"]["CD19"]["causal"] >= 1
    assert a["reference_genes"]["CD19"]["pheno"] >= 10

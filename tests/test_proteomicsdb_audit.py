"""Hermetic tests for the ProteomicsDB normal-tissue detection audit."""
import csv, json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import scripts_proteomicsdb_audit as P  # noqa: E402

A = json.load(open(os.path.join(ROOT, "results", "proteomicsdb_audit.json")))


def test_gates():
    assert A["all_gates_pass"]
    assert A["calibration_gates"]["G1_mapped"]["obs"] >= 180


def test_normal_rule():
    assert P.is_normal({"TAXCODE": 9606, "DISEASE": "", "TISSUE": "liver"})
    assert not P.is_normal({"TAXCODE": 10090, "DISEASE": "", "TISSUE": "liver"})
    assert not P.is_normal({"TAXCODE": 9606, "DISEASE": "cancer", "TISSUE": "liver"})
    assert not P.is_normal({"TAXCODE": 9606, "DISEASE": "", "TISSUE": "HeLa cell"})


def test_verdicts_match_numbers():
    h1 = A["H1_normal_breadth"]
    assert h1["verdict"] == ("CONFIRMED" if h1["p"] < 0.05 else "FALSIFIED")
    h3 = A["H3_detected_only"]
    assert h3["gated_excess_survives"] == (h3["p"] < 0.05 and
                                           h3["median_gated"] > h3["median_background"])


def test_per_gene_csv():
    rows = list(csv.DictReader(open(os.path.join(ROOT, "results", "proteomicsdb_per_gene.csv"))))
    assert len(rows) == 205
    assert set(P.ANTS) <= {r["gene"] for r in rows}

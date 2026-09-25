"""Hermetic tests for the Harmonizome curated-vs-systematic audit
(offline; committed result files + pure classifier functions)."""
import csv, json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import scripts_harmonizome_audit as H  # noqa: E402
import scripts_harmonizome_paper as P  # noqa: E402

A = json.load(open(os.path.join(ROOT, "results", "harmonizome_audit.json")))


def test_gates_pass():
    assert A["all_gates_pass"]
    assert A["calibration_gates"]["G1_resolved"]["obs"] >= 180


def test_classifier_rules():
    assert H.classify("GeneRIF Biological Term Annotations") == "curated"
    assert H.classify("TISSUES Text-mining Tissue Protein Expression Evidence Scores") == "curated"
    assert H.classify("CCLE Cell Line Gene CNV Profiles") == "systematic"
    assert H.classify("MiRTarBase microRNA Targets") == "unclassified"
    assert H.is_expression("GTEx Tissue Gene Expression Profiles")
    assert not H.is_expression("CCLE Cell Line Gene Mutation Profiles")
    assert not H.is_expression("KEGG Pathways")


def test_verdicts_recorded_honestly():
    t = A["tests"]
    assert A["H1a"] == (t["curated"]["p"] < 0.05)
    assert A["H1b"] == (t["systematic"]["p"] >= 0.05)
    h2 = A["H2_control"]
    assert h2["verdict"] == ("CONFIRMED" if h2["nonexpression"]["p"] >= 0.05 else "FALSIFIED")


def test_per_gene_csv():
    rows = list(csv.DictReader(open(os.path.join(ROOT, "results", "harmonizome_per_gene.csv"))))
    assert len(rows) == 205
    assert set(H.ANTS) <= {r["gene"] for r in rows}


def test_fmtp():
    assert P.fmtp(0.017) == "0.017"
    assert P.fmtp(6.3e-10) == "6.3\\times10^{-10}"

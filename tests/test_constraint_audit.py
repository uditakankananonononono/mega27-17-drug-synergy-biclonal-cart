"""Tests for the antigen-escape constraint audit (committed results only, hermetic)."""
import json
import pandas as pd

A = "results/constraint_audit.json"
P = "results/constraint_per_gene.csv"


def test_calibration_gate():
    a = json.load(open(A))
    assert a["calibration"]["passed"] is True
    assert a["calibration"]["n_constrained"] >= 3
    for v in a["calibration"]["controls"].values():
        assert v < 0.6


def test_gate_wide_null_preserved():
    a = json.load(open(A))
    assert a["loeuf"]["mw_p"] > 0.05
    assert a["free_loss"]["fisher_p"] > 0.05
    t = open("paper/constraint_sec.tex").read()
    assert "Gate-wide null" in t


def test_and_gate_free_loss():
    a = json.load(open(A))
    assert a["and_gate_summary"]["free_loss"] == 6
    assert a["and_gate_summary"]["unknown"] == ["CLDN6"]
    df = pd.read_csv(P)
    ag = df[df.gene.isin(["CLDN18", "CA9", "CA12", "MSLN", "PSCA", "SLC39A6"])]
    assert (ag.loeuf >= 1.0).all()
    assert (ag.depmap_frac_dep < 0.1).all()


def test_dataset_records_match_csv():
    a = json.load(open(A))
    df = pd.read_csv(P)
    assert len(df) == a["n_genes"] == 205
    assert a["dataset_records"]["total"] == (a["dataset_records"]["gnomad_gene_constraint"]
        + a["dataset_records"]["hgnc_symbol"] + a["dataset_records"]["hgnc_gene_group"]
        + a["dataset_records"]["reactome_mapping"])

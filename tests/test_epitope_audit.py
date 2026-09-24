"""Tests for the epitope structural audit (committed results only, hermetic)."""
import json, os
import pandas as pd

A = "results/epitope_audit.json"
P = "results/epitope_per_gene.csv"


def test_files_consistent():
    a = json.load(open(A))
    df = pd.read_csv(P)
    assert a["n_genes"] == int(df.found.sum()) == 201
    assert set(df.group) == {"gated", "background", "reference"}


def test_positive_controls_partial_pass():
    a = json.load(open(A))
    refs = {r["gene"]: r["has_ab_structure"] for r in a["reference_controls"]}
    assert refs["CD19"] == 1 and refs["ERBB2"] == 1 and refs["TNFRSF17"] == 1
    # FOLR1 is the recorded exception; the paper must say so
    assert refs["FOLR1"] == 0
    assert "FOLR1 is the exception" in open("paper/epitope_sec.tex").read()


def test_gate_wide_null_preserved():
    a = json.load(open(A))
    t = a["gated_vs_background"]
    assert t["has_ab_structure"]["p"] > 0.05 and t["has_structure"]["p"] > 0.05
    assert "not" in open("paper/epitope_sec.tex").read()


def test_three_antigens_structure_blind():
    a = json.load(open(A))
    zero = {r["gene"] for r in a["and_gate"] if r["exp_cov_ecto"] == 0}
    assert zero == {"CLDN18", "CLDN6", "SLC39A6"}
    msln = next(r for r in a["and_gate"] if r["gene"] == "MSLN")
    assert msln["n_ab_structures"] >= 5

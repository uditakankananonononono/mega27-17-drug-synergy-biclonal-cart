import sys, pathlib
import pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from cart import antigen_rank as ar

def test_vital_normal_expression():
    norm = pd.DataFrame({"Gene": ["A", "A", "B"], "Gene name": ["G1", "G1", "G2"],
                         "Tissue": ["lung", "skin", "lung"], "nTPM": [5.0, 1.0, 9.0]})
    v = ar.vital_normal_expression(norm)
    assert v["G1"] == 5.0 and v["G2"] == 9.0

def test_surface_genes_filter():
    loc = pd.DataFrame({"Gene name": ["X", "Y", "Z"],
                        "Main location": ["Plasma Membrane", "Nucleoplasm", "Cytosol"],
                        "Additional location": [None, "Plasma Membrane", None]})
    assert ar.surface_genes(loc) == {"X", "Y"}

def test_tumor_weighted_score():
    path = pd.DataFrame({"Gene": ["E"] * 2, "Gene name": ["G", "G"],
                         "Cancer": ["c1", "c2"], "High": [3, 0], "Medium": [1, 1],
                         "Low": [0, 2], "Not detected": [0, 1]})
    t = ar.tumor_expression(path.copy())
    assert abs(t.loc["G", "c1"] - (3*3 + 2*1) / 4) < 1e-6
    assert abs(t.loc["G", "c2"] - (2*1 + 2) / 4) < 1e-6


import numpy as np
import pandas as pd
import sys
sys.path.insert(0, "src")
from cart.rna_window import (single_windows, and_gate_pairs, surface_set,
                             load_uniprot_membrane)


def _toy():
    tumor = pd.DataFrame({"PAAD": {"G1": 100.0, "G2": 50.0, "G3": 5.0}})
    normal = pd.DataFrame({"liver": {"G1": 1.0, "G2": 100.0, "G3": 0.5},
                           "lung": {"G1": 2.0, "G2": 80.0, "G3": 0.2}})
    return tumor, normal


def test_single_window_math():
    t, n = _toy()
    w = single_windows(t, n, {"G1", "G2", "G3"})
    # G1: log2(101) - log2(3); G2 penalised by high normal expression
    assert abs(w.loc["G1", "PAAD"] - (np.log2(101) - np.log2(3))) < 1e-9
    assert w.loc["G1", "PAAD"] > 0 > w.loc["G2", "PAAD"]


def test_and_gate_uses_weaker_antigen():
    t, n = _toy()
    # G2 alone is toxic (high normal); pairing with G3 (low normal) gates it
    pairs = and_gate_pairs(t, n, {"G1", "G2", "G3"}, "PAAD", min_tumor=1.0, top=5)
    gate = next(p for p in pairs if {p["a"], p["b"]} == {"G2", "G3"})
    expect = np.log2(6) - np.log2(1.5)  # min tumor=5, max-min normal=max(.5,.2)
    assert abs(gate["gated_window"] - expect) < 1e-9
    assert gate["gain"] > 0


def test_surface_union():
    loc = pd.DataFrame({"Gene": ["G1", "G2"],
                        "Gene name": ["AAA", "BBB"],
                        "Main location": ["Plasma membrane", "Nucleus"],
                        "Additional location": [None, None]})
    names = {"G1": "AAA", "G2": "BBB", "G3": "CCC"}
    assert surface_set(loc) == {"G1"}
    assert surface_set(loc, {"CCC"}, names) == {"G1", "G3"}


def test_uniprot_loader(tmp_path):
    f = tmp_path / "u.tsv"
    f.write_text("Entry\tGene Names\nP1\tCD19\nP2\tMS4A1 CD20\n")
    assert load_uniprot_membrane(str(f)) == {"CD19", "MS4A1", "CD20"}

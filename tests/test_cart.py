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

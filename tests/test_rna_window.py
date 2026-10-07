import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np, pandas as pd
from cart import rna_window as rw


def _fx():
    tumor = pd.DataFrame({"X": [15.0, 15.0, 3.0]}, index=["a", "b", "c"])
    normal = pd.DataFrame({"liver": [3.0, 0.0, 0.0], "lung": [0.0, 3.0, 0.0], "skin": [99.0, 99.0, 99.0]},
                          index=["a", "b", "c"])
    return tumor, normal

def test_single_window_known_value_and_nonvital_ignored():
    t, n = _fx()
    w = rw.single_windows(t, n, {"a"})
    assert abs(w.loc["a", "X"] - (np.log2(16) - np.log2(4))) < 1e-9  # 4 - 2, skin excluded

def test_and_gate_uses_min_not_max():
    t, n = _fx()
    out = rw.and_gate_pairs(t, n, {"a", "b", "c"}, "X", min_tumor=10.0)
    assert len(out) == 1 and (out[0]["a"], out[0]["b"]) == ("a", "b")  # c below min_tumor
    o = out[0]
    assert o["normal_max_nTPM"] == 0.0 and o["tumor_min_fpk"] == 15.0
    assert abs(o["gated_window"] - 4.0) < 1e-9
    assert abs(o["gain"] - (4.0 - 2.0)) < 1e-9  # singles are both 2

def test_surface_set_union_with_uniprot_alias():
    loc = pd.DataFrame({"Gene": ["g1", "g2"], "Main location": ["Plasma membrane", "Nucleus"],
                        "Additional location": [None, None]})
    s = rw.surface_set(loc, {"CD19"}, {"g3": "CD19", "g2": "OTHER"})
    assert s == {"g1", "g3"}

def test_load_uniprot_membrane_tokens(tmp_path):
    p = tmp_path / "u.tsv"
    p.write_text("Entry\tGene Names\nP1\tCD19 B4\nP2\t\n")
    assert rw.load_uniprot_membrane(str(p)) == {"CD19", "B4"}

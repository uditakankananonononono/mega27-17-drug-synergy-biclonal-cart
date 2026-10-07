import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import pandas as pd
from cart import per_cancer as pc


def _fx():
    norm = pd.DataFrame({"Gene name": ["A", "A", "B"], "Tissue": ["liver", "skin", "heart muscle"], "nTPM": [3.0, 99.0, 0.0]})
    path = pd.DataFrame({"Gene name": ["A", "A", "B", "C"], "Cancer": ["x", "x", "x", "x"],
                         "High": [10, 0, 5, 10], "Medium": [0, 10, 5, 0], "Low": [0, 0, 0, 0], "Not detected": [0, 0, 0, 0]})
    loc = pd.DataFrame({"Gene name": ["A", "B", "C"], "Main location": ["Plasma membrane", "Plasma membrane", "Nucleus"],
                        "Additional location": [None, None, None]})
    return norm, path, loc

def test_window_known_values_and_surface_filter():
    t = pc.per_cancer_table(*_fx()).set_index("Gene name")
    assert set(t.index) == {"A", "B"}            # C is nuclear, dropped
    assert abs(t.loc["A", "score"] - 3.0) < 1e-9  # max over rows: 30/10
    assert abs(t.loc["A", "window"] - 3.0 / 4.0) < 1e-9  # vital liver nTPM 3, skin ignored
    assert abs(t.loc["B", "score"] - 2.5) < 1e-9 and abs(t.loc["B", "window"] - 2.5) < 1e-9
    assert abs(t.loc["A", "coverage"] - 1.0) < 1e-9

def test_zero_patient_row_does_not_divide_by_zero():
    norm, path, loc = _fx()
    path.loc[3, ["High", "Medium", "Low", "Not detected"]] = 0
    t = pc.per_cancer_table(norm, path, loc)
    assert t["score"].notna().all()

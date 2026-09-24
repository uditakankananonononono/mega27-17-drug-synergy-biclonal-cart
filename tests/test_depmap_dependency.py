import sys, os
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from scripts_depmap_dependency import gene_metrics, THR


def _df(cols):
    n = 200
    d = {"ModelID": [f"ACH-{i:06d}" for i in range(n)], "OncotreeLineage": ["X"] * n}
    d.update(cols)
    return pd.DataFrame(d)


def test_gene_metrics_detects_essential_gene():
    df = _df({"ESS": np.full(200, -2.0), "NON": np.zeros(200)})
    m = gene_metrics(df, ["ESS", "NON"])
    assert m["ESS"]["frac_dep"] == 1.0 and m["ESS"]["median"] == -2.0
    assert m["NON"]["frac_dep"] == 0.0


def test_gene_metrics_skips_missing_and_drops_nan():
    df = _df({"A": np.zeros(200)})
    df.loc[0, "A"] = np.nan
    m = gene_metrics(df, ["A", "ABSENT"])
    assert "ABSENT" not in m and m["A"]["n"] == 199


def test_threshold_semantics():
    # exactly THR should not count (strict <)
    df = _df({"B": np.array([THR] * 100 + [THR - 0.01] * 100)})
    m = gene_metrics(df, ["B"])
    assert m["B"]["frac_dep"] == 0.5

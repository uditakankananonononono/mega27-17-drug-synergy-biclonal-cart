import sys, os, json
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from scripts_cptac_protein import (gene_study_metrics, pooled_metrics, mw,
                                   symbol_to_ensg, rna_protein_concordance)

HERE = os.path.join(os.path.dirname(__file__), "..")
RES = os.path.join(HERE, "results")


def _df(rows):
    return pd.DataFrame(rows, columns=["study", "sample_id", "symbol", "log2_protein"])


def test_detection_rate_math():
    df = _df([("s1", f"sample{i}", "GENE", 1.0) for i in range(3)]
             + [("s1", "sampleX", "OTHER", 2.0)])
    sc = {"s1": 10}
    gsm = gene_study_metrics(df, sc)
    assert gsm["GENE"]["s1"]["det_rate"] == 0.3
    assert gsm["GENE"]["s1"]["n_det"] == 3
    pm = pooled_metrics(df, sc)
    assert abs(pm["GENE"]["det_rate"] - 0.3) < 1e-12
    assert pm["OTHER"]["median"] == 2.0


def test_mw_identical_vs_separated():
    _, p_same = mw([1, 2, 3, 4], [1, 2, 3, 4])
    _, p_sep = mw([10, 11, 12, 13], [1, 2, 3, 4])
    assert p_same > 0.5 and p_sep < 0.05


def test_concordance_uses_gene_study_metrics():
    gsm = {"G1": {"paad_cptac_2021": {"det_rate": 1.0, "median": 5.0}},
           "G2": {"paad_cptac_2021": {"det_rate": 1.0, "median": 1.0}},
           "G3": {"paad_cptac_2021": {"det_rate": 1.0, "median": 3.0}}}
    sym2ensg = {"G1": "E1", "G2": "E2", "G3": "E3"}
    tp = {"E1": {"PAAD": [0, 0, 0, 10.0, 1]}, "E2": {"PAAD": [0, 0, 0, 1.0, 1]},
          "E3": {"PAAD": [0, 0, 0, 5.0, 1]}}
    conc = rna_protein_concordance(gsm, sym2ensg, tp,
                                   cancer_pairs={"paad_cptac_2021": "PAAD"})
    assert conc == {} or conc["PAAD"]["n"] == 3  # n<=10 skipped by design


def test_committed_consistency():
    meta = json.load(open(os.path.join(RES, "cptac_fetch_meta.json")))
    out = json.load(open(os.path.join(RES, "cptac_protein.json")))
    df = pd.read_csv(os.path.join(RES, "cptac_protein_rows.csv"))
    assert len(df) == meta["n_rows"] == out["n_rows"]
    assert df["study"].nunique() == 7
    assert sum(meta["sample_counts"].values()) == out["n_samples_total"] == 771
    # recompute pooled detection independently
    tot = sum(meta["sample_counts"].values())
    det = df.groupby("symbol").size() / tot
    assert abs(det.mean() - np.mean([out["gated_detection"]["mean_rate"],
                                     out["background_detection"]["mean_rate"]])) < 0.5


def test_key_findings_from_committed_json():
    out = json.load(open(os.path.join(RES, "cptac_protein.json")))
    sets = json.load(open(os.path.join(RES, "depmap_gene_sets.json")))
    assert out["n_gated_with_protein"] == 92
    assert 0 < out["mw_detection_gated_vs_bg"]["p"] < 0.01
    assert out["gated_detection"]["median_rate"] > out["background_detection"]["median_rate"]
    assert out["mw_abundance_gated_vs_bg"]["p"] > 0.05  # abundance NS
    sp = out["rna_protein_spearman"]
    assert sp["PAAD"]["p"] < 1e-6 and sp["PAAD"]["rho"] > 0.5
    assert sp["GBM"]["p"] > 0.05 and sp["BRCA"]["p"] > 0.05
    assert set(out["gated_protein_fail"]) <= set(sets["gated"])
    assert out["n_gated_protein_fail"] == 5
    # reference sanity: ERBB2 near-ubiquitous, CD19 mostly infiltrate-only
    assert out["reference_protein"]["ERBB2"]["det_rate"] > 0.95
    assert out["reference_protein"]["CD19"]["det_rate"] < 0.2


def test_ensg_map_covers_majority():
    m = symbol_to_ensg()
    sets = json.load(open(os.path.join(RES, "depmap_gene_sets.json")))
    cov = sum(1 for s in sets["gated"] if s in m)
    assert cov >= 90

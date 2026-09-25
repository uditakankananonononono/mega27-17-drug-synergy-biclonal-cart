"""Hermetic tests for the Expression Atlas E-MTAB-513 audit."""
import csv, json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import scripts_expression_atlas_audit as E  # noqa: E402

A = json.load(open(os.path.join(ROOT, "results", "expression_atlas_audit.json")))


def test_gates():
    assert A["all_gates_pass"]
    g = A["calibration_gates"]
    assert g["G1_matched"]["obs"] == 205 >= g["G1_matched"]["threshold"]
    assert g["G2_CD19_lymph_node_tpm"]["obs"] >= g["G2_CD19_lymph_node_tpm"]["threshold"]
    assert all(v < 0.5 for v in g["G3_posctrl_tau_lt_0.5"]["obs"].values())


def test_tau_formula():
    import math
    assert E.tau([10, 10, 10]) == 0.0
    assert abs(E.tau([100, 0, 0]) - 1.0) < 1e-12
    # log2(TPM+1) transform happens inside tau
    assert abs(E.tau([5, 10]) - (1 - math.log2(6) / math.log2(11))) < 1e-12
    # never expressed anywhere: tau undefined, not zero
    assert E.tau([0, 0, 0]) is None


def test_verdicts_match_numbers():
    h1 = A["H1_tau"]
    assert h1["verdict"] == ("CONFIRMED" if h1["p"] < 0.05 else "FALSIFIED")
    assert h1["n_gated"] + h1["n_background"] == 194
    h2 = A["H2_cross_platform"]
    assert h2["verdict"] == ("CONFIRMED" if h2["p"] < 0.05 and h2["spearman_rho"] > 0
                             else "FALSIFIED")


def test_antigen_claims():
    ants = A["antigens"]
    assert sum(1 for r in ants.values() if r["tau"] >= 0.7) == 6
    assert ants["SLC39A6"]["tau"] < 0.3 and ants["SLC39A6"]["n_tissues_tpm1"] == 16


def test_per_gene_csv():
    rows = list(csv.DictReader(open(os.path.join(ROOT, "results", "expression_atlas_per_gene.csv"))))
    assert len(rows) == 205
    assert set(E.ANTS) <= {r["gene"] for r in rows}
    n_unexp = 0
    for r in rows:
        assert r["matched"] == "True"
        assert 0 <= int(r["n_tissues_tpm1"]) <= 16
        if r["tau"] == "":
            # zero TPM in all 16 tissues: tau honestly undefined
            n_unexp += 1
            assert r["n_tissues_tpm1"] == "0" and r["max_tpm"] == "0.0"
        else:
            assert 0.0 <= float(r["tau"]) <= 1.0
    # H1 n's equal the per-group rows with defined tau (references and
    # positive controls in the CSV are not part of the H1 comparison)
    h1 = A["H1_tau"]
    ng = sum(1 for r in rows if r["group"] == "gated" and r["tau"] != "")
    nb = sum(1 for r in rows if r["group"] == "background" and r["tau"] != "")
    assert (h1["n_gated"], h1["n_background"]) == (ng, nb)

"""Hermetic checks on the committed tumor-vs-adjacent window audit."""
import csv, json, os

R = os.path.join(os.path.dirname(__file__), "..", "results")
D = json.load(open(os.path.join(R, "firebrowse_window.json")))
G = D["per_gene"]
TARGETS = ["CLDN18", "MSLN", "CA9", "CA12", "CLDN6", "PSCA", "SLC39A6"]


def test_shape_and_sample_records():
    assert D["n_sample_records"] == 2617
    assert set(D["cohorts"]) == {"PAAD", "COAD", "LUAD", "GBM", "BRCA"}
    rows = list(csv.DictReader(open(os.path.join(R, "firebrowse_window_rows.csv"))))
    assert len({(r["cohort"], r["participant"], r["sample_type"]) for r in rows}) == 2617
    assert G["CLDN18"]["cells"]["PAAD"]["n_nt"] == 4  # small-N normals documented


def test_ca9_is_the_only_strong_window():
    assert G["CA9"]["median_delta"] >= 3.0
    assert G["CA9"]["cohorts_delta_ge1_sig"] == G["CA9"]["n_cohorts"]
    for g in TARGETS:
        if g == "CA9":
            continue
        assert abs(G[g]["median_delta"]) <= 1.0, g


def test_noise_floor_reported_and_ca9_exceeds_it():
    floor = D["controls"]["hk_adj_max_abs_delta"]
    assert floor >= 1.0  # housekeeping residual spread = normalization noise floor
    assert G["CA9"]["median_delta"] > 2 * floor


def test_controls():
    assert 0.0 < D["controls"]["erbb2_brca_delta_adj"] < 1.5
    assert len(D["housekeeping_offset_per_cohort"]) == 5

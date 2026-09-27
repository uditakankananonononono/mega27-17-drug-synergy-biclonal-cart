"""Rankscore v1 structural tests (amendment queue #1/#12/#13)."""
import json, subprocess, sys, os

def _out():
    if not os.path.exists("results/rankscore_v1.json"):
        subprocess.run([sys.executable, "scripts_rankscore.py"], check=True)
    return json.load(open("results/rankscore_v1.json"))

def test_panel_and_weights():
    out = _out()
    assert len(out["ranking"]) == 7 == out["panel_size"]
    assert abs(sum(out["weights"].values()) - 1.0) < 1e-9
    assert set(out["ranking"]) == set(out["rows"])

def test_contributions_consistent():
    out = _out()
    for sym, row in out["rows"].items():
        for k, w in out["weights"].items():
            assert abs(row["contributions"][k] - round(w * row["components"][k], 4)) < 1e-9
        assert 0.0 <= row["score"] <= 1.0

def test_validation_field_honest():
    out = _out()
    v = out["validation"]
    assert set(v["clinically_targeted_anchors"]) == {"MSLN", "CLDN18", "PSCA"}
    assert isinstance(v["anchors_in_top_half"], bool)  # may be False; reported, not forced

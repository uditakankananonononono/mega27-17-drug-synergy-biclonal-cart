"""Hermetic checks on the committed DrugCentral engagement audit."""
import csv, json, os

R = os.path.join(os.path.dirname(__file__), "..", "results")
D = json.load(open(os.path.join(R, "drugcentral_engagement.json")))
AND_GATE = ["CLDN18", "MSLN", "CA9", "CA12", "CLDN6", "PSCA", "SLC39A6"]


def test_shape():
    g = D["per_gene"]
    assert sum(1 for d in g.values() if d["group"] == "gated") == 99
    assert sum(1 for d in g.values() if d["group"] == "background") == 98
    assert sum(1 for d in g.values() if d["group"] == "reference") == 4
    assert D["source"]["rows_human"] == 14301
    assert D["source"]["human_genes"] == 1704


def test_controls():
    refs = D["references"]
    assert refs["CD19"]["tclin"] and refs["CD19"]["n_drugs"] == 4
    assert "blinatumomab" in refs["CD19"]["drugs"]
    assert refs["ERBB2"]["tclin"] and refs["ERBB2"]["n_drugs"] >= 30


def test_and_gate_split():
    ag = D["and_gate"]
    for g in ["CLDN18", "MSLN", "CLDN6", "PSCA", "SLC39A6"]:
        assert ag[g]["n_rows"] == 0, g  # invisible to the small-molecule pharmacopeia
    assert ag["CA9"]["n_drugs"] >= 40 and ag["CA9"]["tclin"]
    assert ag["CA12"]["n_drugs"] >= 40 and ag["CA12"]["tclin"]
    assert "acetazolamide" in ag["CA9"]["drugs"]


def test_nulls_and_rows_csv():
    t = D["tests"]
    for key in ("any", "tclin", "moa"):
        assert t[key]["p"] > 0.5, key  # declared null; fail loudly if rerun flips it
        assert 0 < t[key]["odds"] < 2
    assert D["summary"]["gated"]["any"] == 17 and D["summary"]["background"]["any"] == 19
    rows = list(csv.DictReader(open(os.path.join(R, "drugcentral_engagement_rows.csv"))))
    assert len(rows) == 201
    assert sum(1 for r in rows if r["group"] == "gated" and int(r["n_rows"]) > 0) == 17

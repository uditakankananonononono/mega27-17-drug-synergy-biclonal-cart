"""Fix-wave lanes: MyGene/HPA restore, RxNorm, MyChem, openFDA, DailyMed,
Cellosaurus. Hermetic - asserts on committed results only."""
import csv, json
from pathlib import Path

RES = Path(__file__).resolve().parents[1] / "results"


def rows(name):
    return list(csv.DictReader(open(RES / name)))


def test_mygene_hpa_restore():
    s = json.load(open(RES / "mygene_summary.json"))
    assert s["n_resolved"] == 74 and s["n_hpa_tsv_fetched"] == 74
    m = rows("mygene_symbol_ensg.csv")
    assert len(m) == 74 and all(r["ensembl_gene"].startswith("ENSG") for r in m)


def test_rxnorm_map():
    s = json.load(open(RES / "rxnorm_summary.json"))
    assert s["n_drugs"] == 621
    assert 60 <= s["n_rxcui_matched"] <= 300  # investigational codes unmatched
    r = rows("rxnorm_gdsc_map.csv")
    assert all(x["match"] in {"exact", "unmatched", "error"} for x in r)


def test_mychem_annotations():
    s = json.load(open(RES / "mychem_summary.json"))
    assert s["n_annotated"] >= 60
    r = rows("mychem_drug_annotations.csv")
    assert sum(1 for x in r if x["inchikey"]) == s["n_annotated"]


def test_openfda_dailymed():
    s = json.load(open(RES / "openfda_summary.json"))
    assert s["n_with_fda_label"] >= 40 and s["n_boxed_warning"] >= 10
    d = json.load(open(RES / "dailymed_summary.json"))
    assert d["n_with_spl"] >= 40


def test_cellosaurus():
    s = json.load(open(RES / "cellosaurus_summary.json"))
    assert s["n_cell_lines"] == s["n_verified"] >= 8

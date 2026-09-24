"""Hermetic tests for the BioPlex systematic AP-MS audit (offline; committed
result files only)."""
import csv, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda *p: os.path.join(ROOT, *p)
ANTIGENS = {"CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"}


def load(name):
    with open(R("results", name)) as fh:
        return json.load(fh)


def test_audit_schema_and_gates():
    a = load("bioplex_audit.json")
    assert a["n_genes_total"] == 205
    for gname in ("G1_psma1_proteasome_both_networks",
                  "G2_polr2a_polII_partner_293T", "G3_references_present_293T"):
        assert a["gates"][gname]["pass"], f"calibration gate {gname} failed"
    assert a["all_gates_pass"]


def test_antigen_pair_null_and_testability_bounds():
    a = load("bioplex_audit.json")
    for n in ("293T", "HCT116"):
        h = a["H1_per_network"][n]
        assert h["antigen_pair_edges"] == 0
        assert h["testable_antigen_pairs"] <= 21
    # CA9 and CLDN6 undetected in both cell lines -> their pairs untestable
    pres = a["and_gate_presence"]
    assert not pres["CA9"]["293T"] and not pres["CA9"]["HCT116"]
    assert not pres["CLDN6"]["293T"] and not pres["CLDN6"]["HCT116"]


def test_per_gene_covers_universe():
    rows = list(csv.DictReader(open(R("results", "bioplex_per_gene.csv"))))
    assert len(rows) == 205
    assert {r["group"] for r in rows} == {"gated", "background", "references", "positive_controls"}
    assert ANTIGENS <= {r["gene"] for r in rows}


def test_concordance_partitions():
    a = load("bioplex_audit.json")
    h3 = a["H3_reproducibility_concordance"]
    assert h3["bioplex_and_intact"] + h3["bioplex_only"] == h3["edges_either"]
    assert h3["bioplex_and_intact"] + h3["intact_only"] == h3["intact_universe_edges"]
    assert 0.0 <= h3["cross_network_jaccard"] <= 1.0

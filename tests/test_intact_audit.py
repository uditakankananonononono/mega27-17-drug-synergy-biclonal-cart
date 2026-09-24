"""Hermetic tests for the IntAct curated-interaction audit (offline; uses
only committed result files)."""
import csv, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda *p: os.path.join(ROOT, *p)
ANTIGENS = {"CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"}
REFS = {"CD19", "ERBB2", "FOLR1", "TNFRSF17"}
CTRLS = {"POLR2A", "RPS3", "PCNA", "PSMA1"}
PHYS = {"MI:0407", "MI:0914", "MI:0915", "MI:2364"}


def load(name):
    with open(R("results", name)) as fh:
        return json.load(fh)


def test_audit_schema_and_gates():
    a = load("intact_audit.json")
    assert a["n_genes_total"] == 205
    assert a["n_genes_with_accession"] >= 200
    assert set(a["physical_mi_types"]) == PHYS
    for gname in ("G1_erbb2_canonical_partners", "G2_count_mitab_consistency", "G3_controls_nonzero"):
        assert a["gates"][gname]["pass"], f"calibration gate {gname} failed"


def test_erbb2_canonical_edges_present():
    edges = {tuple(sorted((e["gene_a"], e["gene_b"]))) for e in
             csv.DictReader(open(R("results", "intact_edges.csv")))}
    # EGFR is inside the audited universe; ERBB3 is not, so its canonical
    # edge is certified by gate G1 in the audit JSON instead of the edge file
    assert tuple(sorted(("ERBB2", "EGFR"))) in edges
    g1 = load("intact_audit.json")["gates"]["G1_erbb2_canonical_partners"]
    assert g1["egfr"] and g1["erbb3"]


def test_per_gene_covers_universe_and_groups():
    rows = list(csv.DictReader(open(R("results", "intact_per_gene.csv"))))
    assert len(rows) == 205
    groups = {r["group"] for r in rows}
    assert groups == {"gated", "background", "references", "positive_controls"}
    assert ANTIGENS | REFS | CTRLS <= {r["gene"] for r in rows}


def test_h1_pair_counts_consistent():
    a = load("intact_audit.json")
    h1 = a["H1_pairwise_antigen_complex"]
    hits = h1["antigen_pairs_with_edge"]
    assert h1["antigen_pairs_total"] == 21
    assert len(hits) <= 21
    assert abs(h1["antigen_pair_fraction"] - len(hits) / 21) < 1e-3
    for pair in hits:
        assert set(pair) <= ANTIGENS and len(set(pair)) == 2
    assert 0.0 <= h1["perm_uniform_p"] <= 1.0
    if h1["perm_degmatched_p"] is not None:
        assert 0.0 <= h1["perm_degmatched_p"] <= 1.0


def test_edges_only_physical_types():
    for e in csv.DictReader(open(R("results", "intact_edges.csv"))):
        assert set(e["types"].split("|")) <= PHYS


def test_h2_degree_fields():
    a = load("intact_audit.json")
    h2 = a["H2_degree"]
    assert 0.0 <= h2["mw_p"] <= 1.0
    assert "stratified" in h2
    for t, s in h2["stratified"].items():
        assert t in {"single-pass", "multi-pass", "GPI", "unannotated"}
        assert 0.0 <= s["p"] <= 1.0

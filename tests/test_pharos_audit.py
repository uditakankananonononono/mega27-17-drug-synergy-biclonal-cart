"""Hermetic tests for the Pharos/TCRD illumination audit (offline; uses only
committed result files)."""
import csv, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda *p: os.path.join(ROOT, *p)
ANTIGENS = {"CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"}
REFS = {"CD19", "ERBB2", "FOLR1", "TNFRSF17"}
CTRLS = {"POLR2A", "RPS3", "PCNA", "PSMA1"}
TDLS = {"Tclin", "Tchem", "Tbio", "Tdark"}


def load(name):
    with open(R("results", name)) as fh:
        return json.load(fh)


def test_audit_schema_and_gates():
    a = load("pharos_audit.json")
    assert a["n_genes_total"] == 205
    assert a["n_records_ok"] >= 200
    for gname in ("G1_references_tclin", "G2_coverage_ge_200_of_205",
                  "G3_controls_not_tclin", "G4_egfr_tclin_both_sources"):
        assert a["gates"][gname]["pass"], f"calibration gate {gname} failed"
    assert a["all_gates_pass"]


def test_references_and_controls_direction():
    a = load("pharos_audit.json")
    for ref in REFS:
        assert a["gates"]["G1_references_tclin"]["detail"][ref] == "Tclin"
    not_tclin = sum(1 for c in CTRLS
                    if a["gates"]["G3_controls_not_tclin"]["detail"][c] != "Tclin")
    assert not_tclin >= 3
    assert a["egfr_sentinel"]["pharos_tdl"] == "Tclin"


def test_per_gene_covers_universe():
    rows = list(csv.DictReader(open(R("results", "pharos_per_gene.csv"))))
    assert len(rows) == 205
    assert {r["group"] for r in rows} == {"gated", "background", "references", "positive_controls"}
    ok = [r for r in rows if r["status"] == "ok"]
    assert len(ok) >= 200
    for r in ok:
        assert r["tdl"] in TDLS
        assert int(r["n_drugs"]) >= 0 and int(r["n_ligands"]) >= 0
    assert ANTIGENS <= {r["gene"] for r in rows}


def test_h1_table_consistency():
    a = load("pharos_audit.json")
    h1 = a["H1_tclin"]
    assert h1["gated_tclin"] <= h1["gated_n"] and h1["bg_tclin"] <= h1["bg_n"]
    mix_sum_g = sum(h1["tdl_mix_gated"].values())
    assert mix_sum_g == h1["gated_n"]


def test_concordance_partitions():
    a = load("pharos_audit.json")
    h4 = a["H4_concordance"]
    assert h4["both_tclin"] + h4["pharos_only"] + h4["drugcentral_only"] + h4["neither"] == h4["n"]
    assert len(h4["discordant_genes"]) == h4["pharos_only"] + h4["drugcentral_only"]
    assert 0.0 <= h4["agreement"] <= 1.0

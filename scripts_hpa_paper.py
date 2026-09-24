#!/usr/bin/env python3
"""Generate paper/hpa_sec.tex from paper/hpa_tpl.tex using committed JSON."""
import json, re

A = json.load(open("results/hpa_audit.json"))
tpl = open("paper/hpa_tpl.tex").read()

def pct(k, n):
    return f"{k}/{n} ({100*k/n:.0f}\\%)" if n else "0/0"

g = A["calibration_gates"]
pm = A["predicted_membrane"]
cmh = A["topology_stratified_cmh"] or {"common_odds_ratio": float("nan"), "p_value": float("nan"), "strata_used": []}
rows = []
for r in A["and_gate_antigens"] + A["references_detail"]:
    role = "AND-gate" if r in A["and_gate_antigens"] else "reference"
    rows.append("{} & {} & {} & {} & {} \\\\".format(
        r["gene"], role,
        "yes" if r["predicted_membrane"] else "no",
        "yes" if r["if_plasma_membrane"] else ("no" if r["if_has_data"] else "--"),
        (r["if_main_location"] or "--").replace("_", "\\_")))
rep = {
    "@NGENES@": str(A["n_genes"]),
    "@NOK@": str(A["n_with_hpa_record"]),
    "@NNOREC@": str(len(A["no_record"])),
    "@NORECLIST@": ", ".join(A["no_record"]) if A["no_record"] else "none",
    "@GATEONE@": pct(*g["G1_refs_predicted_membrane"]),
    "@GATETWO@": pct(*g["G2_controls_predicted_membrane"]),
    "@GATETHREE@": pct(*g["G3_refs_if_pm_with_data"]),
    "@GATEFOUR@": pct(*g["G4_controls_if_pm"]),
    "@GATEWORD@": "PASSED" if A["calibration_gate_pass"] else "FAILED",
    "@GATEDPM@": pct(*pm["gated"]),
    "@BGPM@": pct(*pm["background"]),
    "@FISHEROR@": f"{pm['fisher_odds_ratio']:.2f}",
    "@FISHERP@": f"{pm['fisher_p']:.3g}",
    "@CMHOR@": f"{cmh['common_odds_ratio']:.2f}",
    "@CMHP@": f"{cmh['p_value']:.3g}",
    "@NSTRATA@": str(len(cmh["strata_used"])),
    "@IFCOVG@": pct(*A["if_coverage"]["gated"]),
    "@IFCOVB@": pct(*A["if_coverage"]["background"]),
    "@IFCOVOR@": f"{A['if_coverage_fisher']['odds_ratio']:.2f}",
    "@IFCOVP@": f"{A['if_coverage_fisher']['p']:.3g}",
    "@ANTIGENROWS@": "\n".join(rows),
}
out = tpl
for k, v in rep.items():
    out = out.replace(k, v)
assert "@" not in re.sub(r"@\\", "", out), "unreplaced token remains"
open("paper/hpa_sec.tex", "w").write(out)
print("wrote paper/hpa_sec.tex; gate", A["calibration_gate_pass"])

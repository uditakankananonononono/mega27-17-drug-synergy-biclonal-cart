"""Generate paper/clintrials_sec.tex from paper/clintrials_tpl.tex + committed
results/clintrials_audit.json (token replace; no hand-typed numbers)."""
import json

A = json.load(open("results/clintrials_audit.json"))
t = open("paper/clintrials_tpl.tex").read()
S = A["stats"]
C = A["curation"]
CAL = A["calibration"]
tab = {r["gene"]: r for r in A["alias_table"]}

def esc(s):
    return s.replace("_", "\\_")

def fmt(r):
    return ("$\\geq$" if r["capped"] else "") + str(r["total"])

rows = []
for r in A["alias_table"]:
    cls = "AND-gate" if r["class"] == "and_gate" else "reference"
    rows.append(f"{esc(r['gene'])} & {cls} & {fmt(r)} & {r['active']} & {r['verified_strict']} \\\\")
alias_rows = "\n".join(rows)

rep = {
    "@NGENES@": str(A["n_genes"]),
    "@CD19VER@": str(CAL["cd19_verified"]),
    "@NEGSUM@": str(CAL["noise_controls_verified_sum"]),
    "@PCNAVER@": str(CAL["noise_controls"]["PCNA"]),
    "@RAWMWP@": f"{S['strict_onc_raw_mw']['p']:.4f}",
    "@NINSPECT@": str(C["n_inspected"]),
    "@NCOLL@": str(C["n_collision_or_biomarker"]),
    "@NGEN@": str(C["n_genuine"]),
    "@VERMWP@": f"{S['verified_mw']['p']:.3f}",
    "@GANY@": S["verified_fisher"]["gated_any"],
    "@BANY@": S["verified_fisher"]["bg_any"],
    "@VFISHP@": f"{S['verified_fisher']['p']:.3f}",
    "@VACTMWP@": f"{S['verified_active_mw']['p']:.3f}",
    "@MSLNALIAS@": str(tab["MSLN"]["total"]),
    "@MSLNACT@": str(tab["MSLN"]["active"]),
    "@MSLNSTRICT@": str(tab["MSLN"]["verified_strict"]),
    "@CLDN18ALIAS@": str(tab["CLDN18"]["total"]),
    "@CLDN6ALIAS@": str(tab["CLDN6"]["total"]),
    "@PSCAALIAS@": str(tab["PSCA"]["total"]),
    "@CA9ALIAS@": str(tab["CA9"]["total"]),
    "@CA12ALIAS@": str(tab["CA12"]["total"]),
    "@SLCALIAS@": str(tab["SLC39A6"]["total"]),
    "@CD19ALIAS@": str(tab["CD19"]["total"]),
    "@ERBB2ALIAS@": str(tab["ERBB2"]["total"]),
    "@BCMAALIAS@": str(tab["TNFRSF17"]["total"]),
    "@FOLR1ALIAS@": str(tab["FOLR1"]["total"]),
    "@NUNIQNCT@": f"{A['datasets']['unique_nct_globally']:,}",
    "@ALIASROWS@": alias_rows,
}
for k, v in rep.items():
    assert k in t, k
    t = t.replace(k, v)
assert "@" not in t, "unreplaced token"
open("paper/clintrials_sec.tex", "w").write(t)
print("wrote paper/clintrials_sec.tex", len(t.splitlines()), "lines")

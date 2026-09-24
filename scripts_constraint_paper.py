#!/usr/bin/env python3
"""Generate paper/constraint_sec.tex from /tmp template + committed results/constraint_audit.json (token replace; no hand-typed numbers)."""
import json

A = json.load(open("results/constraint_audit.json"))
t = open("/tmp/tpl1.tex").read()
ctrl = A["calibration"]["controls"]
lo, hi = min(ctrl.values()), max(ctrl.values())
def f(x, n=2):
    return f"{x:.{n}f}"
rep = {
    "@CTRLRANGE@": f"{lo:.2f}-{hi:.2f}",
    "@NCTRL@": str(A["calibration"]["n_constrained"]),
    "@GMED@": f(A["loeuf"]["gated_median"]),
    "@BMED@": f(A["loeuf"]["background_median"]),
    "@LMWP@": f"{A['loeuf']['mw_p']:.3f}",
    "@GC@": str(A["loeuf"]["gated_constrained"]),
    "@GN@": str(A["loeuf"]["n_gated"]),
    "@BC@": str(A["loeuf"]["background_constrained"]),
    "@BN@": str(A["loeuf"]["n_background"]),
    "@LFP@": f"{A['loeuf']['fisher_p']:.3f}",
    "@GGS@": str(int(A["paralog_groups"]["gated_median_group_size"])),
    "@BGS@": str(int(A["paralog_groups"]["background_median_group_size"])),
    "@PGSP@": f"{A['paralog_groups']['mw_p']:.3f}",
    "@GRX@": str(int(A["reactome"]["gated_median"])),
    "@BRX@": str(int(A["reactome"]["background_median"])),
    "@RMWP@": f"{A['reactome']['mw_p']:.4f}",
    "@CRMWP@": f"{A['protein_coding_control']['reactome_mw_p']:.4f}",
    "@AGFREE@": str(A["and_gate_summary"]["free_loss"]),
    "@AGMIN@": f(A["and_gate_summary"]["min_loeuf"]),
    "@AGFRAC@": f"{100*A['and_gate_summary']['max_frac_dep']:.1f}\\%",
    "@RHO@": f(A["loeuf_vs_depmap"]["spearman_rho"]),
    "@RHOP@": f"{A['loeuf_vs_depmap']['p']:.3f}",
}
for k, v in rep.items():
    assert k in t, k
    t = t.replace(k, v)
assert "@" not in t
open("paper/constraint_sec.tex", "w").write(t)
print("wrote paper/constraint_sec.tex", len(t.splitlines()), "lines")

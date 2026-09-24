#!/usr/bin/env python3
"""Render paper/glygen_sec.tex from paper/glygen_tpl.tex + results/glygen_audit.json (token replace; no hand-typed numbers)."""
import json

A = json.load(open("results/glygen_audit.json"))
S = A["surface_stratum"]
ST = A["topology_strata"]
E = A["ectodomain_density"]
C = A["calibration"]
AG = {r["gene"]: r for r in A["and_gate"]}
RF = {r["gene"]: r for r in A["references"]}


def p(x):
    if x < 1e-4:
        return "%.1e" % x
    return "%.4f" % x if x < 0.01 else "%.3f" % x


t = open("paper/glygen_tpl.tex").read()
tok = {
    "@NGENES@": str(A["n_genes"]),
    "@NRECS@": str(A["n_record_present"]),
    "@NNOREC@": str(A["n_genes"] - A["n_record_present"]),
    "@MSLN@": str(AG["MSLN"]["n_shield_sites"]),
    "@MSLNS@": str(AG["MSLN"]["n_shield_structures"]),
    "@CD19@": str(RF["CD19"]["n_shield_sites"]),
    "@ERBB2@": str(RF["ERBB2"]["n_shield_sites"]),
    "@NEGMAX@": "%.3f" % C["neg_max_density"],
    "@POSMIN@": "%.2f" % C["pos_min_density"],
    "@GATE@": "3/3 + 4/4 + ordering" if C["gate_pass"] else "FAILED",
    "@GANY@": str(S["gated_with_shield"]),
    "@NG@": str(S["n_gated"]),
    "@BANY@": str(S["background_with_shield"]),
    "@NB@": str(S["n_background"]),
    "@FISHP@": p(S["fisher_p"]),
    "@MWP@": p(S["mw_p"]),
    "@SPP@": p(ST["single-pass"]["fisher_p"]),
    "@MPP@": p(ST["multi-pass"]["fisher_p"]),
    "@GPIP@": p(ST["GPI"]["fisher_p"]),
    "@EGMED@": "%.2f" % E["gated_median"],
    "@EBMED@": "%.2f" % E["background_median"],
    "@ECTOP@": p(E["mw_p"]),
    "@SLC@": str(AG["SLC39A6"]["n_shield_sites"]),
    "@SLCS@": str(AG["SLC39A6"]["n_shield_structures"]),
    "@PSCA@": str(AG["PSCA"]["n_shield_sites"]),
    "@CLDN6@": str(AG["CLDN6"]["n_shield_sites"]),
    "@FOLR1S@": str(RF["FOLR1"]["n_shield_structures"]),
}
for k, v in tok.items():
    t = t.replace(k, v)
assert "@" not in t, "unreplaced token"
open("paper/glygen_sec.tex", "w").write(t)
print("glygen_sec.tex written,", len(tok), "tokens")

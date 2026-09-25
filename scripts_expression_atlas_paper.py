#!/usr/bin/env python3
"""Generate paper/eatlas_sec.tex from paper/eatlas_tpl.tex + committed
results/expression_atlas_audit.json and results/proteomicsdb_audit.json.
Token replacement only."""
import json

from scripts_harmonizome_paper import fmtp


def main():
    a = json.load(open("results/expression_atlas_audit.json"))
    p = json.load(open("results/proteomicsdb_audit.json"))
    h1, h2 = a["H1_tau"], a["H2_cross_platform"]
    # template text hard-codes both verdicts; refuse to render otherwise
    assert a["all_gates_pass"]
    assert h1["verdict"] == "FALSIFIED" and h2["verdict"] == "CONFIRMED"
    ants = a["antigens"]
    specific = [s for s, r in ants.items() if r["tau"] >= 0.7]
    assert len(specific) == 6 and ants["SLC39A6"]["tau"] < 0.3
    blk = []
    for s in sorted(ants):
        r = ants[s]
        blk.append("%s $\\tau=%.2f$ (%d tissues $\\ge1$ TPM; max %s, %g TPM)"
                   % (s, r["tau"], r["n_tissues_tpm1"], r["max_tissue"], r["max_tpm"]))
    g = a["calibration_gates"]
    g3 = g["G3_posctrl_tau_lt_0.5"]["obs"]
    pa = p["antigens"]
    rep = {
        "@NGENES@": a["n_genes_total"],
        "@NMATCH@": g["G1_matched"]["obs"],
        "@CD19TPM@": "%g" % g["G2_CD19_lymph_node_tpm"]["obs"],
        "@G3LO@": "%.2f" % min(g3.values()),
        "@G3HI@": "%.2f" % max(g3.values()),
        "@H1GMED@": "%.3f" % h1["median_gated"],
        "@H1BMED@": "%.3f" % h1["median_background"],
        "@H1NG@": h1["n_gated"],
        "@H1NB@": h1["n_background"],
        "@H1P@": fmtp(h1["p"]),
        "@H2N@": h2["n_pairs"],
        "@H2RHO@": "%.2f" % h2["spearman_rho"],
        "@H2P@": fmtp(h2["p"]),
        "@ANTIGENBLOCK@": "; ".join(blk) + ".",
        "@SLCPROT@": pa["SLC39A6"]["n_normal_tissues"],
        "@CA12PROT@": pa["CA12"]["n_normal_tissues"],
        "@MSLNPROT@": pa["MSLN"]["n_normal_tissues"],
    }
    tpl = open("paper/eatlas_tpl.tex").read()
    for k, v in rep.items():
        tpl = tpl.replace(k, str(v))
    assert "@" not in tpl
    with open("paper/eatlas_sec.tex", "w") as f:
        f.write(tpl)
    print("wrote paper/eatlas_sec.tex (%d chars)" % len(tpl))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate paper/monarch_sec.tex from paper/monarch_tpl.tex + committed
results/monarch_audit.json. Token replacement only."""
import json

ANTS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]


def main():
    a = json.load(open("results/monarch_audit.json"))
    h1, h2, h3, h4 = (a["H1_causal_carriage"], a["H2_pheno_breadth"],
                      a["H3_correlated_carriage"], a["H4_cancer_causal_carriage"])
    ps = [h1["fisher_p"], h2["mw_p"], h3["fisher_p"], h4["fisher_p"]]
    if all(p >= 0.05 for p in ps):
        v = ("No Mendelian-causality or phenotypic-breadth axis distinguishes gated "
             "antigens from background: the panel's surface antigens are as free of "
             "germline Mendelian disease entanglement as the matched background, "
             "completing the germline picture alongside the common-variant GWAS null "
             "and the somatic-selection constraint audit.")
    else:
        v = ("At least one Mendelian axis differs between gated and background, so "
             "germline Mendelian-disease entanglement must be examined per antigen "
             "alongside the common-variant and somatic-constraint results.")
    blk = []
    for s in ANTS:
        r = a["and_gate_antigens"][s]
        if r["causal_total"] == 0:
            blk.append("%s: no causal disease associations" % s)
        else:
            blk.append("%s: %d causal (%s)" % (s, r["causal_total"],
                                               "; ".join(r["causal_labels"][:3])))
    rep = {
        "@NGENES@": a["n_genes_total"], "@NRES@": a["n_resolved"],
        "@CDTWO@": a["gates"]["G2_cd19_causal_ge_1"]["cd19_causal"],
        "@CDTHREE@": a["gates"]["G3_cd19_pheno_ge_10"]["cd19_pheno"],
        "@H1G@": h1["gated"], "@H1GN@": h1["gated_n"], "@H1B@": h1["bg"], "@H1BN@": h1["bg_n"],
        "@H1OR@": "%.2f" % h1["fisher_or"], "@H1P@": "%.3g" % h1["fisher_p"],
        "@H2G@": "%g" % h2["gated_median"], "@H2B@": "%g" % h2["bg_median"],
        "@H2P@": "%.3g" % h2["mw_p"],
        "@H3G@": h3["gated"], "@H3GN@": h3["gated_n"], "@H3B@": h3["bg"], "@H3BN@": h3["bg_n"],
        "@H3OR@": "%.2f" % h3["fisher_or"], "@H3P@": "%.3g" % h3["fisher_p"],
        "@H4G@": h4["gated"], "@H4GN@": h4["gated_n"], "@H4B@": h4["bg"], "@H4BN@": h4["bg_n"],
        "@H4P@": "%.3g" % h4["fisher_p"],
        "@VERDICT@": v, "@ANTIGENBLOCK@": "; ".join(blk) + ".",
    }
    tpl = open("paper/monarch_tpl.tex").read()
    for k, vv in rep.items():
        tpl = tpl.replace(k, str(vv))
    assert "@" not in tpl
    with open("paper/monarch_sec.tex", "w") as f:
        f.write(tpl)
    print("wrote paper/monarch_sec.tex (%d chars)" % len(tpl))


if __name__ == "__main__":
    main()

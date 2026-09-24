#!/usr/bin/env python3
"""Generate paper/gwas_sec.tex from paper/gwas_tpl.tex + committed
results/gwas_audit.json. Token replacement only; asserts no '@' remains."""
import json

ANTS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]

def main():
    a = json.load(open("results/gwas_audit.json"))
    h1, h2, h3, g = a["H1_any_sig"], a["H2_cancer_sig"], a["H3_burden"], a["gates"]
    if h1["fisher_p"] < 0.05 and h2["fisher_p"] >= 0.05:
        v = ("Gated loci carry significantly more germline association signal overall, but the "
             "excess vanishes for cancer traits: GWAS mapping is position-based rather than "
             "curation-based, so the any-trait enrichment reads as locus pleiotropy (or LD-block "
             "size), while germline cancer relevance is NOT enriched in the panel. Antigen "
             "selection by tumor-surface expression does not capture population-level cancer risk "
             "loci; PSCA is the lone strong exception.")
    elif h1["fisher_p"] < 0.05 and h2["fisher_p"] < 0.05:
        v = ("Gated loci are enriched for germline association signal including cancer traits, "
             "so the panel partially overlaps population-level cancer risk loci.")
    else:
        v = ("Germline association signal does not distinguish gated from background loci on "
             "either axis; antigen selection by tumor-surface expression is orthogonal to "
             "population-level trait genetics.")
    blk = []
    for s in ANTS:
        r = a["and_gate_antigens"][s]
        blk.append("%s: %d significant (%d cancer, %d SNPs)" %
                   (s, r["n_sig"], r["n_sig_cancer"], r["n_sig_snps"]))
    rep = {
        "@NROWS@": "{:,}".format(a["n_assoc_rows"]),
        "@GATEONE@": "%d associations: %s" % (g["G1_psca_sig_cancer"]["detail"]["n_sig_cancer"],
                     ", ".join(g["G1_psca_sig_cancer"]["detail"]["traits"][:3])),
        "@GATETHREE@": "%d significant cancer associations" % a["egfr_sentinel"]["n_sig_cancer"],
        "@H1G@": h1["gated"], "@H1GN@": h1["gated_n"], "@H1B@": h1["bg"], "@H1BN@": h1["bg_n"],
        "@H1OR@": "%.2f" % h1["fisher_or"], "@H1P@": "%.3g" % h1["fisher_p"],
        "@H2G@": h2["gated"], "@H2GN@": h2["gated_n"], "@H2B@": h2["bg"], "@H2BN@": h2["bg_n"],
        "@H2OR@": "%.2f" % h2["fisher_or"], "@H2P@": "%.3g" % h2["fisher_p"],
        "@H3G@": "%g" % h3["gated_median"], "@H3B@": "%g" % h3["bg_median"],
        "@H3P@": "%.3g" % h3["mw_p"],
        "@VERDICT@": v,
        "@ANTIGENBLOCK@": "; ".join(blk) + ".",
    }
    tpl = open("paper/gwas_tpl.tex").read()
    for k, vv in rep.items():
        tpl = tpl.replace(k, str(vv))
    assert "@" not in tpl
    with open("paper/gwas_sec.tex", "w") as f:
        f.write(tpl)
    print("wrote paper/gwas_sec.tex (%d chars)" % len(tpl))

if __name__ == "__main__":
    main()

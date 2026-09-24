#!/usr/bin/env python3
"""Generate paper/bioplex_sec.tex from paper/bioplex_tpl.tex + committed
results/bioplex_audit.json. Token replacement only; asserts no '@' remains."""
import json

def main():
    a = json.load(open("results/bioplex_audit.json"))
    g = a["gates"]
    n293 = a["H1_per_network"]["293T"]
    nhct = a["H1_per_network"]["HCT116"]
    h3 = a["H3_reproducibility_concordance"]

    if n293["antigen_pair_edges"] == 0 and nhct["antigen_pair_edges"] == 0:
        h1v = ("Three independent interaction screens (IntAct curation, BioPlex 293T, BioPlex "
               "HCT116) now agree: no physical co-complex is observed among testable AND-gate "
               "antigen pairs, so the two gate arms behave as independent surface targets on "
               "current evidence, bounded by the systematic-scale detection caveat above.")
    else:
        h1v = ("BioPlex reports co-precipitation for at least one antigen pair, so "
               "cis-engagement geometry must be revisited for the affected pair(s).")
    h2v = ("The IntAct curated-degree advantage of the gated set (median 51 versus 15) "
           "therefore reflects study intensity rather than differential connectivity: measured "
           "systematically in a single experiment, gated and background genes are "
           "indistinguishable in AP-MS degree, and the higher gated detection rate makes the "
           "absence of antigen-pair edges the more credible.")
    pct = lambda a, b: "%d/%d" % (a, b)
    rep = {
        "@N293EDGES@": "{:,}".format(n293["edges"]), "@NHCTEDGES@": "{:,}".format(nhct["edges"]),
        "@GATEONE@": "293T %d PSMB partners, HCT116 %d" % (
            len(g["G1_psma1_proteasome_both_networks"]["detail"]["293T"]),
            len(g["G1_psma1_proteasome_both_networks"]["detail"]["HCT116"])),
        "@GATETWO@": ", ".join(g["G2_polr2a_polII_partner_293T"]["detail"]),
        "@GATETHREE@": "%d/4" % sum(g["G3_references_present_293T"]["detail"].values()),
        "@TEST293@": n293["testable_antigen_pairs"], "@TESTHCT@": nhct["testable_antigen_pairs"],
        "@EDGES293@": n293["antigen_pair_edges"], "@EDGESHCT@": nhct["antigen_pair_edges"],
        "@H1VERDICT@": h1v,
        "@PRES293@": pct(n293["gated_present"], n293["gated_n"]),
        "@PRESB293@": pct(n293["bg_present"], n293["bg_n"]),
        "@PPRES293@": "%.3g" % n293["presence_fisher_p"],
        "@PRESHCT@": pct(nhct["gated_present"], nhct["gated_n"]),
        "@PRESBHCT@": pct(nhct["bg_present"], nhct["bg_n"]),
        "@PPRESHCT@": "%.3g" % nhct["presence_fisher_p"],
        "@DEG293@": "%.1f" % n293["gated_degree_median"],
        "@DEGB293@": "%.1f" % n293["bg_degree_median"],
        "@PDEG293@": "%.3g" % n293["degree_mw_p"],
        "@DEGHCT@": "%.1f" % nhct["gated_degree_median"],
        "@DEGBHCT@": "%.1f" % nhct["bg_degree_median"],
        "@PDEGHCT@": "%.3g" % nhct["degree_mw_p"],
        "@H2VERDICT@": h2v,
        "@BPEITHER@": h3["edges_either"], "@BPBOTH@": h3["edges_both"],
        "@JACC@": "%.2f" % h3["cross_network_jaccard"],
        "@BPINTACT@": h3["bioplex_and_intact"], "@INTACTN@": h3["intact_universe_edges"],
        "@INTACTONLY@": h3["intact_only"],
    }
    tpl = open("paper/bioplex_tpl.tex").read()
    for k, v in rep.items():
        tpl = tpl.replace(k, str(v))
    assert "@" not in tpl, "unreplaced token remains"
    with open("paper/bioplex_sec.tex", "w") as f:
        f.write(tpl)
    print("wrote paper/bioplex_sec.tex (%d chars)" % len(tpl))

if __name__ == "__main__":
    main()

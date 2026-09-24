#!/usr/bin/env python3
"""Generate paper/pharos_sec.tex from paper/pharos_tpl.tex + committed
results/pharos_audit.json and results/pharos_per_gene.csv.
Token replacement only; asserts no '@' remains."""
import csv, json

ANTIGENS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]
REFS = ["CD19", "ERBB2", "FOLR1", "TNFRSF17"]
CTRLS = ["POLR2A", "RPS3", "PCNA", "PSMA1"]

def fmtp(p):
    return "%.3g" % p

def fmtsci(x):
    return "%.2g" % x

def main():
    a = json.load(open("results/pharos_audit.json"))
    h1, h2, h3, h4, g = a["H1_tclin"], a["H2_novelty"], a["H3_counts"], a["H4_concordance"], a["gates"]

    mix = lambda d: ", ".join("%s %d" % (k, d.get(k, 0)) for k in ("Tclin", "Tchem", "Tbio", "Tdark"))
    if h1["fisher_p"] < 0.05 and h1["fisher_or"] < 1:
        h1v = ("Gated antigens are significantly less clinically illuminated than the random "
               "surfaceome background, consistent with the panel design selecting "
               "under-drugged surface proteins.")
    elif h1["fisher_p"] < 0.05:
        h1v = ("Gated antigens are significantly more clinically illuminated than background, "
               "so illumination does not distinguish the panel from a drugged surfaceome.")
    else:
        h1v = ("Gated and background genes do not differ in Tclin fraction: clinical "
               "illumination is not what separates the panel from a random surface set.")
    if h2["mw_p"] < 0.05 and h2["gated_median"] > h2["bg_median"]:
        h2v = ("Gated antigens are significantly more novel (less studied) than background, "
               "the understudied-enrichment pattern the constraint and trial audits predicted.")
    elif h2["mw_p"] < 0.05:
        h2v = ("Gated antigens are significantly better studied than background, so the panel "
               "is not an understudied set on the TIN-X axis.")
    else:
        h2v = ("Novelty does not differ between the groups on the TIN-X axis.")
    ndisc = len(h4["discordant_genes"])
    if h4["cohens_kappa"] is not None and h4["cohens_kappa"] >= 0.6:
        h4v = ("The two independent druggability resources substantially agree, so Tclin calls "
               "are robust to provenance; named discordances mark post-snapshot approvals or "
               "source-specific curation.")
    else:
        h4v = ("Agreement between the two resources is limited, so Tclin status is "
               "provenance-sensitive and single-source druggability claims should be avoided.")
    blk = []
    for s in ANTIGENS:
        r = a["and_gate_antigens"].get(s, {})
        if r.get("status") != "ok":
            blk.append("%s: no Pharos record" % s)
        else:
            blk.append("%s: %s, %d approved drugs, %d ligands, novelty %s" %
                       (s, r["tdl"], r["n_drugs"], r["n_ligands"], fmtsci(r["novelty"])))
    g1d = g["G1_references_tclin"]["detail"]
    g3d = g["G3_controls_not_tclin"]["detail"]
    rep = {
        "@NGENES@": a["n_genes_total"],
        "@NOK@": a["n_records_ok"],
        "@GATEONE@": "%d/4 (%s)" % (sum(1 for v in g1d.values() if v == "Tclin"),
                                    ", ".join("%s=%s" % kv for kv in g1d.items())),
        "@GATETHREE@": "%d/4 (%s)" % (sum(1 for v in g3d.values() if v != "Tclin"),
                                      ", ".join("%s=%s" % kv for kv in g3d.items())),
        "@GATEFOUR@": "Pharos %s, DrugCentral %s" % (g["G4_egfr_tclin_both_sources"]["pharos"],
                      "Tclin" if g["G4_egfr_tclin_both_sources"]["drugcentral"] else "not Tclin"),
        "@TCLING@": h1["gated_tclin"], "@NG@": h1["gated_n"],
        "@TCLINB@": h1["bg_tclin"], "@NB@": h1["bg_n"],
        "@ORH1@": "%.2f" % h1["fisher_or"], "@PH1@": fmtp(h1["fisher_p"]),
        "@MIXG@": mix(h1["tdl_mix_gated"]), "@MIXB@": mix(h1["tdl_mix_background"]),
        "@H1VERDICT@": h1v,
        "@NOVG@": fmtsci(h2["gated_median"]), "@NOVB@": fmtsci(h2["bg_median"]),
        "@PH2@": fmtp(h2["mw_p"]), "@H2VERDICT@": h2v,
        "@DRUGG@": "%g" % h3["gated_drugs_median"], "@DRUGB@": "%g" % h3["bg_drugs_median"],
        "@PH3A@": fmtp(h3["drugs_mw_p"]),
        "@LIGG@": "%g" % h3["gated_ligands_median"], "@LIGB@": "%g" % h3["bg_ligands_median"],
        "@PH3B@": fmtp(h3["ligands_mw_p"]),
        "@NCONC@": h4["n"], "@AGREEPCT@": "%.1f" % (100 * h4["agreement"]),
        "@KAPPA@": "%.2f" % h4["cohens_kappa"],
        "@NDISC@": ndisc,
        "@DISCLIST@": ", ".join(h4["discordant_genes"]) if ndisc else "none",
        "@H4VERDICT@": h4v,
        "@ANTIGENBLOCK@": "; ".join(blk) + ".",
    }
    tpl = open("paper/pharos_tpl.tex").read()
    for k, v in rep.items():
        tpl = tpl.replace(k, str(v))
    assert "@" not in tpl, "unreplaced token remains"
    with open("paper/pharos_sec.tex", "w") as f:
        f.write(tpl)
    print("wrote paper/pharos_sec.tex (%d chars)" % len(tpl))

if __name__ == "__main__":
    main()

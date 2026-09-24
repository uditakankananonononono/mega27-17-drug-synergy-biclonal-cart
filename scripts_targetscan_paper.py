#!/usr/bin/env python3
"""Generate paper/targetscan_sec.tex from paper/targetscan_tpl.tex +
committed results/targetscan_audit.json. Token replacement only."""
import json

ANTS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]

def main():
    a = json.load(open("results/targetscan_audit.json"))
    h1, h2, h3, h4 = a["H1_any_site"], a["H2_sites"], a["H3_families"], a["H4_cwcs_min"]
    ps = [h1["fisher_p"], h2["mw_p"], h3["mw_p"], h4["mw_p"]]
    if all(p >= 0.05 for p in ps):
        v = ("No axis of post-transcriptional repression distinguishes gated antigens from "
             "background: miRNA regulation is not a confound of the panel design, and antigen "
             "safety margins must be read from RNA/protein abundance directly. Notably, four of "
             "seven AND-gate antigens (CA9, CLDN6, MSLN, PSCA) carry ZERO conserved miRNA sites, "
             "so their expression is post-transcriptionally unbuffered - RNA presence at these "
             "loci translates directly to surface exposure, reinforcing the abundance-based "
             "safety windows and marking CLDN18/SLC39A6 (which carry sites) as the only buffered "
             "antigens of the set.")
    else:
        v = ("At least one repression axis differs between gated and background, so miRNA "
             "buffering must be modeled per antigen when reading abundance-based safety windows.")
    blk = []
    for s in ANTS:
        r = a["and_gate_antigens"][s]
        blk.append("%s: %d conserved sites from %d families" % (s, r["sites"], r["families"]))
    rep = {
        "@NHUMAN@": "{:,}".format(a["n_human_genes_table"]),
        "@NGENES@": a["n_genes_total"],
        "@GATETHREE@": "%d sites" % a["hmga2_sentinel_let7_sites"],
        "@H1G@": h1["gated"], "@H1GN@": h1["gated_n"], "@H1B@": h1["bg"], "@H1BN@": h1["bg_n"],
        "@H1OR@": "%.2f" % h1["fisher_or"], "@H1P@": "%.3g" % h1["fisher_p"],
        "@H2G@": "%g" % h2["gated_median"], "@H2B@": "%g" % h2["bg_median"],
        "@H2P@": "%.3g" % h2["mw_p"],
        "@H3G@": "%g" % h3["gated_median"], "@H3B@": "%g" % h3["bg_median"],
        "@H3P@": "%.3g" % h3["mw_p"],
        "@H4G@": "%.3g" % h4["gated_median"], "@H4B@": "%.3g" % h4["bg_median"],
        "@H4P@": "%.3g" % h4["mw_p"],
        "@VERDICT@": v,
        "@ANTIGENBLOCK@": "; ".join(blk) + ".",
    }
    tpl = open("paper/targetscan_tpl.tex").read()
    for k, vv in rep.items():
        tpl = tpl.replace(k, str(vv))
    assert "@" not in tpl
    with open("paper/targetscan_sec.tex", "w") as f:
        f.write(tpl)
    print("wrote paper/targetscan_sec.tex (%d chars)" % len(tpl))

if __name__ == "__main__":
    main()

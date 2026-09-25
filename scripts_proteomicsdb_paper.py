#!/usr/bin/env python3
"""Generate paper/proteomicsdb_sec.tex from paper/proteomicsdb_tpl.tex +
committed results/proteomicsdb_audit.json. Token replacement only."""
import json
from scripts_harmonizome_paper import fmtp

ANTS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]


def main():
    a = json.load(open("results/proteomicsdb_audit.json"))
    g, h1, h2, h3 = (a["calibration_gates"], a["H1_normal_breadth"],
                     a["H2_detected_any"], a["H3_detected_only"])
    # template hard-codes the reversed H1 and the surviving H3 excess
    assert h1["verdict"] == "FALSIFIED" and h1["median_gated"] > h1["median_background"]
    assert h3["gated_excess_survives"]
    tab = h2["table_gated_bg_detected_not"]
    v = ("Mass spectrometry does not confirm the RNA-level restriction: gated "
         "proteins are seen in more normal human tissues than background, and "
         "the excess survives the detectability control. Protein-level "
         "off-tumour exposure therefore needs per-antigen review rather than "
         "being assumed from RNA windows.")
    blk = ["%s %d tissues" % (s, a["antigens"][s]["n_normal_tissues"])
           for s in ANTS if a["antigens"].get(s, {}).get("n_normal_tissues") is not None]
    rep = {"@NWL@": a["n_normal_tissue_whitelist"], "@NMAP@": g["G1_mapped"]["obs"],
           "@NGENES@": a["n_genes_total"], "@CDB@": g["G2_CD19_blineage"]["obs"],
           "@PCMIN@": min(g["G3_posctrl_normal_tissues"]["obs"].values()),
           "@H1G@": "%g" % h1["median_gated"], "@H1B@": "%g" % h1["median_background"],
           "@H1P@": fmtp(h1["p"]), "@DG@": tab[0][0], "@DGN@": sum(tab[0]),
           "@DB@": tab[1][0], "@DBN@": sum(tab[1]), "@H2P@": fmtp(h2["fisher_p"]),
           "@H3G@": "%g" % h3["median_gated"], "@H3B@": "%g" % h3["median_background"],
           "@H3GN@": h3["n_gated"], "@H3BN@": h3["n_background"], "@H3P@": fmtp(h3["p"]),
           "@VERDICT@": v, "@ANTIGENBLOCK@": "normal tissues with detection: "
           + "; ".join(blk) + "."}
    tpl = open("paper/proteomicsdb_tpl.tex").read()
    for k, vv in rep.items():
        tpl = tpl.replace(k, str(vv))
    assert "@" not in tpl
    with open("paper/proteomicsdb_sec.tex", "w") as f:
        f.write(tpl)
    print("wrote paper/proteomicsdb_sec.tex (%d chars)" % len(tpl))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate paper/harmonizome_sec.tex from paper/harmonizome_tpl.tex + committed
results/harmonizome_audit.json. Token replacement only."""
import json

ANTS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]


def fmtp(p):
    """p-value for LaTeX math mode: 0.017 or 6.3\\times10^{-10}."""
    s = "%.2g" % p
    if "e" in s:
        m, e = s.split("e")
        return "%s\\times10^{%d}" % (m, int(e))
    return s


def main():
    a = json.load(open("results/harmonizome_audit.json"))
    t, h2 = a["tests"], a["H2_control"]
    # template text hard-codes both falsifications; refuse to render otherwise
    assert a["verdict"].startswith("FALSIFIED") and not a["H1b"]
    assert h2["verdict"] == "FALSIFIED"
    ratio_ex = h2["expression"]["median_gated"] / h2["expression"]["median_background"]
    ratio_ne = h2["nonexpression"]["median_gated"] / h2["nonexpression"]["median_background"]
    v = ("The gated panel is more densely annotated than background on every "
         "Harmonizome axis, curated and systematic alike. Expression data carry "
         "most of the excess (median ratio %.2f versus %.2f for non-expression "
         "data), which fits selection for tumour-restricted expression, but a "
         "smaller residual gap remains. The earlier curated-only reading of the "
         "study-bias pattern does not generalise to this integrator."
         % (ratio_ex, ratio_ne))
    blk = []
    for s in ANTS:
        r = a["antigens"].get(s)
        if r and r.get("n_associations"):
            blk.append("%s %d associations (curated share %.2f)"
                       % (s, r["n_associations"], r["curated_frac"]))
    g = a["calibration_gates"]
    rep = {
        "@NGENES@": 205, "@NRES@": g["G1_resolved"]["obs"],
        "@NCUR@": a["dataset_classes"]["curated"],
        "@NSYS@": a["dataset_classes"]["systematic"],
        "@NUNC@": a["dataset_classes"]["unclassified"],
        "@CDT@": "{:,}".format(g["G2_CD19_total"]["obs"]).replace(",", "{,}"),
        "@PCMIN@": "{:,}".format(min(g["G3_posctrl_resolved_ge1000"]["obs"].values())).replace(",", "{,}"),
        "@CURG@": "%g" % t["curated"]["median_gated"],
        "@CURB@": "%g" % t["curated"]["median_background"],
        "@CURP@": fmtp(t["curated"]["p"]),
        "@SYSG@": "%g" % t["systematic"]["median_gated"],
        "@SYSB@": "%g" % t["systematic"]["median_background"],
        "@SYSP@": fmtp(t["systematic"]["p"]),
        "@FRG@": "%.2f" % t["curated_frac"]["median_gated"],
        "@FRB@": "%.2f" % t["curated_frac"]["median_background"],
        "@FRP@": fmtp(t["curated_frac"]["p"]),
        "@NNE@": h2["n_nonexpression_datasets"],
        "@EXG@": "%g" % h2["expression"]["median_gated"],
        "@EXB@": "%g" % h2["expression"]["median_background"],
        "@EXP@": fmtp(h2["expression"]["p"]),
        "@NEG@": "%g" % h2["nonexpression"]["median_gated"],
        "@NEB@": "%g" % h2["nonexpression"]["median_background"],
        "@NEP@": fmtp(h2["nonexpression"]["p"]),
        "@VERDICT@": v, "@ANTIGENBLOCK@": "; ".join(blk) + ".",
    }
    tpl = open("paper/harmonizome_tpl.tex").read()
    for k, vv in rep.items():
        tpl = tpl.replace(k, str(vv))
    assert "@" not in tpl
    with open("paper/harmonizome_sec.tex", "w") as f:
        f.write(tpl)
    print("wrote paper/harmonizome_sec.tex (%d chars)" % len(tpl))


if __name__ == "__main__":
    main()

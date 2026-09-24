#!/usr/bin/env python3
"""Glycan-shielding audit: extracellular-relevant glycosylation burden per gene from GlyGen detail records. Shield sites = reported N-linked + reported O-linked EXCLUDING subtype O-GlcNAcylation (intracellular single-sugar modification, not a steric shield). O-GlcNAc tracked separately; its presence in the intracellular controls calibrates the exclusion. Gene lengths from UniProt caches (epitope audit). Gated vs background compared within the annotated-surface stratum; AND-gate antigens vs approved references tabulated. Reads data/glygen/gg_*.json + data/epitope/uniprot_*.json + results/epitope_per_gene.csv. Writes results/glygen_audit.json + results/glygen_per_gene.csv."""
import csv, json, os
from scipy.stats import fisher_exact, mannwhitneyu

RAW = "data/glygen"
SETS = json.load(open("results/depmap_gene_sets.json"))
GROUPS = {g: grp for grp in ("gated", "background", "references") for g in SETS[grp]}
for g in SETS["positive_controls"]:
    GROUPS.setdefault(g, "positive_controls")
AND_GATE = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]
SURFACE = ("single-pass", "multi-pass", "GPI")

TOPO = {r["gene"]: r["topology"] for r in csv.DictReader(open("results/epitope_per_gene.csv"))}


def uniprot_len(sym):
    p = f"data/epitope/uniprot_{sym}.json"
    if not os.path.exists(p):
        p = f"data/glygen/uniprot_{sym}.json"
    d = json.load(open(p)) if os.path.exists(p) else None
    res = (d or {}).get("results", [])
    for r in res:
        if (r.get("genes") or [{}])[0].get("geneName", {}).get("value") == sym:
            return (r.get("sequence") or {}).get("length")
    return (res[0].get("sequence") or {}).get("length") if res else None


def load_gene(sym):
    p = os.path.join(RAW, f"gg_{sym}.json")
    d = json.load(open(p)) if os.path.exists(p) else None
    present = isinstance(d, dict) and "glycosylation" in d
    shield, oglcnac, pred, strucs = set(), set(), set(), set()
    if present:
        for g in (d.get("glycosylation") or []):
            pos = g.get("start_pos")
            if pos is None:
                continue
            cat = g.get("site_category") or ""
            sub = g.get("subtype") or ""
            key = (g.get("type"), pos)
            if sub == "O-GlcNAcylation":
                if cat.startswith("reported"):
                    oglcnac.add(key)
            elif cat.startswith("reported"):
                shield.add(key)
                if g.get("glytoucan_ac"):
                    strucs.add(g["glytoucan_ac"])
            elif cat == "predicted":
                pred.add(key)
    L = uniprot_len(sym)
    return {"gene": sym, "group": GROUPS[sym], "topology": TOPO.get(sym, "unannotated"),
            "record_present": present, "seq_len": L,
            "n_shield_sites": len(shield), "n_shield_structures": len(strucs),
            "n_oglcnac_sites": len(oglcnac), "n_predicted_sites": len(pred)}


def main():
    rows = [load_gene(s) for s in sorted(GROUPS)]
    for r in rows:
        r["shield_per_100aa"] = round(100.0 * r["n_shield_sites"] / r["seq_len"], 3) if r["seq_len"] else None
    rec = [r for r in rows if r["seq_len"]]
    surface = [r for r in rec if r["topology"] in SURFACE]
    g_s = [r for r in surface if r["group"] == "gated"]
    b_s = [r for r in surface if r["group"] == "background"]
    a = sum(r["n_shield_sites"] > 0 for r in g_s)
    b = sum(r["n_shield_sites"] > 0 for r in b_s)
    odds, p_fisher = fisher_exact([[a, len(g_s) - a], [b, len(b_s) - b]])
    gd = [r["shield_per_100aa"] for r in g_s]
    bd = [r["shield_per_100aa"] for r in b_s]
    U, p_mw = mannwhitneyu(gd, bd, alternative="two-sided")
    # topology-composition control: within-stratum Fisher tests
    strata = {}
    for topo in SURFACE:
        gt = [r for r in g_s if r["topology"] == topo]
        bt = [r for r in b_s if r["topology"] == topo]
        if not gt or not bt:
            continue
        ga = sum(r["n_shield_sites"] > 0 for r in gt)
        ba = sum(r["n_shield_sites"] > 0 for r in bt)
        o, p = fisher_exact([[ga, len(gt) - ga], [ba, len(bt) - ba]])
        strata[topo] = {"n_gated": len(gt), "gated_with": ga, "n_background": len(bt),
                        "background_with": ba, "fisher_p": p}
    # ectodomain-length-normalized density (epitope-audit ectodomains >= 20 aa)
    EP = {r["gene"]: r for r in csv.DictReader(open("results/epitope_per_gene.csv"))}
    egd, ebd = [], []
    for r in surface:
        e = EP.get(r["gene"])
        el = int(e["ecto_len"] or 0) if e else 0
        if el < 20:
            continue
        d = 100.0 * r["n_shield_sites"] / el
        (egd if r["group"] == "gated" else ebd).append(d)
    U2, p_ecto = mannwhitneyu(egd, ebd, alternative="two-sided")
    ecto = {"n_gated": len(egd), "n_background": len(ebd),
            "gated_median": sorted(egd)[len(egd) // 2], "background_median": sorted(ebd)[len(ebd) // 2],
            "mw_U": float(U2), "mw_p": p_ecto}
    pos = {s: next(r for r in rows if r["gene"] == s) for s in ("MSLN", "CD19", "ERBB2")}
    neg = {s: next(r for r in rows if r["gene"] == s) for s in ("POLR2A", "RPS3", "PCNA", "PSMA1")}
    gate_pos = sum(r["n_shield_sites"] >= 2 for r in pos.values())
    gate_og = sum(r["n_oglcnac_sites"] >= 1 for r in neg.values())
    neg_density = [r["shield_per_100aa"] for r in neg.values()]
    pos_density = [r["shield_per_100aa"] for r in pos.values()]
    gate_density = max(neg_density) < min(pos_density)
    audit = {
        "n_genes": len(rows),
        "n_record_present": sum(r["record_present"] for r in rows),
        "genes_without_record": sorted(r["gene"] for r in rows if not r["record_present"]),
        "surface_stratum": {"n_gated": len(g_s), "n_background": len(b_s),
                            "gated_with_shield": a, "background_with_shield": b,
                            "fisher_odds": round(odds, 3), "fisher_p": p_fisher,
                            "mw_U": float(U), "mw_p": p_mw,
                            "gated_median_density": sorted(gd)[len(gd) // 2],
                            "background_median_density": sorted(bd)[len(bd) // 2]},
        "topology_strata": strata,
        "ectodomain_density": ecto,
        "calibration": {"positive_with_2plus_shield": gate_pos, "of": 3,
                        "negative_with_oglcnac": gate_og, "neg_of": 4,
                        "neg_max_density": max(neg_density), "pos_min_density": min(pos_density),
                        "density_ordering": gate_density,
                        "gate_pass": gate_pos == 3 and gate_og == 4 and gate_density},
        "and_gate": [r for r in rows if r["gene"] in AND_GATE],
        "references": [r for r in rows if r["group"] == "references"],
        "intracellular_controls": list(neg.values()),
    }
    json.dump(audit, open("results/glygen_audit.json", "w"), indent=1)
    with open("results/glygen_per_gene.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gene", "group", "topology", "record_present", "seq_len",
                                           "n_shield_sites", "n_shield_structures", "n_oglcnac_sites",
                                           "n_predicted_sites", "shield_per_100aa"])
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["group"], r["gene"])))
    print(json.dumps({"n": audit["n_genes"], "records": audit["n_record_present"],
                      "gate": audit["calibration"]["gate_pass"],
                      "fisher_p": round(p_fisher, 5), "mw_p": p_mw,
                      "gated": f"{a}/{len(g_s)}", "bg": f"{b}/{len(b_s)}",
                      "cal": audit["calibration"]}, indent=1))


if __name__ == "__main__":
    main()

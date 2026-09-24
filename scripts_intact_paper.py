#!/usr/bin/env python3
"""Generate paper/intact_sec.tex from paper/intact_tpl.tex + committed
results/intact_audit.json, intact_edges.csv, intact_per_gene.csv.
Token replacement only; asserts no '@' remains."""
import csv, json

ANTIGENS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]
REFS = ["CD19", "ERBB2", "FOLR1", "TNFRSF17"]
CTRLS = ["POLR2A", "RPS3", "PCNA", "PSMA1"]

def esc(s):
    return str(s).replace("_", r"\_")

def main():
    a = json.load(open("results/intact_audit.json"))
    h1, h2, g = a["H1_pairwise_antigen_complex"], a["H2_degree"], a["gates"]
    edges = list(csv.DictReader(open("results/intact_edges.csv")))
    pg = {r["gene"]: r for r in csv.DictReader(open("results/intact_per_gene.csv"))}

    n_hits = len(h1["antigen_pairs_with_edge"])
    if n_hits == 0:
        pairlist = "none"
        h1v = ("No AND-gate antigen pair has a curated physical interaction: the gate's two "
               "arms are independent surface targets in the curated record, with no cis-complex "
               "to couple or short-circuit the two scFv engagements. Absence here is bounded by "
               "the study-bias audit below: low-degree antigens are simply under-probed.")
    elif h1["perm_degmatched_p"] is not None and h1["perm_degmatched_p"] < 0.05:
        h1v = ("Antigen pairs are significantly co-complexed beyond degree expectation, arguing "
               "that cis-engagement geometry must be considered in biclonal CAR design for these pairs.")
        pairlist = "; ".join("--".join(p) for p in h1["antigen_pairs_with_edge"])
    else:
        h1v = ("Curated antigen-pair edges do not exceed the degree-matched background rate, so "
               "there is no evidence that the gated antigens preferentially co-complex; observed "
               "edges are consistent with study intensity rather than biology of the gate.")
        pairlist = "; ".join("--".join(p) for p in h1["antigen_pairs_with_edge"])
    strat_parts = []
    for t in ("single-pass", "multi-pass", "GPI", "unannotated"):
        if t in h2["stratified"]:
            s = h2["stratified"][t]
            strat_parts.append(f"{t} $p={s['p']}$")
    sig_strata = [t for t in ("single-pass", "multi-pass", "GPI", "unannotated")
                  if t in h2["stratified"] and h2["stratified"][t]["p"] < 0.05]
    stratsum = ("non-significant in every topology stratum (" + ", ".join(strat_parts) + ")"
                if not sig_strata else
                "significant in " + ", ".join(sig_strata) + " (" + ", ".join(strat_parts) + ")")
    if h2["mw_p"] < 0.05:
        h2v = ("Gated and background genes differ in curated interaction degree, so the "
               "interaction record is study-biased across the design split; pairwise conclusions "
               "must be read against the degree-matched null, which controls this bias.")
    else:
        h2v = ("Gated and background genes do not differ in curated interaction degree, so the "
               "pairwise null is not confounded by differential study intensity across the split.")

    # per-antigen block: degree + partners within the audited focus sets
    focus = set(ANTIGENS + REFS + CTRLS)
    by_gene = {}
    for e in edges:
        by_gene.setdefault(e["gene_a"], []).append(e["gene_b"])
        by_gene.setdefault(e["gene_b"], []).append(e["gene_a"])
    blocks = []
    for ag in ANTIGENS:
        deg = pg.get(ag, {}).get("count", "0") or "0"
        fpart = sorted(set(by_gene.get(ag, [])) & focus)
        fp = ("; audited-set partners: " + ", ".join(esc(x) for x in fpart)) if fpart else ""
        blocks.append(f"\\textbf{{{esc(ag)}}} (curated degree {deg}{fp})")
    antigen_block = "; ".join(blocks) + "."

    g2ok = sum(1 for v in g["G2_count_mitab_consistency"]["per_gene"].values() if v["consistent"])
    repl = {
        "NGENES": a["n_genes_total"], "NOK": a["n_genes_with_accession"],
        "GATEONE": "present" if g["G1_erbb2_canonical_partners"]["pass"] else "MISSING",
        "GATETWO": f"{g2ok}/15", "GATETHREE": "4/4" if g["G3_controls_nonzero"]["pass"] else "FAILED",
        "NEDGES": n_hits, "PAIRLIST": pairlist, "PAIRFRAC": f"{h1['antigen_pair_fraction']:.3f}",
        "NPERM": "20{,}000", "PUNIF": h1["perm_uniform_p"],
        "PMATCH": h1["perm_degmatched_p"], "NMATCHED": f"{h1['perm_n_matched']:,}".replace(",", "{,}"),
        "H1VERDICT": h1v, "MEDG": h2["median_gated"], "MEDB": h2["median_bg"], "PMW": h2["mw_p"],
        "COVG": h2["coverage_gated"], "COVB": h2["coverage_bg"],
        "COVOR": h2["coverage_fisher_or"], "COVP": h2["coverage_fisher_p"],
        "STRATSUM": stratsum, "H2VERDICT": h2v, "ANTIGENBLOCK": antigen_block,
    }
    tex = open("paper/intact_tpl.tex").read()
    for k, v in repl.items():
        tex = tex.replace(f"@{k}@", str(v))
    assert "@" not in tex, "unreplaced token remains"
    open("paper/intact_sec.tex", "w").write(tex)
    print("wrote paper/intact_sec.tex")

if __name__ == "__main__":
    main()

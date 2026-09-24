#!/usr/bin/env python3
"""Antigen-escape constraint audit (item 17).
Question: can tumours lose gated CAR-T antigens cheaply? Three orthogonal axes per gene:
 (1) germline LoF intolerance, gnomAD v4 LOEUF (oe_lof_upper; constrained if < 0.6);
 (2) paralog redundancy, HGNC gene-group size;
 (3) pathway embedding, number of Reactome pathways (lowest-level mappings).
Calibration gate: essential controls (POLR2A, RPS3, PCNA, PSMA1) must read LoF-constrained (>=3/4) before any claim.
Escape-risk quadrant combines LOEUF with committed DepMap CRISPR dependency (results/depmap_dependency.json).
Writes results/constraint_audit.json + results/constraint_per_gene.csv."""
import csv, json, os, time, urllib.request
import numpy as np
from scipy.stats import mannwhitneyu, fisher_exact, spearmanr

RAW = "data/constraint"
AND_GATE = ["CLDN18", "MSLN", "CA9", "CA12", "CLDN6", "PSCA", "SLC39A6"]
LOEUF_CUT = 0.6


def load(p):
    p = os.path.join(RAW, p)
    return json.load(open(p)) if os.path.exists(p) else None


def group_size(gid):
    p = os.path.join(RAW, f"hgncgroup_{gid}.json")
    if not os.path.exists(p):
        url = f"https://rest.genenames.org/search/gene_group_id:{gid}"
        for a in range(4):
            try:
                req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "Mozilla/5.0"})
                with urllib.request.urlopen(req, timeout=60) as r:
                    open(p, "wb").write(r.read())
                break
            except Exception:
                time.sleep(3 * (a + 1))
    if not os.path.exists(p):
        return None
    return json.load(open(p))["response"]["numFound"]


def gene_record(g, grp, dep):
    rec = {"gene": g, "group": grp}
    gn = load(f"gnomad_{g}.json") or {}
    gene = ((gn.get("data") or {}).get("gene")) or {}
    c = gene.get("gnomad_constraint") or {}
    rec["ensembl"] = gene.get("gene_id")
    rec["loeuf"] = c.get("oe_lof_upper")
    rec["pli"] = c.get("pli")
    rec["mis_z"] = c.get("mis_z")
    rec["obs_lof"] = c.get("obs_lof")
    rec["exp_lof"] = c.get("exp_lof")
    h = load(f"hgnc_{g}.json") or {}
    docs = ((h.get("response") or {}).get("docs")) or []
    d0 = docs[0] if docs else {}
    rec["hgnc_id"] = d0.get("hgnc_id")
    rec["locus_type"] = d0.get("locus_type")
    gids = d0.get("gene_group_id") or []
    rec["gene_groups"] = ";".join(d0.get("gene_group") or [])
    sizes = [group_size(x) for x in gids]
    sizes = [s for s in sizes if s]
    rec["max_group_size"] = max(sizes) if sizes else 0
    meta = load(f"meta_{g}.json") or {}
    rec["accession"] = meta.get("accession")
    r = load(f"reactome_{g}.json")
    rec["n_reactome"] = len(r) if isinstance(r, list) else 0
    dd = dep.get(g)
    rec["depmap_median"] = dd["median"] if dd else None
    rec["depmap_frac_dep"] = dd["frac_dep"] if dd else None
    return rec


def mw(a, b):
    a = [x for x in a if x is not None]; b = [x for x in b if x is not None]
    if len(a) < 3 or len(b) < 3:
        return None
    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)


def main():
    s = json.load(open("results/depmap_gene_sets.json"))
    depj = json.load(open("results/depmap_dependency.json"))
    dep = {}
    for grp in depj["per_gene"].values():
        dep.update(grp)
    dep.update(depj["references"]); dep.update(depj["positive_controls"])
    genes = {}
    for grp in ("gated", "background", "references", "positive_controls"):
        for g in s[grp]:
            genes.setdefault(g, grp)
    rows = [gene_record(g, genes[g], dep) for g in sorted(genes)]
    by = {k: [r for r in rows if r["group"] == k] for k in ("gated", "background", "references", "positive_controls")}
    out = {"design": __doc__.strip().splitlines()[0], "loeuf_cut": LOEUF_CUT,
           "sources": {"gnomad": "gnomAD v4 GraphQL gene.gnomad_constraint (GRCh38)",
                       "hgnc": "HGNC REST fetch/symbol + search/gene_group_id",
                       "reactome": "Reactome ContentService data/mapping/UniProt/<acc>/pathways?species=9606"}}
    ctrl = by["positive_controls"]
    ok = [r["gene"] for r in ctrl if r["loeuf"] is not None and r["loeuf"] < LOEUF_CUT]
    out["calibration"] = {"controls": {r["gene"]: r["loeuf"] for r in ctrl}, "n_constrained": len(ok),
                          "passed": len(ok) >= 3}
    G, B = by["gated"], by["background"]
    def col(rs, k):
        return [r[k] for r in rs if r[k] is not None]
    gl, bl = col(G, "loeuf"), col(B, "loeuf")
    gc = sum(x < LOEUF_CUT for x in gl); bc = sum(x < LOEUF_CUT for x in bl)
    out["loeuf"] = {"n_gated": len(gl), "n_background": len(bl),
                    "gated_median": float(np.median(gl)), "background_median": float(np.median(bl)),
                    "mw_p": mw(gl, bl), "gated_constrained": gc, "background_constrained": bc,
                    "fisher_p": float(fisher_exact([[gc, len(gl) - gc], [bc, len(bl) - bc]])[1])}
    gg, bg = col(G, "max_group_size"), col(B, "max_group_size")
    gin = sum(x > 1 for x in gg); bin_ = sum(x > 1 for x in bg)
    out["paralog_groups"] = {"gated_median_group_size": float(np.median(gg)), "background_median_group_size": float(np.median(bg)),
                             "mw_p": mw(gg, bg), "gated_in_multigene_group": gin, "background_in_multigene_group": bin_,
                             "n_gated": len(gg), "n_background": len(bg),
                             "fisher_p": float(fisher_exact([[gin, len(gg) - gin], [bin_, len(bg) - bin_]])[1])}
    gr, br = col(G, "n_reactome"), col(B, "n_reactome")
    gz = sum(x == 0 for x in gr); bz = sum(x == 0 for x in br)
    out["reactome"] = {"gated_median": float(np.median(gr)), "background_median": float(np.median(br)), "mw_p": mw(gr, br),
                       "gated_pathway_dark": gz, "background_pathway_dark": bz, "n_gated": len(gr), "n_background": len(br),
                       "fisher_p": float(fisher_exact([[gz, len(gr) - gz], [bz, len(br) - bz]])[1])}
    both = [r for r in G + B if r["loeuf"] is not None and r["depmap_median"] is not None]
    rho, p = spearmanr([r["loeuf"] for r in both], [r["depmap_median"] for r in both])
    out["loeuf_vs_depmap"] = {"n": len(both), "spearman_rho": float(rho), "p": float(p)}
    # escape-risk quadrant: tolerant germline AND non-essential in cancer lines (frac_dep < 0.1)
    def quad(r):
        if r["loeuf"] is None or r["depmap_frac_dep"] is None:
            return "unknown"
        t = r["loeuf"] >= LOEUF_CUT; ne = r["depmap_frac_dep"] < 0.1
        return {(True, True): "free_loss", (True, False): "tumour_dependent", (False, True): "germline_constrained",
                (False, False): "doubly_constrained"}[(t, ne)]
    for r in rows:
        r["escape_quadrant"] = quad(r)
    def qc(rs):
        c = {}
        for r in rs:
            c[r["escape_quadrant"]] = c.get(r["escape_quadrant"], 0) + 1
        return c
    out["quadrants"] = {"gated": qc(G), "background": qc(B)}
    kg = [r for r in G if r["escape_quadrant"] != "unknown"]; kb = [r for r in B if r["escape_quadrant"] != "unknown"]
    fg = sum(r["escape_quadrant"] == "free_loss" for r in kg); fb = sum(r["escape_quadrant"] == "free_loss" for r in kb)
    out["free_loss"] = {"gated": fg, "gated_n": len(kg), "background": fb, "background_n": len(kb),
                        "fisher_p": float(fisher_exact([[fg, len(kg) - fg], [fb, len(kb) - fb]])[1])}
    # confound control: restrict to protein-coding loci outside olfactory-receptor / Ig / TCR groups
    def clean(rs):
        return [r for r in rs if r["locus_type"] == "gene with protein product"
                and not any(k in (r["gene_groups"] or "") for k in ("Olfactory receptor", "Immunoglobulin kappa", "Immunoglobulin lambda", "T cell receptor"))]
    cG, cB = clean(G), clean(B)
    cgr, cbr = col(cG, "n_reactome"), col(cB, "n_reactome")
    cgz = sum(x == 0 for x in cgr); cbz = sum(x == 0 for x in cbr)
    cgl, cbl = col(cG, "loeuf"), col(cB, "loeuf")
    out["protein_coding_control"] = {"n_gated": len(cG), "n_background": len(cB),
        "reactome_gated_median": float(np.median(cgr)), "reactome_background_median": float(np.median(cbr)),
        "reactome_mw_p": mw(cgr, cbr), "gated_pathway_dark": cgz, "background_pathway_dark": cbz,
        "reactome_dark_fisher_p": float(fisher_exact([[cgz, len(cgr) - cgz], [cbz, len(cbr) - cbz]])[1]),
        "loeuf_gated_median": float(np.median(cgl)), "loeuf_background_median": float(np.median(cbl)),
        "loeuf_mw_p": mw(cgl, cbl), "n_loeuf_gated": len(cgl), "n_loeuf_background": len(cbl)}
    ag = [r for r in rows if r["gene"] in AND_GATE]
    out["and_gate_summary"] = {"n": len(ag), "free_loss": sum(r["escape_quadrant"] == "free_loss" for r in ag),
        "unknown": [r["gene"] for r in ag if r["escape_quadrant"] == "unknown"],
        "min_loeuf": min(r["loeuf"] for r in ag if r["loeuf"] is not None),
        "max_frac_dep": max(r["depmap_frac_dep"] for r in ag if r["depmap_frac_dep"] is not None)}
    refs = [r for r in by["references"]]
    out["reference_summary"] = {r["gene"]: r["escape_quadrant"] for r in refs}
    keys = ["gene", "loeuf", "pli", "mis_z", "obs_lof", "exp_lof", "gene_groups", "max_group_size", "n_reactome",
            "depmap_median", "depmap_frac_dep", "escape_quadrant"]
    out["references"] = [{k: r[k] for k in keys} for r in by["references"]]
    out["and_gate"] = [{k: r[k] for k in keys} for r in rows if r["gene"] in AND_GATE]
    out["n_genes"] = len(rows)
    out["n_with_loeuf"] = sum(r["loeuf"] is not None for r in rows)
    out["n_with_hgnc"] = sum(r["hgnc_id"] is not None for r in rows)
    out["n_with_reactome_mapping"] = sum(os.path.exists(os.path.join(RAW, f"reactome_{r['gene']}.json")) for r in rows)
    out["n_hgnc_groups_fetched"] = len([f for f in os.listdir(RAW) if f.startswith("hgncgroup_")])
    import glob
    def _ok(f, test):
        try:
            return test(json.load(open(f)))
        except Exception:
            return False
    out["dataset_records"] = {
        "gnomad_gene_constraint": sum(_ok(f, lambda d: bool(((d or {}).get("data") or {}).get("gene"))) for f in glob.glob(RAW + "/gnomad_*.json")),
        "hgnc_symbol": sum(_ok(f, lambda d: ((d or {}).get("response") or {}).get("numFound", 0) > 0) for f in glob.glob(RAW + "/hgnc_*.json")),
        "hgnc_gene_group": len(glob.glob(RAW + "/hgncgroup_*.json")),
        "reactome_mapping": sum(_ok(f, lambda d: isinstance(d, list)) for f in glob.glob(RAW + "/reactome_*.json"))}
    out["dataset_records"]["total"] = sum(out["dataset_records"].values())
    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/constraint_audit.json", "w"), indent=1)
    with open("results/constraint_per_gene.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(json.dumps({k: out[k] for k in ("calibration", "loeuf", "paralog_groups", "reactome", "loeuf_vs_depmap",
                                           "quadrants", "free_loss", "n_with_loeuf", "protein_coding_control", "and_gate_summary", "reference_summary")}, indent=1))
    for r in out["and_gate"] + out["references"]:
        print(r)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Pharos (TCRD/IDG) target-illumination audit for the 205-gene study set.
Question: what is the illumination state (Target Development Level, TIN-X
novelty, approved-drug and ligand counts) of gated antigens vs the random
surfaceome background? And does Pharos Tclin (current) agree with the
DrugCentral 2021_09_01 Tclin snapshot (cross-source calibration)?
Pre-registered gates:
  G1: all 4 references (CD19, ERBB2, FOLR1, TNFRSF17) are Tclin in Pharos.
  G2: >= 200/205 gene-set records resolve to a Pharos target record.
  G3: >= 3/4 essential positive controls (POLR2A, RPS3, PCNA, PSMA1) are
      NOT Tclin (negative-direction sanity: housekeeping is undrugged).
  G4: EGFR sentinel is Tclin in BOTH Pharos and DrugCentral (anchor).
Analyses:
  H1: Tclin fraction gated vs background (Fisher exact); full TDL mix.
  H2: TIN-X novelty gated vs background (Mann-Whitney).
  H3: approved-drug and ligand counts gated vs background (MW); AND-gate
      antigen table.
  H4: Pharos-Tclin vs DrugCentral-Tclin 2x2 concordance over the 197
      universe genes + 4 references; discordant genes named."""
import csv, json, os
from scipy.stats import mannwhitneyu, fisher_exact

DATA = "data/pharos"
ANTIGENS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]
REFS = ["CD19", "ERBB2", "FOLR1", "TNFRSF17"]
CTRLS = ["POLR2A", "RPS3", "PCNA", "PSMA1"]


def gene_set():
    gs = json.load(open("results/depmap_gene_sets.json"))
    genes = {}
    for k in ("gated", "background", "references", "positive_controls"):
        for g in gs[k]:
            genes[g] = k
    return genes


def load_record(sym):
    path = os.path.join(DATA, sym + ".json")
    if not os.path.exists(path):
        return None
    rec = json.load(open(path))
    if "_error" in rec:
        return None
    tgt = (rec.get("data") or {}).get("target")
    return tgt


def count_list(tgt, key, name):
    for e in (tgt.get(key) or []):
        if e.get("name") == name:
            return e.get("value") or 0
    return 0


def main():
    genes = gene_set()
    dc = json.load(open("results/drugcentral_engagement.json"))["per_gene"]
    rows = {}
    n_ok = 0
    for sym, grp in genes.items():
        tgt = load_record(sym)
        if tgt is None:
            rows[sym] = {"group": grp, "status": "no_record"}
            continue
        n_ok += 1
        rows[sym] = {
            "group": grp, "status": "ok",
            "tdl": tgt.get("tdl"),
            "novelty": tgt.get("novelty"),
            "n_drugs": count_list(tgt, "ligandCounts", "drug"),
            "n_ligands": count_list(tgt, "ligandCounts", "ligand"),
            "publication_count": tgt.get("publicationCount") or 0,
            "generif_count": tgt.get("generifCount") or 0,
            "gwas_total": sum((e.get("value") or 0) for e in (tgt.get("gwasCounts") or [])),
            "ppi_total": sum((e.get("value") or 0) for e in (tgt.get("ppiCounts") or [])),
        }
    # EGFR sentinel
    egfr = load_record("EGFR")
    egfr_tdl = egfr.get("tdl") if egfr else None

    # gates
    g1 = all(rows.get(r, {}).get("tdl") == "Tclin" for r in REFS)
    g1_detail = {r: rows.get(r, {}).get("tdl") for r in REFS}
    g2 = n_ok >= 200
    g3_detail = {c: rows.get(c, {}).get("tdl") for c in CTRLS}
    g3 = sum(1 for c in CTRLS if rows.get(c, {}).get("tdl") not in (None, "Tclin")) >= 3
    g4 = (egfr_tdl == "Tclin") and bool(dc.get("EGFR", {}).get("tclin"))
    gates = {
        "G1_references_tclin": {"pass": g1, "detail": g1_detail},
        "G2_coverage_ge_200_of_205": {"pass": g2, "n_ok": n_ok},
        "G3_controls_not_tclin": {"pass": g3, "detail": g3_detail},
        "G4_egfr_tclin_both_sources": {"pass": g4, "pharos": egfr_tdl,
                                       "drugcentral": bool(dc.get("EGFR", {}).get("tclin"))},
    }

    def grp_rows(g):
        return {s: r for s, r in rows.items() if r.get("status") == "ok" and r["group"] == g}

    gated, bg = grp_rows("gated"), grp_rows("background")

    # H1 Tclin fraction
    a = sum(1 for r in gated.values() if r["tdl"] == "Tclin")
    b = len(gated) - a
    c = sum(1 for r in bg.values() if r["tdl"] == "Tclin")
    d = len(bg) - c
    or_, p_tclin = fisher_exact([[a, b], [c, d]])
    from collections import Counter
    tdl_mix_gated = Counter(r["tdl"] for r in gated.values())
    tdl_mix_bg = Counter(r["tdl"] for r in bg.values())

    # H2 novelty
    ng = [r["novelty"] for r in gated.values() if r["novelty"] is not None]
    nb = [r["novelty"] for r in bg.values() if r["novelty"] is not None]
    u2, p_nov = mannwhitneyu(ng, nb)
    import statistics as st
    # H3 drug / ligand counts
    dg = [r["n_drugs"] for r in gated.values()]
    db = [r["n_drugs"] for r in bg.values()]
    _, p_drug = mannwhitneyu(dg, db)
    lg = [r["n_ligands"] for r in gated.values()]
    lb = [r["n_ligands"] for r in bg.values()]
    _, p_lig = mannwhitneyu(lg, lb)

    antigen_table = {s: rows.get(s, {"status": "no_record"}) for s in ANTIGENS}

    # H4 concordance over universe + references
    both = []
    for sym, r in rows.items():
        if r.get("status") != "ok" or sym not in dc:
            continue
        both.append((sym, r["tdl"] == "Tclin", bool(dc[sym].get("tclin"))))
    tt = sum(1 for _, p, d in both if p and d)
    tf = sum(1 for _, p, d in both if p and not d)
    ft = sum(1 for _, p, d in both if not p and d)
    ff = sum(1 for _, p, d in both if not p and not d)
    n_b = len(both)
    agree = (tt + ff) / n_b if n_b else None
    pe = ((tt + tf) * (tt + ft) + (ft + ff) * (tf + ff)) / (n_b * n_b) if n_b else None
    kappa = (agree - pe) / (1 - pe) if pe is not None and pe < 1 else None
    discordant = sorted(s for s, p, d in both if p != d)

    out = {
        "tool": "Pharos GraphQL API (pharos-api.ncats.io; TCRD/IDG illumination)",
        "n_genes_total": len(genes), "n_records_ok": n_ok,
        "gates": gates, "all_gates_pass": all(g["pass"] for g in gates.values()),
        "H1_tclin": {"gated_tclin": a, "gated_n": len(gated),
                     "bg_tclin": c, "bg_n": len(bg),
                     "fisher_or": or_, "fisher_p": p_tclin,
                     "tdl_mix_gated": dict(tdl_mix_gated),
                     "tdl_mix_background": dict(tdl_mix_bg)},
        "H2_novelty": {"gated_median": st.median(ng), "bg_median": st.median(nb),
                       "mw_p": p_nov},
        "H3_counts": {"gated_drugs_median": st.median(dg), "bg_drugs_median": st.median(db),
                      "drugs_mw_p": p_drug,
                      "gated_ligands_median": st.median(lg), "bg_ligands_median": st.median(lb),
                      "ligands_mw_p": p_lig},
        "H4_concordance": {"n": n_b, "both_tclin": tt, "pharos_only": tf,
                           "drugcentral_only": ft, "neither": ff,
                           "agreement": agree, "cohens_kappa": kappa,
                           "discordant_genes": discordant},
        "and_gate_antigens": antigen_table,
        "egfr_sentinel": {"pharos_tdl": egfr_tdl},
    }
    with open("results/pharos_audit.json", "w") as f:
        json.dump(out, f, indent=1)

    with open("results/pharos_per_gene.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "group", "status", "tdl", "novelty", "n_drugs",
                    "n_ligands", "publication_count", "generif_count",
                    "gwas_total", "ppi_total"])
        for sym in sorted(rows):
            r = rows[sym]
            if r.get("status") != "ok":
                w.writerow([sym, r["group"], r.get("status"), "", "", "", "", "", "", "", ""])
            else:
                w.writerow([sym, r["group"], "ok", r["tdl"], r["novelty"], r["n_drugs"],
                            r["n_ligands"], r["publication_count"], r["generif_count"],
                            r["gwas_total"], r["ppi_total"]])

    print(json.dumps({"gates": {k: v["pass"] for k, v in gates.items()},
                      "n_ok": n_ok,
                      "H1": out["H1_tclin"],
                      "H2": out["H2_novelty"],
                      "H3": out["H3_counts"],
                      "H4": {k: v for k, v in out["H4_concordance"].items() if k != "discordant_genes"},
                      "discordant": discordant}, indent=1))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""IntAct curated-interaction audit for the 205-gene DepMap universe.
Question: do AND-gate antigen pairs physically interact (cis-complex),
which would affect biclonal CAR avidity / epitope access?
Pre-registered gates:
  G1: ERBB2 partners include EGFR and ERBB3 (canonical heterodimers).
  G2: for the 15 focus genes, paginated MITAB rows == PSICQUIC count.
  G3: all 4 intracellular positive controls have count > 0.
Analyses:
  H1 pairwise: fraction of 21 antigen-antigen pairs with a human physical
     curated edge; permutation null over background 7-subsets (uniform and
     degree-matched).
  H2 degree: curated interaction count gated vs background (Mann-Whitney),
     stratified by epitope topology stratum; any-coverage Fisher.
Physical evidence types (PSI-MI): MI:0407 direct interaction,
  MI:0914 association, MI:0915 physical association, MI:2364 proximity.
Both interactors required taxid:9606. Self edges excluded from pairwise."""
import csv, glob, json, os, random, re
from collections import Counter, defaultdict
from scipy.stats import mannwhitneyu, fisher_exact

PHYS = {"MI:0407", "MI:0914", "MI:0915", "MI:2364"}
MD = "data/intact/mitab"
ANTIGENS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]
REFS = ["CD19", "ERBB2", "FOLR1", "TNFRSF17"]
CTRLS = ["POLR2A", "RPS3", "PCNA", "PSMA1"]

def base_acc(tok):
    m = re.search(r"uniprotkb:([A-Z0-9]+(?:-\d+)?)", tok)
    if not m:
        return None
    return re.sub(r"-\d+$", "", m.group(1))

def load_inputs():
    gs = json.load(open("results/depmap_gene_sets.json"))
    genes = {}
    for k in ("gated", "background", "references", "positive_controls"):
        for g in gs[k]:
            genes[g] = k
    man = {r["symbol"]: r for r in json.load(open("data/intact/intact_manifest.json"))}
    topo = {}
    for r in csv.DictReader(open("results/epitope_per_gene.csv")):
        topo[r["gene"]] = r["topology"]
    return genes, man, topo

def parse_mitab(sym, self_acc):
    """Return list of partner dicts for one gene's MITAB cache."""
    f = f"{MD}/{sym}.mitab"
    out = []
    if not os.path.exists(f):
        return out
    for line in open(f):
        c = line.rstrip("\n").split("\t")
        if len(c) < 15:
            continue
        a, b = base_acc(c[0]), base_acc(c[1])
        if not a or not b:
            continue  # non-protein interactor (e.g. intact-only id)
        if "taxid:9606" not in c[9] or "taxid:9606" not in c[10]:
            continue
        typ = re.search(r'MI:\d{4}', c[11])
        typ = typ.group(0) if typ else "unknown"
        det = re.search(r'MI:\d{4}', c[6])
        det = det.group(0) if det else "unknown"
        pm = re.findall(r"pubmed:(\d+)", c[8])
        partner = b if a == self_acc else a
        out.append({"partner": partner, "type": typ, "det": det, "pubmeds": pm})
    return out

def main():
    genes, man, topo = load_inputs()
    sym2acc = {s: r["accession"] for s, r in man.items() if r.get("accession")}
    acc2sym = {a: s for s, a in sym2acc.items()}
    # ---- per-gene edges within universe (physical only) ----
    per_gene_edges = {}
    type_counter = Counter()
    for sym, acc in sym2acc.items():
        rows = parse_mitab(sym, acc)
        for r in rows:
            type_counter[r["type"]] += 1
        phys = [r for r in rows if r["type"] in PHYS]
        keep = {}
        for r in phys:
            p = r["partner"]
            if p == acc:
                continue
            ent = keep.setdefault(p, {"types": set(), "det": set(), "pubmeds": set()})
            ent["types"].add(r["type"]); ent["det"].add(r["det"]); ent["pubmeds"].update(r["pubmeds"])
        per_gene_edges[sym] = keep
    # ---- union adjacency between universe genes ----
    edges = {}
    for s1, k in per_gene_edges.items():
        for pacc, ev in k.items():
            s2 = acc2sym.get(pacc)
            if not s2 or s2 == s1:
                continue
            key = tuple(sorted((s1, s2)))
            e = edges.setdefault(key, {"types": set(), "det": set(), "pubmeds": set()})
            e["types"] |= ev["types"]; e["det"] |= ev["det"]; e["pubmeds"] |= ev["pubmeds"]
    # ---- gates ----
    erbb2_partners = {acc2sym.get(p, p) for p in per_gene_edges.get("ERBB2", {})}
    # partners hold symbols where resolvable into the universe, raw accessions otherwise
    g1_egfr = bool({"EGFR", "P00533"} & erbb2_partners)
    g1_erbb3 = bool({"ERBB3", "P21860"} & erbb2_partners)
    g1 = g1_egfr and g1_erbb3
    g2_rows = {}
    focus = ANTIGENS + REFS + CTRLS
    g2 = True
    for s in focus:
        r = man.get(s, {})
        exp, got = r.get("count"), r.get("mitab_rows")
        ok = (exp is not None and got == exp)
        g2 &= ok
        g2_rows[s] = {"count": exp, "mitab_rows": got, "consistent": ok}
    g3_rows = {s: man.get(s, {}).get("count") for s in CTRLS}
    g3 = all((v or 0) > 0 for v in g3_rows.values())
    gates = {
        "G1_erbb2_canonical_partners": {"pass": g1, "egfr": g1_egfr, "erbb3": g1_erbb3,
            "detail": "ERBB2 partners include EGFR+ERBB3" if g1 else f"MISSING: egfr={g1_egfr} erbb3={g1_erbb3}"},
        "G2_count_mitab_consistency": {"pass": g2, "per_gene": g2_rows},
        "G3_controls_nonzero": {"pass": g3, "counts": g3_rows},
    }
    # ---- H1: pairwise antigen edges ----
    def pair_frac(subset, ed):
        ks = [tuple(sorted(x)) for i, x in enumerate(subset) for y in subset[i+1:]]
        if not ks:
            return 0.0, []
        hit = [k for k in ks if k in ed]
        return len(hit) / len(ks), sorted(hit)
    obs_frac, obs_hits = pair_frac(ANTIGENS, edges)
    rng = random.Random(20260925)
    bg = [g for g, grp in genes.items() if grp == "background" and g in sym2acc]
    cnt = {s: man[s].get("count") or 0 for s in sym2acc}
    N = 20000
    null_uniform, null_matched = 0, 0
    matched_sets = []
    ag_deg = sorted(cnt[a] for a in ANTIGENS)
    for _ in range(N):
        f1, _ = pair_frac(rng.sample(bg, 7), edges)
        if f1 >= obs_frac:
            null_uniform += 1
        # degree-matched: for each antigen degree pick random bg gene within +/-25%
        picks = []
        pool = set(bg)
        okm = True
        for d in ag_deg:
            cand = [g for g in pool if abs(cnt[g] - d) <= max(1, 0.25 * d)]
            if not cand:
                okm = False; break
            c = rng.choice(cand); picks.append(c); pool.discard(c)
        if okm:
            f2, _ = pair_frac(picks, edges)
            matched_sets.append(f2 >= obs_frac)
            if f2 >= obs_frac:
                null_matched += 1
    n_matched = len(matched_sets)
    p_unif = (1 + null_uniform) / (1 + N)
    p_match = (1 + null_matched) / (1 + n_matched) if n_matched else None
    h1 = {"antigen_pairs_total": 21, "antigen_pairs_with_edge": sorted(map(list, obs_hits)),
          "antigen_pair_fraction": round(obs_frac, 4),
          "edges_detail": {"|".join(k): {"types": sorted(edges[k]["types"]), "det": sorted(edges[k]["det"]),
                              "pubmeds": sorted(edges[k]["pubmeds"])} for k in obs_hits},
          "perm_uniform_p": round(p_unif, 5), "perm_degmatched_p": (round(p_match, 5) if p_match is not None else None),
          "perm_n_matched": n_matched}
    # ---- H2: degree ----
    def grp_deg(grp):
        return [cnt[g] for g, gr in genes.items() if gr == grp and g in sym2acc]
    gd, bd = grp_deg("gated"), grp_deg("background")
    U, pdeg = mannwhitneyu(gd, bd, alternative="two-sided")
    g_cov = sum(1 for x in gd if x > 0); b_cov = sum(1 for x in bd if x > 0)
    orc, pcov = fisher_exact([[g_cov, len(gd) - g_cov], [b_cov, len(bd) - b_cov]])
    strat = {}
    for t in ("single-pass", "multi-pass", "GPI", "unannotated"):
        xs = [cnt[g] for g, gr in genes.items() if gr == "gated" and topo.get(g) == t and g in sym2acc]
        ys = [cnt[g] for g, gr in genes.items() if gr == "background" and topo.get(g) == t and g in sym2acc]
        if xs and ys:
            _, pt = mannwhitneyu(xs, ys, alternative="two-sided")
            strat[t] = {"n_gated": len(xs), "n_bg": len(ys), "med_gated": sorted(xs)[len(xs)//2],
                        "med_bg": sorted(ys)[len(ys)//2], "p": round(float(pt), 4)}
    h2 = {"median_gated": sorted(gd)[len(gd)//2], "median_bg": sorted(bd)[len(bd)//2],
          "mean_gated": round(sum(gd)/len(gd), 2), "mean_bg": round(sum(bd)/len(bd), 2),
          "mw_p": round(float(pdeg), 6),
          "coverage_gated": f"{g_cov}/{len(gd)}", "coverage_bg": f"{b_cov}/{len(bd)}",
          "coverage_fisher_or": round(float(orc), 3), "coverage_fisher_p": round(float(pcov), 6),
          "stratified": strat}
    # ---- outputs ----
    with open("results/intact_per_gene.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene", "group", "accession", "count", "mitab_rows", "topology", "n_universe_partners_physical"])
        for s, grp in sorted(genes.items()):
            r = man.get(s, {})
            w.writerow([s, grp, r.get("accession"), r.get("count"), r.get("mitab_rows"),
                        topo.get(s, "control"), len(per_gene_edges.get(s, {}))])
    with open("results/intact_edges.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene_a", "gene_b", "types", "detection", "pubmeds"])
        for (a, b), ev in sorted(edges.items()):
            w.writerow([a, b, "|".join(sorted(ev["types"])), "|".join(sorted(ev["det"])),
                        "|".join(sorted(ev["pubmeds"]))])
    audit = {"tool": "IntAct (EBI PSICQUIC REST, mitab tab25 + count)",
             "n_genes_with_accession": len(sym2acc), "n_genes_total": len(genes),
             "physical_mi_types": sorted(PHYS),
             "type_distribution_raw_rows": dict(type_counter.most_common()),
             "gates": gates, "all_gates_pass": all(g["pass"] for g in gates.values()),
             "H1_pairwise_antigen_complex": h1, "H2_degree": h2,
             "n_universe_edges": len(edges)}
    json.dump(audit, open("results/intact_audit.json", "w"), indent=1)
    print(json.dumps({"gates_pass": audit["all_gates_pass"], "H1": {k: v for k, v in h1.items() if k != "edges_detail"},
                      "H2_medians": [h2["median_gated"], h2["median_bg"]], "H2_mw_p": h2["mw_p"],
                      "n_edges": len(edges)}, indent=1))

if __name__ == "__main__":
    main()

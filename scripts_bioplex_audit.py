#!/usr/bin/env python3
"""BioPlex systematic AP-MS co-complex audit for the 205-gene study set.
Question: do AND-gate antigen pairs co-precipitate in the two large systematic
AP-MS networks (BioPlex HEK293T 10K and HCT116 5.5K, Dec 2019 releases)? This
is the orthogonal experimental check on the IntAct curated 0/21 null: IntAct
could miss complexes by under-curation; BioPlex is systematic but only covers
proteins expressed/detected in those two cell lines, so only pairs with BOTH
antigens present in a network are testable there.
Pre-registered gates:
  G1: PSMA1 has >=1 proteasome beta-subunit (PSMB*) partner in BOTH networks.
  G2: POLR2A has >=1 Pol II subunit (POLR2x) partner in the 293T network.
  G3: >=3/4 CAR-T reference targets are present in the 293T network.
Analyses:
  H1: antigen-pair edges among testable pairs (both antigens present).
  H2: AP-MS degree gated vs background among present genes (MW), per network;
      presence (detection) gated vs background (Fisher), per network.
  H3: cross-network edge reproducibility (293T vs HCT116, Jaccard over
      universe-restricted edges) and IntAct concordance (BioPlex-either vs
      curated human physical edges over universe pairs)."""
import csv, json, os
from collections import Counter
from scipy.stats import mannwhitneyu, fisher_exact

ANTS = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]
REFS = ["CD19", "ERBB2", "FOLR1", "TNFRSF17"]
NETS = {"293T": "data/bioplex/BioPlex_293T_Network_10K_Dec_2019.tsv",
        "HCT116": "data/bioplex/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv"}


def gene_set():
    gs = json.load(open("results/depmap_gene_sets.json"))
    genes = {}
    for k in ("gated", "background", "references", "positive_controls"):
        for g in gs[k]:
            genes[g] = k
    return genes


def load_net(path):
    nodes, edges = set(), set()
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            a, b = r["SymbolA"], r["SymbolB"]
            if a == b:
                continue
            nodes.add(a); nodes.add(b)
            edges.add(tuple(sorted((a, b))))
    return nodes, edges


def main():
    genes = gene_set()
    nets = {k: load_net(p) for k, p in NETS.items()}
    universe = set(genes)

    # gates
    n293_nodes, n293_edges = nets["293T"]
    nh_nodes, nh_edges = nets["HCT116"]
    psma1_b = {n: sorted({b for a, b in nets[n][1] if a == "PSMA1" and b.startswith("PSMB")} |
                         {a for a, b in nets[n][1] if b == "PSMA1" and a.startswith("PSMB")})
               for n in nets}
    polr2_partners_293 = sorted({b for a, b in n293_edges if a == "POLR2A" and b.startswith("POLR2")} |
                                {a for a, b in n293_edges if b == "POLR2A" and a.startswith("POLR2")})
    ref_presence_293 = {r: r in n293_nodes for r in REFS}
    gates = {
        "G1_psma1_proteasome_both_networks": {
            "pass": all(len(v) >= 1 for v in psma1_b.values()), "detail": psma1_b},
        "G2_polr2a_polII_partner_293T": {
            "pass": len(polr2_partners_293) >= 1, "detail": polr2_partners_293},
        "G3_references_present_293T": {
            "pass": sum(ref_presence_293.values()) >= 3, "detail": ref_presence_293},
    }

    out_nets = {}
    per_gene = {}
    for name, (nodes, edges) in nets.items():
        testable = [(a, b) for i, a in enumerate(ANTS) for b in ANTS[i + 1:]
                    if a in nodes and b in nodes]
        hits = [p for p in testable if p in edges]
        # degree among universe genes present
        deg = Counter()
        for a, b in edges:
            if a in universe: deg[a] += 1
            if b in universe: deg[b] += 1
        gd = [deg[g] for g in universe if genes[g] == "gated" and g in nodes]
        bd = [deg[g] for g in universe if genes[g] == "background" and g in nodes]
        _, p_deg = mannwhitneyu(gd, bd) if gd and bd else (None, None)
        import statistics as st
        # presence Fisher gated vs background
        gp = sum(1 for g in universe if genes[g] == "gated" and g in nodes)
        ga = sum(1 for g in universe if genes[g] == "gated")
        bp = sum(1 for g in universe if genes[g] == "background" and g in nodes)
        ba = sum(1 for g in universe if genes[g] == "background")
        or_p, p_pres = fisher_exact([[gp, ga - gp], [bp, ba - bp]])
        out_nets[name] = {
            "nodes": len(nodes), "edges": len(edges),
            "testable_antigen_pairs": len(testable),
            "antigen_pair_edges": len(hits), "edges_detail": hits,
            "gated_present": gp, "gated_n": ga, "bg_present": bp, "bg_n": ba,
            "presence_fisher_or": or_p, "presence_fisher_p": p_pres,
            "gated_degree_median": st.median(gd) if gd else None,
            "bg_degree_median": st.median(bd) if bd else None,
            "degree_mw_p": p_deg,
        }
        for g in universe:
            per_gene.setdefault(g, {"group": genes[g]})
            per_gene[g][name + "_present"] = g in nodes
            per_gene[g][name + "_degree"] = deg.get(g, 0) if g in nodes else None

    # H3 reproducibility + IntAct concordance over universe pairs
    e293 = {e for e in n293_edges if e[0] in universe and e[1] in universe}
    ehct = {e for e in nh_edges if e[0] in universe and e[1] in universe}
    both_bp = e293 & ehct
    union_bp = e293 | ehct
    jacc = len(both_bp) / len(union_bp) if union_bp else None
    intact_pairs = set()
    with open("results/intact_edges.csv") as fh:
        for r in csv.DictReader(fh):
            intact_pairs.add(tuple(sorted((r["gene_a"], r["gene_b"]))))
    bp_either = union_bp
    inter = bp_either & intact_pairs
    h3 = {"edges_293T_universe": len(e293), "edges_HCT116_universe": len(ehct),
          "edges_both": len(both_bp), "edges_either": len(union_bp),
          "cross_network_jaccard": jacc,
          "intact_universe_edges": len(intact_pairs),
          "bioplex_and_intact": len(inter),
          "bioplex_only": len(bp_either - intact_pairs),
          "intact_only": len(intact_pairs - bp_either)}

    out = {
        "tool": "BioPlex AP-MS networks (bioplex.hms.harvard.edu; HEK293T 10K + HCT116 5.5K, Dec 2019)",
        "n_genes_total": len(genes),
        "gates": gates, "all_gates_pass": all(g["pass"] for g in gates.values()),
        "H1_per_network": out_nets,
        "H3_reproducibility_concordance": h3,
        "and_gate_presence": {a: {n: a in nets[n][0] for n in nets} for a in ANTS},
    }
    with open("results/bioplex_audit.json", "w") as f:
        json.dump(out, f, indent=1)
    with open("results/bioplex_per_gene.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "group", "293T_present", "293T_degree", "HCT116_present", "HCT116_degree"])
        for g in sorted(per_gene):
            r = per_gene[g]
            w.writerow([g, r["group"], r["293T_present"], r["293T_degree"],
                        r["HCT116_present"], r["HCT116_degree"]])
    print(json.dumps({"gates": {k: v["pass"] for k, v in gates.items()},
                      "H1": {n: {k: v[k] for k in ("testable_antigen_pairs", "antigen_pair_edges",
                                                   "gated_present", "gated_n", "bg_present", "bg_n",
                                                   "presence_fisher_p", "gated_degree_median",
                                                   "bg_degree_median", "degree_mw_p")}
                             for n, v in out_nets.items()},
                      "H3": h3}, indent=1))


if __name__ == "__main__":
    main()

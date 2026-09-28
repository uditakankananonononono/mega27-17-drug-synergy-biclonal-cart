"""Amendment-queue Tier-2 #4: unified 10k label-permutation empirical p for the
pre-specified primary endpoints (results/stats_prespec.json, declared 2026-09-28
BEFORE this rerun) + BH FDR across the family (#5 part 2). Each endpoint recomputed
from its committed per-gene CSV with 10,000 group-label permutations, seed
20260928. intact keeps its existing degree-matched 20k null (already >= 10k;
pair-level, not gene-label permutation) - restated in the unified table with that
provenance noted. Output: results/unified_permutation_p.json."""
import json, csv
import numpy as np

N_PERM = 10000
rng = np.random.default_rng(20260928)
prespec = json.load(open("results/stats_prespec.json"))

def load(path):
    rows = list(csv.DictReader(open(path)))
    g = [r for r in rows if r["group"] == "gated"]
    b = [r for r in rows if r["group"] == "background"]
    return g, b

def num_or_nan(x):
    if x in ("", "NA", "nan"): return float("nan")
    return num(x)

def num(x):
    if x in ("True", "TRUE", "true"): return 1.0
    if x in ("False", "FALSE", "false", "", "NA"): return 0.0
    return float(x)

def carriage_perm(g_vals, b_vals, direction):
    """carriage = fraction > 0; stat = gated_rate - bg_rate."""
    g = np.array(g_vals, float); b = np.array(b_vals, float)
    obs = (g > 0).mean() - (b > 0).mean()
    pooled = np.concatenate([(g > 0).astype(float), (b > 0).astype(float)])
    n = len(g)
    null = np.array([rng.permutation(pooled)[:n].mean() - rng.permutation(pooled)[n:].mean()
                     for _ in range(0)])  # placeholder replaced below
    # vectorized permutation
    perms = np.array([rng.permutation(pooled) for _ in range(N_PERM)])
    null = perms[:, :n].mean(1) - perms[:, n:].mean(1)
    return emp_p(obs, null, direction)

def median_perm(g_vals, b_vals, direction):
    g = np.array(g_vals, float); b = np.array(b_vals, float)
    obs = np.median(g) - np.median(b)
    pooled = np.concatenate([g, b])
    n = len(g)
    perms = np.array([rng.permutation(pooled) for _ in range(N_PERM)])
    null = np.median(perms[:, :n], axis=1) - np.median(perms[:, n:], axis=1)
    return emp_p(obs, null, direction)

def emp_p(obs, null, direction):
    if direction == "greater":
        k = (null >= obs).sum()
    elif direction == "less":
        k = (null <= obs).sum()
    else:
        k = (np.abs(null) >= abs(obs)).sum()
    return float((1 + k) / (1 + len(null))), float(obs), float(np.median(null))

results = []
g, b = load("results/bioplex_per_gene.csv")
p, obs, nmed = carriage_perm([num(r["293T_present"]) for r in g], [num(r["293T_present"]) for r in b], "two-sided")
results.append({"audit": "bioplex", "endpoint": "H1_per_network presence (293T)", "observed_stat": obs, "null_median": nmed, "empirical_p": p})

g, b = load("results/expression_atlas_per_gene.csv")
p, obs, nmed = median_perm([v for v in (num_or_nan(r["tau"]) for r in g) if v == v], [v for v in (num_or_nan(r["tau"]) for r in b) if v == v], "greater")
results.append({"audit": "expression_atlas", "endpoint": "H1_tau", "observed_stat": obs, "null_median": nmed, "empirical_p": p})

g, b = load("results/gwas_per_gene.csv")
p, obs, nmed = carriage_perm([num(r["n_sig"]) for r in g], [num(r["n_sig"]) for r in b], "two-sided")
results.append({"audit": "gwas", "endpoint": "H1_any_sig", "observed_stat": obs, "null_median": nmed, "empirical_p": p})

g, b = load("results/harmonizome_per_gene.csv")
p, obs, nmed = median_perm([v for v in (num_or_nan(r["systematic"]) for r in g) if v == v], [v for v in (num_or_nan(r["systematic"]) for r in b) if v == v], "two-sided")
results.append({"audit": "harmonizome", "endpoint": "H2_control nonexpression (systematic)", "observed_stat": obs, "null_median": nmed, "empirical_p": p})

ia = json.load(open("results/intact_audit.json"))["H1_pairwise_antigen_complex"]
results.append({"audit": "intact", "endpoint": "H1_pairwise_antigen_complex",
                "observed_stat": ia["antigen_pair_fraction"], "null_median": None,
                "empirical_p": float(ia["perm_degmatched_p"]),
                "provenance": "existing degree-matched 20,000-draw null from intact_audit.json (>= 10k; pair-level design, not a gene-label permutation)"})

g, b = load("results/monarch_per_gene.csv")
p, obs, nmed = carriage_perm([num(r["causal_total"]) for r in g], [num(r["causal_total"]) for r in b], "two-sided")
results.append({"audit": "monarch", "endpoint": "H1_causal_carriage", "observed_stat": obs, "null_median": nmed, "empirical_p": p})

g, b = load("results/pharos_per_gene.csv")
p, obs, nmed = carriage_perm([1 if r["tdl"] == "Tclin" else 0 for r in g], [1 if r["tdl"] == "Tclin" else 0 for r in b], "two-sided")
results.append({"audit": "pharos", "endpoint": "H1_tclin", "observed_stat": obs, "null_median": nmed, "empirical_p": p})

g, b = load("results/proteomicsdb_per_gene.csv")
p, obs, nmed = median_perm([v for v in (num_or_nan(r["n_normal_tissues"]) for r in g) if v == v], [v for v in (num_or_nan(r["n_normal_tissues"]) for r in b) if v == v], "less")
results.append({"audit": "proteomicsdb", "endpoint": "H1_normal_breadth", "observed_stat": obs, "null_median": nmed, "empirical_p": p})

g, b = load("results/targetscan_per_gene.csv")
p, obs, nmed = carriage_perm([num(r["conserved_sites"]) for r in g], [num(r["conserved_sites"]) for r in b], "two-sided")
results.append({"audit": "targetscan", "endpoint": "H1_any_site", "observed_stat": obs, "null_median": nmed, "empirical_p": p})

ps = np.array([r["empirical_p"] for r in results])
order = np.argsort(ps)
m = len(ps)
bh = np.empty(m)
prev = 1.0
for rank, idx in enumerate(order[::-1]):
    adj = ps[idx] * m / (m - rank)
    prev = min(prev, adj)
    bh[idx] = prev
for i, r in enumerate(results):
    r["bh_q"] = float(bh[i])
    r["fdr_significant_q05"] = bool(bh[i] < 0.05)

out = {"n_perm": N_PERM, "seed": 20260928, "family": "9 pre-specified primary endpoints (results/stats_prespec.json)",
       "fdr_method": "BH q=0.05", "n_significant": int(sum(r["fdr_significant_q05"] for r in results)),
       "endpoints": results}
print(json.dumps(out, indent=1))
json.dump(out, open("results/unified_permutation_p.json", "w"), indent=1)

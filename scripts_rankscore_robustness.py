"""Tier-2 #17 (rankscore leg): robustness sweep of the locked panel ranking.
Locked weights from results/rankscore_v1.json. Sweep: 1,000 Dirichlet weight draws
(concentration = locked weights x 200, seed 20260929) + one-at-a-time zero-out.
Reported as measured: median Spearman vs locked ranking, P(CA9 stays #1),
P(top-3 set holds), P(top-5 set holds), zero-out podium table.
Reads ONLY committed component values; writes results/rankscore_robustness.json."""
import json
import numpy as np
from scipy import stats

d = json.load(open("results/rankscore_v1.json"))
genes = d["ranking"]
comps = list(d["weights"].keys())
W0 = np.array([d["weights"][c] for c in comps])
C = np.array([[d["rows"][g]["components"][c] for c in comps] for g in genes])
locked_order = list(np.argsort(-(C @ W0)))

rng = np.random.default_rng(20260929)
draws = rng.dirichlet(W0 * 200.0, size=1000)
spears, ca9_top, top3_hold, top5_hold = [], 0, 0, 0
locked_top3, locked_top5 = set(np.argsort(-(C @ W0))[:3]), set(np.argsort(-(C @ W0))[:5])
for w in draws:
    order = np.argsort(-(C @ w))
    spears.append(stats.spearmanr(order, locked_order).statistic)
    ca9_top += int(order[0] == locked_order[0])
    top3_hold += int(set(order[:3]) == locked_top3)
    top5_hold += int(set(order[:5]) == locked_top5)
zero_out = {}
for i, c in enumerate(comps):
    w = W0.copy(); w[i] = 0.0; w = w / w.sum()
    order = np.argsort(-(C @ w))
    zero_out[c] = {"new_ranking": [genes[j] for j in order], "podium_changed": bool(set(order[:3]) != locked_top3)}
spears = np.array(spears)
out = {
 "design": "1,000 Dirichlet draws, concentration = locked weights x 200, seed 20260929; one-at-a-time zero-out renormalized",
 "locked_weights": d["weights"], "locked_ranking": genes,
 "spearman_vs_locked": {"median": float(np.median(spears)), "iqr": [float(np.quantile(spears, 0.25)), float(np.quantile(spears, 0.75))], "min": float(spears.min())},
 "p_ca9_stays_first": ca9_top / 1000.0,
 "p_top3_set_holds": top3_hold / 1000.0,
 "p_top5_set_holds": top5_hold / 1000.0,
 "zero_out_sensitivity": zero_out,
 "boundary": "robustness of the ranking WITHIN the locked 7-gene panel; displacement by non-panel candidates needs the full candidate component matrix - logged as remaining #17 gap",
}
json.dump(out, open("results/rankscore_robustness.json", "w"), indent=1)
print(json.dumps({k: out[k] for k in ["spearman_vs_locked", "p_ca9_stays_first", "p_top3_set_holds", "p_top5_set_holds"]}, indent=1))
for c, z in zero_out.items():
    print(f"zero {c}: podium_changed={z['podium_changed']} top3={z['new_ranking'][:3]}")

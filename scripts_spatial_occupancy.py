"""Amendment-queue Tier-4 #10: spatial occupancy of candidate antigen pairs.
Prespec: results/spatial_occupancy_prespec.json (committed 9ad8274 BEFORE this run).
Jaccard of positive-spot sets per pair per HER2ST section vs 1,000 toroidal-shift
nulls (autocorrelation-preserving). Verdicts per the pre-declared rule; positive
controls gate the method. Cross-lane data: Lane B features_cache (public, same owner).
Run: python3 scripts_spatial_occupancy.py -> results/spatial_occupancy.json"""
import json, os
import numpy as np
from scipy.spatial import cKDTree

SEED = 20260929
N_SHIFT = 1000
MIN_POS = 10
LANEB = "/home/sandbox/mega27-17s-spatial-morphoscan/features_cache"
ENSG = {"CA9": "ENSG00000107159", "ENPP3": "ENSG00000154269", "SLC17A3": "ENSG00000124564",
        "MSLN": "ENSG00000102854", "PSCA": "ENSG00000167653", "SLC39A6": "ENSG00000141424",
        "CA12": "ENSG00000074410", "ERBB2": "ENSG00000141736", "TACSTD2": "ENSG00000184292",
        "MUC1": "ENSG00000185499"}
PAIRS = {"primary": [("CA9", "ENPP3"), ("CA9", "SLC17A3")],
         "secondary": [("MSLN", "CA9"), ("MSLN", "PSCA")],
         "positive_controls": [("ERBB2", "TACSTD2"), ("ERBB2", "MUC1"), ("SLC39A6", "CA12")]}
rng = np.random.default_rng(SEED)

def jacc(a, b):
    u = np.logical_or(a, b).sum()
    return float(np.logical_and(a, b).sum() / u) if u else 0.0

sections = sorted(f for f in os.listdir(LANEB) if f.endswith(".npz"))
pair_sec = {}
for fn in sections:
    z = np.load(os.path.join(LANEB, fn), allow_pickle=True)
    genes = list(z["genes"])
    M = z["counts_sparse"].item()
    px = z["px"]
    sec = fn.replace(".npz", "")
    cols, pos = {}, {}
    for g, e in ENSG.items():
        if e in genes:
            c = M[:, genes.index(e)].toarray().ravel()
            cols[g] = c
            pos[g] = c > 0
    # per-section torus extents
    mins, maxs = px.min(0), px.max(0)
    spans = maxs - mins
    tree = cKDTree(px)
    for cls, pairs in PAIRS.items():
        for g1, g2 in pairs:
            if g1 not in pos or g2 not in pos:
                continue
            if pos[g1].sum() < MIN_POS or pos[g2].sum() < MIN_POS:
                continue
            obs = jacc(pos[g1], pos[g2])
            # shift null on gene2's COUNT vector (keeps its marginal intensity map)
            v2 = cols[g2]
            nulls = np.empty(N_SHIFT)
            for s in range(N_SHIFT):
                off = rng.random(2) * spans
                q = (px + off - mins) % spans + mins
                _, nn = tree.query(q)
                nulls[s] = jacc(pos[g1], v2[nn] > 0)
            p_hi = float((np.sum(nulls >= obs) + 1) / (N_SHIFT + 1))  # observed ABOVE null -> co-occupancy
            p_lo = float((np.sum(nulls <= obs) + 1) / (N_SHIFT + 1))  # observed BELOW null -> segregation
            pair_sec.setdefault(f"{g1}x{g2}", {"class": cls, "sections": []})
            pair_sec[f"{g1}x{g2}"]["sections"].append(
                {"section": sec, "n_pos_g1": int(pos[g1].sum()), "n_pos_g2": int(pos[g2].sum()),
                 "observed_jaccard": obs, "null_median": float(np.median(nulls)),
                 "p_co": p_hi, "p_seg": p_lo})
    print(f"{sec} done", flush=True)

out = {"seed": SEED, "prespec": "results/spatial_occupancy_prespec.json", "n_shift": N_SHIFT,
       "min_pos": MIN_POS, "data": "Lane B features_cache HER2ST (cross-lane, public, same owner)",
       "pairs": {}}
method_failed = []
for pair, rec in pair_sec.items():
    secs = rec["sections"]
    if not secs:
        out["pairs"][pair] = {"class": rec["class"], "n_sections": 0, "verdict": "NO QUALIFYING SECTIONS"}
        continue
    obs_med = float(np.median([s["observed_jaccard"] for s in secs]))
    null_med = float(np.median([s["null_median"] for s in secs]))
    f_co = float(np.mean([s["p_co"] < 0.05 for s in secs]))
    f_seg = float(np.mean([s["p_seg"] < 0.05 for s in secs]))
    if f_co >= 0.6 and obs_med > null_med:
        verdict = "CO-OCCUPANCY"
    elif f_seg >= 0.6 and obs_med < null_med:
        verdict = "SEGREGATION"
    else:
        verdict = "INDETERMINATE"
    out["pairs"][pair] = {"class": rec["class"], "n_sections": len(secs),
                          "median_observed_jaccard": obs_med, "median_null_jaccard": null_med,
                          "frac_sections_co_p05": f_co, "frac_sections_seg_p05": f_seg,
                          "verdict": verdict, "sections": secs}
    if rec["class"] == "positive_controls" and verdict != "CO-OCCUPANCY":
        method_failed.append(pair)
    print(f"{pair} [{rec['class']}]: n={len(secs)} obs {obs_med:.3f} vs null {null_med:.3f} co% {f_co:.2f} seg% {f_seg:.2f} -> {verdict}")
out["method_check"] = {"failed_controls": method_failed,
                       "status": "METHOD OK - pair verdicts stand" if not method_failed else "METHOD FAILED on listed controls - per prespec, pair verdicts WITHHELD"}
json.dump(out, open("results/spatial_occupancy.json", "w"), indent=1)
print("wrote results/spatial_occupancy.json |", out["method_check"]["status"])

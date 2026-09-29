"""Tier-1 #17 gap closure: cross-candidate displacement from the committed full
battery (results/rankscore_battery_full.csv). PRE-DECLARED here: additive
diagnostic, no gate, weights stay frozen. For each axis, drop it and renormalize
the remaining frozen weights; report every gene's rank in each leave-one-out
ranking plus the locked ranking, and the top-20 displacement table (rank delta).
Run: python3 src/tools/rankscore_displacement.py -> results/rankscore_displacement.json"""
import json
import pandas as pd
import numpy as np

W = {"tumor_signal": 0.30, "normal_safety": 0.25, "protein_breadth": 0.15,
     "tractability": 0.15, "nonessentiality": 0.15}
df = pd.read_csv("results/rankscore_battery_full.csv")
df = df[df["rankscore"].notna()].copy()

def score(frame, drop=None):
    axes = [a for a in W if a != drop]
    s = np.zeros(len(frame))
    wsum = np.zeros(len(frame))
    for a in axes:
        v = frame[a].values.astype(float)
        ok = ~np.isnan(v)
        s[ok] += W[a] * v[ok]
        wsum[ok] += W[a]
    out = np.full(len(frame), np.nan)
    valid = wsum > 0
    out[valid] = s[valid] / wsum[valid]
    return out

df["rank_locked"] = df["rankscore"].rank(ascending=False).astype(int)
for a in W:
    df[f"rank_no_{a}"] = pd.Series(score(df, drop=a), index=df.index).rank(ascending=False)

out = {"item": "Tier-1 #17 cross-candidate displacement (leave-one-axis-out, frozen weights)",
       "weights": W, "n_scored": int(len(df))}
top = df.nsmallest(20, "rank_locked")
rows = []
for _, r in top.iterrows():
    row = {"gene": r["gene"], "rank_locked": int(r["rank_locked"])}
    for a in W:
        row[f"no_{a}"] = int(r[f"rank_no_{a}"])
        row[f"d_{a}"] = int(r[f"rank_no_{a}"]) - int(r["rank_locked"])
    rows.append(row)
out["top20_displacement"] = rows
# podium stability per leave-one-out (top-3 set)
def top3(col):
    return set(df.nsmallest(3, col)["gene"])
t0 = top3("rank_locked")
out["podium_leave_one_out"] = {f"no_{a}": {"top3": sorted(top3(f"rank_no_{a}")),
                                           "podium_holds": top3(f"rank_no_{a}") == t0} for a in W}
# largest displacements anywhere in the corpus
disp = []
for _, r in df.iterrows():
    for a in W:
        disp.append({"gene": r["gene"], "axis_dropped": a,
                     "delta": int(r[f"rank_no_{a}"]) - int(r["rank_locked"])})
disp.sort(key=lambda x: -abs(x["delta"]))
out["largest_displacements"] = disp[:15]
json.dump(out, open("results/rankscore_displacement.json", "w"), indent=1)
print(json.dumps(out["podium_leave_one_out"], indent=1))
print("top5 largest:", disp[:5])

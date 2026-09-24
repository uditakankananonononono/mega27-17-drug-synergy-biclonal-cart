"""STRING network proximity of ALMANAC drug targets vs synergy (item 17).

Hypothesis (falsifiable, motivated by the GDSC same-pathway finding p=2.1e-35):
drug pairs whose top ChEMBL targets directly interact in the STRING human
network synergize LESS often than pairs with non-interacting targets.

Targets: per-drug top-5 ChEMBL targets with max pChEMBL>=6 (same definition as
scripts_chembl_analyze.py). Non-protein ChEMBL target records (cell lines,
organisms, 'Unchecked', 'NON-PROTEIN TARGET') are dropped by a documented name
heuristic: keep names that contain a space and none of the excluded tokens.
Names resolved via STRING API get_string_ids (species 9606, best hit); edges from
STRING API network (required_score 400 = medium confidence) among the resolved set.
Pair label: max MAXSCORE over cell lines >= 50 (caveat: single-test maxima).
Output: results/string_proximity.json; API responses cached in data/string/.
"""
import json, pathlib, time, urllib.request, urllib.parse
from collections import Counter
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

API = "https://string-db.org/api/json/"
C = pathlib.Path("data/string")
EXCL = ("cell line", "unchecked", "unknown", "non-protein", "plasmodium", "cells",
        "virus", "enterococcus", "neurone", "carbonate dehydratase")
# Manual corrections after inspecting every STRING best hit (text match picked the wrong
# gene for these ChEMBL target names); the correct HGNC symbol is queried instead.
OVERRIDE = {"Protein smoothened": "SMO", "Dihydrofolate reductase": "DHFR",
            "Alpha-1A adrenergic receptor": "ADRA1A"}


def post(method, params, cache):
    f = C / cache
    if f.exists():
        return json.loads(f.read_text())
    data = urllib.parse.urlencode({**params, "caller_identity": "mega27_laneF"}).encode()
    for att in range(4):
        try:
            r = json.load(urllib.request.urlopen(urllib.request.Request(API + method, data=data), timeout=90))
            f.write_text(json.dumps(r)); return r
        except Exception:
            time.sleep(3 * (att + 1))
    raise SystemExit("STRING failed: " + method)


drug_targets = {}
for f in pathlib.Path("data/chembl").glob("*.json"):
    rec = json.loads(f.read_text()); tg = Counter()
    for a in rec.get("assays", []):
        if a.get("target") and a.get("pchembl_max", 0) >= 6.0:
            tg[a["target"]] = max(tg[a["target"]], a["pchembl_max"])
    tops = [t for t, _ in tg.most_common(5) if " " in t and not any(x in t.lower() for x in EXCL)]
    if tops:
        drug_targets[int(rec["nsc"])] = tops
names = sorted({t for v in drug_targets.values() for t in v})
q = [OVERRIDE.get(n, n) for n in names]
res = post("get_string_ids", {"identifiers": "\r".join(q), "species": 9606, "limit": 1, "echo_query": 1}, "ids.json")
back = {OVERRIDE.get(n, n): n for n in names}
name2sid = {back[r["queryItem"]]: r["stringId"] for r in res if r["queryItem"] in back}
sid2pref = {r["stringId"]: r["preferredName"] for r in res}
sids = sorted(set(name2sid.values()))
net = post("network", {"identifiers": "\r".join(sids), "species": 9606, "required_score": 400}, "network.json")
edges = {}
for e in net:
    a, b = e["stringId_A"], e["stringId_B"]
    edges[frozenset((a, b))] = max(edges.get(frozenset((a, b)), 0), e["score"])
dsid = {d: {name2sid[t] for t in ts if t in name2sid} for d, ts in drug_targets.items()}
dsid = {d: s for d, s in dsid.items() if s}

df = pd.read_csv("data/almanac_synergy.tsv", sep="\t")
pair = df.groupby(["NSC1", "NSC2"])["MAXSCORE"].max().reset_index()
rows = []
for r in pair.itertuples():
    A, B = dsid.get(r.NSC1), dsid.get(r.NSC2)
    if not A or not B:
        continue
    shared = bool(A & B)
    inter = any(frozenset((a, b)) in edges for a in A for b in B if a != b)
    rows.append((r.NSC1, r.NSC2, r.MAXSCORE >= 50, shared, inter))
R = pd.DataFrame(rows, columns=["a", "b", "syn", "shared", "interact"])
out = {"n_drugs_with_protein_targets": len(drug_targets), "n_target_names": len(names),
       "n_resolved_string": len(name2sid), "n_string_proteins": len(sids),
       "n_edges_score400": len(edges), "n_drugs_resolved": len(dsid), "n_pairs": len(R),
       "resolved_map": {k: sid2pref[v] for k, v in name2sid.items()}, "manual_overrides": OVERRIDE}
for label, mask in [("all_pairs", np.ones(len(R), bool)), ("excluding_shared_target", ~R.shared.values)]:
    S = R[mask]
    a = int((S.interact & S.syn).sum()); b = int((S.interact & ~S.syn).sum())
    c = int((~S.interact & S.syn).sum()); d = int((~S.interact & ~S.syn).sum())
    odds, p = fisher_exact([[a, b], [c, d]])
    out[label] = {"interacting": {"n": a + b, "synergistic": a, "rate": round(a / max(a + b, 1), 4)},
                  "non_interacting": {"n": c + d, "synergistic": c, "rate": round(c / max(c + d, 1), 4)},
                  "fisher_odds": float(odds), "p": float(p)}
json.dump(out, open("results/string_proximity.json", "w"), indent=1)
print(json.dumps(out, indent=1))

"""DGIdb v5 (GraphQL) druggability annotation for AND-gate CAR-T targets.
Real API calls to https://dgidb.org/api/graphql; results committed."""
import json, time, urllib.request

TARGETS = ["CLDN18", "MSLN", "CA9", "CA12", "CLDN6", "PSCA", "SLC39A6", "ERBB2"]
URL = "https://dgidb.org/api/graphql"
Q = """{ genes(names: %s) { nodes { name interactions { drug { name } interactionTypes { type directionality } interactionScore sources { sourceDbName } } } } }"""

out = {}
for gene in TARGETS:
    payload = json.dumps({"query": Q % json.dumps([gene])}).encode()
    req = urllib.request.Request(URL, data=payload, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        data = json.load(r)
    nodes = data.get("data", {}).get("genes", {}).get("nodes", [])
    inters = nodes[0]["interactions"] if nodes else []
    drugs = {}
    srcs = set()
    typed = {"inhibitor": 0, "antibody": 0, "agonist": 0, "other": 0}
    for it in inters:
        d = it["drug"]["name"]
        drugs[d] = max(drugs.get(d, 0.0), it.get("interactionScore") or 0.0)
        for s in it.get("sources", []):
            srcs.add(s["sourceDbName"])
        kinds = {t.get("type") for t in it.get("interactionTypes", []) if t.get("type")}
        if "antibody" in kinds: typed["antibody"] += 1
        elif "inhibitor" in kinds: typed["inhibitor"] += 1
        elif "agonist" in kinds: typed["agonist"] += 1
        else: typed["other"] += 1
    top = sorted(drugs.items(), key=lambda kv: -kv[1])[:10]
    out[gene] = {
        "n_interactions": len(inters),
        "n_unique_drugs": len(drugs),
        "n_sources": len(srcs),
        "sources": sorted(srcs),
        "typed_counts": typed,
        "top_drugs_by_score": [{"drug": d, "score": round(s, 4)} for d, s in top],
    }
    print(f"{gene}: {len(inters)} interactions, {len(drugs)} unique drugs, {len(srcs)} sources")
    time.sleep(1.5)

res = {
    "source": "DGIdb v5 GraphQL API (https://dgidb.org/api/graphql)",
    "date": "2026-09-24",
    "genes": out,
    "note": "interactionScore = DGIdb normalized score (publications x sources). Typed counts classify by interaction type with antibody > inhibitor > agonist precedence.",
}
with open("results/dgidb_interactions.json", "w") as f:
    json.dump(res, f, indent=2)
print("WROTE results/dgidb_interactions.json")

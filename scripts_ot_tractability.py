"""Open Targets Platform tractability of AND-gate CAR-T antigen pairs (item 17).

Question: does the expression-gated AND-gate pipeline surface antigens that are
antibody-tractable per Open Targets (Brown et al. tractability assessment)?
Gated set = unique genes across all per-cancer top AND-gate pairs
(results/andgate_p75_igfiltered.json). Background = 100 random UniProt KW-1003
cell-membrane genes (same surfaceome universe the gate draws from).
clinical_AB = any Open Targets AB-modality label in {Approved Drug,
Advanced Clinical, Phase 1 Clinical}. Fisher exact one-sided (greater).
CAR-T reference antigens (CD19, BCMA, ERBB2, FOLR1) annotated as controls.
Responses cached in data/opentargets/*.json. Output: results/ot_tractability.json
"""
import csv, json, random, time, urllib.request
from pathlib import Path
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "data" / "opentargets"; CACHE.mkdir(exist_ok=True)
API = "https://api.platform.opentargets.org/api/v4/graphql"
CLIN = {"Approved Drug", "Advanced Clinical", "Phase 1 Clinical"}
REFS = {"CD19": "ENSG00000177455", "TNFRSF17": "ENSG00000048462",
        "ERBB2": "ENSG00000141736", "FOLR1": "ENSG00000110195"}

pairs = json.load(open(ROOT / "results" / "andgate_p75_igfiltered.json"))["and_gate_pairs_p75_igfiltered"]
gated = {}
for cancer, ps in pairs.items():
    for p in ps:
        gated[p["a"]] = p["a_name"]; gated[p["b"]] = p["b_name"]

surf = []
with open(ROOT / "data" / "uniprot_cellmembrane.tsv") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        g = (row["Gene Names"] or "").split()
        if g:
            surf.append(g[0])
random.seed(0)
bg_syms = random.sample(sorted(set(surf)), 100)

def gql(query):
    req = urllib.request.Request(API, data=json.dumps({"query": query}).encode(),
                                 headers={"Content-Type": "application/json"})
    for _ in range(4):
        try:
            return json.load(urllib.request.urlopen(req, timeout=60))
        except Exception as e:
            print("retry", e, flush=True); time.sleep(5)
    raise RuntimeError("api failed")

def fetch_batch(ids):
    """ids: list of (alias, ensemblId); returns alias -> tractability list"""
    key = "batch_" + "_".join(a for a, _ in ids)[:60] + f"_{len(ids)}.json"
    f = CACHE / key
    if f.exists():
        return json.load(open(f))
    q = "query { " + " ".join(
        f'{a}: target(ensemblId: "{e}") {{ approvedSymbol tractability {{ modality label value }} }}'
        for a, e in ids) + " }"
    d = gql(q)["data"]
    f.write_text(json.dumps(d))
    time.sleep(1)
    return d

def resolve_symbol(sym):
    f = CACHE / f"search_{sym}.json"
    if f.exists():
        d = json.load(open(f))
    else:
        d = gql(f'{{ search(queryString: "{sym}", entityNames: ["target"], page: {{index:0, size:1}}) {{ hits {{ id name entity }} }} }}')
        f.write_text(json.dumps(d)); time.sleep(1)
    for h in d.get("data", {}).get("search", {}).get("hits", []):
        if h.get("entity") == "target" and h.get("name", "").upper() == sym.upper():
            return h["id"]
    return None

def clin_ab(tract):
    return any(t["modality"] == "AB" and t["label"] in CLIN and t["value"] for t in tract)

def ab_labels(tract):
    return sorted({t["label"] for t in tract if t["modality"] == "AB" and t["value"]})

# gated set
gated_rows, gid_list = {}, sorted(gated)
for i in range(0, len(gid_list), 25):
    batch = [(f"t{j}", e) for j, e in enumerate(gid_list[i:i + 25])]
    d = fetch_batch(batch)
    for a, e in batch:
        t = (d.get(a) or {})
        tr = t.get("tractability") or []
        gated_rows[e] = {"symbol": gated[e], "clinical_ab": clin_ab(tr), "ab_labels": ab_labels(tr)}
g_pos = sum(1 for r in gated_rows.values() if r["clinical_ab"])

# background: resolve symbols to ensembl, then tractability
bg_rows, unresolved = {}, []
ens = []
for s in bg_syms:
    e = resolve_symbol(s)
    (ens.append((s, e)) if e else unresolved.append(s))
for i in range(0, len(ens), 25):
    chunk = ens[i:i + 25]
    d = fetch_batch([(f"t{j}", e) for j, (s, e) in enumerate(chunk)])
    for j, (s, e) in enumerate(chunk):
        tr = (d.get(f"t{j}") or {}).get("tractability") or []
        bg_rows[s] = {"ensembl": e, "clinical_ab": clin_ab(tr)}
b_pos = sum(1 for r in bg_rows.values() if r["clinical_ab"])

# references
ref_rows = {}
d = fetch_batch([(s, e) for s, e in REFS.items()])
for s, e in REFS.items():
    tr = (d.get(s) or {}).get("tractability") or []
    ref_rows[s] = {"clinical_ab": clin_ab(tr), "ab_labels": ab_labels(tr)}

table = [[g_pos, len(gated_rows) - g_pos], [b_pos, len(bg_rows) - b_pos]]
odds, p = fisher_exact(table, alternative="greater")
out = {"gated": {"n": len(gated_rows), "clinical_ab": g_pos, "rows": gated_rows},
       "background": {"n": len(bg_rows), "clinical_ab": b_pos, "unresolved_symbols": unresolved,
                      "rows": bg_rows},
       "references": ref_rows,
       "fisher": {"table": table, "odds": round(odds, 3), "p_one_sided": p}}
json.dump(out, open(ROOT / "results" / "ot_tractability.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k in ("references", "fisher")}, indent=1))
print("gated", g_pos, "/", len(gated_rows), "bg", b_pos, "/", len(bg_rows), flush=True)

#!/usr/bin/env python3
"""GTEx normal-tissue median expression of top AND-gate CAR-T targets (off-tumor safety check).
Real data: GTEx Portal API v2, dataset gtex_v8, per-sample TPM vectors -> median per tissue."""
import json, statistics, time, urllib.request

TARGETS = ["CLDN18","MSLN","CA9","CA12","CLDN6","PSCA","SLC39A6"]
BASE = "https://gtexportal.org/api/v2"

def get(url):
    req = urllib.request.Request(url, headers={"User-Agent":"mega27/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())

out = {}
for sym in TARGETS:
    ref = get(f"{BASE}/reference/gene?geneId={sym}")["data"]
    gid = ref[0]["gencodeId"]
    expr = get(f"{BASE}/expression/geneExpression?gencodeId={gid}&datasetId=gtex_v8")["data"]
    med = {}
    nsamp = 0
    for d in expr:
        vals = d["data"]
        med[d["tissueSiteDetailId"]] = statistics.median(vals)
        nsamp += len(vals)
    top = sorted(med.items(), key=lambda kv: -kv[1])[:6]
    vals_sorted = sorted(med.values())
    out[sym] = {"gencode": gid, "n_tissues": len(med), "n_samples": nsamp,
                "max_median_tpm": top[0][1] if top else None,
                "top_tissues": top,
                "global_median": vals_sorted[len(vals_sorted)//2] if vals_sorted else None}
    print(sym, gid, "tissues", len(med), "samples", nsamp, "max", top[0] if top else None)
    time.sleep(0.5)

json.dump(out, open("results/gtex_safety.json","w"), indent=1)
print("wrote results/gtex_safety.json")

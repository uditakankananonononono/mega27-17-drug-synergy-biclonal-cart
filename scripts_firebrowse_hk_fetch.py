#!/usr/bin/env python3
"""Fetch extra housekeeping genes (ACTB, RPLP0, TBP) for the tumor-vs-normal
normalization offset (same FireBrowse endpoint/protocol as the main panel).
Cached in data/firebrowse/hk_<cohort>.json."""
import json, os, time, urllib.parse, urllib.request

GENES = ["ACTB", "RPLP0", "TBP"]
COHORTS = ["PAAD", "COAD", "LUAD", "GBM", "BRCA"]
BASE = "https://firebrowse.org/api/v1/Samples/mRNASeq"

for cohort in COHORTS:
    path = f"data/firebrowse/hk_{cohort}.json"
    if os.path.exists(path):
        print(cohort, "cached"); continue
    recs, page = [], 1
    while True:
        q = urllib.parse.urlencode({"gene": ",".join(GENES), "cohort": cohort,
                                    "sample_type": "TP,NT", "format": "json",
                                    "page_size": 2000, "page": page})
        req = urllib.request.Request(f"{BASE}?{q}", headers={"User-Agent": "mega27/1.0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            batch = json.loads(r.read().decode())["mRNASeq"]
        recs.extend(b for b in batch if b.get("expression_log2") is not None)
        if len(batch) < 2000:
            break
        page += 1
        time.sleep(0.3)
    json.dump(recs, open(path, "w"))
    print(cohort, len(recs))
    time.sleep(0.5)
print("done")

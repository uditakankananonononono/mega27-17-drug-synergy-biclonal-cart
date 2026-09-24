#!/usr/bin/env python3
"""Fetch TCGA per-participant tumor vs adjacent-normal expression for the
AND-gate CAR-T antigen panel via the FireBrowse (Broad Firehose) REST API.

Genes: 7 AND-gate targets (CLDN18, MSLN, CA9, CA12, CLDN6, PSCA, SLC39A6),
references (MS4A1/CD19, TNFRSF17/BCMA, ERBB2), control GAPDH.
Cohorts: PAAD, COAD, LUAD, GBM, BRCA (same five as the CPTAC protein audit).
Endpoint: GET /api/v1/Samples/mRNASeq (RSEM log2), comma-separated
genes/cohorts/sample_types, paginated 2,000 records/page.
Raw pages cached in data/firebrowse/<cohort>.json (untracked; refetch if wiped)."""
import json, os, time, urllib.parse, urllib.request

GENES = ["CLDN18", "MSLN", "CA9", "CA12", "CLDN6", "PSCA", "SLC39A6",
         "MS4A1", "TNFRSF17", "ERBB2", "GAPDH"]
COHORTS = ["PAAD", "COAD", "LUAD", "GBM", "BRCA"]
BASE = "https://firebrowse.org/api/v1/Samples/mRNASeq"

def get(cohort, page):
    q = urllib.parse.urlencode({"gene": ",".join(GENES), "cohort": cohort,
                                "sample_type": "TP,NT", "format": "json",
                                "page_size": 2000, "page": page})
    req = urllib.request.Request(f"{BASE}?{q}", headers={"User-Agent": "mega27/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read().decode())["mRNASeq"]

for cohort in COHORTS:
    path = f"data/firebrowse/{cohort}.json"
    if os.path.exists(path):
        print(cohort, "cached"); continue
    recs, page = [], 1
    while True:
        batch = get(cohort, page)
        recs.extend(batch)
        if len(batch) < 2000:
            break
        page += 1
        time.sleep(0.3)
    json.dump(recs, open(path, "w"))
    print(cohort, len(recs), "records,", page, "pages")
    time.sleep(0.5)
print("done")

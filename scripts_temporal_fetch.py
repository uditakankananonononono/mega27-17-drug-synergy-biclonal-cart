#!/usr/bin/env python3
"""#14 temporal: re-fetch strict_onc tier for the 205-gene corpus (start dates needed).
Emits results/temporal_clintrials_trials.csv (gene, nct, start, status, phase).
Same query and fields as scripts_clintrials_fetch.py strict_onc tier."""
import csv, json, time, urllib.request, urllib.parse
API = "https://clinicaltrials.gov/api/v2/studies"
FIELDS = "NCTId,BriefTitle,OverallStatus,Phase,StartDate"
ONC = "(cancer OR neoplasm OR tumor OR tumour OR carcinoma OR lymphoma OR leukemia OR leukaemia OR melanoma OR myeloma OR sarcoma OR glioma OR blastoma OR adenocarcinoma)"

def query_all(term, cap=1000):
    recs, token = [], None
    while True:
        q = {"query.term": term, "pageSize": 100, "fields": FIELDS}
        if token: q["pageToken"] = token
        url = API + "?" + urllib.parse.urlencode(q)
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        for attempt in range(4):
            try:
                with urllib.request.urlopen(req, timeout=30) as r:
                    d = json.load(r)
                break
            except Exception:
                if attempt == 3: raise
                time.sleep(2 * (attempt + 1))
        for s in d.get("studies", []):
            ps = s.get("protocolSection", {})
            idm = ps.get("identificationModule", {})
            st = ps.get("statusModule", {})
            dm = ps.get("designModule", {})
            recs.append({"nct": idm.get("nctId"), "status": st.get("overallStatus", ""),
                         "phase": (dm.get("phases") or [""])[0] if dm.get("phases") else "",
                         "start": (st.get("startDateStruct") or {}).get("date", "")})
        token = d.get("nextPageToken")
        if not token or len(recs) >= cap: break
        time.sleep(0.12)
    return recs

sets = json.load(open("results/depmap_gene_sets.json"))
genes = sets["gated"] + sets["background"] + sets["references"] + sets["positive_controls"]
rows = []
for i, g in enumerate(genes, 1):
    recs = query_all(f"AREA[InterventionName] {g} AND {ONC}")
    for r in recs:
        rows.append([g, r["nct"], r["start"], r["status"], r["phase"]])
    if i % 20 == 0: print(f"{i}/{len(genes)} rows={len(rows)}", flush=True)
    time.sleep(0.12)
with open("results/temporal_clintrials_trials.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["gene", "nct", "start", "status", "phase"])
    w.writerows(rows)
meta = {"fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "api": API, "fields": FIELDS, "n_genes": len(genes), "n_rows": len(rows)}
json.dump(meta, open("results/temporal_clintrials_meta.json", "w"), indent=1)
print("done", len(rows))

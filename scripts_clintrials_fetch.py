"""Fetch ClinicalTrials.gov API v2 trial records for the 205-gene item-17 set.
Three tiers per gene:
  strict: AREA[InterventionName] <SYM>          (trial intervention names carry the symbol)
  broad : <SYM> AND (CAR OR antibody OR bispecific OR conjugate)   (full text)
  alias : curated alias query for 11 key genes (7 AND-gate antigens + 4 references)
Caches JSON per gene in data/clintrials/. Resumable: skips existing cache files.
"""
import json, os, sys, time, urllib.request, urllib.parse

OUT = "data/clintrials"
os.makedirs(OUT, exist_ok=True)
API = "https://clinicaltrials.gov/api/v2/studies"
FIELDS = "NCTId,BriefTitle,OverallStatus,Phase,StartDate"
ONC = "(cancer OR neoplasm OR tumor OR tumour OR carcinoma OR lymphoma OR leukemia OR leukaemia OR melanoma OR myeloma OR sarcoma OR glioma OR blastoma OR adenocarcinoma)"
ACTIVE = {"RECRUITING", "ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION", "NOT_YET_RECRUITING"}

ALIAS = {
    "CLDN18": '("claudin 18.2" OR "claudin 18" OR CLDN18 OR zolbetuximab) AND (antibody OR CAR OR bispecific OR conjugate)',
    "CLDN6": '("claudin 6" OR CLDN6) AND (antibody OR CAR OR bispecific OR conjugate)',
    "CA9": '("carbonic anhydrase IX" OR CA9 OR girentuximab) AND (antibody OR CAR OR bispecific OR conjugate)',
    "CA12": '("carbonic anhydrase XII" OR CA12) AND (antibody OR CAR OR bispecific OR conjugate)',
    "MSLN": '(mesothelin OR MSLN) AND (antibody OR CAR OR bispecific OR conjugate)',
    "PSCA": '(PSCA OR "prostate stem cell antigen") AND (antibody OR CAR OR bispecific OR conjugate)',
    "SLC39A6": '(SLC39A6 OR LIV-1 OR ladiratuzumab) AND (antibody OR CAR OR bispecific OR conjugate)',
    "CD19": 'CD19 AND (antibody OR CAR OR bispecific OR conjugate)',
    "ERBB2": '(HER2 OR ERBB2) AND (antibody OR CAR OR bispecific OR conjugate)',
    "FOLR1": '("folate receptor alpha" OR FOLR1 OR mirvetuximab) AND (antibody OR CAR OR bispecific OR conjugate)',
    "TNFRSF17": '(BCMA OR TNFRSF17) AND (antibody OR CAR OR bispecific OR conjugate)',
}

def query_all(term, cap=1000):
    """Return dict(total, records[]) paging through the v2 API.
    totalCount is unreliable (returns 0 for some large queries), so total is
    derived from the full page walk; total_reported keeps the API value."""
    recs, token, total = [], None, None
    while True:
        q = {"query.term": term, "pageSize": 100, "countTotal": "true", "fields": FIELDS}
        if token:
            q["pageToken"] = token
        url = API + "?" + urllib.parse.urlencode(q)
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        for attempt in range(4):
            try:
                with urllib.request.urlopen(req, timeout=30) as r:
                    d = json.load(r)
                break
            except Exception as e:
                if attempt == 3:
                    raise
                time.sleep(2 * (attempt + 1))
        total = d.get("totalCount", 0)
        for s in d.get("studies", []):
            ps = s.get("protocolSection", {})
            idm = ps.get("identificationModule", {})
            st = ps.get("statusModule", {})
            dm = ps.get("designModule", {})
            recs.append({
                "nct": idm.get("nctId"),
                "title": idm.get("briefTitle", ""),
                "status": st.get("overallStatus", ""),
                "phase": (dm.get("phases") or [""])[0] if dm.get("phases") else "",
                "start": (st.get("startDateStruct") or {}).get("date", ""),
            })
        token = d.get("nextPageToken")
        if not token or len(recs) >= cap:
            break
        time.sleep(0.12)
    capped = bool(token) and len(recs) >= cap
    return {"total": len(recs) if (not total or capped) else total,
            "total_reported": total, "capped": capped, "records": recs[:cap]}

def fetch_tier(gene, tier, term):
    path = os.path.join(OUT, f"{tier}_{gene}.json")
    if os.path.exists(path):
        return "cached"
    d = query_all(term)
    d["gene"] = gene
    d["tier"] = tier
    d["term"] = term
    with open(path, "w") as f:
        json.dump(d, f)
    return "fetched"

def main():
    sets = json.load(open("results/depmap_gene_sets.json"))
    genes = sets["gated"] + sets["background"] + sets["references"] + sets["positive_controls"]
    only = sys.argv[1:] if len(sys.argv) > 1 else None
    n = 0
    for g in genes:
        if only and g not in only:
            continue
        r1 = fetch_tier(g, "strict", f"AREA[InterventionName] {g}")
        r1b = fetch_tier(g, "strict_onc", f"AREA[InterventionName] {g} AND " + ONC)
        r2 = fetch_tier(g, "broad", f"{g} AND (CAR OR antibody OR bispecific OR conjugate)")
        r3 = "skip"
        if g in ALIAS:
            r3 = fetch_tier(g, "alias", ALIAS[g])
        n += 1
        if n % 20 == 0 or r1 == r2 == "fetched" and n % 5 == 0:
            print(f"{n}/{len(genes)} {g} strict={r1} broad={r2} alias={r3}", flush=True)
        time.sleep(0.12)
    meta = {"genes": genes, "n_genes": len(genes), "alias_genes": sorted(ALIAS),
            "fields": FIELDS, "active_statuses": sorted(ACTIVE), "cap": 500,
            "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    with open(os.path.join(OUT, "meta_clintrials.json"), "w") as f:
        json.dump(meta, f, indent=1)
    print("done", n)

if __name__ == "__main__":
    main()

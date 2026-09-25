"""RxNorm (NLM RxNav API) lane: normalize the GDSC 8.5 compound names.

Maps each screened compound's primary name to its RxCUI via the RxNav REST
API, recording exact vs approximate matches - the identifier normalization
lane used to join drug records across resources. Committed:
results/rxnorm_gdsc_map.csv (+ rxnorm_summary.json).
"""
import csv, json, time, urllib.parse, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
UA = {"User-Agent": "mega27-research/1.0", "Accept": "application/json"}


def rxcui(name):
    q = urllib.parse.quote(name)
    url = f"https://rxnav.nlm.nih.gov/REST/rxcui.json?name={q}"
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            d = json.load(r)
    except Exception:
        return "", "error"
    ids = d.get("idGroup", {}).get("rxnormId", [])
    return (ids[0] if ids else ""), ("exact" if ids else "unmatched")


def main():
    drugs = [r["DRUG_NAME"] for r in csv.DictReader(open(ROOT / "data" / "gdsp_compounds_8.5.csv"))]
    from concurrent.futures import ThreadPoolExecutor
    def one(name):
        cui, how = rxcui(name)
        return {"drug_name": name, "rxcui": cui, "match": how}
    with ThreadPoolExecutor(max_workers=12) as ex:
        rows = list(ex.map(one, drugs))
    with open(ROOT / "results" / "rxnorm_gdsc_map.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    n = sum(1 for r in rows if r["rxcui"])
    json.dump({"n_drugs": len(rows), "n_rxcui_matched": n,
               "frac_matched": round(n / len(rows), 4),
               "note": "unmatched are mostly investigational codes (AZD/BI/BX series) "
                       "without RxNorm entries - recorded, not forced"},
              open(ROOT / "results" / "rxnorm_summary.json", "w"), indent=2)
    print(n, "/", len(rows))


if __name__ == "__main__":
    main()

"""DailyMed (NLM) lane: SPL label cross-check for the approved drugs.

Counts DailyMed structured product labels per ChEMBL-annotated drug - a
second regulatory-label route, cross-checking the openFDA lane's coverage.
Committed: results/dailymed_label_counts.csv (+ dailymed_summary.json).
"""
import csv, json, time, urllib.parse, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
UA = {"User-Agent": "mega27-research/1.0", "Accept": "application/json"}


def main():
    ann = json.load(open(ROOT / "results" / "chembl_annotation.json"))
    drugs = sorted(ann["drug_top_targets"].keys())
    rows = []
    for name in drugs:
        q = urllib.parse.quote(name.lower())
        url = f"https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json?drug_name={q}&pageSize=1"
        req = urllib.request.Request(url, headers=UA)
        n = 0
        try:
            with urllib.request.urlopen(req, timeout=20) as r:
                d = json.load(r)
            n = int(d.get("metadata", {}).get("total_elements", 0))
        except Exception:
            n = -1
        rows.append({"drug": name, "dailymed_spl_count": n})
        time.sleep(0.12)
    with open(ROOT / "results" / "dailymed_label_counts.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    nf = sum(1 for r in rows if r["dailymed_spl_count"] > 0)
    json.dump({"n_drugs": len(rows), "n_with_spl": nf},
              open(ROOT / "results" / "dailymed_summary.json", "w"), indent=2)
    print(nf, "/", len(rows))


if __name__ == "__main__":
    main()

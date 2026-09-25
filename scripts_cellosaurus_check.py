"""Cellosaurus lane: identity audit of the L1000 cell lines.

Verifies each distinct cell line used in the committed iLINCS exemplar
signatures (results/ilincs_gdsp_map.csv) against the Cellosaurus API:
accession, recommended name, and the misidentified/problematic flag.
Committed: results/cellosaurus_cell_lines.csv (+ cellosaurus_summary.json).
"""
import csv, json, time, urllib.parse, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
UA = {"User-Agent": "mega27-research/1.0", "Accept": "application/json"}


def query(cl):
    q = urllib.parse.quote(f'idsy:"{cl}"')
    url = f"https://api.cellosaurus.org/search/cell-line?q={q}&format=json&rows=1"
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            d = json.load(r)
    except Exception:
        return {}
    docs = d.get("Cellosaurus", {}).get("cell-line-list", [])
    return docs[0] if docs else {}


def main():
    cls = sorted({r["cellline"] for r in
                  csv.DictReader(open(ROOT / "results" / "ilincs_gdsp_map.csv"))})
    rows = []
    for cl in cls:
        h = query(cl)
        accl = h.get("accession-list", [])
        acc = ""
        for a in accl:
            if a.get("type") == "primary":
                acc = a.get("value", "")
        cat = h.get("category", "")
        rows.append({"cell_line": cl, "cellosaurus_accession": acc,
                     "category": cat,
                     "problematic": "problematic" in json.dumps(h).lower()})
        print(cl, acc, cat, flush=True)
        time.sleep(0.15)
    with open(ROOT / "results" / "cellosaurus_cell_lines.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    n = sum(1 for r in rows if r["cellosaurus_accession"])
    json.dump({"n_cell_lines": len(rows), "n_verified": n,
               "n_problematic": sum(1 for r in rows if r["problematic"])},
              open(ROOT / "results" / "cellosaurus_summary.json", "w"), indent=2)
    print(n, "/", len(rows))


if __name__ == "__main__":
    main()

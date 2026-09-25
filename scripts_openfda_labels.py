"""openFDA lane: label safety check for the approved synergy drugs.

For each ChEMBL-annotated ALMANAC drug, queries the openFDA drug/label
endpoint and records whether an FDA label exists and whether it carries a
boxed warning - the regulatory-safety lane for combo prioritization.
Committed: results/openfda_label_flags.csv (+ openfda_summary.json).
"""
import csv, json, time, urllib.parse, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
UA = {"User-Agent": "mega27-research/1.0", "Accept": "application/json"}


def label(name):
    q = urllib.parse.quote(f'openfda.generic_name:"{name}"')
    url = (f"https://api.fda.gov/drug/label.json?search={q}"
           "&limit=1&count=boxed_warning.exact")
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            d = json.load(r)
    except Exception:
        return None
    res = d.get("results", [])
    return bool(res)


def main():
    ann = json.load(open(ROOT / "results" / "chembl_annotation.json"))
    drugs = sorted(ann["drug_top_targets"].keys())
    rows = []
    for name in drugs:
        q = urllib.parse.quote(f'openfda.generic_name:"{name}"')
        url = f"https://api.fda.gov/drug/label.json?search={q}&limit=1"
        req = urllib.request.Request(url, headers=UA)
        found, boxed = False, False
        try:
            with urllib.request.urlopen(req, timeout=20) as r:
                d = json.load(r)
            res = d.get("results", [])
            found = bool(res)
            boxed = bool(res and res[0].get("boxed_warning"))
        except Exception:
            pass
        rows.append({"drug": name, "fda_label_found": found, "boxed_warning": boxed})
        time.sleep(0.15)
    with open(ROOT / "results" / "openfda_label_flags.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    nf = sum(1 for r in rows if r["fda_label_found"])
    nb = sum(1 for r in rows if r["boxed_warning"])
    json.dump({"n_drugs": len(rows), "n_with_fda_label": nf, "n_boxed_warning": nb,
               "note": "unfound = investigational/never-FDA-labeled - recorded"},
              open(ROOT / "results" / "openfda_summary.json", "w"), indent=2)
    print(nf, nb)


if __name__ == "__main__":
    main()

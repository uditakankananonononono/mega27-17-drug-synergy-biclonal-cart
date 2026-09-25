"""MyChem.info lane: chemical annotation of the synergy drugs.

Annotates the ChEMBL-annotated ALMANAC drugs (the 70 drugs with target
annotations in results/chembl_annotation.json) through the MyChem.info API:
InChIKey, formula, ChEBI/DrugBank cross-references. Used to verify that the
synergy pairs join on the same chemical entities across resources.
Committed: results/mychem_drug_annotations.csv (+ mychem_summary.json).
"""
import csv, json, time, urllib.parse, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
UA = {"User-Agent": "mega27-research/1.0", "Accept": "application/json"}


def find_key(obj, names):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in names and isinstance(v, str):
                return v
            got = find_key(v, names)
            if got:
                return got
    elif isinstance(obj, list):
        for v in obj:
            got = find_key(v, names)
            if got:
                return got
    return ""


def mychem(name):
    q = urllib.parse.quote(f'"{name}"')
    url = f"https://mychem.info/v1/query?q={q}&size=1"
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            d = json.load(r)
    except Exception:
        return {}
    hits = d.get("hits", [])
    return hits[0] if hits else {}


def main():
    ann = json.load(open(ROOT / "results" / "chembl_annotation.json"))
    drugs = sorted(ann["drug_top_targets"].keys())
    rows = []
    for name in drugs:
        h = mychem(name)
        rows.append({"drug": name,
                     "mychem_id": h.get("_id", ""),
                     "inchikey": find_key(h, {"inchikey", "inchi_key"}),
                     "chebi_id": find_key(h.get("chebi"), {"id"}),
                     "chebi_name": find_key(h.get("chebi"), {"name"}),
                     "drugbank_id": find_key(h.get("drugbank"), {"id"})})
        print(name, rows[-1]["inchikey"][:14], flush=True)
        time.sleep(0.12)
    with open(ROOT / "results" / "mychem_drug_annotations.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    n = sum(1 for r in rows if r["inchikey"])
    json.dump({"n_drugs": len(rows), "n_annotated": n},
              open(ROOT / "results" / "mychem_summary.json", "w"), indent=2)
    print(n, "/", len(rows))


if __name__ == "__main__":
    main()

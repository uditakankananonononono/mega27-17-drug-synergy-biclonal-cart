"""Fetch the DrugCentral drug-target interaction snapshot.

Source: https://unmtid-dbs.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz
(curated quantitative drug-target activities with TDL target-development levels).
One bulk snapshot file = one dataset under the program counting rule.
Output: data/drugcentral/drug.target.interaction.tsv.gz
"""
import gzip, pathlib, urllib.request

URL = "https://unmtid-dbs.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz"
OUT = pathlib.Path("data/drugcentral"); OUT.mkdir(exist_ok=True)
F = OUT / "drug.target.interaction.tsv.gz"

if not F.exists():
    req = urllib.request.Request(URL, headers={"User-Agent": "mega27-research"})
    with urllib.request.urlopen(req, timeout=180) as r, open(F, "wb") as fh:
        fh.write(r.read())
    print("downloaded", F, F.stat().st_size, "bytes")
else:
    print("cached", F, F.stat().st_size, "bytes")

with gzip.open(F, "rt") as fh:
    n = sum(1 for _ in fh) - 1
assert n > 19000, f"row count sanity failed: {n}"
print("rows:", n)

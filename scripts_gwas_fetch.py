#!/usr/bin/env python3
"""Download the GWAS Catalog full ontology-annotated associations file
(bulk TSV, latest release) used by the germline-association audit.
Lands in data/gwas/ (untracked cache; re-run to refetch)."""
import os, urllib.request, zipfile

URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
       "gwas-catalog-associations_ontology-annotated-full.zip")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "gwas")

def main():
    os.makedirs(OUT, exist_ok=True)
    if any(fn.endswith(".tsv") for fn in os.listdir(OUT)):
        print("cached: associations TSV present")
        return
    zp = os.path.join(OUT, "gwas_assoc.zip")
    print("downloading:", URL)
    urllib.request.urlretrieve(URL, zp)
    print("unzipping")
    with zipfile.ZipFile(zp) as z:
        z.extractall(OUT)
    print("done")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Download the two BioPlex AP-MS network releases (bulk TSVs) used by the
co-complex audit:
  BioPlex 3.0 HEK293T 10K network (Dec 2019)  118,162 edges
  BioPlex HCT116 5.5K network (Dec 2019)       70,966 edges
Files land in data/bioplex/ (untracked cache; re-run this script to refetch)."""
import os, urllib.request

BASE = "https://bioplex.hms.harvard.edu/data/"
FILES = ["BioPlex_293T_Network_10K_Dec_2019.tsv",
         "BioPlex_HCT116_Network_5.5K_Dec_2019.tsv"]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "bioplex")

def main():
    os.makedirs(OUT, exist_ok=True)
    for fn in FILES:
        dst = os.path.join(OUT, fn)
        if os.path.exists(dst) and os.path.getsize(dst) > 10 ** 6:
            print("cached:", fn)
            continue
        print("downloading:", fn)
        urllib.request.urlretrieve(BASE + fn, dst)
        print("  ->", os.path.getsize(dst), "bytes")

if __name__ == "__main__":
    main()

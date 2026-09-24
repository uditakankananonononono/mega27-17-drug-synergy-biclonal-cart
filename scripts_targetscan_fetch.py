#!/usr/bin/env python3
"""Download TargetScan 8.0 human prediction summaries used by the miRNA
regulatory-burden audit:
  Summary_Counts.default_predictions.txt.zip (per gene x miRNA family)
  Conserved_Site_Context_Scores.txt.zip     (per conserved site)
Land in data/targetscan/ (untracked cache; re-run to refetch)."""
import os, urllib.request, zipfile

BASE = "https://www.targetscan.org/vert_80/vert_80_data_download/"
FILES = ["Summary_Counts.default_predictions.txt.zip",
         "Conserved_Site_Context_Scores.txt.zip"]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "targetscan")

def main():
    os.makedirs(OUT, exist_ok=True)
    for fn in FILES:
        plain = fn[:-4]
        if os.path.exists(os.path.join(OUT, plain)):
            print("cached:", plain)
            continue
        print("downloading:", fn)
        zp = os.path.join(OUT, fn)
        urllib.request.urlretrieve(BASE + fn, zp)
        with zipfile.ZipFile(zp) as z:
            z.extractall(OUT)
        print("  done")
if __name__ == "__main__":
    main()

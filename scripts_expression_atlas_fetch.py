#!/usr/bin/env python3
"""Download EBI Expression Atlas baseline experiment E-MTAB-513 (Illumina
Human Body Map 2.0; 16 normal tissues, RNA-seq) gene TPMs + configuration
(assay group -> tissue) + condensed SDRF into data/expression_atlas/."""
import os, urllib.request

BASE = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-513/"
FILES = ["E-MTAB-513-tpms.tsv", "E-MTAB-513-configuration.xml", "E-MTAB-513.condensed-sdrf.tsv"]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "expression_atlas")


def main():
    os.makedirs(OUT, exist_ok=True)
    for f in FILES:
        p = os.path.join(OUT, f)
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            urllib.request.urlretrieve(BASE + f, p)
        print(f, os.path.getsize(p))


if __name__ == "__main__":
    main()

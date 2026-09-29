#!/bin/bash
# HPA per-sample TCGA tumor expression (rna_cancer_sample.tsv.zip), needed by src/cart/cli.py.
# www.proteinatlas.org/download/ retired this file (404 as of 2026-09-29; v25 redirects to the
# retired path). v23 is the last versioned host serving it - same HPA era as the study's
# committed analyses. Live-verified HTTP 200 on 2026-09-29.
set -e
curl -L --retry 3 -o data/rna_cancer_sample.tsv.zip \
  "https://v23.proteinatlas.org/download/rna_cancer_sample.tsv.zip"   # 1233429065 bytes

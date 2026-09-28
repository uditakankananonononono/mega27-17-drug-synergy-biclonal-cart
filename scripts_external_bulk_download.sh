#!/bin/bash
# Tier-1 #15: independent GEO bulk validation of the rankscore panel.
# 6 series, all public no-auth, URLs + byte sizes live-verified HTTP 200 on 2026-09-28.
set -e
mkdir -p /home/sandbox/geoext && cd /home/sandbox/geoext
B='https://ftp.ncbi.nlm.nih.gov/geo/series'
curl -O "$B/GSE53nnn/GSE53757/matrix/GSE53757_series_matrix.txt.gz"  #  32,826,277 ccRCC vs matched normal kidney, GPL570, 145 samples (72 pairs) - CA9, CA12
curl -O "$B/GSE62nnn/GSE62165/matrix/GSE62165_series_matrix.txt.gz"  #  51,179,939 PDAC, GPL13667, 132 samples - MSLN
curl -O "$B/GSE42nnn/GSE42568/matrix/GSE42568_series_matrix.txt.gz"  #  22,565,943 breast (17 normal + 104 tumor), GPL570 - SLC39A6
curl -O "$B/GSE46nnn/GSE46602/matrix/GSE46602_series_matrix.txt.gz"  #  11,663,959 prostate (14 benign + 36 cancer), GPL570 - PSCA
curl -O "$B/GSE13nnn/GSE13861/matrix/GSE13861_series_matrix.txt.gz"  #  14,543,125 gastric, GPL6884, 91 samples - CLDN18
curl -O "$B/GSE26nnn/GSE26712/matrix/GSE26712_series_matrix.txt.gz"  #  23,678,195 ovarian (10 normal HOSE + 185 tumor), GPL96 - CLDN6

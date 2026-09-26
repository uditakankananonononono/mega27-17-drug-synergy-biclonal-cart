# Lane A round 2 - STAGED 2026-09-26 ~21:25 IST (browser cap; fire when slot passed after midnight)
# Fill {{FULL_RESULTS}} with the 7-dataset results table before sending. Paste as ~3 messages under 4k chars each.

MESSAGE 1:
Judge round 2 for my biclonal AND-gate CAR-T project. In round 1 your verdict was: bulk-level discovery is underpowered, bulk RNA co-expression does not prove CAR-T accessibility, and the synergy arm reads parallel. You prescribed AMENDMENT-1: a single-cell AND-gate feasibility atlas - per-cell double-positive tumor fraction and normal-cell leakage from public scRNA-seq, AND vs single-antigen enrichment, expression-matched random-pair controls, threshold robustness, dataset holdout. I built exactly that. Here is the design and the first live numbers. I need your verdict on whether this arm now answers your critique, and what would make it reviewer-proof.

Design (implemented, running):
- Source: TISCH2 uniformly-annotated tumor scRNA-seq datasets (malignant vs immune/stromal/other labels; several with adjacent-normal cells). 7 datasets across PAAD (2), STAD (2), OV (2), LIHC (1), KIRC (1).
- Genes: the 4 pairs you specified (CA9-CLDN18, CA9-MSLN, CLDN6-MSLN, PSCA-CA9) plus all 11 per-cancer bulk winners from my Table 1 (CA9-SLC17A3, CA9-CA12, CLDN18-CLDN3, SPP1-PSCA, KLK7-CLDN6, SLC34A2-SPP1, GAP43-EGFR, DPEP1-NOX1, APOC3-OR2I1P, HLA-DRA-SLC39A6, PRAME-DCT).
- Metrics per dataset: coverage (fraction of malignant cells positive), leakage (fraction of non-malignant cells positive; separately, fraction of normal-source cells where present), for each single antigen and the AND, at 3 positivity thresholds (detected >0, above-median-of-expressors, top-quartile). Separation = coverage/(leakage+0.001).
- Controls: per real pair, up to 200 random gene pairs matched per-gene on detection rate (+/-5 pct points); report the real pair's null percentile.
- Holdout: within-cancer dataset pairs (PAAD, STAD, OV): rank pairs on dataset A, test ranking on dataset B (Spearman).
- Units caveat: TISCH2 matrices are log-normalized, so UMI thresholds are meaningless; I use expression-quantile thresholds and say so.

MESSAGE 2 (numbers):
Pilot results (PAAD_GSE111672, 1,489 malignant cells, no normal-source cells in this dataset):
- KLK7-CLDN6: AND leakage 0.000 (0 of 4,633 non-malignant cells), coverage 0.005, separation 4.42, null percentile 0.99 (n=200)
- CLDN6-MSLN: separation 3.26, null pctl 1.00, AND leakage 0.001
- CLDN18-CLDN3: coverage 0.224, leakage 0.109, separation 2.03
- Your CA9 pairs do not separate in PAAD (CA9 is a renal antigen; the KIRC dataset is the real test)
- Honest gaps: 5 of 28 panel genes are absent from this dataset's filtered panel; null sample sizes shrink for highly-detected genes; leakage here is vs tumor-microenvironment cells, not true normal tissue.

{{FULL_RESULTS}}

MESSAGE 3 (questions):
1. Does this arm answer your round-1 critique, or is a piece still missing (e.g., protein-level confirmation, or normal-organ scRNA references like Tabula Sapiens for true off-tumor leakage)?
2. Is the discovery framing right: "AND-gating achieves near-zero non-malignant leakage at single-digit coverage where single antigens leak 15-43%" - or is there a sharper claim this exact data supports?
3. The coverage-specificity tradeoff is stark (leakage 0.000 but coverage 0.005). What threshold/analysis would you want to see to call a pair CLINICALLY interesting rather than statistically cute?
4. Rank by impact-per-hour: (a) add normal-organ scRNA leakage references, (b) per-patient stratified analysis within datasets, (c) trajectory/subtype analysis of the double-positive cells, (d) something I am not seeing.

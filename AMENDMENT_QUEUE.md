# AMENDMENT_QUEUE.md - Lane A (mega27-17), LOCKED from ChatGPT judge verdict 2026-09-27

Locked from JUDGE_VERDICT_CHATGPT_20260927.md (conversation 6ab8bdde-c768-83e8-b0fd-5e9144ebbd9f,
gate MET entry in JUDGE_ROUNDS.md, commit 02235bf). Protocol: queue locked BEFORE execution;
each item closes only with committed, live-verified evidence (script + output + hash).
Verdict overall: computational strength 8/10; highest-impact improvement = bias-corrected
antigen-ranking framework with matched controls and external validation (items 1+2+13+12+15).

Status: DONE = evidence committed. PARTIAL = exists, gap named. TODO = not started.

## Tier 1 - the transform package (verdict's highest-impact improvement)
| # | Addition | Status | Evidence / Action |
|---|----------|--------|-------------------|
| 1 | Formal multi-objective target-ranking algorithm | PARTIAL | src/cart/antigen_rank.py (tumor/vital-normal window, benchmarked on CD19/BCMA/HER2/MSLN/GPC3). Gap: objectives limited to expression window; add druggability, novelty, interaction-independence terms with explicit weights + optimization objective. -> scripts_rankscore.py |
| 2 | Propensity-score matched backgrounds | TODO | Match background genes on expression, length, protein abundance, cancer association, literature count; rerun key enrichments vs matched nulls. |
| 12 | Explainable feature contributions per antigen | TODO | Ship with #1: per-candidate percentage contributions (specificity / normal penalty / protein evidence / novelty). |
| 13 | Unified antigen safety score Safety=f(RNA, Protein, GWAS, Mendelian, Essentiality) | TODO | Inputs already committed: results/gtex_safety*, cptac_protein.json, gwas_audit.json, depmap_dependency.json. Validate: known successful targets must rank appropriately. |
| 15 | External cohorts beyond TCGA-like data | PARTIAL | CPTAC protein (results/cptac_protein.json) + 8 external scRNA cohorts (atlas). Gap: ICGC / independent GEO bulk validation of the ranking. |

## Tier 2 - statistical unification (cheap, high rigor value)
| # | Addition | Status | Evidence / Action |
|---|----------|--------|-------------------|
| 4 | Permutation-based null testing, 10k perms, empirical p | PARTIAL | scAtlas null-pctl = 1.0 on all pairs (results/scatlas_per_dataset.json). Gap: unified 10k-permutation empirical-p across the enrichment audits. |
| 5 | FDR-controlled global statistics + pre-specified primary endpoints | TODO | Apply BH across all audit batteries; declare primary endpoints in a pre-spec file before reruns (no PREREG gate edits - additive stats file). |
| 17 | Robustness analysis (threshold/cohort/cutoff sweeps) | PARTIAL | p75/p90 windows (results/rna_window_percentile.json, firebrowse_window.json). Gap: systematic sweep table across all headline claims. |

## Tier 3 - the bias-correction core (technical basis of the #20 reframe)
| # | Addition | Status | Evidence / Action |
|---|----------|--------|-------------------|
| 3 | Annotation-bias correction model (observed - expected annotation) | TODO | Train expected-annotation-density model; report corrected enrichment. |
| 6 | Missingness model for databases | TODO | P(observation) = f(expression, popularity, assay coverage, protein properties); test conclusion survival. |
| 7 | Protein detectability prediction | TODO | MS detectability from peptide/abundance/sequence features; observed vs expected breadth. |
| 8 | RNA->protein translation model | PARTIAL | CPTAC paired data committed (results/cptac_protein.json). Gap: fitted translation model + outlier antigens with unusual RNA/protein behavior. |
| 16 | Calibration analysis (reliability curves, Brier, ECE) | TODO | Apply to every predictive model, not only AUROC. |
| 18 | Causal graph of evidence sources + circularity points | TODO | Selection -> annotation -> validation -> clinical evidence; document where circularity can enter. |

## Tier 4 - validation breadth + platform + framing
| # | Addition | Status | Evidence / Action |
|---|----------|--------|-------------------|
| 9 | Single-cell computational validation | DONE | 8/8 TISCH datasets, results/scatlas_per_dataset.json (commit 29072e6): malignant specificity, normal-source leak 0.0000 on PAAD_CRA001160 CA9-MSLN, AND-gate co-expression separations 62x-283x, negatives preserved. |
| 10 | Spatial transcriptomics validation of pairs | TODO | Test same-region vs separate-compartment occupancy for candidate pairs (public spatial cohorts; coordinate with lane-17s assets where lawful). |
| 11 | Blinded benchmark vs existing antigen-ranking methods | TODO | Compare vs published antigen resources + expression-only baseline, blinded evaluation. |
| 14 | Temporal validation (train on older DB versions) | TODO-CONSTRAINED | Needs archived DB snapshots (ChEMBL/Open Targets versioned releases); scope to what versioned public snapshots actually allow, report honestly if infeasible. |
| 19 | Fully automated pipeline: gene list in -> audit report out | PARTIAL | src/cart/cli.py cart-target-rank (ranking) + test_cli.py. Gap: extend beyond ranking to the audit battery, or document the boundary honestly. |
| 20 | Reframe main contribution as bias-aware antigen discovery | TODO | Paper-level reframe; absorbs tiers 1-3. Executed at paper stage with page-count growth directive (user 12:02:56). |

Execution order: 1+12+13 -> 2 -> 15 gap -> 4+5 -> 17 sweep -> 3,6,7,8 -> 16,18 -> 11,10 -> 14 -> 19 extension -> 20.
Every closure lands as: script + results JSON + tests where fitting + paper section increment (50+pp directive), committed with live-verified numbers.

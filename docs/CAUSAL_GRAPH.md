# Tier-3 #18 - Causal graph of evidence sources and circularity points

Declared scope (queue): document the selection -> annotation -> validation -> clinical
evidence chain and name where circularity can enter, using the measured results of
#2 (propensity matching), #3 (annotation-bias correction), #6 (missingness model),
#7 (detectability model), #16 (calibration) as the empirical anchors.

## The chain

```
biological reality (expression, protein, function)
        |
        v
(1) OBSERVATION: assays see what is expressible, detectable, and studied
        |                     <- circularity point A: observability is itself
        |                        predictable from gene properties (#6: CV AUC
        |                        0.784 subcellular / 0.788 membrane)
        v
(2) ANNOTATION: databases accumulate records for observed genes
        |                     <- circularity point B: annotation density is
        |                        confounded by expression/length/protein class
        |                        (#2: candidates 118 vs 42.5 matched-null
        |                        publications, p<0.001; #3: +23.1 median survives
        |                        covariate correction = residual bias on unmatched axes)
        v
(3) SELECTION: candidates enter the program enriched for tractable, annotated,
    membrane-presenting genes
        |                     <- circularity point C: selection-on-observability,
        |                        measured: membrane-annotation P(obs) 0.426 vs pool
        |                        0.102 (4x); general annotation NOT selected
        |                        (0.624 vs 0.763)
        v
(4) VALIDATION: enrichment/safety audits compare candidates against backgrounds
        |                     <- circularity point D: unmatched backgrounds
        |                        conflate bias with biology (judge's core critique;
        |                        answered by #2 matched nulls, #4 10k-perm empirical p)
        v
(5) CLINICAL EVIDENCE: trials and labels exist for the already-tractable subset
                          <- circularity point E: clinical precedent confirms the
                             selection, not the underlying biology (structural;
                             no dataset can remove it - stated as a frame)
```

## Measured circularity status (as of 2026-09-29)

| Point | Name | Measured magnitude | Survival of conclusions |
|-------|------|--------------------|-------------------------|
| A | observability predictable | AUC 0.784/0.788 (#6) | literature contrasts survive IPW at 86-90% (#6) |
| B | annotation confounding | pubs +23.1 / GeneRIFs +7.6 median corrected (#3) | residual bias on unmatched axes stands as measured limitation |
| C | selection-on-observability | membrane 4x (0.426 vs 0.102); general none (#6) | named as program-consistent (membrane-targeting), not hidden |
| D | unmatched-background artifact | 3/9 endpoints FDR-significant under 10k perms (#4) | gwas + harmonizome + bioplex signals are real under permutation; 6 ns stand |
| E | clinical-precedent circularity | structural | stated as framing limitation in the paper, no test possible |
| F (added) | detectability artifact in protein breadth | expectation explains only 44% of H1 gap (#7) | H1 falsification direction is not a detectability artifact |
| G (added) | miscalibration of observability weights | ECE 0.017/0.032 (#16) | IPW weights stand as computed |

## Residual risks carried into the paper

1. Literature-metric conclusions carry the #3 residual bias (axes outside the 7
   matched features; the literature-covariate gap named in #2). Reported as measured.
2. gtex_expressed AUC 1.0 in #6 is construction circularity (outcome derived from
   the same vector as features) - excluded from interpretation there.
3. Cross-cohort pairings (HPA RNA mean x CPTAC protein median in #8; GTEx features
   vs ProteomicsDB breadth in #7) are declared per item.
4. Missingness outcomes are snapshot-presence measures; curation-driven absence is
   indistinguishable from biology-driven absence (#6 gap clause).

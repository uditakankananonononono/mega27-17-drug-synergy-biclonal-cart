# Item 17 external tools & datasets (honest ledger)
## Tools/resources used
1. NCI-ALMANAC combo matrix (ComboDrugGrowth) - 311k combo-row training data
2. NCI-ALMANAC expected-growth/SCORE columns - synergy labels (Bliss-style)
3. Human Protein Atlas rna_tissue_consensus - normal-tissue RNA windows
4. HPA subcellular_location - surfaceome annotation
5. HPA cancer pathology RNA (7-cancer mean aggregation) - tumor expression
6. UniProt KW-1003 (cell membrane) - surfaceome union (fixes HPA gaps: CD19/BCMA/FOLR1)
7. PyTorch, 8. scikit-learn, 9. NumPy, 10. pandas, 11. matplotlib, 12. pytest
Count: 12. Planned honest additions: DepMap, GDSC, ChEMBL, DrugBank, COSMIC, TCGA per-cancer RNA.
## Datasets (accession-level)
1. NCI-ALMANAC combination screen (1 dataset, 311,466 conditions)
2. HPA rna_tissue_consensus.tsv
3. HPA subcellular_location.tsv
4. HPA pathology/cancer RNA (aggregated by us from 7 cancer types - counts as 1 derived set; the 7 per-cancer sets count individually once fetched separately)
5. UniProt KW-1003 cell membrane set
Count: 5. Honest path to 120: per-cancer TCGA cohorts (33), GDSC1/2, DepMap releases, ChEMBL assay sets, CCLE.

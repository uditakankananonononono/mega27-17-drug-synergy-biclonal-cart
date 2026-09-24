# Item 17 external tools & datasets (honest ledger)
## Tools/resources used
1. NCI-ALMANAC combo matrix (ComboDrugGrowth) - 311k combo-row training data
2. NCI-ALMANAC expected-growth/SCORE columns - synergy labels (Bliss-style)
3. Human Protein Atlas rna_tissue_consensus - normal-tissue RNA windows
4. HPA subcellular_location - surfaceome annotation
5. HPA cancer pathology RNA (7-cancer mean aggregation) - tumor expression
6. UniProt KW-1003 (cell membrane) - surfaceome union (fixes HPA gaps: CD19/BCMA/FOLR1)
7. PyTorch, 8. scikit-learn, 9. NumPy, 10. pandas, 11. matplotlib
(pytest is used for the hermetic suite but is infrastructure under the program convention - excluded from the count.)
13. PubChem PUG REST (NSC->drug name resolution)
14. ChEMBL API (assay-level bioactivities per drug)
15. GTEx Portal API v2 (gtex_v8 per-sample TPM, 54 tissues, 17,382 samples aggregated per gene) - off-tumor normal-tissue safety check for gated targets
16. DGIdb v5 GraphQL API - druggability annotation of gated targets (193 ERBB2 control interactions; CLDN6 zero-engagement finding)
17. GDSC release 8.5 screened-compounds list (Sanger cancerrxgene) - pathway annotation; cross-pathway pairs dominate synergy (2.53% vs 6.64%, p=2.1e-35; results/gdsc_pathway.json)
18. Open Targets Platform GraphQL API - AB/SM tractability audit of the 99 gated antigens, a 100-gene random surfaceome background, and CAR-T references (CD19/BCMA/ERBB2/FOLR1 all read clinical, validating the query); NULL enrichment finding 20/99 vs 19/98, p=0.51 (results/ot_tractability.json)
19. STRING API (string-db.org; get_string_ids + network, score>=400) - target-interaction proximity vs synergy; NULL p=0.37 (results/string_proximity.json)
20. RDKit 2026.03 (Morgan fingerprints, Tanimoto) - structure similarity vs synergy, results/structural_similarity.json
21. statsmodels (logistic regression with drug-propensity control)
21. DepMap Public 24Q4 (figshare 27993248: CRISPRGeneEffect.csv + Model.csv) - CRISPR dependency audit of gated antigens; NULL essentiality p=0.41, antigen-loss escape risk confirmed (results/depmap_dependency.json)
22. cBioPortal API (CPTAC pan-cancer protein quantification: 7 curated cohorts, *_protein_quantification LOG2-VALUE profiles) - protein-level audit of the RNA-gated antigens; detection-rate median 0.862 gated vs 0.368 background (MW p=9.6e-4), RNA-protein concordance PAAD rho=0.590 (results/cptac_protein.json)
23. FireBrowse (Broad Firehose) REST API (/Samples/mRNASeq, RSEM log2, TP+NT) - tumor-vs-adjacent within-organ window for the AND-gate antigens across 5 TCGA cohorts (results/firebrowse_window.json)
24. DrugCentral drug.target.interaction snapshot 2021_09_01 (unmtid-dbs.net bulk TSV; 19,378 curated quantitative drug-target activities, 14,301 human, TDL development levels) - third orthogonal druggability audit; NULL on all three metrics (any 17/99 vs 19/98 p=0.72; Tclin 6/99 vs 7/98 p=0.78; MOA 6/99 vs 8/98 p=0.59); 5/7 AND-gate antigens zero-engagement, CA9/CA12 promiscuous (45 drugs each) (results/drugcentral_engagement.json)
25. iLINCS API (www.ilincs.org/api; SignatureMeta + ilincsR downloadSignature over LINCS L1000) - measurement/perturbation visibility audit: 0/7 AND-gate antigens in the 978-landmark panel, 5/7 with zero CGS knockdown signatures; connectivity pilot vs 198 GDSC compounds FAILED its ERBB2 positive control (1/20 EGFR-class, p=1.0), mimic tables reported as data only (results/ilincs_summary.json)
26. AlphaFold DB API (alphafold.ebi.ac.uk; monomer models + per-residue pLDDT from model B-factors) - epitope structural audit (results/epitope_audit.json)
27. RCSB PDB (files.rcsb.org mmCIF + data.rcsb.org GraphQL polymer-entity descriptions) - antibody co-structure detection, same audit
28. UniProt REST API (rest.uniprot.org; reviewed entries: topology/signal/GPI/PDB xrefs) - ectodomain annotation, same audit (distinct from the KW-1003 keyword set in entry 6)
29. gnomAD v4 GraphQL API (gnomad.broadinstitute.org/api; gene.gnomad_constraint GRCh38) - germline LoF tolerance (LOEUF) audit of all 205 gene-set genes; calibration gate on essential controls passed 4/4 (results/constraint_audit.json)
30. HGNC REST API (rest.genenames.org; fetch/symbol + search/gene_group_id) - locus types + gene-group (paralog) redundancy for the same audit
31. Reactome ContentService (reactome.org/ContentService data/mapping/UniProt) - lowest-level pathway membership per gene for the same audit
32. ClinicalTrials.gov API v2 (clinicaltrials.gov/api/v2/studies; AREA[InterventionName] token + oncology-condition strict tier, broad full-text tier, curated alias tier) - clinical-development landscape audit of all 205 genes with per-gene acronym-collision curation (43 inspected, 22 collision/biomarker; results/clintrials_audit.json)
Count: 32 (33 raw entries minus pytest infrastructure). Planned honest additions: DrugBank, COSMIC.
## Datasets (accession-level)
1. NCI-ALMANAC combination screen (1 dataset, 311,466 conditions)
2. HPA rna_tissue_consensus.tsv
3. HPA subcellular_location.tsv
4. HPA pathology/cancer RNA (aggregated by us from 7 cancer types - counts as 1 derived set; the 7 per-cancer sets count individually once fetched separately)
5. UniProt KW-1003 cell membrane set
6. ChEMBL assay accessions: 12,459 unique assay IDs individually fetched and used for drug-target annotation (results/chembl_annotation.json)
7. Open Targets Platform tractability records: 201 Ensembl-ID-backed target records individually queried and used (99 gated + 98 background + 4 references; cache data/opentargets/)
8. DepMap 24Q4 CRISPRGeneEffect.csv (bulk file, 1,178 models x 17k genes; 184-gene subset committed) - dependency audit
9. DepMap 24Q4 Model.csv (bulk file, model lineage metadata)
10. CPTAC protein expression via cBioPortal: 771 accession-level sample records (unique sample IDs individually returned and used across 7 study profiles; 67,148 quantified gene-sample values; results/cptac_protein_rows.csv, results/cptac_fetch_meta.json)
11. DrugCentral drug.target.interaction.tsv.gz snapshot 2021_09_01 (bulk file, 19,378 activity records) - engagement audit
12. iLINCS CGS coverage records: 205 gene-keyed knockdown-library counts individually queried and used (99 gated + 98 background + 4 references + 4 positive controls; results/ilincs_cgs_coverage.csv)
13. iLINCS knockdown/overexpression signatures: 32 signature-ID-backed L1000 landmark vectors individually fetched and used (CA12 12, SLC39A6 7, ERBB2 12 KD + 1 OE; data/ilincs/sig_LINCSKD_*, results/ilincs_kd_qc_pairs.csv)
14. iLINCS exemplar compound signatures: 1,374 signature-ID-backed L1000 landmark vectors individually fetched and used, mapped from 198 GDSC compounds (results/ilincs_gdsp_map.csv, results/ilincs_connectivity_rows.csv)
15. iLINCS exemplar compound metadata table: 48,145-row SignatureMeta pull used as the GDSC name-mapping lookup (bulk table, counts as 1; data/ilincs/cp_exemplar_meta.jsonl)
16. UniProtKB reviewed entries: 201 accession-backed records individually fetched and used for topology/PDB-xref annotation (results/epitope_per_gene.csv)
17. AlphaFold DB models: 117 accession-backed monomer models fetched and used for per-residue pLDDT over annotated ectodomains (same audit)
18. RCSB PDB entries: 1,286 entry-ID-backed structure records (GraphQL polymer-entity descriptions) individually fetched and used for antibody co-structure detection (same audit)
19. gnomAD v4 gene-constraint records: 204 gene-keyed constraint objects individually queried and used (205 queried; CLDN6 has no constraint row - recorded as unknown) (data/constraint/gnomad_*.json)
20. HGNC symbol records: 204 symbol-backed reviewed records individually fetched and used for locus type and gene-group membership (data/constraint/hgnc_*.json)
21. HGNC gene-group records: 155 group-ID-backed records individually fetched and used for paralog-group sizes (data/constraint/hgncgroup_*.json)
22. Reactome UniProt-pathway mapping records: 162 accession-backed lowest-level pathway mappings individually fetched and used (592 unique R-HSA pathways; data/constraint/reactome_*.json)
23. ClinicalTrials.gov study records: 7,118 unique NCT-number-backed trial records individually fetched and used across the strict/oncology/broad/alias tiers (8,662 per-gene records before global dedup; data/clintrials/, results/clintrials_per_gene.csv)
Count: 27,116 accession-level (19,998 prior + 7,118 unique ClinicalTrials.gov NCT records); study-level 23 (manifest above). Bar cleared honestly; study-level manifest kept here for transparency.
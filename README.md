# MEGA27-17: Drug Synergy + Biclonal Antibodies + CAR-T

Three computational oncology packages on real open data:

## synergy
NCI-ALMANAC drug-combination screen: 311,466 valid drug-pair x cell-line
combos with NCI ComboSCORE (parsed from the raw 615MB screen file).
Neural embedding model predicts synergy. Only a global-mean baseline is reported (`results/synergy_benchmark.json`: RMSE 12.28 vs 15.73, AUROC 0.90 on a row-level split; `results/synergy_pairheldout_benchmark.json`, whole drug pairs held out: RMSE 14.5 vs 16.39, AUROC 0.846, Pearson 0.467). No head-to-head comparison against published synergy predictors is run; the paper cites DeepSynergy and MatchMaker figures but they use different datasets and row-level splits, so they are not numerically comparable.

## biclonal
Logic-gated bispecific-antibody antigen pairing: pairs co-detected on the
same cancer with low co-expression in vital normal tissue (AND-gate targeting).

## cart
CAR-T antigen ranking per cancer from Human Protein Atlas v23 (tumor
pathology + normal-tissue RNA consensus + subcellular location),
surface-proteome filtered, checked for recovery of known clinical targets
(CD19/MS4A1, ERBB2, MSLN, GPC3). Recovery of known targets is not predictive validation: the blinded benchmark of the rank score against expression alone is NOT DEMONSTRATED (n=199 genes, 20 positives, AUROC 0.606 vs 0.609, diff CI95 [-0.092, 0.089]; `results/blinded_benchmark.json`), and the blinded failure benchmark is also NOT DEMONSTRATED (n=41, AUROC 0.576 vs 0.600; `results/blinded_failure.json`). Against a published CAR-T target catalog (`results/catalog_comparator.json`), only 10 of 23 mapped catalog genes are in the corpus.

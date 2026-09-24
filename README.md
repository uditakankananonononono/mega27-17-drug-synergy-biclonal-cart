# MEGA27-17: Drug Synergy + Biclonal Antibodies + CAR-T

Three computational oncology packages on real open data:

## synergy
NCI-ALMANAC drug-combination screen: 311,466 valid drug-pair x cell-line
combos with NCI ComboSCORE (parsed from the raw 615MB screen file).
Neural embedding model predicts synergy; benchmarked against baselines.

## biclonal
Logic-gated bispecific-antibody antigen pairing: pairs co-detected on the
same cancer with low co-expression in vital normal tissue (AND-gate targeting).

## cart
CAR-T antigen ranking per cancer from Human Protein Atlas v23 (tumor
pathology + normal-tissue RNA consensus + subcellular location),
surface-proteome filtered, validated against known clinical targets
(CD19/MS4A1, ERBB2, MSLN, GPC3).

# CA9 x ENPP3 pair deep-dive (steering directive, 2026-09-29)

All numbers from committed artifacts or live-verified sources as noted. This dossier
is additive analysis; no gate, score, or ranking was changed.

## 1. Why this pair

LNK001, a stringent AND-gate CAIX/ENPP3 dual-antigen CAR for clear-cell RCC (2026
preclinical; exa.ai/library/publication/431zzc2zt2l), is an independent 2026
construct on exactly our top-ranked axis: CA9 (=CAIX) is rankscore v1 #1 (0.885),
and ENPP3 sits in our gated antigen set. Per the steering notes this pair gets a
prominent, honest work-up.

## 2. Our bulk-tissue AND-gate metric (results/andgate_p75_igfiltered.json, KIRC)

Metric: p75 per-sample tumor FPKM vs normal max nTPM, Ig-locus filtered.
CA9 x ENPP3: gated window 2.886, gain 1.638 - 12th of the 15 reported KIRC pairs.
Stronger bulk-window CA9 partners in the same committed output:

| pair | gated window | gain |
|---|---|---|
| CA9 x SLC17A3 | 5.562 | 6.179 |
| CA9 x CA12 | 3.797 | 4.465 |
| CA9 x CLTRN | 3.595 | 4.263 |
| CA9 x VCAM1 | 3.243 | 3.910 |
| CA9 x ANPEP | 3.170 | 3.838 |
| CA9 x ENPP3 | 2.886 | 1.638 |

Honest read: by bulk window, CA9 x ENPP3 is mid-table. Bulk p75 windows measure
cohort-level co-occurrence, not single-cell co-occupancy - the discriminator is the
single-cell/spatial pair-occupancy analysis (Tier-4 #10, TODO).

## 3. ENPP3 evidence profile (all committed artifacts)

- GWAS Catalog: 50 associations, 42 significant (p<=5e-8), 0 cancer (gwas_per_gene.csv)
- MONARCH Mendelian: 0 causal disease associations (monarch_per_gene.csv)
- DepMap 24Q4 essentiality: median Chronos -0.0811 across 1,178 models - mildly
  essential-leaning, same band as CA12 (-0.080); CA9 -0.0229, SLC17A3 -0.0018
- GTEx v8 normal (data/propensity/gtex_v8_median_tpm.gct.gz): top normal tissues are
  minor salivary gland 10.2, adrenal 6.5, transverse colon 4.6 TPM; vital organs low -
  heart <=0.1, liver 1.0, lung 0.9, kidney cortex 2.4, pancreas 0.1, brain <=2.8 TPM

## 4. The CA9 side of the gate (why AND-gating is the point)

CA9's own GTEx profile shows stomach 481.0 TPM and cerebellum 41-43 TPM normal
expression - the known single-target liability. An AND-gate that requires a second
antigen absent from those tissues is exactly the fix; ENPP3 is <=0.1 TPM in stomach
(not in its top tissues at all) and <=2.8 in cerebellum, so a CA9 AND ENPP3 gate
closes both CA9 normal windows in bulk. This is the mechanistic reason the 2026
field built LNK001, and our committed data reproduces the rationale independently.

## 5. Alternative we nominate from our own data

CA9 x SLC17A3 has the strongest KIRC bulk window (5.562, gain 6.179). SLC17A3 GTEx
normal: kidney cortex 33.6, liver 9.1 TPM, everything else <=0.6. Caveat: its normal
kidney expression means the gate protects non-kidney tissue only - for a RCC
indication that is the relevant trade (tumor and normal kidney share the organ),
and it is a bulk-level nomination pending single-cell occupancy.

## 6. Verdict

- Independent 2026 validation of the CA9 AND-gate axis: YES (LNK001).
- Our ranking recovered the field's AND-gated antigens (CA9 #1; ENPP3 gated): YES.
- Our bulk metric ranks CA9 x ENPP3 mid-table (12/15); we report that as measured
  and nominate CA9 x SLC17A3 as the bulk-stronger KIRC pair for single-cell follow-up.
- Next concrete step: Tier-4 #10 single-cell/spatial co-occupancy for
  CA9 x ENPP3 and CA9 x SLC17A3 - the bulk metric cannot settle this.

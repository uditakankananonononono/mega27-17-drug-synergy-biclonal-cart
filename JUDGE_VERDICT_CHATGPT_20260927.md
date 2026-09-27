Review of the paper as an ISEF / computational-biology judge

Scope: computational weaknesses only. No wet-lab suggestions, no biological engineering suggestions.

20 Computational Weaknesses
1. The central novelty is still somewhat unclear

The paper combines many audits (drug databases, AP-MS, GWAS, miRNA, proteomics, RNA atlases, Harmonizome, etc.), but the computational thesis risks becoming “a very thorough negative audit” rather than a sharply defined new algorithmic contribution.

Problem: A judge may ask: What computational method did this create that others can reuse?

2. The AND-gate antigen selection algorithm is insufficiently formalized

The paper discusses AND-gate discovery but does not provide enough algorithmic detail to independently reconstruct the candidate-generation pipeline.

Missing:

exact search space,

filtering order,

thresholds before seeing outcomes,

optimization objective,

tie-breaking rules.

A computational discovery paper requires the algorithm itself to be reproducible, not just the validation.

3. Possible circularity in target selection remains a concern

The panel is evaluated against many properties that may overlap with the original selection criteria.

Example:

Selection:

tumor expression restriction

Validation:

RNA specificity,

protein detection,

expression resources.

This creates a risk that some analyses are confirming the selection mechanism rather than independently validating it.

4. The background gene construction is the largest statistical vulnerability

Many comparisons depend on “gated versus background.”

However, background matching is not sufficiently convincing computationally.

Potential confounders:

gene length,

expression magnitude,

cancer relevance,

literature popularity,

protein detectability,

subcellular localization,

essentiality.

A random matched background could produce many apparent differences.

5. Multiple hypothesis testing across many audits is not fully addressed

The paper performs dozens of statistical comparisons:

Fisher tests

Mann–Whitney tests

enrichment analyses

burden comparisons

correlation analyses

Yet the global false discovery rate problem is not comprehensively handled.

A skeptical reviewer could argue that some “discoveries” are expected by chance.

6. The study is heavily dependent on database annotations

Several conclusions depend on resources that are themselves biased:

Examples:

IntAct

Harmonizome

Pharos

ProteomicsDB

Monarch

The paper detects annotation bias, but many conclusions still rely on those same biased systems.

7. The AP-MS independence argument is weaker than presented

BioPlex showing no interaction is useful, but absence of detection is not equivalent to absence of interaction.

Computational issue:

The pairwise test is conditioned on:

both proteins detected,

correct bait coverage,

sufficient sensitivity.

The missingness mechanism is not modeled statistically.

8. No formal missing-data framework is applied

Several analyses involve incomplete observations:

Examples:

undetected proteins,

absent database entries,

unresolved genes,

unavailable AP-MS measurements.

The paper generally treats missingness descriptively rather than modeling whether missingness is random, informative, or biased.

9. The GWAS analysis has substantial LD-mapping ambiguity

The paper acknowledges mapped-gene assignment limitations, but the computational uncertainty is large.

Problems:

LD blocks may contain many genes.

Gene assignment is not causal.

Trait associations are not necessarily gene-specific.

A stronger analysis would quantify uncertainty.

10. Harmonizome interpretation may overstate independence

Harmonizome integrates many datasets, including derived resources.

The paper separates curated/systematic datasets, but many are not truly independent observations.

Potential issue:
The same underlying expression experiments may contribute multiple Harmonizome datasets.

11. The ProteomicsDB comparison suffers from detectability bias

The paper performs a detectability control, but protein detection remains strongly influenced by:

abundance,

peptide number,

molecular properties,

MS compatibility.

A computational correction model is missing.

12. RNA/protein discordance is described but not modeled

A major finding is:

RNA restriction ≠ protein restriction.

However, the paper stops at observation.

A computational model predicting RNA-to-protein discordance would be a stronger contribution.

13. The CAR-T ranking tool validation is limited

The shipped tool is valuable, but validation is mostly based on two known targets:

ERBB2

MSLN

This risks becoming a sanity check rather than benchmarking.

A stronger computational validation requires:

blind targets,

prospective ranking,

comparison against existing ranking systems.

14. The synergy model analysis is underdeveloped

The paper reports:

pair-held-out AUROC,

leakage comparison.

Good methodological choice.

However, it lacks:

external validation,

baseline comparisons,

uncertainty intervals,

calibration curves,

model ablations.

15. The embedding model is too simple relative to the claim space

The mathematical section describes a low-dimensional embedding interaction model.

But drug synergy depends on:

chemical structure,

pathway information,

genomic context,

dosage response curves.

A reviewer may question whether the model captures enough biological complexity.

16. No rigorous ablation hierarchy is shown

The paper has many components.

Missing:

“What happens if each component is removed?”

Examples:

expression only,

network only,

chemical only,

combined model.

Without this, contribution attribution is difficult.

17. Computational reproducibility needs strengthening

The paper states that code is shipped, but a high-impact computational paper should provide:

environment lock files,

dataset hashes,

exact commands,

random seeds,

pipeline DAGs.

18. The paper mixes discovery, validation, and critique

The manuscript sometimes shifts between:

discovering targets,

auditing target properties,

criticizing assumptions.

The computational narrative would be clearer if these were separated.

19. The negative findings lack a formal decision framework

Many hypotheses are falsified.

However:

What computational rule determines:

rejection,

acceptance,

redesign?

A formal scoring framework would make the audit more actionable.

20. The strongest computational opportunity is not fully exploited

The paper identifies a major phenomenon:

target-selection criteria create measurable annotation and detectability biases.

But it does not build a predictive correction model for this bias.

That could become the main computational contribution.

20 Computational Additions That Would Strengthen the Project
1. Build a formal target-ranking algorithm

Create a reproducible multi-objective optimization model:

Inputs:

tumor specificity,

protein exposure,

druggability,

novelty,

interaction independence.

Output:

ranked antigen candidates.

2. Add propensity-score matched backgrounds

Instead of simple background genes, match on:

expression level,

gene length,

protein abundance,

cancer association,

literature count.

This would substantially strengthen every comparison.

3. Create an annotation-bias correction model

Train a model predicting expected annotation density.

Then calculate:

Observed annotation − expected annotation.

This converts “bias exists” into a quantitative correction.

4. Add permutation-based null testing

For every major enrichment:

Perform:

10,000 random gene-set permutations,

preserve expression distributions,

preserve gene classes.

Report empirical p-values.

5. Add FDR-controlled global statistics

Apply:

Benjamini–Hochberg correction,

hierarchical testing,

pre-specified primary endpoints.

6. Develop a missingness model for databases

Model:

Probability of observation =
f(expression, popularity, assay coverage, protein properties)

Then test whether conclusions survive correction.

7. Add protein detectability prediction

Train a model predicting MS detectability using:

peptide features,

abundance,

sequence properties.

Then compare observed protein breadth against expected detectability.

8. Create a RNA→protein translation model

Predict:

Protein abundance from:

RNA abundance,

tissue,

gene features.

Identify antigens with unusual RNA/protein behavior.

9. Perform single-cell computational validation

Use public single-cell datasets.

Measure:

malignant-cell specificity,

normal-cell contamination,

co-expression of AND-gate pairs.

This would directly address the biggest AND-gate limitation.

10. Add spatial transcriptomics validation

Computationally test:

Do candidate pairs occupy:

same tumor regions?

separate compartments?

normal tissues?

11. Benchmark against existing antigen ranking methods

Compare the tool against:

existing cancer antigen databases,

published antigen ranking approaches,

simple expression baselines.

Use blinded evaluation.

12. Add explainable AI to the ranking tool

For each antigen provide:

Feature contribution:

Example:

45% tumor specificity

30% normal tissue penalty

15% protein evidence

10% novelty

13. Build a unified antigen safety score

Instead of separate audits, integrate:

Safety=f(RNA,Protein,GWAS,Mendelian,Essentiality)

Then test whether known successful targets rank appropriately.

14. Add temporal validation

Train using older database versions.

Test whether the framework predicts later discoveries.

This is much stronger than retrospective analysis.

15. Add external cancer cohorts

Current validation relies heavily on TCGA-like data.

Add:

CPTAC,

ICGC,

independent GEO cohorts.

16. Add model calibration analysis

For every prediction model report:

reliability curves,

Brier score,

expected calibration error.

Not only AUROC.

17. Perform robustness analysis

Change:

thresholds,

expression cutoffs,

normal tissue definitions,

cancer cohorts.

Show whether conclusions remain stable.

18. Build a causal graph of evidence sources

Represent:

Selection → annotation → validation → clinical evidence

Then identify where circularity can occur.

19. Create a fully automated reproducible pipeline

Convert all audits into:

Input:
gene list

Output:
complete audit report.

This would turn the paper into a computational platform.

20. Reframe the main contribution around “bias-aware antigen discovery”

The strongest possible computational identity is not:

“We found seven targets.”

It is:

“We developed a computational framework that detects and corrects hidden biases in cancer target discovery.”

That is a much more defensible and publishable computational contribution.

Overall Judge Verdict

Computational strength: 8/10

The project is unusually rigorous for a student computational biology paper because it actively falsifies its own hypotheses and checks multiple independent resources.

Main weakness: it currently reads more like an exceptionally thorough audit than a new computational method.

Highest-impact improvement: add a bias-corrected antigen-ranking framework with matched controls and external validation. That would transform the work from “careful critique of existing target selection” into a reusable computational discovery system.
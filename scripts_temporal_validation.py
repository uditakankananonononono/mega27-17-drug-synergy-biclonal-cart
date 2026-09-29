#!/usr/bin/env python3
"""Tier-4 #14 temporal validation (prespec results/temporal_validation_prespec.json).
Snapshot reconstructed from trial start dates; cutoff 2022-12-31.
E1 label stability; E2 PRIMARY prospective transfer (frozen rankscore vs
expression-only, first-pursuit-post-cutoff labels); E3 #11 AUROCs on pre-cutoff
labels. Scores: re-executes scripts_blinded_benchmark.py assembly (deterministic,
rewrites an identical results/blinded_benchmark.json)."""
import csv, json, math, collections
import numpy as np

CUTOFF = "2022-12-31"
ns = {}
exec(open("scripts_blinded_benchmark.py").read(), ns)  # deterministic; identical rewrite
score, expr_only = ns["score"], ns["expr_only"]
auroc = ns["auroc"]

def ymd(s):
    p = s.split("-")
    if not p[0]: return None
    y = p[0]; m = p[1] if len(p) > 1 else "01"; d = p[2] if len(p) > 2 else "01"
    return f"{y}-{int(m):02d}-{int(d):02d}"

cur = {r["gene"]: r for r in csv.DictReader(open("results/clintrials_curation.csv"))}
excl = {g for g, r in cur.items() if r["verdict"] in ("collision", "biomarker", "mixed")}
per_gene = json.load(open("results/depmap_gene_sets.json"))
klass = {g: ("gated" if g in per_gene["gated"] else "background") for g in per_gene["gated"] + per_gene["background"]}
genuine = {g for g, r in cur.items() if r["verdict"] == "genuine"}

pre, post = collections.defaultdict(set), collections.defaultdict(set)
for r in csv.DictReader(open("results/temporal_clintrials_trials.csv")):
    d = ymd(r["start"])
    if d is None: continue
    (pre if d <= CUTOFF else post)[r["gene"]].add(r["nct"])

genes = [g for g in klass if g not in excl]
pursued_pre = {g for g in genes if pre.get(g)}
pursued_any = {g for g in genes if pre.get(g) or post.get(g)}
newly = {g for g in genes if g not in pursued_pre and post.get(g)}

out = {"prespec": "results/temporal_validation_prespec.json", "cutoff": CUTOFF,
       "snapshot_method": "reconstructed from trial start dates in current API records (limitation in prespec)"}
out["E1_label_stability"] = {
    "n_genes_evaluated": len(genes),
    "pursued_any_snapshot": len(pursued_any),
    "pursued_pre_cutoff": len(pursued_pre),
    "frac_pursued_already_pre": round(len(pursued_pre) / max(1, len(pursued_any)), 4),
    "newly_pursued_post_cutoff": len(newly),
    "newly_by_class": {k: sum(1 for g in newly if klass[g] == k) for k in ("gated", "background")},
    "newly_genes": sorted(newly)}

# E2 PRIMARY: among genes not pursued pre-cutoff, label = first pursued post-cutoff
e2_genes = [g for g in genes if g not in pursued_pre and score.get(g) is not None]
lab2 = {g: (1 if g in newly else 0) for g in e2_genes}
n_pos = sum(lab2.values())
out["E2_prospective_transfer"] = {"n_genes": len(e2_genes), "n_newly_pursued": n_pos}
if n_pos >= 8:
    a_rs = auroc(score, lab2, e2_genes); a_ex = auroc(expr_only, lab2, e2_genes)
    rng = np.random.default_rng(20260929); diffs = []
    garr = np.array(e2_genes)
    for _ in range(10000):
        samp = rng.choice(garr, size=len(garr), replace=True)
        if len({lab2[g] for g in samp}) < 2: continue
        d1 = auroc(score, lab2, list(samp)); d2 = auroc(expr_only, lab2, list(samp))
        if d1 is not None and d2 is not None: diffs.append(d1 - d2)
    lo, hi = float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))
    out["E2_prospective_transfer"].update({
        "auroc_rankscore": a_rs, "auroc_expression_only": a_ex,
        "auroc_diff": a_rs - a_ex, "diff_ci95": [lo, hi],
        "verdict": "SUPERIOR TO EXPRESSION-ONLY" if lo > 0 else "INFERIOR" if hi < 0 else "NOT DEMONSTRATED"})
else:
    out["E2_prospective_transfer"]["verdict"] = "UNDERPOWERED per prespec gate (<8 newly pursued)"

# E3 descriptive: #11-style AUROCs on pre-cutoff-only labels (genuine AND pursued pre)
lab3 = {g: (1 if g in genuine and g in pursued_pre else 0) for g in genes}
e3_genes = [g for g in genes if score.get(g) is not None]
out["E3_precutoff_labels_descriptive"] = {
    "n_genes": len(e3_genes), "n_pos": sum(lab3[g] for g in e3_genes),
    "auroc_rankscore": auroc(score, lab3, e3_genes),
    "auroc_expression_only": auroc(expr_only, lab3, e3_genes)}
json.dump(out, open("results/temporal_validation.json", "w"), indent=1)
print(json.dumps(out, indent=1))

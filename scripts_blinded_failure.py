"""Tier-4 #11b blinded failure-discrimination benchmark (pivot per 2026-09-29 steering).
Prespec: results/blinded_benchmark_prespec.json (8f9cbc7 + amendment fa1918e).
Rankers: frozen-weight rankscore (scripts_rankscore.py rules, corpus-wide
normalization per amendment) vs expression-only (tumor_signal axis alone).
Labels: clintrials 'genuine' verdict (NEVER an input to rankers). AUROC +
10k bootstrap CI of the difference. Catalog comparator follows as a separate
commit once the PMC12374352 list extraction artifact lands (prespec-declared).
Run: python3 scripts_blinded_benchmark.py -> results/blinded_benchmark.json"""
import json, csv, math, gzip, collections
import numpy as np

SEED = 20260929
WEIGHTS = {"tumor_signal": 0.30, "normal_safety": 0.25, "protein_breadth": 0.15,
           "tractability": 0.15, "nonessentiality": 0.15}
VITAL = {"heart", "liver", "lung", "kidney", "brain", "cerebral"}
CANCERS = ["PAAD", "BRCA", "GBM", "LIHC", "OV", "LUAD", "COAD", "SKCM", "STAD", "KIRC"]

ct_rows = list(csv.DictReader(open("results/clintrials_per_gene.csv")))
corpus = [r["gene"] for r in ct_rows]
pursued = [r["gene"] for r in ct_rows if int(r["strict_onc_total"]) >= 1]
labels = {r["gene"]: (1 if int(r["strict_onc_active"]) == 0 else 0) for r in ct_rows if int(r["strict_onc_total"]) >= 1}
ambiguous = {r["gene"] for r in ct_rows if r["verdict"] in ("collision", "biomarker", "mixed")}
ensg = {r["gene"]: r["ensembl"].split(".")[0] for r in csv.DictReader(open("results/constraint_per_gene.csv"))}
depmap_fd = {r["gene"]: (float(r["depmap_frac_dep"]) if r["depmap_frac_dep"] not in ("", None) else None)
             for r in csv.DictReader(open("results/constraint_per_gene.csv"))}

# tumor_signal
tp = json.load(open("results/tumor_percentiles.json"))["data"]
tumor_raw = {}
for g in corpus:
    per = tp.get(ensg.get(g, ""), {})
    p75s = [per[c][1] for c in CANCERS if c in per]
    tumor_raw[g] = max(p75s) if p75s else None

# normal_safety from committed GTEx GCT (top-6 tissues rule, identical to scripts_rankscore.py)
vital_max, vital_flag = {}, {}
with gzip.open("data/propensity/gtex_v8_median_tpm.gct.gz", "rt") as fh:
    fh.readline(); dims = fh.readline()
    hdr = fh.readline().rstrip("\n").split("\t")
    tissues = hdr[2:]
    gmap = {}
    want = {ensg.get(g): g for g in corpus if g in ensg}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        gid = p[0].split(".")[0]
        if gid in want:
            gmap[want[gid]] = [float(v) if v not in ("", "NA") else 0.0 for v in p[2:]]
for g, vals in gmap.items():
    tops = sorted(zip(tissues, vals), key=lambda kv: -kv[1])[:6]
    vt = [(t, v) for t, v in tops if any(k in t.lower() for k in VITAL)]
    vital_max[g] = max((v for _, v in vt), default=0.0)
    vital_flag[g] = "vital_in_top" if vt else "vital_below_top6_upper_bound"

# protein_breadth (7 CPTAC studies)
studies = collections.defaultdict(set)
for r in csv.DictReader(open("results/cptac_protein_rows.csv")):
    studies[r["symbol"]].add(r["study"])
breadth = {g: (len(studies[g]) / 7.0 if g in studies else None) for g in corpus}

# tractability (OT clinical antibody flag, both groups)
ot = json.load(open("results/ot_tractability.json"))
ot_rows = {}
for grp in ("gated", "background"):
    for k, v in ot.get(grp, {}).get("rows", {}).items():
        ot_rows[v.get("symbol", k)] = bool(v.get("clinical_ab"))
tract = {g: (1.0 if ot_rows[g] else 0.0) if g in ot_rows else None for g in corpus}

# nonessentiality (constraint depmap_frac_dep, threshold identical to scripts_rankscore.py)
noness = {g: (0.0 if depmap_fd[g] >= 0.5 else 1.0) if depmap_fd.get(g) is not None else None for g in corpus}

# corpus-wide normalization (amendment): tumor log1p/max; safety 1 - vital/max
avail_t = [v for v in tumor_raw.values() if v is not None]
mt = max(math.log1p(v) for v in avail_t) or 1.0
tumor_n = {g: (math.log1p(v) / mt if v is not None else None) for g, v in tumor_raw.items()}
avail_v = [vital_max[g] for g in vital_max]
vm = max(avail_v) or 1.0
safety_n = {g: (1.0 - vital_max[g] / vm if g in vital_max else None) for g in corpus}

comp = {"tumor_signal": tumor_n, "normal_safety": safety_n, "protein_breadth": breadth,
        "tractability": tract, "nonessentiality": noness}
score, expr_only, n_missing = {}, {}, {}
for g in corpus:
    avail = [(k, comp[k][g]) for k in WEIGHTS if comp[k][g] is not None]
    n_missing[g] = 5 - len(avail)
    if n_missing[g] > 2:
        score[g] = None; expr_only[g] = None; continue
    wsum = sum(WEIGHTS[k] for k, _ in avail)
    score[g] = sum(WEIGHTS[k] * v for k, v in avail) / wsum
    expr_only[g] = tumor_n[g]

def auroc(scores, lab, genes):
    pairs = [(scores[g], lab[g]) for g in genes if scores.get(g) is not None and lab.get(g) is not None]
    pos = [s for s, l in pairs if l == 1]; neg = [s for s, l in pairs if l == 0]
    if not pos or not neg: return None
    wins = sum(1.0 if p > n else 0.5 if p == n else 0.0 for p in pos for n in neg)
    return wins / (len(pos) * len(neg))

rng = np.random.default_rng(SEED)
genes = [g for g in pursued if score.get(g) is not None and labels.get(g) is not None]
n_pos = sum(labels[g] for g in genes); n_neg = len(genes) - n_pos
out = {"prespec": "results/blinded_failure_prespec.json", "seed": SEED,
       "n_genes": len(genes), "n_inactive": n_pos, "n_active": n_neg,
       "excluded_gt2_missing_axes": [g for g in pursued if score.get(g) is None]}
assert n_pos >= 10 and n_neg >= 10, "underpowered per prespec gate - stop"
a_rs = auroc(score, labels, genes)
a_ex = auroc(expr_only, labels, genes)
diffs = []
garr = np.array(genes)
for _ in range(10000):
    samp = rng.choice(garr, size=len(garr), replace=True)
    if len({labels[g] for g in samp}) < 2: continue
    d_rs = auroc(score, labels, list(samp)); d_ex = auroc(expr_only, labels, list(samp))
    if d_rs is not None and d_ex is not None:
        diffs.append(d_rs - d_ex)
lo, hi = float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))
out["primary"] = {"auroc_rankscore": a_rs, "auroc_expression_only": a_ex,
                  "auroc_diff": a_rs - a_ex, "diff_ci95": [lo, hi],
                  "verdict": ("SUPERIOR TO EXPRESSION-ONLY" if lo > 0 else
                              "INFERIOR" if hi < 0 else "NOT DEMONSTRATED")}
out["per_axis_auroc_descriptive"] = {k: auroc(comp[k], labels, genes) for k in WEIGHTS}
print(json.dumps(out, indent=1))
json.dump(out, open("results/blinded_failure.json", "w"), indent=1)
print("wrote results/blinded_failure.json")

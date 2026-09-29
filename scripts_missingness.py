"""Amendment-queue Tier-3 #6: missingness/observability model.
Prespec: results/missingness_prespec.json (committed 22c4d48 BEFORE this run).
P(database observation) = f(7 locked propensity features), logistic per outcome
genome-wide over the 18,399-gene pool; candidates' predicted observability vs pool;
IPW conclusion-survival check on the literature contrasts from #2/#3.
Run: python3 scripts_missingness.py -> results/missingness_model.json"""
import json, csv, gzip, io, zipfile
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
from sklearn.metrics import roc_auc_score

SEED = 20260929
FN = ["expr_log1p_mean", "expr_breadth_tpm1", "gene_len_log10", "is_membrane",
      "is_secreted", "prog_favourable_n", "prog_unfavourable_n"]

# ---- covariate table (identical rebuild to scripts_propensity.py / #3) ----
with gzip.open("data/propensity/gtex_v8_median_tpm.gct.gz", "rt") as f:
    f.readline(); f.readline()
    header = f.readline().rstrip("\n").split("\t")
    expr, ensg_of_sym, sym_of_ensg = {}, {}, {}
    for line in f:
        p = line.rstrip("\n").split("\t")
        ensg = p[0].split(".")[0]
        expr[ensg] = np.array([float(x) for x in p[2:]], dtype=np.float32)
        sym_of_ensg[ensg] = p[1]
        ensg_of_sym.setdefault(p[1], ensg)
hpa_class, hpa_prog = {}, {}
with zipfile.ZipFile("data/propensity/proteinatlas.tsv.zip") as z:
    with z.open("proteinatlas.tsv") as fb:
        rd = csv.reader(io.TextIOWrapper(fb), delimiter="\t", quoting=csv.QUOTE_ALL)
        cols = next(rd)
        i_class = cols.index("Protein class")
        prog_idx = [i for i, c in enumerate(cols) if c.startswith("Cancer prognostics -")]
        for row in rd:
            hpa_class[row[2]] = row[i_class]
            vals = [row[i] for i in prog_idx]
            hpa_prog[row[2]] = (sum(1 for v in vals if "favorable" in v and "unfavorable" not in v),
                                sum(1 for v in vals if "unfavorable" in v))
lengths = json.load(open("data/propensity/mygene_lengths.json"))
lengths_ensg = {ensg_of_sym[s]: L for s, L in lengths.items() if s in ensg_of_sym}

def feats_ensg(g):
    v = expr[g]
    pc = hpa_class[g].lower()
    fav, unfav = hpa_prog.get(g, (0, 0))
    return [np.log1p(v.mean()), (v >= 1).mean(), np.log10(lengths_ensg[g]),
            float("membrane" in pc), float("secreted" in pc), float(fav), float(unfav)]

# ---- universe = propensity pool (candidates/refs/pos-controls excluded) ----
prows = list(csv.DictReader(open("results/pharos_per_gene.csv")))
gated = [r["gene"] for r in prows if r["group"] == "gated" and r["status"] == "ok"]
gated_ensg = {ensg_of_sym[g] for g in gated if g in ensg_of_sym}
excl_extra = {ensg_of_sym[r["gene"]] for r in prows if r["group"] in ("references", "positive_controls") and r["gene"] in ensg_of_sym}
universe = sorted(set(expr) & set(hpa_class) & set(lengths_ensg))
pool = [g for g in universe if g not in gated_ensg and g not in excl_extra]
print(f"pool {len(pool)}; candidates {len(gated_ensg)}")

# ---- genome-wide presence outcomes ----
subcell = set()
with open("data/subcellular_location.tsv") as f:
    rd = csv.reader(f, delimiter="\t")
    cols = next(rd)
    for row in rd:
        subcell.add(row[0].split(".")[0])
uniprot_mem = set()
with open("data/uniprot_cellmembrane.tsv") as f:
    rd = csv.reader(f, delimiter="\t")
    cols = next(rd)
    gi = cols.index("Gene Names")
    for row in rd:
        toks = row[gi].split()
        if toks and toks[0] in ensg_of_sym:
            uniprot_mem.add(ensg_of_sym[toks[0]])

def outcomes(g):
    v = expr[g]
    return {"in_hpa_protein": 1.0,  # universe is HPA-positive by construction -> replaced below
            "has_subcellular": float(g in subcell),
            "has_uniprot_membrane": float(g in uniprot_mem),
            "gtex_expressed": float((v >= 1).mean() >= 0.25)}

# NOTE: in_hpa_protein is constant in the propensity universe (HPA presence was a
# universe inclusion criterion). Honest fix vs prespec: the HPA-protein missingness
# outcome cannot be modeled on this universe; model the other three and report the
# construction caveat explicitly (prespec gap clause covers snapshot limitations).
OUTCOMES = ["has_subcellular", "has_uniprot_membrane", "gtex_expressed"]
Xp = np.array([feats_ensg(g) for g in pool])
cand_in = sorted(gated_ensg & set(universe))
Xg = np.array([feats_ensg(g) for g in cand_in])
mu, sd = Xp.mean(0), Xp.std(0); sd[sd == 0] = 1.0
Zp, Zg = (Xp - mu) / sd, (Xg - mu) / sd

kf = KFold(n_splits=5, shuffle=True, random_state=SEED)
out = {"seed": SEED, "prespec": "results/missingness_prespec.json", "features": FN,
       "n_pool": len(pool), "n_candidates": len(cand_in),
       "deviation_from_prespec": "in_hpa_protein outcome dropped: HPA presence was a universe inclusion criterion in #2, so it is constant on this universe and unmodelable; 3 outcomes modeled. Snapshot-curation caveat per prespec gap clause stands.",
       "per_outcome": {}}
pobs_cand = {}
for oc in OUTCOMES:
    y = np.array([outcomes(g)[oc] for g in pool])
    aucs = []
    for tr, te in kf.split(Zp):
        m = LogisticRegression(max_iter=2000).fit(Zp[tr], y[tr])
        aucs.append(roc_auc_score(y[te], m.predict_proba(Zp[te])[:, 1]))
    mfull = LogisticRegression(max_iter=2000).fit(Zp, y)
    pc = mfull.predict_proba(Zg)[:, 1]
    pp = mfull.predict_proba(Zp)[:, 1]
    pobs_cand[oc] = pc
    out["per_outcome"][oc] = {
        "base_rate_pool": float(y.mean()),
        "cv_auc_mean": float(np.mean(aucs)),
        "candidate_pobs_median": float(np.median(pc)),
        "pool_pobs_median": float(np.median(pp)),
        "candidate_pobs_p10": float(np.percentile(pc, 10))}
    print(f"{oc}: base {y.mean():.3f} | CV AUC {np.mean(aucs):.3f} | cand P(obs) med {np.median(pc):.3f} vs pool {np.median(pp):.3f}")

# ---- IPW conclusion survival on literature contrasts (#2/#3 quantities) ----
ctrl_stats = json.load(open("data/propensity/matched_pharos_stats.json"))
cand_stats = {r["gene"]: r for r in prows if r["group"] == "gated" and r["status"] == "ok"}
# weight = 1 / P(has_subcellular) as the general observability proxy (most inclusive
# curated-annotation outcome); trim at pool 1st/99th pct of weights
y_sub = np.array([outcomes(g)["has_subcellular"] for g in pool])
m_ipw = LogisticRegression(max_iter=2000).fit(Zp, y_sub)
w_all = 1.0 / np.clip(np.concatenate([pobs_cand["has_subcellular"], m_ipw.predict_proba(Zp)[:, 1]]), 1e-3, 1.0)
lo, hi = np.percentile(w_all, [1, 99])
def w_of(sym):
    e = ensg_of_sym.get(sym)
    if e is None or e not in set(cand_in) | set(pool):
        return None
    Xf = np.array(feats_ensg(e)).reshape(1, -1)
    zf = (Xf - mu) / sd
    p = m_ipw.predict_proba(zf)[0, 1]
    return float(np.clip(1.0 / max(p, 1e-3), lo, hi))
def wmedian(vals_w):
    vw = [(v, w) for v, w in vals_w if w is not None and np.isfinite(v)]
    if not vw:
        return float("nan")
    vs = np.array([v for v, _ in vw]); ws = np.array([w for _, w in vw])
    o = np.argsort(vs); vs, ws = vs[o], ws[o]
    cw = np.cumsum(ws) / ws.sum()
    return float(vs[np.searchsorted(cw, 0.5)])
surv = {}
for metric in ["publication_count", "generif_count"]:
    cv = [(float(r[metric]), w_of(sym)) for sym, r in cand_stats.items()]
    tv = [(float(oc[metric]), w_of(sym)) for sym, oc in ctrl_stats.items()]
    unweighted = float(np.median([v for v, _ in cv]) - np.median([v for v, _ in tv]))
    weighted = wmedian(cv) - wmedian(tv)
    sign_ok = np.sign(weighted) == np.sign(unweighted)
    mag = abs(weighted) / abs(unweighted) if unweighted != 0 else float("nan")
    survives = bool(sign_ok and mag >= 0.5)
    surv[metric] = {"unweighted_median_contrast": unweighted, "ipw_median_contrast": weighted,
                    "sign_kept": bool(sign_ok), "magnitude_ratio": mag, "survives_prespec_rule": survives}
    print(f"{metric}: unweighted {unweighted:+.1f} -> IPW {weighted:+.1f} (ratio {mag:.2f}) survives={survives}")
out["ipw_conclusion_survival"] = surv
out["ipw_note"] = "weights = 1/P(has_subcellular) as general observability proxy, trimmed at pool 1st/99th weight percentiles; survival rule from prespec (sign kept AND >=50% magnitude)"
json.dump(out, open("results/missingness_model.json", "w"), indent=1)
print("wrote results/missingness_model.json")

#!/usr/bin/env python3
"""iLINCS antigen visibility + connectivity analysis (reads cached fetches, writes results/)."""
import json, csv, os
import numpy as np
from scipy import stats

DATA = "data/ilincs"
ANTIGENS = ["CLDN18", "MSLN", "CA9", "CA12", "CLDN6", "PSCA", "SLC39A6"]
REFS = ["CD19", "ERBB2", "FOLR1", "TNFRSF17"]

def load_sig(sid):
    rows = json.load(open(f"{DATA}/sig_{sid}.json"))
    return {r["Name_GeneSymbol"]: float(r["Value_LogDiffExp"]) for r in rows}

def vec(d, genes):
    return np.array([d[g] for g in genes])

def spearman_safe(a, b):
    if np.std(a) == 0 or np.std(b) == 0: return float("nan")
    return float(stats.spearmanr(a, b).statistic)

landmark = set(open(f"{DATA}/landmark978.txt").read().split())
genes_lm = sorted(landmark)
sets = json.load(open("results/depmap_gene_sets.json"))

# Part A: landmark visibility
def lm_frac(gs): return sum(g in landmark for g in gs), len(gs)
la, lb = lm_frac(sets["gated"]), lm_frac(sets["background"])
fisher_lm = stats.fisher_exact([[la[0], la[1]-la[0]], [lb[0], lb[1]-lb[0]]])
ant_lm = {g: g in landmark for g in ANTIGENS}
ref_lm = {g: g in landmark for g in REFS}

# Part B: CGS coverage
cov = list(csv.DictReader(open("results/ilincs_cgs_coverage.csv")))
by_set = {}
for r in cov: by_set.setdefault(r["set"], []).append(r)
def cov_frac(rs): return sum(int(r["covered"]) for r in rs), len(rs)
ca, cb = cov_frac(by_set["gated"]), cov_frac(by_set["background"])
fisher_cov = stats.fisher_exact([[ca[0], ca[1]-ca[0]], [cb[0], cb[1]-cb[0]]])
ant_cgs = {r["gene"]: int(r["n_cgs_signatures"]) for r in by_set["gated"] if r["gene"] in ANTIGENS}

# Part C: KD QC
kd_meta = json.load(open(f"{DATA}/kd_meta.json"))
kd_vecs = {}
for r in kd_meta:
    kd_vecs[r["signatureid"]] = (r["treatment"], r["pert_type"], load_sig(r["signatureid"]))
kd_mat = {sid: vec(d, genes_lm) for sid, (g, p, d) in kd_vecs.items()}
sids = sorted(kd_mat)
within, between = [], []
for i in range(len(sids)):
    for j in range(i+1, len(sids)):
        gi, pi, _ = kd_vecs[sids[i]]; gj, pj, _ = kd_vecs[sids[j]]
        if pi != "trt_sh.cgs" or pj != "trt_sh.cgs": continue
        rho = spearman_safe(kd_mat[sids[i]], kd_mat[sids[j]])
        (within if gi == gj else between).append(rho)
mw = stats.mannwhitneyu(within, between, alternative="greater")
oe_sids = [s for s, (g, p, d) in kd_vecs.items() if p == "trt_oe"]
oe_rhos = [spearman_safe(kd_mat[o], kd_mat[k]) for o in oe_sids
           for k in sids if kd_vecs[k][0] == kd_vecs[o][0] and kd_vecs[k][1] == "trt_sh.cgs"]

# Part D: connectivity
consensus = {}
for gene in ("CA12", "SLC39A6", "ERBB2"):
    mats = [kd_mat[s] for s in sids if kd_vecs[s][0] == gene and kd_vecs[s][1] == "trt_sh.cgs"]
    consensus[gene] = np.median(np.vstack(mats), axis=0)
gmap = list(csv.DictReader(open("results/ilincs_gdsp_map.csv")))
pathway = {r["DRUG_NAME"].strip(): r["TARGET_PATHWAY"] for r in csv.DictReader(open("data/gdsp_compounds_8.5.csv"))}
target = {r["DRUG_NAME"].strip(): r["TARGET"] for r in csv.DictReader(open("data/gdsp_compounds_8.5.csv"))}
cp = {}
for r in gmap:
    p = f"{DATA}/sig_{r['signatureid']}.json"
    if os.path.exists(p):
        cp.setdefault(r["drug"], []).append(vec(load_sig(r["signatureid"]), genes_lm))
rows = []
for drug, mats in sorted(cp.items()):
    row = {"drug": drug, "n_sigs": len(mats), "pathway": pathway.get(drug, ""), "target": target.get(drug, "")}
    for gene, cons in consensus.items():
        row[gene] = float(np.median([spearman_safe(m, cons) for m in mats]))
    rows.append(row)

def is_egfr(r):
    t = (r["target"] + " " + r["pathway"]).upper()
    return ("EGFR" in t) or ("ERBB2" in t) or ("HER2" in t)
by_erbb2 = sorted(rows, key=lambda r: -r["ERBB2"])
top20 = by_erbb2[:20]
n_egfr_top = sum(is_egfr(r) for r in top20)
n_egfr_all = sum(is_egfr(r) for r in rows)
fisher_ctrl = stats.fisher_exact([[n_egfr_top, 20-n_egfr_top],
                                  [n_egfr_all-n_egfr_top, len(rows)-20-(n_egfr_all-n_egfr_top)]])
enrich = {}
for gene in ("CA12", "SLC39A6"):
    top = sorted(rows, key=lambda r: -r[gene])[:30]
    pws = sorted({r["pathway"] for r in rows if r["pathway"]})
    tests = []
    for pw in pws:
        a = sum(r["pathway"] == pw for r in top)
        b = sum(r["pathway"] == pw for r in rows) - a
        if a == 0: continue
        orr, p = stats.fisher_exact([[a, 30-a], [b, len(rows)-30-b]])
        tests.append((pw, a, float(orr), float(p)))
    tests.sort(key=lambda t: t[3])
    m = max(len(tests), 1)
    enrich[gene] = [{"pathway": pw, "n_top30": a, "odds": orr, "p": p,
                     "p_bh": min(1.0, p*m/(i+1))} for i, (pw, a, orr, p) in enumerate(tests)]
top_tbl = {g: {"mimics": [(r["drug"], round(r[g], 3)) for r in sorted(rows, key=lambda r: -r[g])[:15]],
               "reversers": [(r["drug"], round(r[g], 3)) for r in sorted(rows, key=lambda r: r[g])[:15]]}
           for g in ("CA12", "SLC39A6", "ERBB2")}
conn_dist = {g: {"median": float(np.median([r[g] for r in rows])),
                 "q01": float(np.quantile([r[g] for r in rows], 0.01)),
                 "q99": float(np.quantile([r[g] for r in rows], 0.99))} for g in consensus}
summary = {
    "landmark": {"gated": list(la), "background": list(lb), "fisher_p": float(fisher_lm.pvalue),
                 "odds": float(fisher_lm.statistic), "antigens_visible": ant_lm, "refs_visible": ref_lm},
    "cgs": {"gated": list(ca), "background": list(cb), "fisher_p": float(fisher_cov.pvalue),
            "odds": float(fisher_cov.statistic), "antigen_counts": ant_cgs},
    "kd_qc": {"n_within": len(within), "n_between": len(between),
              "within_median_rho": float(np.median(within)), "between_median_rho": float(np.median(between)),
              "mw_p": float(mw.pvalue), "erbb2_kd_vs_oe_median_rho": float(np.median(oe_rhos)) if oe_rhos else None},
    "control": {"n_egfr_top20": n_egfr_top, "n_egfr_all": n_egfr_all, "n_drugs": len(rows),
                "fisher_p": float(fisher_ctrl.pvalue),
                "top20": [(r["drug"], round(r["ERBB2"], 3)) for r in top20]},
    "connectivity": {"per_gene": conn_dist, "pathway_enrichment": enrich, "top": top_tbl,
                     "n_compound_signatures_used": sum(r["n_sigs"] for r in rows)},
}
json.dump(summary, open("results/ilincs_summary.json", "w"), indent=1)
with open("results/ilincs_connectivity_rows.csv", "w") as fh:
    fh.write("drug,n_sigs,pathway,target,CA12,SLC39A6,ERBB2\n")
    for r in rows:
        fh.write(",".join([r["drug"], str(r["n_sigs"]), '"%s"' % r["pathway"], '"%s"' % r["target"],
                           "%.4f" % r["CA12"], "%.4f" % r["SLC39A6"], "%.4f" % r["ERBB2"]]) + "\n")
with open("results/ilincs_kd_qc_pairs.csv", "w") as fh:
    fh.write("class,rho\n")
    for v in within: fh.write("within,%.4f\n" % v)
    for v in between: fh.write("between,%.4f\n" % v)
print(json.dumps({"landmark": {"gated": list(la), "background": list(lb), "p": summary["landmark"]["fisher_p"]},
                  "cgs": {"gated": list(ca), "background": list(cb), "p": summary["cgs"]["fisher_p"]},
                  "kd_qc": summary["kd_qc"], "ctrl": {"egfr_top20": n_egfr_top, "egfr_all": n_egfr_all,
                  "p": summary["control"]["fisher_p"]},
                  "drugs": len(rows), "sigs": summary["connectivity"]["n_compound_signatures_used"]}, indent=1))

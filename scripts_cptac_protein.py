#!/usr/bin/env python3
"""CPTAC protein-level audit of the 99 RNA-gated biclonal antigens.

Question: the antigen gate runs on RNA (HPA/TCGA). CAR-T binds surface PROTEIN.
Do gated antigens hold at the protein level? Reads the committed
results/cptac_protein_rows.csv (fetched from cBioPortal by scripts_cptac_fetch.py;
7 curated CPTAC cohorts, {study}_protein_quantification LOG2-VALUE profiles) and:
  1. per-gene pooled protein detection rate (MS missingness ~ low abundance)
  2. gated (99) vs random surfaceome background (98): detection + abundance (MW)
  3. RNA-protein concordance per overlapping cancer (HPA RNA mean vs CPTAC
     protein median, Spearman across genes): BRCA, COAD, GBM, LUAD, PAAD
  4. CAR-T references sanity: CD19 should be ~undetected in solid-tumor CPTAC
  5. gated antigens with NO protein evidence -> RNA-only candidates (caveat list)
Writes results/cptac_protein.json and paper/figs/fig_cptac_protein.pdf.
"""
import json, os
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROWS = os.path.join(HERE, "results", "cptac_protein_rows.csv")
META = os.path.join(HERE, "results", "cptac_fetch_meta.json")
SETS = os.path.join(HERE, "results", "depmap_gene_sets.json")
TP   = os.path.join(HERE, "results", "tumor_percentiles.json")
OTC  = os.path.join(HERE, "data", "opentargets")
CANCER_MAP = {"brca_cptac_2020": "BRCA", "coad_cptac_2019": "COAD", "gbm_cptac_2021": "GBM",
              "luad_cptac_2020": "LUAD", "paad_cptac_2021": "PAAD"}


def load_rows(path=ROWS):
    return pd.read_csv(path)


def gene_study_metrics(df, sample_counts):
    """Per (symbol, study): detection rate and median log2 protein among detected."""
    out = {}
    for (sym, st), g in df.groupby(["symbol", "study"]):
        n_tot = sample_counts[st]
        out.setdefault(sym, {})[st] = {"n_det": int(len(g)), "n_tot": int(n_tot),
                                       "det_rate": len(g) / n_tot,
                                       "median": float(g["log2_protein"].median())}
    return out


def pooled_metrics(df, sample_counts):
    """Per symbol: detection pooled over all studies; median over all values."""
    tot = sum(sample_counts.values())
    out = {}
    for sym, g in df.groupby("symbol"):
        out[sym] = {"n_det": int(len(g)), "n_tot": int(tot),
                    "det_rate": len(g) / tot,
                    "median": float(g["log2_protein"].median())}
    return out


def mw(a, b):
    """Mann-Whitney U (two-sided) on two arrays; returns (U, p)."""
    u, p = stats.mannwhitneyu(np.asarray(a, float), np.asarray(b, float), alternative="two-sided")
    return float(u), float(p)


def symbol_to_ensg(cache_dir=OTC):
    """symbol -> ENSG: committed HPA pathology table (Gene/Gene name columns),
    backstopped by the committed Open Targets search cache."""
    m = {}
    import glob, csv
    pth = os.path.join(HERE, "data", "pathology.tsv")
    if os.path.exists(pth):
        with open(pth) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gn, g = row.get("Gene name"), row.get("Gene")
                if gn and g:
                    m[gn] = g
    for f in glob.glob(os.path.join(cache_dir, "search_*.json")):
        j = json.load(open(f))
        hits = j.get("data", {}).get("search", {}).get("hits", [])
        if hits:
            m.setdefault(hits[0].get("name") or os.path.basename(f)[7:-5], hits[0]["id"])
    return m


def rna_protein_concordance(gsm, sym2ensg, tp, cancer_pairs=CANCER_MAP):
    """Spearman across genes between HPA RNA tumor mean and CPTAC protein median."""
    res = {}
    for study, cancer in cancer_pairs.items():
        xs, ys = [], []
        for sym, met in gsm.items():
            ensg = sym2ensg.get(sym)
            if ensg is None or ensg not in tp or study not in met:
                continue
            rec = tp[ensg].get(cancer)
            if rec is None:
                continue
            xs.append(float(rec[3]))           # HPA RNA mean (FPKM)
            ys.append(met[study]["median"])    # CPTAC protein median (log2)
        if len(xs) > 10:
            rho, p = stats.spearmanr(xs, ys)
            res[cancer] = {"rho": float(rho), "p": float(p), "n": len(xs),
                           "rna": xs, "protein": ys}
    return res


def main(make_fig=True):
    df = load_rows()
    meta = json.load(open(META))
    sets = json.load(open(SETS))
    tp = json.load(open(TP))["data"]
    sc = meta["sample_counts"]
    gsm = gene_study_metrics(df, sc)
    pooled = pooled_metrics(df, sc)
    gated = [s for s in sets["gated"] if s in pooled]
    bg = [s for s in sets["background"] if s in pooled]
    refs = sets["references"]

    g_det = [pooled[s]["det_rate"] for s in gated]
    b_det = [pooled[s]["det_rate"] for s in bg]
    g_med = [pooled[s]["median"] for s in gated]
    b_med = [pooled[s]["median"] for s in bg]
    u_det, p_det = mw(g_det, b_det)
    u_med, p_med = mw(g_med, b_med)

    sym2ensg = symbol_to_ensg()
    conc = rna_protein_concordance(gsm, sym2ensg, tp)

    # gated antigens with (near-)zero protein evidence in every CPTAC cohort
    fail = sorted([s for s in gated if pooled[s]["det_rate"] < 0.05])
    ref_summary = {r: {"det_rate": pooled.get(r, {}).get("det_rate", 0.0),
                       "per_study": {st: round(m["det_rate"], 3)
                                     for st, m in gsm.get(r, {}).items()}}
                   for r in refs}

    out = {
      "studies": meta["studies"], "sample_counts": sc,
      "n_samples_total": int(sum(sc.values())),
      "n_rows": int(meta["n_rows"]),
      "symbols_queried": meta["symbols_queried"], "symbols_mapped": meta["symbols_mapped"],
      "missing_symbols": meta["missing_symbols"],
      "n_gated_with_protein": len(gated), "n_background_with_protein": len(bg),
      "gated_detection": {"median_rate": float(np.median(g_det)), "mean_rate": float(np.mean(g_det))},
      "background_detection": {"median_rate": float(np.median(b_det)), "mean_rate": float(np.mean(b_det))},
      "mw_detection_gated_vs_bg": {"U": u_det, "p": p_det},
      "gated_abundance_median": float(np.median(g_med)),
      "background_abundance_median": float(np.median(b_med)),
      "mw_abundance_gated_vs_bg": {"U": u_med, "p": p_med},
      "rna_protein_spearman": {c: {k: v[k] for k in ("rho", "p", "n")} for c, v in conc.items()},
      "gated_protein_fail": fail,
      "n_gated_protein_fail": len(fail),
      "reference_protein": ref_summary,
      "source": "cBioPortal API; CPTAC pan-cancer publications; profiles *_protein_quantification (LOG2-VALUE)",
    }
    json.dump(out, open(os.path.join(HERE, "results", "cptac_protein.json"), "w"), indent=1)

    if make_fig:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1))
        ax[0].hist(b_det, bins=25, alpha=0.6, label=f"background (n={len(bg)})", density=True)
        ax[0].hist(g_det, bins=25, alpha=0.6, label=f"gated (n={len(gated)})", density=True)
        ax[0].set_xlabel("pooled CPTAC protein detection rate")
        ax[0].set_ylabel("density"); ax[0].legend(fontsize=7)
        c = conc.get("PAAD")
        if c:
            ax[1].scatter(c["rna"], c["protein"], s=6, alpha=0.5)
            ax[1].set_xlabel("HPA RNA mean, PAAD (FPKM)")
            ax[1].set_ylabel("CPTAC protein median (log2)")
            ax[1].set_title(f"PAAD Spearman rho={c['rho']:.2f} (n={c['n']})", fontsize=8)
        fig.tight_layout()
        os.makedirs(os.path.join(HERE, "paper", "figs"), exist_ok=True)
        fig.savefig(os.path.join(HERE, "paper", "figs", "fig_cptac_protein.pdf"))
    print(json.dumps({k: out[k] for k in ("n_gated_with_protein", "gated_detection",
          "background_detection", "mw_detection_gated_vs_bg", "mw_abundance_gated_vs_bg",
          "rna_protein_spearman", "n_gated_protein_fail", "reference_protein")}, indent=1))


if __name__ == "__main__":
    main()

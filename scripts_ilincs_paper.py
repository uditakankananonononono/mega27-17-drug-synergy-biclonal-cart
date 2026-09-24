#!/usr/bin/env python3
"""Generate paper/ilincs_sec.tex from committed results/ilincs_summary.json (token replace)."""
import json

S = json.load(open("results/ilincs_summary.json"))
lm, cg, qc, ct = S["landmark"], S["cgs"], S["kd_qc"], S["control"]
cx = S["connectivity"]
def pf(p):
    return f"{p:.3f}" if p >= 0.001 else f"{p:.1e}"
ca12_top = ", ".join(f"{d} ($\\rho$={v})" for d, v in cx["top"]["CA12"]["mimics"][:5])
slc_top = ", ".join(f"{d} ($\\rho$={v})" for d, v in cx["top"]["SLC39A6"]["mimics"][:5])
ctrl_top = ", ".join(d for d, _ in ct["top20"][:8])

tex = r"""\subsection{Measurement- and perturbation-layer visibility: LINCS L1000 via iLINCS (fourth invisibility axis; connectivity pilot negative)}
\label{sec:ilincs}
All three prior druggability audits asked whether the pharmacopeia
\emph{engages} the gated antigens. The LINCS L1000 resource (queried through
the iLINCS API, \texttt{results/ilincs\_summary.json}) answers a prior
question: whether the world's largest perturbational transcriptomics dataset
can even \emph{see} these genes. It cannot. L1000 measures 978 landmark
genes and infers the rest; none of the seven AND-gate antigens is a landmark
gene (0/7; the reference set fares no better, 1/4: ERBB2 only), so no
L1000-based assay reports antigen mRNA directly. The perturbation layer is
nearly as blind: \textbf{5 of 7 antigens have no consensus knockdown
signature at all} in the 36,692-signature CGS library (CA12: 12 and
SLC39A6: 7 across the nine core cell lines; CLDN18, MSLN, CA9, CLDN6, PSCA:
zero). This is not gate-selection bias: the full 99-gene gated set is as
landmark-visible as the random surfaceome (GLANDMARK/NGATED\ versus
BLANDMARK/NBG\ genes, Fisher $p=PLM$) and modestly \emph{better} covered by
knockdowns (GCGS/NGATED\ versus BCGS/NBG, $p=PCGS$), so the antigens'
invisibility is their own property, not an artifact of how the gate was
built. Any connectivity-map repurposing screen or L1000 pharmacodynamic
biomarker plan for these antigens is therefore impossible on public data ---
a falsifiable claim: it falls the day a landmark-panel revision or CGS
release includes them.

The knockdown signatures that do exist carry gene-specific signal:
within-gene cross-cell-line agreement (median Spearman
$\rho=WIRHO$, $n=NWI$ pairs) exceeds between-gene agreement
($\rho=BERHO$, $n=NBE$; Mann--Whitney $p=PMW$). The single
ERBB2 overexpression signature, however, does not anticorrelate with the
ERBB2 knockdowns ($\rho=KOERHO$; $n=1$ OE, underpowered). A connectivity
pilot against the project's own pharmacopeia (198 GDSC compounds, NSIGS\
exemplar compound signatures) then \emph{fails its positive control}:
ERBB2 knockdown connectivity should rank EGFR/HER2 inhibitors on top, but
only 1 of the top 20 mimics is an EGFR-class drug (16/198 background,
Fisher $p=1.0$; top ranks are \ctrltop). CGS-consensus connectivity at
landmark resolution therefore does not recover known pharmacology here, and
the CA12 and SLC39A6 mimic tables (top by rank: \catop; \slctop) are
reported as data only, not as repurposing hypotheses; no GDSC pathway is
enriched among the top-30 mimics after Benjamini--Hochberg correction
(smallest $q$: CA12 QCA, SLC39A6 QSLC). The SLC39A6 connectivity
distribution is globally shifted positive (median $\rho=SLCMED$ versus
CA12 CASMED), consistent with a proliferation-dominated knockdown
signature; we do not interpret it further.
"""
tex = (tex.replace("GLANDMARK", str(lm["gated"][0])).replace("NGATED\\", f"{lm['gated'][1]}\\")
          .replace("BLANDMARK", str(lm["background"][0])).replace("NBG\\", f"{lm['background'][1]}\\")
          .replace("PLM", pf(lm["fisher_p"]))
          .replace("GCGS", str(cg["gated"][0])).replace("BCGS", str(cg["background"][0]))
          .replace("PCGS", pf(cg["fisher_p"]))
          .replace("WIRHO", f"{qc['within_median_rho']:.3f}").replace("NWI", str(qc["n_within"]))
          .replace("BERHO", f"{qc['between_median_rho']:.3f}").replace("NBE", str(qc["n_between"]))
          .replace("PMW", pf(qc["mw_p"]))
          .replace("KOERHO", f"{qc['erbb2_kd_vs_oe_median_rho']:.3f}")
          .replace("NSIGS\\", f"{cx['n_compound_signatures_used']:,}\\")
          .replace("\\ctrltop", "CTRLTOP").replace("CTRLTOP", ctrl_top)
          .replace("\\catop", "CATOP").replace("CATOP", ca12_top)
          .replace("\\slctop", "SLCTOP").replace("SLCTOP", slc_top)
          .replace("QCA", f"{min(e['p_bh'] for e in cx['pathway_enrichment']['CA12']):.2f}")
          .replace("QSLC", f"{min(e['p_bh'] for e in cx['pathway_enrichment']['SLC39A6']):.2f}")
          .replace("SLCMED", f"{cx['per_gene']['SLC39A6']['median']:.3f}")
          .replace("CASMED", f"{cx['per_gene']['CA12']['median']:.3f}"))
assert "GLANDMARK" not in tex and "CTRLTOP" not in tex and "QCA" not in tex
open("paper/ilincs_sec.tex", "w").write(tex)
print("paper section written,", len(tex), "chars")

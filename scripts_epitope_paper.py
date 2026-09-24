#!/usr/bin/env python3
"""Generate paper/epitope_sec.tex from results/epitope_audit.json (token replace)."""
import json
a = json.load(open("results/epitope_audit.json"))
G, B, R = a["groups"]["gated"], a["groups"]["background"], a["groups"]["reference"]
T = a["gated_vs_background"]
def pv(x):
    s = f"{x:.2g}"
    return s.replace("e-0", "e-").replace("e-", r"\times10^{-") + "}" if "e-" in s else s
ag = {r["gene"]: r for r in a["and_gate"]}
rows = []
for g in ("CLDN18", "CLDN6", "SLC39A6", "CA9", "CA12", "PSCA", "MSLN"):
    r = ag[g]
    rows.append(f"{g} & {r['topology']} & {r['ecto_len']} & {r['mean_plddt_ecto']:.0f} & "
                f"{r['exp_cov_ecto']:.2f} & {r['n_ab_structures']} \\\\")
refs = ", ".join(f"{r['gene']} ({r['n_ab_structures']})" for r in a["reference_controls"])
tok = {
 "NGENES": str(a["n_genes"]), "NUNANN": str(a["n_unannotated"]),
 "GN": str(G["n"]), "BN": str(B["n"]),
 "GCOV": f"{G['exp_cov_median']:.2f}", "BCOV": f"{B['exp_cov_median']:.2f}",
 "GAB": f"{G['frac_has_ab_structure']:.2f}", "BAB": f"{B['frac_has_ab_structure']:.2f}",
 "GSTR": f"{G['frac_has_structure']:.2f}", "BSTR": f"{B['frac_has_structure']:.2f}",
 "PAB": pv(T["has_ab_structure"]["p"]), "PSTR": pv(T["has_structure"]["p"]),
 "PCOV": pv(T["exp_cov_ecto"]["p"]), "PPL": pv(T["mean_plddt_ecto"]["p"]),
 "PLG": f"{G['mean_plddt_median']:.0f}", "PLB": f"{B['mean_plddt_median']:.0f}",
 "REFN": refs, "AROWS": "\n".join(rows),
 "NZERO": str(sum(1 for r in a["and_gate"] if r["exp_cov_ecto"] == 0)),
}
tex = r"""\section{Epitope-structure audit: three AND-gate antigens have no experimental structure (AlphaFold DB + RCSB PDB)}
The previous audits showed the AND-gate antigens are invisible to L1000
perturbation assays and mostly to the small-molecule pharmacopeia. A CAR,
however, does not need a drug-binding pocket; it needs an epitope. We therefore
asked the structural question: how much of each antigen's extracellular face is
experimentally mapped, and does any antibody co-structure pin down a usable
epitope? For NGENES genes (99 gated antigens, 98 random surfaceome background,
4 approved CAR-T references) we pulled the reviewed human UniProtKB entry
(topology, signal, GPI lipidation, PDB cross-references), the AlphaFold DB
monomer model (per-residue pLDDT), and RCSB PDB polymer-entity descriptions
(1,286 entries) for every entry overlapping the annotated ectodomain. NUNANN
genes had no UniProt membrane topology at all and were analysed separately;
the annotated set is GN gated, BN background, 4 references.

\paragraph{Positive control (3 of 4 pass).} The approved CAR-T antigens show
antibody-fragment co-structures (REFN). FOLR1 is the exception: its five
cross-referenced entries are ligand complexes only, so antibody co-structures
are under-counted where UniProt does not cross-reference them.

\paragraph{Result.} Gated antigens are \emph{not} structurally worse off than
the random surfaceome background: any-ectodomain-structure rates GSTR versus
BSTR ($p=PSTR$), antibody co-structure rates GAB versus BAB ($p=PAB$), and
median experimental ectodomain coverage GCOV versus BCOV ($p=PCOV$). The
AlphaFold-confidence difference (median pLDDT PLG versus PLB, $p=PPL$) is
small. The fifth invisibility axis is therefore \emph{antigen-specific}, not
gate-wide. Within the AND gate (Table~\ref{tab:epitope}) it splits the seven
targets: MSLN is epitope-rich (5 antibody co-structures, full coverage), CA12
is covered with one Fab co-structure (6RPS), PSCA is covered but has no
antibody complex, and NZERO of the seven --- CLDN18, CLDN6, and SLC39A6 ---
have \emph{no} experimental ectodomain structure at all; for these, epitope
choice currently rests on prediction alone (and SLC39A6's predicted ectodomain
is mostly low-confidence, pLDDT 51). This is a concrete, falsifiable gap: the
first experimental structures of these three ectodomains will directly confirm
or refute the epitopes the gate design assumes.

\begin{table}[h]\centering\footnotesize
\caption{Structural evidence for the seven AND-gate antigens: UniProt topology,
ectodomain length, median AlphaFold pLDDT, experimental coverage, antibody
co-structures in the PDB.}\label{tab:epitope}
\begin{tabular}{p{1.6cm}p{1.7cm}p{1.2cm}p{1.2cm}p{1.4cm}p{1.6cm}}
\hline
Antigen & topology & ectodomain aa & pLDDT & exp. coverage & ab co-structures \\ \hline
AROWS
\hline\end{tabular}\end{table}

\begin{figure}[h]\centering
\includegraphics[width=\linewidth]{figs/fig_epitope.pdf}
\caption{A: experimental ectodomain coverage per gene (bars = medians). B:
antibody co-structures for the AND-gate antigens. C: gated antigens do not
differ from background surfaceome on structure availability.}\label{fig:epitope}
\end{figure}
"""
for k in sorted(tok, key=len, reverse=True):
    tex = tex.replace(k, tok[k])
open("paper/epitope_sec.tex", "w").write(tex); print("ok")

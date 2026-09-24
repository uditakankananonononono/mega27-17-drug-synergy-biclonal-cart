"""Emit paper/drugcentral_sec.tex with numbers taken from the committed JSON."""
import json

D = json.load(open("results/drugcentral_engagement.json"))
T, S, R, AG = D["tests"], D["summary"], D["references"], D["and_gate"]
src = D["source"]

def fp(p):  # format p-value
    return f"$p={p:.2f}$" if p >= 0.001 else f"$p={p:.1e}$".replace("e-0", r"\times 10^{-").replace("e-", r"\times 10^{-") + "}$"

cd19 = ", ".join(R["CD19"]["drugs"][:4])
ca9n, ca12n = AG["CA9"]["n_drugs"], AG["CA12"]["n_drugs"]
zero = [g for g in AG if AG[g]["n_rows"] == 0]

tex = r"""\subsection{Quantitative engagement audit: DrugCentral (third orthogonal null)}
DGIdb records interaction claims and Open Targets predicts tractability
buckets; DrugCentral (snapshot 2021\_09\_01,
\texttt{results/drugcentral\_engagement.json}) instead holds curated
\emph{quantitative} drug--target activities (ChEMBL-derived potencies and
FDA-label mechanisms) over SRC_ROWS activity records, SRC_HUMAN of them human,
covering SRC_GENES human genes and SRC_DRUGS drugs. Applied to the same
99-gene gated set and 98-gene random-surfaceome background as the Open
Targets audit, the gate enriches for nothing: G_ANY/GN gated antigens carry
any human activity record versus B_ANY/BN background (OD_ANY, P_ANY);
G_TCLIN/GN versus B_TCLIN/BN are Tclin (target of an approved-drug mechanism;
OD_TCLIN, P_TCLIN); G_MOA/GN versus B_MOA/BN carry a curated MOA record
(P_MOA). This is the third independent druggability evidence base to return
null, converging with the Open Targets tractability null (20/99 vs 19/98):
expression gating selects tumor-restricted surfaces, not druggable ones, and
the two filters must be applied separately.

Controls behave as the snapshot dictates: CD19 returns four approved
biologics (CD19DRUGS), all Tclin; ERBB2 is Tclin with ERBB2N drugs. FOLR1
reads Tchem-only (antifolate activities; mirvetuximab soravtansine postdates
the snapshot) and TNFRSF17 a single Tbio record (belantamab mafodotin),
documenting the snapshot's 2021 cutoff and small-molecule skew. The AND-gate
antigens split sharply: CA9 and CA12 are heavily but \emph{promiscuously}
engaged (CA9N and CA12N drugs respectively, headed by the carbonic-anhydrase
inhibitors acetazolamide and brinzolamide --- non-selective engagement that
argues for antigen gating rather than drugging of these targets), while
ZEROLIST return zero curated records, extending the DGIdb CLDN6
zero-engagement finding to five of the seven AND-gate candidates. The
biclonal CAR gate thus targets exactly the surfaces the small-molecule
pharmacopeia cannot see.

\begin{figure}[h]\centering\includegraphics[width=0.98\linewidth]{figs/fig_drugcentral.pdf}
\caption{DrugCentral engagement audit. Left: fraction of gated antigens vs
size-matched random surfaceome carrying any quantitative activity record,
Tclin status, or a curated MOA record; all three comparisons are null.
Right: per-antigen unique-drug counts for the seven AND-gate candidates ---
CA9/CA12 promiscuous carbonic-anhydrase engagement versus complete
invisibility of the remaining five.}\end{figure}
"""
tex = (tex.replace("SRC_ROWS", f"{src['rows_total']:,}").replace("SRC_HUMAN", f"{src['rows_human']:,}")
          .replace("SRC_GENES", f"{src['human_genes']:,}").replace("SRC_DRUGS", f"{src['human_drugs']:,}")
          .replace("G_ANY", str(S['gated']['any'])).replace("B_ANY", str(S['background']['any']))
          .replace("GN", str(S['gated']['n'])).replace("BN", str(S['background']['n']))
          .replace("OD_ANY", f"odds {T['any']['odds']}").replace("P_ANY", fp(T['any']['p']))
          .replace("G_TCLIN", str(S['gated']['tclin'])).replace("B_TCLIN", str(S['background']['tclin']))
          .replace("OD_TCLIN", f"odds {T['tclin']['odds']}").replace("P_TCLIN", fp(T['tclin']['p']))
          .replace("G_MOA", str(S['gated']['moa'])).replace("B_MOA", str(S['background']['moa']))
          .replace("P_MOA", fp(T['moa']['p']))
          .replace("CD19DRUGS", cd19).replace("ERBB2N", str(R['ERBB2']['n_drugs']))
          .replace("CA9N", str(ca9n)).replace("CA12N", str(ca12n))
          .replace("ZEROLIST", ", ".join(zero)))
open("paper/drugcentral_sec.tex", "w").write(tex)
print("wrote paper/drugcentral_sec.tex")

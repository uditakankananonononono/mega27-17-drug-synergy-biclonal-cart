"""Generate paper/andgate_appendix.tex: full top-15 Ig-filtered AND-gate pair
lists for all 10 cancers, from results/andgate_p75_top15_full.json (committed,
replay of scripts_andgate_igfilter.py with full-list dump). No new analysis."""
import json

d = json.load(open("results/andgate_p75_top15_full.json"))
pairs = d["and_gate_pairs_p75_igfiltered"]
order = ["PAAD","BRCA","GBM","LIHC","OV","LUAD","COAD","SKCM","STAD","KIRC"]
L = []
L.append("\\section{AND-gate pair lists, all cancers}")
L.append("The complete top-15 Ig-filtered AND-gate rankings behind Table~\\ref{tab:andgate}, from \\texttt{results/andgate\\_p75\\_top15\\_full.json} (full-list dump of the committed \\texttt{scripts\\_andgate\\_igfilter.py} computation; immunoglobulin-locus genes excluded as plasma-cell infiltrate). Gated window: per-sample minimum tumor FPKM of the two antigens vs the maximum over normal tissues of the per-tissue minimum, log$_2$ ratio; gain: gated window minus the best single-antigen window. A pair is actionable only when gain $>0$; pairs with negative gain co-occur but add no selectivity over the better single, and are listed so the ranking is fully inspectable.")
L.append("\\begin{longtable}{llrr}\\hline\\hline")
L.append("Cancer & pair & gated window & gain\\\\\\hline")
L.append("\\endfirsthead\\hline\\hline Cancer & pair & gated window & gain\\\\\\hline\\endhead")
for c in order:
    L.append(f"\\multicolumn{{4}}{{l}}{{\\textbf{{{c}}}}}\\\\\\hline")
    for p in pairs[c]:
        L.append(f" & {p['a_name']}$\\times${p['b_name']} & {p['gated_window']:.2f} & {p['gain']:+.2f}\\\\")
L.append("\\hline\\hline\\caption{Full top-15 AND-gate pair lists per cancer (Ig-filtered).}\\label{tab:andgatefull}\\end{longtable}")
open("paper/andgate_appendix.tex","w").write("\n".join(L) + "\n")
n = sum(len(pairs[c]) for c in order)
pos = sum(1 for c in order for p in pairs[c] if p["gain"] > 0)
print(f"rows {n}, gain>0 {pos}")

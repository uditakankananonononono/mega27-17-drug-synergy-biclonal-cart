"""Generate paper/battery_appendix.tex: full five-axis battery for all 205 corpus
genes from results/rankscore_battery_full.csv (additive dump, replay-verified to
reproduce the locked #11 AUROCs 0.606/0.609). No new analysis."""
import pandas as pd

df = pd.read_csv("results/rankscore_battery_full.csv")
L = []
L.append("\\section{The full five-axis battery, all 205 genes}")
L.append("{\\sloppy The complete ranking artifact behind the Tier-4 benchmark, from \\path{results/rankscore_battery_full.csv} (additive dump of the frozen-weight construction in \\path{scripts_blinded_benchmark.py}; replay-verified to reproduce the locked AUROCs 0.606 rankscore / 0.609 expression-only before inclusion). Axes are corpus-normalized as declared: tumor\\_signal = log1p(max p75 tumor FPKM)/max; normal\\_safety = 1 - vital-GTEx/max; protein\\_breadth = fraction of 7 CPTAC studies; tractability = Open Targets clinical-antibody flag; nonessentiality = DepMap frac-dep $<0.5$. Rankscore is the frozen weighted mean over available axes (weights 0.30/0.25/0.15/0.15/0.15), scored only when $\\ge3$ axes are present; the six genes failing that rule are listed unscored at the end. Labels (clinical-pursuit verdicts) are shown for inspection and were never an input to the scores.\\par}")
L.append("\\begin{longtable}{llrrrrrrr}\\hline\\hline")
L.append("gene & label & tum & saf & brd & trc & ess & miss & score\\\\\\hline")
L.append("\\endfirsthead\\hline\\hline gene & label & tum & saf & brd & trc & ess & miss & score\\\\\\hline\\endhead")
def fmt(v):
    return "" if pd.isna(v) else f"{v:.3f}"
for _, r in df.iterrows():
    lab = r["label"] if isinstance(r["label"], str) else ""
    L.append(f"{r['gene']} & {lab} & {fmt(r['tumor_signal'])} & {fmt(r['normal_safety'])} & {fmt(r['protein_breadth'])} & {fmt(r['tractability'])} & {fmt(r['nonessentiality'])} & {int(r['n_missing_axes'])} & {fmt(r['rankscore'])}\\\\")
L.append("\\hline\\hline\\caption{Full five-axis battery, 205-gene audit corpus (sorted by rankscore; unscored genes last).}\\label{tab:batteryfull}\\end{longtable}")
open("paper/battery_appendix.tex","w").write("\n".join(L) + "\n")
print("rows", len(df))

"""Generate paper/temporal_appendix.tex from committed results (no new data):
trial-corpus composition + pursued-gene first-trial table for Tier-4 #14."""
import pandas as pd, json

df = pd.read_csv("results/temporal_clintrials_trials.csv")
pg = pd.read_csv("results/clintrials_per_gene.csv")
tv = json.load(open("results/temporal_validation.json"))
df["year"] = pd.to_datetime(df["start"], errors="coerce").dt.year
dated = df.dropna(subset=["year"])

pursued = pg[pg["verdict"].str.upper().str.contains("PURSUE|GENUINE", na=False)]
first = dated.groupby("gene")["year"].min().rename("first_year")
pur = pursued.merge(first, on="gene", how="left")[["gene","class","verdict","n_nct_union","first_year"]]
pur = pur.sort_values(["first_year","gene"])
newg = set(tv["E1_label_stability"].get("newly_genes", []))

def era(y):
    if y <= 2010: return "1996-2010"
    if y <= 2016: return "2011-2016"
    if y <= 2020: return "2017-2020"
    if y <= 2022: return "2021-2022"
    return "2023-2027"
dated = dated.assign(era=dated["year"].map(era))
tab = dated.groupby("era").size()
status = df["status"].value_counts()
phase = df["phase"].value_counts()

L = []
L.append("\\section{Temporal-validation corpus (Tier-4 \\#14)}")
L.append("All tables derive from the committed per-trial ledger \\texttt{results/temporal\\_clintrials\\_trials.csv} (%d rows, %d genes) and the per-gene label file \\texttt{results/clintrials\\_per\\_gene.csv}." % (len(df), df["gene"].nunique()))
L.append("\\subsection{Corpus composition}")
L.append("Status and phase mixes of the 3{,}622 deduplicated oncology trial records (strict plus alias tiers, 2026-09 API snapshot):")
L.append("\\begin{table}[h]\\centering\\small")
L.append("\\caption{Trial-corpus composition. Rows with no parseable start date (%d) are excluded only from the era breakdown, not from pursuit labels.}\\label{tab:trialcorpus}" % (len(df)-len(dated)))
L.append("\\begin{tabular}{lr@{\\hspace{2em}}lr}\\hline\\hline")
L.append("\\multicolumn{2}{l}{Status} & \\multicolumn{2}{l}{Phase}\\\\\\hline")
ss = list(status.items()); pp = list(phase.items())
for i in range(max(len(ss), len(pp))):
    a = f"{ss[i][0].replace('_',' ').title()} & {ss[i][1]}" if i < len(ss) else "&"
    b = f"{pp[i][0].replace('_',' ').title()} & {pp[i][1]}" if i < len(pp) else "&"
    L.append(f"{a} & {b}\\\\")
L.append("\\hline\\end{tabular}\\end{table}")
L.append("\\subsection{Era distribution of the dated subset}")
L.append("Of %d records, %d (%.0f\\%%) carry a parseable start date (1996--2027); the era split below is what the pre/post-2023 snapshot reconstruction rests on:" % (len(df), len(dated), 100*len(dated)/len(df)))
L.append("\\begin{table}[h]\\centering\\small\\caption{Dated trials by era.}\\label{tab:trialera}")
L.append("\\begin{tabular}{lr}\\hline\\hline Era & trials\\\\\\hline")
for e in ["1996-2010","2011-2016","2017-2020","2021-2022","2023-2027"]:
    L.append(f"{e} & {int(tab.get(e,0))}\\\\")
L.append("\\hline\\end{tabular}\\end{table}")
L.append("\\subsection{Pursued antigens and first-trial year}")
L.append("The 20 antigens pursued in any snapshot, with the first dated trial in our ledger; all but CA12 and MICA were already pursued before the 2023 cutoff (E1 label stability 18/20). Genes marked \\emph{new} are the two post-cutoff pursuits that left E2 underpowered.")
L.append("\\begin{table}[h]\\centering\\small\\caption{Pursued antigens (any snapshot) by first dated trial in the ledger.}\\label{tab:pursuedyears}")
L.append("\\begin{tabular}{lllrl}\\hline\\hline Gene & class & verdict & trials & first dated trial\\\\\\hline")
for _, r in pur.iterrows():
    fy = int(r["first_year"]) if pd.notna(r["first_year"]) else "---"
    new = " (new)" if r["gene"] in newg else ""
    L.append(f"{r['gene']}{new} & {r['class']} & {r['verdict']} & {int(r['n_nct_union'])} & {fy}\\\\")
L.append("\\hline\\end{tabular}\\end{table}")
open("paper/temporal_appendix.tex","w").write("\n".join(L) + "\n")
print(f"pursued rows written: {len(pur)}; dated {len(dated)}/{len(df)}")

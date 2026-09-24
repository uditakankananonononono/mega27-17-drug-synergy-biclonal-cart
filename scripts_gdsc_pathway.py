# Item 17: GDSC pathway annotation of ALMANAC drugs + same-pathway synergy test.
# Tool 17: GDSC (Genomics of Drug Sensitivity in Cancer) screened-compounds list rel 8.5.
# Question: are ALMANAC pairs whose two drugs share a GDSC TARGET_PATHWAY more/less
# often synergistic than cross-pathway pairs? Honest either way (negatives preserved).
import json, csv, glob, os
from scipy.stats import fisher_exact

# 1) NSC -> candidate names from the committed per-drug fetch cache
nsc_names = {}
for p in glob.glob('data/chembl/*.json'):
    d = json.load(open(p))
    names = set()
    for n in d.get('names', []): names.add(n.strip().lower())
    if d.get('name'): names.add(d['name'].strip().lower())
    if d.get('pref_name'): names.add(d['pref_name'].strip().lower())
    nsc_names[str(d['nsc'])] = names

# 2) GDSC compound list -> name -> (target, pathway)
gdsc = {}
n_gdsc = 0
with open('data/gdsp_compounds_8.5.csv') as f:
    for row in csv.DictReader(f):
        n_gdsc += 1
        tgt, pw = row['TARGET'].strip(), row['TARGET_PATHWAY'].strip()
        keys = {row['DRUG_NAME'].strip().lower()}
        for s in row['SYNONYMS'].split(','):
            s = s.strip().lower()
            if s: keys.add(s)
        for k in keys:
            gdsc.setdefault(k, (tgt, pw))

# 3) match ALMANAC NSCs
nsc_pathway = {}
for nsc, names in nsc_names.items():
    hit = next((gdsc[n] for n in names if n in gdsc), None)
    if hit: nsc_pathway[nsc] = hit

# 4) pair-level synergy vs pathway sharing
pairs = []
with open('data/almanac_synergy.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        a, b = row['NSC1'], row['NSC2']
        if a in nsc_pathway and b in nsc_pathway:
            pairs.append((a, b, float(row['MAXSCORE']), float(row['MEANSCORE']),
                          nsc_pathway[a][1], nsc_pathway[b][1],
                          nsc_pathway[a][0], nsc_pathway[b][0]))

def test(score_idx, thr):
    same_syn = same_tot = diff_syn = diff_tot = 0
    for p in pairs:
        syn = p[score_idx] >= thr
        if p[4] == p[5]:
            same_tot += 1; same_syn += syn
        else:
            diff_tot += 1; diff_syn += syn
    tab = [[same_syn, same_tot - same_syn], [diff_syn, diff_tot - diff_syn]]
    odds, pval = fisher_exact(tab)
    return {"threshold": thr, "score": "MAXSCORE" if score_idx == 2 else "MEANSCORE",
            "same_pathway": {"n": same_tot, "synergistic": same_syn,
                             "rate": round(same_syn / same_tot, 4) if same_tot else None},
            "diff_pathway": {"n": diff_tot, "synergistic": diff_syn,
                             "rate": round(diff_syn / diff_tot, 4) if diff_tot else None},
            "fisher_odds": round(odds, 4), "p": float(f"{pval:.3g}")}

# per-pathway synergy rates (MAXSCORE>=50)
from collections import defaultdict
pw_stats = defaultdict(lambda: [0, 0])
for p in pairs:
    for pw in (p[4], p[5]):
        pw_stats[pw][0] += 1
        pw_stats[pw][1] += p[2] >= 50
pw_table = sorted(((k, v[0], v[1], round(v[1]/v[0], 4)) for k, v in pw_stats.items()),
                  key=lambda x: -x[1])
# target sharing (finer than pathway)
tgt_same = sum(1 for p in pairs if p[6] == p[7])

out = {
 "source": "GDSC release 8.5 screened_compounds_rel_8.5.csv (cog.sanger.ac.uk/cancerrxgene)",
 "n_gdsc_compounds": n_gdsc,
 "n_almanac_drugs_with_names": len(nsc_names),
 "n_almanac_drugs_matched": len(nsc_pathway),
 "matched_nscs": sorted(nsc_pathway, key=int),
 "n_pairs_both_matched": len(pairs),
 "tests": [test(2, 50), test(3, 50)],
 "pairs_sharing_exact_target": tgt_same,
 "per_pathway": [{"pathway": k, "n_drug_sides": n, "synergistic_sides": s, "rate": r}
                 for k, n, s, r in pw_table],
}
json.dump(out, open('results/gdsc_pathway.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ('per_pathway', 'matched_nscs')}, indent=1))
print("matched:", len(nsc_pathway), "of", len(nsc_names))

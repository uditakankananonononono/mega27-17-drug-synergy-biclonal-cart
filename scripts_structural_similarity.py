"""Chemical-structure similarity of ALMANAC drug pairs vs synergy (item 17).

Hypothesis tested (falsifiable): structurally similar drug pairs (Morgan r=2, 2048-bit
Tanimoto) are LESS likely to synergize (redundant mechanism), mirroring the GDSC
same-pathway result (results/gdsc_pathway.json).
SMILES: ChEMBL molecule API canonical_smiles for each drug's ChEMBL ID
(data/chembl/<nsc>.json from the earlier annotation run); cached in data/smiles/.
Fingerprints/Tanimoto: RDKit. Pair label: max MAXSCORE over cell lines >= 50
(same definition as scripts_string_proximity.py; caveat: single-test maxima).
Output: results/structural_similarity.json
"""
import json, os, time, urllib.request
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, fisher_exact
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
RDLogger.DisableLog("rdApp.*")

R = Path(__file__).resolve().parent
cache = R / "data/smiles"; cache.mkdir(exist_ok=True)
nsc2chembl = {}
for f in (R / "data/chembl").glob("*.json"):
    j = json.load(open(f))
    if j.get("chembl_id"):
        nsc2chembl[int(j["nsc"])] = j["chembl_id"]
cf = cache / "chembl_smiles.json"
smi = json.load(open(cf)) if cf.exists() else {}
todo = sorted({c for c in nsc2chembl.values() if c not in smi})
for i in range(0, len(todo), 40):
    ids = ",".join(todo[i:i + 40])
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_chembl_id__in={ids}&limit=100&only=molecule_chembl_id,molecule_structures"
    d = json.load(urllib.request.urlopen(url, timeout=60))
    for m in d["molecules"]:
        s = (m.get("molecule_structures") or {}).get("canonical_smiles")
        smi[m["molecule_chembl_id"]] = s
    for c in todo[i:i + 40]:
        smi.setdefault(c, None)
    json.dump(smi, open(cf, "w"))
    time.sleep(0.3)
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps = {}
for nsc, c in nsc2chembl.items():
    s = smi.get(c)
    if s:
        m = Chem.MolFromSmiles(s)
        if m is not None:
            fps[nsc] = gen.GetFingerprint(m)
df = pd.read_csv(R / "data/almanac_synergy.tsv", sep="\t")
pair = df.groupby(["NSC1", "NSC2"])["MAXSCORE"].max().reset_index()
rows = []
for r in pair.itertuples():
    if r.NSC1 in fps and r.NSC2 in fps and r.NSC1 != r.NSC2:
        rows.append((r.NSC1, r.NSC2, DataStructs.TanimotoSimilarity(fps[r.NSC1], fps[r.NSC2]), r.MAXSCORE, r.MAXSCORE >= 50))
T = pd.DataFrame(rows, columns=["a", "b", "tanimoto", "maxscore", "syn"])
q = T.tanimoto.quantile([.25, .5, .75]).tolist()
T["quartile"] = pd.qcut(T.tanimoto.rank(method="first"), 4, labels=[1, 2, 3, 4])
hi = T.tanimoto >= 0.4
tab = [[int((hi & T.syn).sum()), int((hi & ~T.syn).sum())], [int((~hi & T.syn).sum()), int((~hi & ~T.syn).sum())]]
orr, fp = fisher_exact(tab)
mw = mannwhitneyu(T.tanimoto[T.syn], T.tanimoto[~T.syn])
rho, rp = spearmanr(T.tanimoto, T.maxscore)
out = {"n_drugs_chembl": len(nsc2chembl), "n_drugs_fingerprinted": len(fps),
       "n_pairs": len(T), "synergy_base_rate": round(float(T.syn.mean()), 4),
       "tanimoto_quartile_cuts": [round(x, 4) for x in q],
       "synergy_rate_by_quartile": {int(k): round(float(v), 4) for k, v in T.groupby("quartile", observed=True).syn.mean().items()},
       "median_tanimoto_syn": round(float(T.tanimoto[T.syn].median()), 4),
       "median_tanimoto_nonsyn": round(float(T.tanimoto[~T.syn].median()), 4),
       "mannwhitney_p": float(mw.pvalue),
       "spearman_tanimoto_vs_maxscore": {"rho": round(float(rho), 4), "p": float(rp)},
       "high_similarity_ge_0.4": {"n": int(hi.sum()), "syn_rate": round(float(T.syn[hi].mean()), 4) if hi.sum() else None,
                                   "low_syn_rate": round(float(T.syn[~hi].mean()), 4), "fisher_or": float(orr), "fisher_p": float(fp)}}
T.sort_values("tanimoto", ascending=False).head(15).to_csv(R / "results/structural_similarity_top_pairs.csv", index=False)
json.dump(out, open(R / "results/structural_similarity.json", "w"), indent=1)
print(json.dumps(out, indent=1))

# ---- confound controls: per-drug synergy propensity (leave-pair-out), test volume, molecule size
import statsmodels.api as sm
tests = df.groupby(["NSC1", "NSC2"])["NTESTS"].sum()
T["log_tests"] = np.log([tests.loc[(a, b)] for a, b in zip(T.a, T.b)])
na = {n: Chem.MolFromSmiles(smi[nsc2chembl[n]]).GetNumHeavyAtoms() for n in fps}
T["log_atoms_min"] = np.log([min(na[a], na[b]) for a, b in zip(T.a, T.b)])
deg = pd.concat([T[["a", "syn"]].rename(columns={"a": "d"}), T[["b", "syn"]].rename(columns={"b": "d"})])
s_sum = deg.groupby("d").syn.sum(); s_n = deg.groupby("d").syn.count()
def lpo(d, y):
    return (s_sum[d] - y + 1) / (s_n[d] - 1 + 2)
lg = lambda p: np.log(p / (1 - p))
T["prop_a"] = [lg(lpo(a, y)) for a, y in zip(T.a, T.syn)]
T["prop_b"] = [lg(lpo(b, y)) for b, y in zip(T.b, T.syn)]
X = sm.add_constant(T[["tanimoto", "prop_a", "prop_b", "log_tests", "log_atoms_min"]])
fit = sm.Logit(T.syn.astype(int), X).fit(disp=0)
# drug-label permutation null for the adjusted Tanimoto coefficient: permute fingerprints across drugs
rng = np.random.default_rng(0)
drugs = sorted(fps); perm_coefs = []
for _ in range(200):
    p = dict(zip(drugs, rng.permutation(drugs)))
    tp = [DataStructs.TanimotoSimilarity(fps[p[a]], fps[p[b]]) for a, b in zip(T.a, T.b)]
    Xp = X.copy(); Xp["tanimoto"] = tp
    Xp["log_atoms_min"] = np.log([min(na[p[a]], na[p[b]]) for a, b in zip(T.a, T.b)])
    perm_coefs.append(sm.Logit(T.syn.astype(int), Xp).fit(disp=0).params["tanimoto"])
perm_coefs = np.array(perm_coefs)
c0 = fit.params["tanimoto"]
out["adjusted_logit"] = {"covariates": ["prop_a_lpo", "prop_b_lpo", "log_tests", "log_atoms_min"],
    "coef": {k: round(float(v), 4) for k, v in fit.params.items()},
    "p": {k: float(v) for k, v in fit.pvalues.items()},
    "tanimoto_coef_perm_null_mean": round(float(perm_coefs.mean()), 4),
    "tanimoto_coef_perm_null_sd": round(float(perm_coefs.std()), 4),
    "perm_p_two_sided": float((np.abs(perm_coefs - perm_coefs.mean()) >= abs(c0 - perm_coefs.mean())).mean() + 1/201),
    "n_perm": 200}
json.dump(out, open(R / "results/structural_similarity.json", "w"), indent=1)
print(json.dumps(out["adjusted_logit"], indent=1))

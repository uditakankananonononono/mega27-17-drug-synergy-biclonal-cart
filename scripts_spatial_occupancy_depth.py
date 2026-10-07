"""P6-SO2: depth-stratified permutation null for pooled spatial co-occupancy.
Prespec: results/spatial_occupancy_depth_prespec.json (commit b3326a0, before this script ran).
Run: python3 scripts_spatial_occupancy_depth.py -> results/spatial_occupancy_depth.json"""
import json, os
import numpy as np
from scipy.stats import wilcoxon, chi2
SEED=20261008; N=1000; MIN_POS=10; NBIN=5
LANEB="/home/sandbox/mega27-17s-spatial-morphoscan/features_cache"
ENSG={"CA9":"ENSG00000107159","ENPP3":"ENSG00000154269","MSLN":"ENSG00000102854","PSCA":"ENSG00000167653",
"SLC39A6":"ENSG00000141424","CA12":"ENSG00000074410","ERBB2":"ENSG00000141736","TACSTD2":"ENSG00000184292","MUC1":"ENSG00000185499"}
PAIRS={"target":[("MSLN","CA9"),("MSLN","PSCA")],
"positive_controls":[("ERBB2","TACSTD2"),("ERBB2","MUC1"),("SLC39A6","CA12")]}
rng=np.random.default_rng(SEED)
def jac(a,b):
    u=np.logical_or(a,b).sum(); return float(np.logical_and(a,b).sum()/u) if u else 0.0
res={}
for fn in sorted(f for f in os.listdir(LANEB) if f.endswith('.npz')):
    z=np.load(os.path.join(LANEB,fn),allow_pickle=True); genes=list(z['genes']); M=z['counts_sparse'].item().tocsc()
    depth=np.asarray(M.sum(1)).ravel()
    edges=np.quantile(depth,np.linspace(0,1,NBIN+1)[1:-1]); bins=np.searchsorted(edges,depth,side='right')
    pos={g:(M[:,genes.index(e)].toarray().ravel()>0) for g,e in ENSG.items() if e in genes}
    sec=fn[:-4]
    for cls,prs in PAIRS.items():
        for g1,g2 in prs:
            if g1 not in pos or g2 not in pos or pos[g1].sum()<MIN_POS or pos[g2].sum()<MIN_POS: continue
            obs=jac(pos[g1],pos[g2]); idx=[np.where(bins==b)[0] for b in range(NBIN)]
            nul=np.empty(N)
            for s in range(N):
                p2=pos[g2].copy()
                for ix in idx: p2[ix]=rng.permutation(pos[g2][ix])
                nul[s]=jac(pos[g1],p2)
            res.setdefault(f"{g1}x{g2}",{"class":cls,"sections":[]})["sections"].append(
              {"section":sec,"obs":obs,"null_median":float(np.median(nul)),"p_co":float((np.sum(nul>=obs)+1)/(N+1))})
out={"prespec":"results/spatial_occupancy_depth_prespec.json","seed":SEED,"pairs":{}}
brng=np.random.default_rng(SEED+1)
for pair,r in res.items():
    S=r["sections"]; eff=np.array([s["obs"]-s["null_median"] for s in S]); p=np.array([s["p_co"] for s in S])
    w=float(wilcoxon(eff).pvalue) if len(eff)>=8 and np.any(eff!=0) else None
    f=float(chi2.sf(-2*np.log(p).sum(),2*len(p)))
    bs=[np.median(brng.choice(eff,len(eff))) for _ in range(10000)]
    ok=(w is not None and w<0.05 and f<0.05 and np.median(eff)>0)
    out["pairs"][pair]={"class":r["class"],"n_sections":len(S),"median_effect":float(np.median(eff)),
      "effect_ci95":[float(np.percentile(bs,2.5)),float(np.percentile(bs,97.5))],"wilcoxon_p":w,"fisher_p":f,"passes_rule":bool(ok),"sections":S}
ctrl=all(v["passes_rule"] for v in out["pairs"].values() if v["class"]=="positive_controls")
out["controls_all_pass"]=bool(ctrl)
for k,v in out["pairs"].items():
    if v["class"]=="target": v["verdict"]="SURVIVES" if (v["passes_rule"] and ctrl) else "DOES NOT SURVIVE"
json.dump(out,open("results/spatial_occupancy_depth.json","w"),indent=1)
for k,v in out["pairs"].items(): print(k,v["class"],v["n_sections"],round(v["median_effect"],4),v["effect_ci95"],v["wilcoxon_p"],v["fisher_p"],v.get("verdict",v["passes_rule"]))
print("controls_all_pass",ctrl)

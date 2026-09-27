"""Amendment-queue Tier-1 (#1 + #12 + #13): multi-objective antigen ranking with
explainable feature contributions and a unified safety score, over the 7 gated
panel antigens, computed ONLY from committed results JSONs/CSVs.

Components (each normalized 0-1 across the panel, higher = better target):
  tumor_signal   max p75 tumor nTPM across 10 cancers (results/tumor_percentiles.json)
  normal_safety  1 - scaled max median TPM in vital tissues (results/gtex_safety.json;
                 vital = heart/liver/lung/kidney/brain; flag if no vital tissue in top list)
  protein_breadth fraction of 7 CPTAC studies with detection (results/cptac_protein_rows.csv)
  tractability   Open Targets clinical antibody flag (results/ot_tractability.json)
  nonessentiality absent from DepMap top dependencies with frac_dep>=0.5
                 (results/depmap_dependency.json; WEAK: only top-list available)
Weights declared before computation: tumor .30, normal_safety .25, protein .15,
tractability .15, nonessentiality .15. Contributions (#12) = weight x component.
Validation anchor (#13): clinically targeted MSLN/CLDN18/PSCA should rank top-half.
Novelty reported as a SEPARATE axis (1 - normalized strict clinical-trial count,
results/clintrials_per_gene.csv) - not mixed into the score (novelty vs validation
tradeoff is a judgment, not a sum).
Output: results/rankscore_v1.json
"""
import json, csv, math, collections

WEIGHTS = {"tumor_signal": 0.30, "normal_safety": 0.25, "protein_breadth": 0.15,
           "tractability": 0.15, "nonessentiality": 0.15}
VITAL = {"heart", "liver", "lung", "kidney", "brain", "cerebral"}
CANCERS = ["PAAD", "BRCA", "GBM", "LIHC", "OV", "LUAD", "COAD", "SKCM", "STAD", "KIRC"]

gtex = json.load(open("results/gtex_safety.json"))
panel = {sym: g["gencode"].split(".")[0] for sym, g in gtex.items()}

# 1. tumor signal: max p75 across cancers
tp = json.load(open("results/tumor_percentiles.json"))["data"]
tumor_raw = {}
for sym, gid in panel.items():
    per = tp.get(gid, {})
    p75s = [per[c][1] for c in CANCERS if c in per]
    tumor_raw[sym] = max(p75s) if p75s else 0.0

# 2. normal safety from GTEx vital tissues
vital_max, overall_max, vital_flag = {}, {}, {}
for sym, g in gtex.items():
    tops = g["top_tissues"]
    overall_max[sym] = tops[0][1] if tops else 0.0
    vt = [(t, v) for t, v in tops if any(k in t.lower() for k in VITAL)]
    vital_max[sym] = max((v for _, v in vt), default=0.0)
    vital_flag[sym] = "vital_in_top" if vt else "vital_below_top6_upper_bound"

# 3. CPTAC protein breadth + abundance
studies_seen = collections.defaultdict(set)
abund = collections.defaultdict(list)
with open("results/cptac_protein_rows.csv") as fh:
    for row in csv.DictReader(fh):
        s = row["symbol"]
        if s in panel:
            studies_seen[s].add(row["study"])
            abund[s].append(float(row["log2_protein"]))
breadth = {s: len(studies_seen.get(s, ())) / 7.0 for s in panel}
med_abund = {s: (sorted(abund[s])[len(abund[s]) // 2] if abund.get(s) else None) for s in panel}

# 4. OT tractability
ot = json.load(open("results/ot_tractability.json"))["gated"]["rows"]
gid2sym = {v["symbol"]: k for k, v in ot.items()}
tract = {s: (1.0 if ot.get(gid, {}).get("clinical_ab") else 0.0) for s, gid in panel.items() if gid in ot}
for s in panel:
    tract.setdefault(s, 0.0)

# 5. DepMap nonessentiality (weak: top-dependencies list only)
dd = json.load(open("results/depmap_dependency.json"))
top_dep = {d["gene"]: d["frac_dep"] for d in dd["gated_top_dependencies"]}
noness = {s: (0.0 if top_dep.get(s, 0) >= 0.5 else 1.0) for s in panel}

# 6. novelty axis: inverse strict trial count
ct = {}
with open("results/clintrials_per_gene.csv") as fh:
    for row in csv.DictReader(fh):
        if row["gene"] in panel:
            ct[row["gene"]] = int(row["strict_total"])

# normalize tumor (log1p scaled to panel max), vital (scale to panel max, invert)
def norm_log(d):
    m = max((math.log1p(v) for v in d.values()), default=1.0) or 1.0
    return {k: math.log1p(v) / m for k, v in d.items()}

tumor_n = norm_log(tumor_raw)
vm = max(vital_max.values()) or 1.0
safety_n = {s: 1.0 - vital_max[s] / vm for s in panel}

rows = {}
for s in panel:
    comp = {"tumor_signal": round(tumor_n[s], 4), "normal_safety": round(safety_n[s], 4),
            "protein_breadth": round(breadth[s], 4), "tractability": tract[s],
            "nonessentiality": noness[s]}
    score = sum(WEIGHTS[k] * comp[k] for k in WEIGHTS)
    contrib = {k: round(WEIGHTS[k] * comp[k], 4) for k in WEIGHTS}
    rows[s] = {"score": round(score, 4), "components": comp, "contributions": contrib,
               "raw": {"tumor_p75_max_nTPM": tumor_raw[s], "gtex_vital_max_tpm": vital_max[s],
                       "gtex_overall_max_tpm": overall_max[s], "vital_flag": vital_flag[s],
                       "cptac_studies_detected": len(studies_seen.get(s, ())),
                       "cptac_median_log2": med_abund[s], "strict_trials": ct.get(s, 0),
                       "novelty_axis": None}}
ranked = sorted(rows, key=lambda s: -rows[s]["score"])
# novelty axis after normalization
ctm = max(ct.values()) if ct else 1
for s in panel:
    rows[s]["raw"]["novelty_axis"] = round(1.0 - (ct.get(s, 0) / (ctm or 1)), 4)

anchors = [s for s in ("MSLN", "CLDN18", "PSCA") if s in rows]
anchor_pos = {s: ranked.index(s) + 1 for s in anchors}
out = {"design": __doc__.splitlines()[0], "weights": WEIGHTS, "panel_size": len(panel),
       "ranking": ranked, "rows": rows,
       "validation": {"clinically_targeted_anchors": anchor_pos,
                      "anchors_in_top_half": all(p <= math.ceil(len(panel) / 2) for p in anchor_pos.values())},
       "limitations": ["nonessentiality from DepMap top-dependencies list only (absence != measured non-essentiality)",
                        "vital tissue max is an upper bound when no vital tissue appears in GTEx top-6",
                        "CPTAC detection = presence in committed rows, not a detection-limit model"]}
json.dump(out, open("results/rankscore_v1.json", "w"), indent=1)
print(f"{'rank':<5}{'gene':<10}{'score':<8}{'tumor':<8}{'safety':<8}{'prot':<6}{'tract':<7}{'noness':<7}")
for i, s in enumerate(ranked, 1):
    c = rows[s]["components"]
    print(f"{i:<5}{s:<10}{rows[s]['score']:<8}{c['tumor_signal']:<8}{c['normal_safety']:<8}{c['protein_breadth']:<6}{c['tractability']:<7}{c['nonessentiality']:<7}")
print("anchors:", anchor_pos, "top-half:", out["validation"]["anchors_in_top_half"])

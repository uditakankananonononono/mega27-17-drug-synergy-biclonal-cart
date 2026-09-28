"""Tier-1 #15: independent GEO bulk validation of rankscore panel.
Pre-spec: results/external_bulk_prespec.json (locked 2026-09-28 BEFORE analysis).
Primary: >=6/7 antigens tumor-enriched (median tumor > median normal, one-sided
Mann-Whitney p, BH<0.05 across 7 tests). GSE53757 is paired; pairing NOT exploited
(same rank-sum test everywhere, as pre-specified). Probe->symbol via GEO platform
annotation; multi-mapping probes (a /// b) match if ANY symbol equals the target;
collapse = max-median probe over ALL samples, chosen before any group comparison.
Writes results/external_bulk_validation.json."""
import gzip, json, io, os
import numpy as np
from scipy import stats

GE = "/home/sandbox/geoext"
PANEL = {"GSE53757": ["CA9","CA12"], "GSE62165": ["MSLN"], "GSE42568": ["SLC39A6"],
         "GSE46602": ["PSCA"], "GSE13861": ["CLDN18"], "GSE26712": ["CLDN6"]}
PLAT = {"GSE53757":"GPL570","GSE42568":"GPL570","GSE46602":"GPL570",
        "GSE26712":"GPL96","GSE13861":"GPL6884","GSE62165":"GPL13667"}

def load_map_annot(path):
    m = {}
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("ID\t"):
                break
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) > 2 and p[2]:
                m[p[0]] = p[2]
    return m

def load_map_soft(path):
    m = {}
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                hdr = next(f).rstrip("\n").split("\t")
                si = hdr.index("Gene Symbol")
                for line2 in f:
                    if line2.startswith("!platform_table_end"): break
                    p = line2.rstrip("\n").split("\t")
                    if len(p) > si and p[si] and p[si] != "---":
                        m[p[0]] = p[si]
                break
    return m

def load_series(gse):
    path = os.path.join(GE, f"{gse}_series_matrix.txt.gz")
    with gzip.open(path, "rt", errors="replace") as f:
        titles = chars1 = None
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and chars1 is None:
                chars1 = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                hdr = next(f).rstrip("\n").split("\t")
                ids, rows = [], []
                for line2 in f:
                    if line2.startswith("!series_matrix_table_end"): break
                    p = line2.rstrip("\n").split("\t")
                    ids.append(p[0].strip('"'))
                    rows.append([float(x) if x not in ("","null") else np.nan for x in p[1:]])
                break
    X = np.array(rows)
    if np.nanmax(X) > 100:  # not log scale
        X = np.log2(X + 1.0)
    return titles, chars1, ids, X

def groups(gse, titles, chars1):
    t, n = [], []
    for i, ti in enumerate(titles):
        if gse == "GSE53757":
            (t if ti.endswith("T") else n).append(i)
        elif gse == "GSE42568":
            (n if ti.startswith("Normal") else t).append(i)
        elif gse == "GSE46602":
            (n if ti.startswith("Benign") else t).append(i)
        elif gse == "GSE26712":
            (n if ti.startswith("Normal") else t).append(i)
        elif gse == "GSE13861":
            c = chars1[i]
            if c == "gastric adenocarcinoma": t.append(i)
            elif c == "normal surrounding gastric tissue": n.append(i)
        elif gse == "GSE62165":
            c = chars1[i]
            if "non-tumoral" in c: n.append(i)
            elif "tumor" in c: t.append(i)
    return np.array(t), np.array(n)

def symbols_for(raw, gene):
    out = set()
    for part in raw.split("///"):
        s = part.strip()
        if s == gene: out.add(s)
    return out

results = {"prespec": "results/external_bulk_prespec.json", "per_antigen": {}}
tests = []
for gse, genes in PANEL.items():
    plat = PLAT[gse]
    smap = load_map_soft(os.path.join(GE, "GPL13667.soft")) if plat == "GPL13667" \
           else load_map_annot(os.path.join(GE, f"{plat}.annot.gz"))
    titles, chars1, ids, X = load_series(gse)
    ti, ni = groups(gse, titles, chars1)
    for gene in genes:
        probes = [i for i, pid in enumerate(ids) if gene in symbols_for(smap.get(pid, ""), gene)]
        if not probes:
            results["per_antigen"][gene] = {"gse": gse, "error": "no probe mapped"}
            continue
        meds = np.nanmedian(X[probes, :], axis=1)
        best = probes[int(np.argmax(meds))]
        xt, xn = X[best, ti], X[best, ni]
        u = stats.mannwhitneyu(xt, xn, alternative="greater")
        rec = {"gse": gse, "platform": plat, "probe": ids[best], "n_tumor": int(len(ti)),
               "n_normal": int(len(ni)), "median_tumor": float(np.median(xt)),
               "median_normal": float(np.median(xn)), "p_one_sided": float(u.pvalue),
               "n_probes_mapped": len(probes)}
        results["per_antigen"][gene] = rec
        tests.append((gene, u.pvalue))
# BH
order = sorted(range(len(tests)), key=lambda k: tests[k][1])
m = len(tests); q = [None]*m; prev = 1.0
for rank in range(m-1, -1, -1):
    k = order[rank]; val = min(prev, tests[k][1]*m/(rank+1)); q[k] = val; prev = val
n_pass = 0
for k, (gene, p) in enumerate(tests):
    results["per_antigen"][gene]["q_bh"] = float(q[k])
    ok = bool(q[k] < 0.05 and results["per_antigen"][gene]["median_tumor"] > results["per_antigen"][gene]["median_normal"])
    results["per_antigen"][gene]["enriched"] = ok
    n_pass += ok
results["primary_endpoint"] = {"rule": ">=6/7 enriched at BH<0.05", "n_enriched": n_pass,
                               "n_tested": m, "PASS": bool(n_pass >= 6)}
# secondary: Spearman external tumor medians vs rankscore tumor_signal ordering
rk = json.load(open("results/rankscore_v1.json"))
pairs = [(results["per_antigen"][g]["median_tumor"], rk["rows"][g]["raw"]["tumor_p75_max_nTPM"])
         for g in rk["ranking"] if g in results["per_antigen"] and "median_tumor" in results["per_antigen"][g]]
if len(pairs) >= 5:
    ext = [p[0] for p in pairs]; intn = [p[1] for p in pairs]
    rho, sp = stats.spearmanr(ext, intn)
    results["secondary_endpoint"] = {"genes": [g for g in rk["ranking"] if g in results["per_antigen"] and "median_tumor" in results["per_antigen"][g]],
        "external_tumor_medians": ext, "rankscore_tumor_signal_raw": intn,
        "spearman_rho": float(rho), "p": float(sp), "note": "reported as measured, no threshold (prespec)"}
json.dump(results, open("results/external_bulk_validation.json", "w"), indent=1)
print(json.dumps({"n_enriched": n_pass, "PASS": results["primary_endpoint"]["PASS"],
                  "per_gene": {g: r.get("q_bh") for g, r in results["per_antigen"].items()},
                  "spearman": results.get("secondary_endpoint", {}).get("spearman_rho")}, indent=1))

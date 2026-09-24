#!/usr/bin/env python3
"""iLINCS (LINCS L1000) fetcher for the antigen visibility/connectivity audit.

Subcommands (all resumable, polite-rate-limited):
  coverage  - per-gene trt_sh.cgs consensus-signature counts for the 205 audit genes
  kd        - knockdown/overexpression signature metadata + landmark vectors (CA12, SLC39A6, ERBB2)
  cpmeta    - page all exemplar trt_cp signature metadata
  cpmap     - map GDSC compounds -> iLINCS exemplar signatures (core cell lines), write id list
  cpsigs    - download mapped compound signature landmark vectors
"""
import json, os, sys, time, urllib.request, urllib.parse

BASE = "https://www.ilincs.org/api"
DATA = "data/ilincs"
os.makedirs(DATA, exist_ok=True)
SLEEP = 0.12

def get(path, params=None, tries=4):
    url = BASE + path + ("?" + urllib.parse.urlencode(params) if params else "")
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36"})
            with urllib.request.urlopen(req, timeout=40) as r:
                time.sleep(SLEEP)
                return json.loads(r.read().decode())
        except Exception as e:
            if i == tries - 1: raise
            time.sleep(1.5 * (i + 1))

def filt(where, **kw):
    f = {"where": where}; f.update(kw)
    return {"filter": json.dumps(f)}

def cmd_coverage(maxn):
    sets = json.load(open("results/depmap_gene_sets.json"))
    genes = ([(g, "gated") for g in sets["gated"]] +
             [(g, "background") for g in sets["background"]] +
             [(g, "reference") for g in sets["references"]] +
             [(g, "positive_control") for g in sets["positive_controls"]])
    out = "results/ilincs_cgs_coverage.csv"
    done = set()
    if os.path.exists(out):
        for ln in open(out).read().splitlines()[1:]:
            if ln.strip(): done.add(ln.split(",")[0])
    else:
        open(out, "w").write("gene,set,n_cgs_signatures,covered\n")
    n = 0
    with open(out, "a") as fh:
        for g, s in genes:
            if g in done: continue
            c = get("/SignatureMeta/count", {"where": json.dumps({"pert_type": "trt_sh.cgs", "treatment": g})})["count"]
            fh.write(f"{g},{s},{c},{1 if c > 0 else 0}\n"); fh.flush()
            n += 1
            if n >= maxn: break
    print(f"coverage: +{n} (total {len(done)+n}/{len(genes)})")

def dl_sig(sid):
    p = f"{DATA}/sig_{sid}.json"
    if os.path.exists(p): return False
    d = get("/ilincsR/downloadSignature", {"sigID": sid})
    sig = d["data"]["signature"]
    json.dump(sig, open(p, "w"))
    return True

def cmd_kd():
    meta_p = f"{DATA}/kd_meta.json"
    if not os.path.exists(meta_p):
        recs = []
        for g in ("CA12", "SLC39A6", "ERBB2"):
            recs += get("/SignatureMeta", filt({"pert_type": "trt_sh.cgs", "treatment": g}))
            recs += get("/SignatureMeta", filt({"pert_type": "trt_oe", "treatment": g}))
        json.dump(recs, open(meta_p, "w"))
    recs = json.load(open(meta_p))
    n = sum(dl_sig(r["signatureid"]) for r in recs)
    lm = sorted({row["Name_GeneSymbol"] for row in json.load(open(f'{DATA}/sig_{recs[0]["signatureid"]}.json'))})
    open(f"{DATA}/landmark978.txt", "w").write("\n".join(lm) + "\n")
    print(f"kd: {len(recs)} meta records, +{n} signatures, landmark={len(lm)}")

def cmd_cpmeta():
    p = f"{DATA}/cp_exemplar_meta.jsonl"
    have = sum(1 for _ in open(p)) if os.path.exists(p) else 0
    total = get("/SignatureMeta/count", {"where": json.dumps({"pert_type": "trt_cp", "is_exemplar": 1})})["count"]
    page = 5000
    with open(p, "a") as fh:
        while have < total:
            recs = get("/SignatureMeta", filt({"pert_type": "trt_cp", "is_exemplar": 1},
                                              limit=page, skip=have,
                                              fields=["signatureid", "compound", "clueIoCompound",
                                                       "cellline", "integratedMoas", "treatment", "lincsSigID"]))
            if not recs: break
            for r in recs: fh.write(json.dumps(r) + "\n")
            have += len(recs); fh.flush()
            print(f"cpmeta: {have}/{total}")
            if len(recs) < page: break
    print(f"cpmeta done: {have}")

CORE = {"A375", "A549", "HA1E", "HCC515", "HEPG2", "HT29", "MCF7", "PC3", "VCAP"}

def cmd_cpmap(cap):
    import csv
    names = {}
    for row in csv.DictReader(open("data/gdsp_compounds_8.5.csv")):
        keys = {row["DRUG_NAME"].strip().lower()}
        for s in (row.get("SYNONYMS") or "").split(","):
            if s.strip(): keys.add(s.strip().lower())
        for k in keys: names.setdefault(k, row["DRUG_NAME"].strip())
    matched, seen = {}, set()
    for ln in open(f"{DATA}/cp_exemplar_meta.jsonl"):
        r = json.loads(ln)
        if r.get("cellline") not in CORE: continue
        for cand in (r.get("compound"), r.get("clueIoCompound"), r.get("treatment")):
            k = (cand or "").strip().lower()
            if k in names:
                matched.setdefault(names[k], []).append(r); break
    rows, ids = [], []
    for drug in sorted(matched):
        sigs = sorted({r["signatureid"]: r for r in matched[drug]}.values(),
                      key=lambda r: (r["cellline"], r["signatureid"]))
        for r in sigs:
            ids.append(r["signatureid"])
            rows.append({"drug": drug, "signatureid": r["signatureid"], "cellline": r["cellline"],
                         "moa": r.get("integratedMoas") or "", "compound": r.get("compound") or ""})
    rows = rows[:cap]
    with open("results/ilincs_gdsp_map.csv", "w") as fh:
        fh.write("drug,signatureid,cellline,moa,compound\n")
        for r in rows:
            fh.write(",".join('"%s"' % str(v).replace('"', "'") if "," in str(v) else str(v)
                              for v in (r["drug"], r["signatureid"], r["cellline"], r["moa"], r["compound"])) + "\n")
    open(f"{DATA}/cp_sig_ids.txt", "w").write("\n".join(r["signatureid"] for r in rows) + "\n")
    print(f"cpmap: {len(matched)} GDSC drugs matched -> {len(rows)} signatures (cap {cap})")

def cmd_cpsigs(maxn):
    from concurrent.futures import ThreadPoolExecutor
    ids = [l.strip() for l in open(f"{DATA}/cp_sig_ids.txt") if l.strip()]
    todo = [sid for sid in ids if not os.path.exists(f"{DATA}/sig_{sid}.json")][:maxn]
    n = 0
    with ThreadPoolExecutor(max_workers=6) as ex:
        for ok in ex.map(dl_sig, todo):
            n += bool(ok)
    have = sum(1 for sid in ids if os.path.exists(f"{DATA}/sig_{sid}.json"))
    print(f"cpsigs: +{n} ({have}/{len(ids)})")

if __name__ == "__main__":
    cmd = sys.argv[1]
    maxn = int(sys.argv[2]) if len(sys.argv) > 2 else 100000
    {"coverage": cmd_coverage, "kd": lambda m: cmd_kd(), "cpmeta": lambda m: cmd_cpmeta(),
     "cpmap": lambda m: cmd_cpmap(m), "cpsigs": cmd_cpsigs}[cmd](maxn)

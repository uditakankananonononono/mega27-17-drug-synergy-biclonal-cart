#!/usr/bin/env python3
"""Fetch IntAct (EBI PSICQUIC) curated interaction records for the 205-gene
DepMap universe. Per gene: format=count + full paginated MITAB tab25.
Caches under data/intact/{counts,mitab}/. Resumable: skips existing outputs.
Positive-control accessions resolved live via UniProt REST (cached)."""
import json, glob, os, re, sys, time, urllib.request, urllib.parse

BASE = "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
CD = "data/intact/counts"; MD = "data/intact/mitab"
UA = {"User-Agent": "mega27-item17-intact-audit (research)"}

def get(url, tries=3, timeout=40):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace"), r.status
        except Exception as e:
            if i == tries - 1:
                return f"__ERROR__ {e}", 0
            time.sleep(2 * (i + 1))
    return "__ERROR__", 0

def accession_map():
    m = {}
    for f in glob.glob("data/epitope/uniprot_*.json"):
        sym = os.path.basename(f)[8:-5]
        try:
            d = json.load(open(f))
            res = d.get("results", [])
            if res:
                m[sym] = res[0]["primaryAccession"]
        except Exception:
            pass
    # 4 positive controls: resolve live via UniProt REST, cache
    for sym in ["POLR2A", "RPS3", "PCNA", "PSMA1"]:
        cf = f"data/intact/uniprot_ctrl_{sym}.json"
        if os.path.exists(cf):
            m[sym] = json.load(open(cf))["accession"]
            continue
        q = urllib.parse.quote(f"gene_exact:{sym} AND organism_id:9606 AND reviewed:true")
        txt, st = get(f"{UNIPROT}?query={q}&fields=accession,gene_names&size=1&format=json")
        if st == 200 and not txt.startswith("__ERROR__"):
            res = json.loads(txt).get("results", [])
            if res:
                m[sym] = res[0]["primaryAccession"]
                json.dump({"symbol": sym, "accession": m[sym]}, open(cf, "w"), indent=1)
        time.sleep(0.3)
    return m

def fetch_count(sym, acc):
    f = f"{CD}/{sym}.txt"
    if os.path.exists(f) and os.path.getsize(f) > 0:
        return None
    txt, st = get(f"{BASE}/uniprotkb:{acc}?format=count")
    open(f, "w").write(txt.strip())
    return st

def fetch_mitab(sym, acc, expected):
    f = f"{MD}/{sym}.mitab"
    if os.path.exists(f):
        n = sum(1 for _ in open(f))
        if n >= (expected or 1):
            return None
    rows, first = [], 0
    while True:
        url = f"{BASE}/uniprotkb:{acc}?format=tab25&firstResult={first}&maxResults=1000"
        txt, st = get(url)
        if st != 200 or txt.startswith("__ERROR__"):
            break
        page = [l for l in txt.splitlines() if l.strip()]
        rows.extend(page)
        if len(page) < 1000:
            break
        first += 1000
        time.sleep(0.2)
    open(f, "w").write("\n".join(rows) + ("\n" if rows else ""))
    return len(rows)

def main():
    gs = json.load(open("results/depmap_gene_sets.json"))
    genes = []
    for k in ("gated", "background", "references", "positive_controls"):
        genes += [(g, k) for g in gs[k]]
    amap = accession_map()
    man = []
    t0 = time.time()
    for i, (sym, grp) in enumerate(genes):
        acc = amap.get(sym)
        rec = {"symbol": sym, "group": grp, "accession": acc}
        if not acc:
            rec["status"] = "no_uniprot_record"; man.append(rec); continue
        st = fetch_count(sym, acc)
        cf = f"{CD}/{sym}.txt"
        ctxt = open(cf).read().strip() if os.path.exists(cf) else ""
        rec["count"] = int(ctxt) if re.fullmatch(r"\d+", ctxt) else None
        rec["count_http"] = st
        n = fetch_mitab(sym, acc, rec["count"])
        mf = f"{MD}/{sym}.mitab"
        rec["mitab_rows"] = sum(1 for _ in open(mf)) if os.path.exists(mf) else 0
        rec["status"] = "ok" if rec["count"] is not None else "count_error"
        man.append(rec)
        if (i + 1) % 20 == 0:
            print(f"{i+1}/205 fetched, {time.time()-t0:.0f}s", flush=True)
        time.sleep(0.15)
    json.dump(man, open("data/intact/intact_manifest.json", "w"), indent=1)
    ok = sum(1 for r in man if r.get("status") == "ok")
    print(f"DONE ok={ok}/{len(man)} elapsed={time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()

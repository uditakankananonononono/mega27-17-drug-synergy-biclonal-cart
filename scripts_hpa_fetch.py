#!/usr/bin/env python3
"""Fetch Human Protein Atlas per-gene TSV records for the 17 gene universe.

Genes: results/depmap_gene_sets.json (gated/background/references/positive_controls).
Ensembl IDs resolved via HGNC REST (cached). HPA per-gene TSV at
https://www.proteinatlas.org/<ENSEMBL>.tsv (HPA is the new tool under audit).
Resumable: skips genes whose TSV already exists. Writes data/hpa/hpa_manifest.json.
"""
import json, os, time, urllib.request, urllib.error

UA = {"User-Agent": "mega27-laneF-research/1.0 (academic audit)"}
HGNC_DIR = "data/hpa/hgnc"
TSV_DIR = "data/hpa/tsv"
os.makedirs(HGNC_DIR, exist_ok=True)
os.makedirs(TSV_DIR, exist_ok=True)

def get(url, binary=False):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.read() if binary else r.read().decode()

def hgnc_ensembl(sym):
    p = os.path.join(HGNC_DIR, f"{sym}.json")
    if os.path.exists(p):
        d = json.load(open(p))
    else:
        d = json.loads(get(f"https://rest.genenames.org/fetch/symbol/{sym}",
                           ).replace("\\", "\\\\") if False else
                urllib.request.urlopen(urllib.request.Request(
                    f"https://rest.genenames.org/fetch/symbol/{sym}",
                    headers={**UA, "Accept": "application/json"}), timeout=30).read().decode())
        json.dump(d, open(p, "w"))
        time.sleep(0.35)
    docs = d.get("response", {}).get("docs", [])
    for doc in docs:
        if doc.get("symbol", "").upper() == sym.upper() and doc.get("ensembl_gene_id"):
            return doc["ensembl_gene_id"]
    return None

def main():
    sets = json.load(open("results/depmap_gene_sets.json"))
    genes = []
    for grp in ["gated", "background", "references", "positive_controls"]:
        for g in sets[grp]:
            genes.append((g, grp))
    manifest = []
    done = skipped = 0
    for i, (sym, grp) in enumerate(genes):
        out = os.path.join(TSV_DIR, f"{sym}.tsv")
        rec = {"gene": sym, "group": grp, "ensembl": None, "hpa_tsv": None, "status": None}
        if os.path.exists(out) and os.path.getsize(out) > 100:
            rec["status"] = "cached"; rec["hpa_tsv"] = out; done += 1
            manifest.append(rec); continue
        try:
            ens = hgnc_ensembl(sym)
        except Exception as e:
            rec["status"] = f"hgnc_error:{e}"; manifest.append(rec); continue
        rec["ensembl"] = ens
        if not ens:
            rec["status"] = "no_ensembl"; manifest.append(rec); continue
        try:
            tsv = get(f"https://www.proteinatlas.org/{ens}.tsv")
            if "\t" not in tsv.splitlines()[0]:
                rec["status"] = "no_hpa_record"
            else:
                open(out, "w").write(tsv)
                rec["status"] = "ok"; rec["hpa_tsv"] = out; done += 1
        except urllib.error.HTTPError as e:
            rec["status"] = f"hpa_http_{e.code}"
        except Exception as e:
            rec["status"] = f"hpa_error:{e}"
        manifest.append(rec)
        time.sleep(0.35)
        if (i + 1) % 25 == 0:
            print(f"{i+1}/{len(genes)} done={done}", flush=True)
            json.dump(manifest, open("data/hpa/hpa_manifest.json", "w"), indent=1)
    json.dump(manifest, open("data/hpa/hpa_manifest.json", "w"), indent=1)
    from collections import Counter
    print("DONE", Counter(r["status"] for r in manifest))

if __name__ == "__main__":
    main()

"""MyGene.info lane + restore committed HPA per-gene TSVs.

Resolves each HPA-audit gene symbol to its Ensembl gene id through the
MyGene.info v3 API, then re-fetches the per-gene HPA TSV
(proteinatlas.org/<ENSG>.tsv) into data/hpa/tsv/<SYMBOL>.tsv for every gene
whose committed audit row has if_plasma_membrane=1 (the raw grounding the
test suite asserts on). Committed: results/mygene_symbol_ensg.csv,
results/mygene_summary.json, data/hpa/tsv/*.tsv.
"""
import csv, json, time, urllib.parse, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
UA = {"User-Agent": "mega27-research/1.0", "Accept": "application/json"}


def get(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.read()


def mygene(symbol):
    q = urllib.parse.quote(f"symbol:{symbol}")
    d = json.loads(get(f"https://mygene.info/v3/query?q={q}&species=human&fields=ensembl.gene,symbol,type_of_gene&size=5"))
    hits = d.get("hits", [])
    if not hits:
        return ""
    def ens_list(h):
        ens = h.get("ensembl", {})
        if isinstance(ens, dict):
            ens = [ens]
        return [e.get("gene", "") for e in ens if e.get("gene")]
    # Prefer exact-symbol protein-coding hits: MyGene's first ensembl entry can
    # be a readthrough/novel locus (e.g. ENSG0000028xxxx) that HPA does not
    # carry. Return ALL candidate ids in preference order; the caller tries
    # each against HPA and keeps the first that resolves.
    exact_pc = [h for h in hits if h.get("symbol") == symbol and h.get("type_of_gene") == "protein-coding"]
    exact = [h for h in hits if h.get("symbol") == symbol]
    out = []
    for pool in (exact_pc, exact, hits):
        for h in pool:
            out.extend(ens_list(h))
    seen, uniq = set(), []
    for e in out:
        if e not in seen:
            seen.add(e); uniq.append(e)
    return uniq


def main():
    rows = list(csv.DictReader(open(ROOT / "results" / "hpa_per_gene.csv")))
    need = [r["gene"] for r in rows if r["status"] == "ok" and r["if_plasma_membrane"] == "1"]
    outdir = ROOT / "data" / "hpa" / "tsv"
    outdir.mkdir(parents=True, exist_ok=True)
    from concurrent.futures import ThreadPoolExecutor
    def one(g):
        ens = ""
        got = False
        for cand in mygene(g):
            try:
                tsv = get(f"https://www.proteinatlas.org/{cand}.tsv").decode()
                (outdir / f"{g}.tsv").write_text(tsv)
                ens, got = cand, True
                break
            except Exception:
                continue
        return {"symbol": g, "ensembl_gene": ens}, got
    with ThreadPoolExecutor(max_workers=8) as ex:
        out = list(ex.map(one, sorted(need)))
    maps = [m for m, _ in out]
    fetched = sum(1 for _, g in out if g)
    with open(ROOT / "results" / "mygene_symbol_ensg.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["symbol", "ensembl_gene"])
        w.writeheader(); w.writerows(maps)
    json.dump({"n_symbols": len(maps), "n_resolved": sum(1 for m in maps if m["ensembl_gene"]),
               "n_hpa_tsv_fetched": fetched},
              open(ROOT / "results" / "mygene_summary.json", "w"), indent=2)
    print(fetched, "/", len(maps))


if __name__ == "__main__":
    main()

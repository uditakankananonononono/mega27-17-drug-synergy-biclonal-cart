#!/usr/bin/env python3
"""Fetch ProteomicsDB (TUM, api_v2 OData) protein-level expression for the
205-gene study set.

Per gene: UniProt accession (HGNC cache) -> ProteomicsDB PROTEIN_ID (target,
not decoy) -> every ProteinExpression row (SAMPLE_ID, EXPRESSION,
CALCULATION_METHOD, PEPTIDES). Sample annotations (tissue, disease) are
fetched once as a paged table.

Caches: data/proteomicsdb/genes/<SYM>.json (resumable), data/proteomicsdb/samples.json
Question this feeds: mass-spectrometry protein detection breadth in normal
tissues - the protein-level off-tumour check on gated antigens.
"""
import json, os, sys, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor

API = "https://www.proteomicsdb.org/proteomicsdb/logic/api_v2/api.xsodata"
BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "data", "proteomicsdb")
GDIR = os.path.join(OUT, "genes")


def get(path, params, tries=4, timeout=100):
    q = urllib.parse.urlencode(dict(params, **{"$format": "json"}))
    url = "%s/%s?%s" % (API, path, q)
    for t in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "mega27-audit/1.0"})
            return json.load(urllib.request.urlopen(req, timeout=timeout))["d"]["results"]
        except Exception as e:
            if t == tries - 1:
                return {"_error": str(e)}
            time.sleep(2 * (t + 1))


def gene_set():
    gs = json.load(open(os.path.join(BASE, "results", "depmap_gene_sets.json")))
    out = []
    for grp in ("gated", "background", "references", "positive_controls"):
        out.extend(gs[grp])
    return out


def uniprot_for(sym):
    p = os.path.join(BASE, "data", "constraint", "hgnc_%s.json" % sym)
    try:
        docs = json.load(open(p))["response"]["docs"]
        return (docs[0].get("uniprot_ids") or [None])[0]
    except Exception:
        return None


def fetch_samples(page=2000, n_total=None):
    """Paged Sample table, pages fetched in parallel (each page ~40 s server-side)."""
    p = os.path.join(OUT, "samples.json")
    if os.path.exists(p):
        return
    if n_total is None:
        n_total = int(urllib.request.urlopen(API + "/Sample/$count", timeout=60).read())
    skips = list(range(0, n_total, page))

    def one(skip):
        return get("Sample", {"$select": "SAMPLE_ID,EXPERIMENT_ID,TAXCODE,TISSUE_ID,TISSUE,DISEASE",
                              "$top": page, "$skip": skip, "$orderby": "SAMPLE_ID"})
    with ThreadPoolExecutor(len(skips)) as ex:
        pages = list(ex.map(one, skips))
    rows = []
    for r in pages:
        if isinstance(r, dict):
            raise SystemExit("sample page failed: %s" % r["_error"])
        rows += [{k: x.get(k) for k in ("SAMPLE_ID", "EXPERIMENT_ID", "TAXCODE", "TISSUE_ID",
                                         "TISSUE", "DISEASE")} for x in r]
    assert len({r["SAMPLE_ID"] for r in rows}) == len(rows) == n_total
    json.dump(rows, open(p, "w"))


def fetch_gene(sym):
    p = os.path.join(GDIR, sym + ".json")
    if os.path.exists(p):
        return "skip"
    acc = uniprot_for(sym)
    rec = {"gene": sym, "uniprot": acc}
    if not acc:
        rec["_error"] = "no_uniprot"
    else:
        prot = get("Protein", {"$filter": "UNIQUE_IDENTIFIER eq '%s'" % acc,
                               "$select": "PROTEIN_ID,GENE_NAME,DECOY"})
        if isinstance(prot, dict):
            return "err"  # transient: leave uncached so it retries
        ids = [x["PROTEIN_ID"] for x in prot if not x.get("DECOY")]
        rec["protein_ids"] = ids
        rec["expression"] = []
        for pid in ids[:1]:
            ex = get("ProteinExpression", {"$filter": "PROTEIN_ID eq %d" % pid,
                                           "$select": "SAMPLE_ID,EXPRESSION,CALCULATION_METHOD,PEPTIDES"})
            if isinstance(ex, dict):
                return "err"
            rec["expression"] = [[x["SAMPLE_ID"], float(x["EXPRESSION"]),
                                  x["CALCULATION_METHOD"], x["PEPTIDES"]] for x in ex]
        if not ids:
            rec["_error"] = "not_in_proteomicsdb"
    json.dump(rec, open(p, "w"))
    return "ok"


def main(budget_s=90, workers=32):
    os.makedirs(GDIR, exist_ok=True)
    if "--samples" in sys.argv:
        fetch_samples()
        return
    t0 = time.time()
    todo = [g for g in gene_set() if not os.path.exists(os.path.join(GDIR, g + ".json"))]
    done = []
    with ThreadPoolExecutor(workers) as ex:
        for i in range(0, len(todo), workers):
            if time.time() - t0 > budget_s:
                break
            done += list(ex.map(fetch_gene, todo[i:i + workers]))
    print("ok", done.count("ok"), "err", done.count("err"),
          "remaining", len(todo) - done.count("ok"))


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1][0].isdigit() else 90)

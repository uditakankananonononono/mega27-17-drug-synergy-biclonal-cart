#!/usr/bin/env python3
"""Fetch Pharos (TCRD / Illuminating the Druggable Genome) target records for the
205-gene study set plus an EGFR sentinel via the Pharos GraphQL API.

Per-gene cache: data/pharos/<SYM>.json  (resumable: existing ok files skipped)
Manifest:       data/pharos/manifest.json

Question this feeds: target-illumination audit (TDL, TIN-X novelty, drug/ligand
counts) of gated antigens vs background, and cross-source concordance of Pharos
Tclin (current) vs DrugCentral Tclin (2021_09_01 snapshot).
"""
import json, os, sys, time, urllib.request

API = "https://pharos-api.ncats.io/graphql"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "pharos")

QUERY = """{ target(q:{sym:\"%s\"}) { sym tdl novelty publicationCount generifCount
  ligandCounts { name value } ppiCounts { name value } gwasCounts { name value } } }"""
QUERY_ACC = """{ target(q:{uniprot:\"%s\"}) { sym tdl novelty publicationCount generifCount
  ligandCounts { name value } ppiCounts { name value } gwasCounts { name value } } }"""


def accession_for(sym):
    """UniProt accession from the epitope audit cache (symbol-rename fallback)."""
    import json as _json, os as _os
    p = _os.path.join("data", "epitope", "uniprot_%s.json" % sym)
    if not _os.path.exists(p):
        return None
    try:
        return _json.load(open(p))["results"][0]["primaryAccession"]
    except Exception:
        return None


def gene_set():
    with open("results/depmap_gene_sets.json") as f:
        gs = json.load(f)
    genes = []
    for grp in ("gated", "background", "references", "positive_controls"):
        for g in gs[grp]:
            genes.append((g, grp))
    genes.append(("EGFR", "sentinel"))
    return genes


def fetch(sym, tries=4):
    q = QUERY % sym
    acc = accession_for(sym)
    if acc:
        q = QUERY % acc + "__OR__" + QUERY % sym
    return _fetch_multi(sym, q, tries)


def _fetch_multi(sym, q, tries=4):
    """Symbol query first; on null target, fall back to UniProt accession."""
    if "__OR__" in q:
        acc_q, sym_q = q.split("__OR__")
        rec = _post(sym_q, tries)
        tgt = (rec.get("data") or {}).get("target") if "_error" not in rec else None
        if tgt:
            return rec
        rec2 = _post(acc_q, tries)
        tgt2 = (rec2.get("data") or {}).get("target") if "_error" not in rec2 else None
        if tgt2:
            rec2["_fallback_accession"] = True
        return rec2 if tgt2 else rec
    return _post(q, tries)


def _post(q, tries=4):
    body = json.dumps({"query": q}).encode()
    req = urllib.request.Request(API, data=body,
                                 headers={"Content-Type": "application/json"})
    for t in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except Exception as e:
            if t == tries - 1:
                return {"_error": "%s: %s" % (type(e).__name__, e)}
            time.sleep(2 ** t)


def main():
    os.makedirs(OUT, exist_ok=True)
    genes = gene_set()
    manifest = []
    done = 0
    for sym, grp in genes:
        path = os.path.join(OUT, sym + ".json")
        if os.path.exists(path):
            try:
                with open(path) as f:
                    rec = json.load(f)
                if "_error" not in rec:
                    manifest.append({"sym": sym, "group": grp, "status": "cached"})
                    done += 1
                    continue
            except Exception:
                pass
        rec = fetch(sym)
        with open(path, "w") as f:
            json.dump(rec, f)
        tgt = (rec.get("data") or {}).get("target") if "_error" not in rec else None
        status = "error" if "_error" in rec else ("ok" if tgt else "no_record")
        manifest.append({"sym": sym, "group": grp, "status": status})
        done += 1
        print("[%d/%d] %s %s" % (done, len(genes), sym, status), flush=True)
        time.sleep(0.3)
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump({"fetched": time.strftime("%Y-%m-%d %H:%M:%S"),
                   "n": len(manifest), "records": manifest}, f, indent=1)
    n_err = sum(1 for m in manifest if m["status"] == "error")
    print("done: %d records, %d errors" % (len(manifest), n_err))


if __name__ == "__main__":
    sys.exit(main())

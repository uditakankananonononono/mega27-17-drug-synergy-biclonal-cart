#!/usr/bin/env python3
"""Fetch Monarch Initiative v3 gene-disease (causal + correlated) and gene-phenotype
(HPO) associations for the 205-gene study set via the Monarch v3 API.

Per-gene cache: data/monarch/<SYM>.json  (resumable: good files skipped)
Manifest:       data/monarch/manifest.json

Question this feeds: Mendelian-disease causality (OMIM/Orphanet causal and
correlated gene-disease associations) and HPO phenotypic breadth of gated
antigens vs background - the rare/Mendelian complement to the GWAS common-
variant germline audit (results/gwas_audit.json).
"""
import json, os, time, urllib.parse, urllib.request

API = "https://api.monarchinitiative.org/v3/api"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "monarch")
CATS = ["CausalGeneToDiseaseAssociation", "CorrelatedGeneToDiseaseAssociation",
        "GeneToPhenotypicFeatureAssociation"]


def gene_set():
    gs = json.load(open("results/depmap_gene_sets.json"))
    genes = []
    for grp in ("gated", "background", "references", "positive_controls"):
        genes.extend(gs[grp])
    return genes


def hgnc_id_for(sym):
    """HGNC id from the constraint-audit cache; Monarch search fallback."""
    p = os.path.join("data", "constraint", "hgnc_%s.json" % sym)
    if os.path.exists(p):
        try:
            docs = json.load(open(p))["response"]["docs"]
            for d in docs:
                if d.get("symbol") == sym and d.get("hgnc_id"):
                    return d["hgnc_id"]
            if docs and docs[0].get("hgnc_id"):
                return docs[0]["hgnc_id"]
        except Exception:
            pass
    d = get("%s/search?term=%s&category=biolink:Gene&limit=10"
            % (API, urllib.parse.quote(sym)))
    for it in (d.get("items") or []):
        if it.get("symbol") == sym or it.get("name") == sym:
            return it.get("id")
    return None


def get(url, tries=4):
    for t in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "mega27-audit/1.0"})
            return json.load(urllib.request.urlopen(req, timeout=40))
        except Exception as e:
            if t == tries - 1:
                return {"_error": str(e)}
            time.sleep(1.5 * (t + 1))


def fetch_gene(sym):
    rec = {"gene": sym}
    hg = hgnc_id_for(sym)
    rec["hgnc_id"] = hg
    if not hg:
        rec["_error"] = "no_hgnc_id"
        return rec
    for cat in CATS:
        # phenotypes: only the total is needed (limit=1); diseases: keep labels
        lim = 1 if "Phenotypic" in cat else 200
        u = "%s/association?subject=%s&category=biolink:%s&limit=%d" % (
            API, urllib.parse.quote(hg), cat, lim)
        d = get(u)
        if "_error" in d:
            rec["_error_" + cat] = d["_error"]
            continue
        rec[cat + "_total"] = d.get("total", 0)
        if "Disease" in cat:
            rec[cat + "_items"] = [
                {"id": i.get("object"), "label": i.get("object_label"),
                 "source": i.get("primary_knowledge_source")}
                for i in d.get("items", [])]
        time.sleep(0.12)
    return rec


def main():
    os.makedirs(OUT, exist_ok=True)
    genes = gene_set()
    for i, sym in enumerate(genes):
        p = os.path.join(OUT, "%s.json" % sym)
        if os.path.exists(p):
            try:
                r = json.load(open(p))
                if "_error" not in r and not any(k.startswith("_error_") for k in r):
                    continue
            except Exception:
                pass
        json.dump(fetch_gene(sym), open(p, "w"))
        if (i + 1) % 25 == 0:
            print("%d/%d" % (i + 1, len(genes)), flush=True)
    done, errors = [], []
    for f in os.listdir(OUT):
        if not f.endswith(".json") or f == "manifest.json":
            continue
        r = json.load(open(os.path.join(OUT, f)))
        if "_error" in r or any(k.startswith("_error_") for k in r):
            errors.append(r.get("gene", f[:-5]))
        else:
            done.append(r.get("gene", f[:-5]))
    json.dump({"n_genes": len(genes), "done": sorted(done), "errors": sorted(errors)},
              open(os.path.join(OUT, "manifest.json"), "w"), indent=1)
    print("done=%d errors=%d" % (len(done), len(errors)))


if __name__ == "__main__":
    main()

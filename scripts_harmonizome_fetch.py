#!/usr/bin/env python3
"""Fetch Harmonizome (Ma'ayan Lab) per-gene association profiles for the
205-gene study set, plus dataset-level metadata (datasetGroup, measurement,
attributeGroup) for every Harmonizome dataset.

Harmonizome integrates ~114 omics / curated / text-mined datasets into
thresholded gene-attribute associations. Per gene we keep a compact summary:
{dataset: [n_up, n_down]} (thresholdValue +1 / -1), not the raw 10k-row lists.

Per-gene cache: data/harmonizome/genes/<SYM>.json (resumable; good files skipped)
Dataset meta:   data/harmonizome/dataset_meta.json
Question this feeds: is the gated panel's elevated annotation density (IntAct,
HPA-IF, BioPlex detection) confined to curated/literature-derived resources,
while systematic high-throughput resources show no difference?
"""
import json, os, sys, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor

API = "https://maayanlab.cloud/Harmonizome/api/1.0"
BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "data", "harmonizome")
GDIR = os.path.join(OUT, "genes")
META_KEYS = ("name", "association", "datasetGroup", "measurement",
             "attributeGroup", "attributeType")


def get(url, tries=4, timeout=90):
    for t in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "mega27-audit/1.0"})
            return json.load(urllib.request.urlopen(req, timeout=timeout))
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


def dataset_names():
    names, url = [], API + "/dataset"
    while url:
        d = get(url)
        names += [e["name"] for e in d.get("entities", [])]
        nxt = d.get("next")
        url = ("https://maayanlab.cloud" + nxt) if nxt else None
    return names


def fetch_meta():
    p = os.path.join(OUT, "dataset_meta.json")
    if os.path.exists(p):
        return json.load(open(p))
    names = dataset_names()

    def one(n):
        d = get(API + "/dataset/" + urllib.parse.quote_plus(n))
        m = {k: d.get(k) for k in META_KEYS}
        m["name"] = n
        m["n_gene_sets"] = len(d.get("geneSets") or [])
        if "_error" in d:
            m["_error"] = d["_error"]
        return m
    with ThreadPoolExecutor(8) as ex:
        meta = list(ex.map(one, names))
    json.dump(meta, open(p, "w"), indent=1)
    return meta


def summarize(sym, d):
    rec = {"gene": sym, "symbol_returned": d.get("symbol"),
           "entrez": d.get("ncbiEntrezGeneId"), "per_dataset": {}}
    for a in d.get("associations") or []:
        gsn = a["geneSet"]["name"]
        ds = gsn.rsplit("/", 1)[1] if "/" in gsn else gsn
        slot = rec["per_dataset"].setdefault(ds, [0, 0])
        if (a.get("thresholdValue") or 0) < 0:
            slot[1] += 1
        else:
            slot[0] += 1
    rec["n_associations"] = len(d.get("associations") or [])
    return rec


def fetch_gene(sym):
    p = os.path.join(GDIR, sym + ".json")
    if os.path.exists(p):
        try:
            if "_error" not in json.load(open(p)):
                return "skip"
        except Exception:
            pass
    d = get("%s/gene/%s?showAssociations=true" % (API, urllib.parse.quote(sym)))
    if "_error" in d or not d.get("symbol"):
        rec = {"gene": sym, "_error": d.get("_error", "no_record")}
    else:
        rec = summarize(sym, d)
    json.dump(rec, open(p, "w"))
    return "err" if "_error" in rec else "ok"


def main(budget_s=100, workers=32):
    os.makedirs(GDIR, exist_ok=True)
    fetch_meta()
    t0 = time.time()
    todo = [g for g in gene_set()
            if not (os.path.exists(os.path.join(GDIR, g + ".json")) and
                    "_error" not in open(os.path.join(GDIR, g + ".json")).read())]
    done = []
    with ThreadPoolExecutor(workers) as ex:
        for i in range(0, len(todo), workers):
            if time.time() - t0 > budget_s:
                break
            done += list(ex.map(fetch_gene, todo[i:i + workers]))
    left = len(todo) - len(done)
    json.dump({"source": API, "n_genes_target": len(gene_set()),
               "remaining": left, "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%S")},
              open(os.path.join(OUT, "manifest.json"), "w"), indent=1)
    print("fetched", done.count("ok"), "err", done.count("err"), "remaining", left)


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 100)

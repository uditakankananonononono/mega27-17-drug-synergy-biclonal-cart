#!/usr/bin/env python3
"""Fetch GlyGen protein-detail records (glycosylation sites: reported, reported_with_glycan, predicted) for the 205-gene set + 4 intracellular calibration controls. Accessions reused from data/epitope/uniprot_<SYM>.json caches (UniProt REST, epitope audit). Raw -> data/glygen/ (untracked). Resumable."""
import json, os, time, urllib.request, urllib.error

RAW = "data/glygen"
UA = {"User-Agent": "Mozilla/5.0 mega27-laneF"}
SETS = json.load(open("results/depmap_gene_sets.json"))
GENES = {g: grp for grp in ("gated", "background", "references") for g in SETS[grp]}
for g in SETS["positive_controls"]:
    GENES.setdefault(g, "positive_controls")


def get(url, dest):
    p = os.path.join(RAW, dest)
    if not os.path.exists(p):
        for a in range(4):
            try:
                req = urllib.request.Request(url, headers=UA)
                with urllib.request.urlopen(req, timeout=90) as r:
                    open(p, "wb").write(r.read())
                break
            except urllib.error.HTTPError as e:
                if e.code in (400, 404):
                    open(p, "w").write("null")
                    break
                time.sleep(3 * (a + 1))
            except Exception:
                time.sleep(3 * (a + 1))
    return p


def accession(sym):
    p = f"data/epitope/uniprot_{sym}.json"
    if not os.path.exists(p):
        q = f"gene_exact:{sym}+AND+organism_id:9606+AND+reviewed:true"
        f = "accession,gene_primary,length"
        url = f"https://rest.uniprot.org/uniprotkb/search?query={q}&fields={f}&format=json&size=5"
        d = json.load(open(get(url, f"uniprot_{sym}.json")))
    else:
        d = json.load(open(p))
    res = (d or {}).get("results", [])
    for r in res:
        g = (r.get("genes") or [{}])[0].get("geneName", {}).get("value")
        if g == sym:
            return r.get("primaryAccession")
    return res[0].get("primaryAccession") if res else None


def main():
    os.makedirs(RAW, exist_ok=True)
    n_done = n_new = n_null = 0
    for i, (sym, grp) in enumerate(sorted(GENES.items())):
        dest = f"gg_{sym}.json"
        if os.path.exists(os.path.join(RAW, dest)):
            n_done += 1
            continue
        acc = accession(sym)
        if not acc:
            open(os.path.join(RAW, dest), "w").write(json.dumps({"error": "no_uniprot"}))
            n_null += 1
            continue
        get(f"https://api.glygen.org/protein/detail/{acc}/", dest)
        n_new += 1
        if (i + 1) % 25 == 0:
            print(f"{i+1}/{len(GENES)} ({n_new} new)", flush=True)
    print(f"done: {len(GENES)} genes, {n_new} fetched this run, {n_done} cached, {n_null} no accession")


if __name__ == "__main__":
    main()

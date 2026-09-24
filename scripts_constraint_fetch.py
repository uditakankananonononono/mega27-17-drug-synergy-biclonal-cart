#!/usr/bin/env python3
"""Fetch antigen-escape evidence for gated antigens, background surfaceome, CAR-T references and essential controls:
gnomAD v4 GraphQL LoF constraint, HGNC REST gene groups, Reactome ContentService pathway mappings.
Raw -> data/constraint/ (untracked). Resumable."""
import json, os, time, urllib.request, urllib.error
RAW = "data/constraint"
UA = {"User-Agent": "Mozilla/5.0 mega27-laneF", "Accept": "application/json"}
GQ = ("{ gene(gene_symbol:\"%s\", reference_genome: GRCh38){ gene_id symbol gnomad_constraint"
      "{ pli oe_lof oe_lof_upper oe_mis mis_z lof_z exp_lof obs_lof } } }")


def fetch(url, dest, data=None):
    p = os.path.join(RAW, dest)
    if os.path.exists(p):
        return json.load(open(p))
    for a in range(5):
        try:
            h = dict(UA)
            if data is not None:
                h["Content-Type"] = "application/json"
            req = urllib.request.Request(url, data=data, headers=h)
            with urllib.request.urlopen(req, timeout=60) as r:
                raw = r.read()
            obj = json.loads(raw) if raw.strip() else None
            json.dump(obj, open(p, "w"))
            return obj
        except urllib.error.HTTPError as e:
            if e.code == 404:
                json.dump(None, open(p, "w"))
                return None
            time.sleep(4 * (a + 1))
        except Exception:
            time.sleep(4 * (a + 1))
    return None


def accession(sym):
    p = f"data/epitope/uniprot_{sym}.json"
    if not os.path.exists(p):
        return None
    d = json.load(open(p)) or {}
    for r in d.get("results", []):
        g = (r.get("genes") or [{}])[0].get("geneName", {}).get("value")
        if g == sym:
            return r["primaryAccession"]
    res = d.get("results", [])
    return res[0]["primaryAccession"] if res else None


def main():
    os.makedirs(RAW, exist_ok=True)
    s = json.load(open("results/depmap_gene_sets.json"))
    genes = {}
    for grp in ("gated", "background", "references", "positive_controls"):
        for g in s[grp]:
            genes.setdefault(g, grp)
    order = sorted(genes, reverse=bool(os.environ.get("REV")))
    for i, g in enumerate(order):
        fetch("https://gnomad.broadinstitute.org/api", f"gnomad_{g}.json",
              json.dumps({"query": GQ % g}).encode())
        fetch(f"https://rest.genenames.org/fetch/symbol/{g}", f"hgnc_{g}.json")
        acc = accession(g)
        if acc is None:
            h = fetch(f"https://rest.genenames.org/fetch/symbol/{g}", f"hgnc_{g}.json")
            docs = ((h or {}).get("response") or {}).get("docs") or []
            ids = docs[0].get("uniprot_ids") if docs else None
            acc = ids[0] if ids else None
        if acc:
            fetch(f"https://reactome.org/ContentService/data/mapping/UniProt/{acc}/pathways?species=9606",
                  f"reactome_{g}.json")
        json.dump({"gene": g, "group": genes[g], "accession": acc}, open(os.path.join(RAW, f"meta_{g}.json"), "w"))
        if i % 25 == 0:
            print(i, g, flush=True)
    print("done", len(genes))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Fetch CAR-epitope structural evidence (UniProt REST topology + PDB xrefs, AlphaFold DB pLDDT, RCSB PDB GraphQL entity descriptions) for gated antigens, background surfaceome and CAR-T references. Raw -> data/epitope/ (untracked). Resumable."""
import json, os, re, time, urllib.request, urllib.error
RAW = "data/epitope"
UA = {"User-Agent": "Mozilla/5.0 mega27-laneF"}


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


def uniprot(sym):
    q = f"gene_exact:{sym}+AND+organism_id:9606+AND+reviewed:true"
    f = "accession,gene_primary,length,ft_topo_dom,ft_transmem,ft_signal,ft_chain,ft_lipid,xref_pdb"
    url = f"https://rest.uniprot.org/uniprotkb/search?query={q}&fields={f}&format=json&size=5"
    d = json.load(open(get(url, f"uniprot_{sym}.json")))
    res = (d or {}).get("results", [])
    for r in res:
        g = (r.get("genes") or [{}])[0].get("geneName", {}).get("value")
        if g == sym:
            return r
    return res[0] if res else None


def ectodomain(r):
    """Set of extracellular residue positions (1-based) and topology class."""
    feats = r.get("features", [])
    ext = set()
    for f in feats:
        if f["type"] == "Topological domain" and "Extracellular" in f.get("description", ""):
            s, e = f["location"]["start"]["value"], f["location"]["end"]["value"]
            if isinstance(s, int) and isinstance(e, int):
                ext |= set(range(s, e + 1))
    ntm = sum(1 for f in feats if f["type"] == "Transmembrane")
    if ext:
        return ext, ("multi-pass" if ntm > 1 else "single-pass"), ntm
    gpi = [f for f in feats if f["type"] == "Lipidation" and "GPI" in f.get("description", "")]
    chains = [f for f in feats if f["type"] == "Chain"]
    if gpi and chains:
        e = gpi[0]["location"]["start"]["value"]
        # use the mature chain that carries the GPI site (e.g. MSLN: mesothelin, not MPF)
        hit = [c for c in chains if isinstance(c["location"]["start"]["value"], int)
               and isinstance(c["location"]["end"]["value"], int) and isinstance(e, int)
               and c["location"]["start"]["value"] <= e <= c["location"]["end"]["value"]]
        hit.sort(key=lambda c: c["location"]["end"]["value"] - c["location"]["start"]["value"])  # shortest = mature form
        s = (hit or chains)[0]["location"]["start"]["value"]
        if isinstance(s, int) and isinstance(e, int):
            return set(range(s, e + 1)), "GPI", ntm
    return set(), "unannotated", ntm


def alphafold_plddt(acc):
    meta = json.load(open(get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", f"af_{acc}.json")))
    if not meta:
        return None
    m = meta[0]
    p = get(m["pdbUrl"], f"af_{acc}.pdb")
    from Bio.PDB import PDBParser
    try:
        s = PDBParser(QUIET=True).get_structure(acc, p)
    except Exception:
        return None
    pl = {}
    for res in s[0]["A"]:
        ca = res["CA"] if "CA" in res else None
        if ca is not None:
            pl[res.id[1]] = float(ca.get_bfactor())
    return pl


ABRE = re.compile(r"antibod|\bfab\b|fab fragment|scfv|nanobod|heavy chain|light chain|immunoglobulin|vhh|igg", re.I)


def pdb_entries(r, ext):
    out = []
    for x in r.get("uniProtKBCrossReferences", []):
        if x.get("database") != "PDB":
            continue
        props = {p["key"]: p["value"] for p in x.get("properties", [])}
        cov = set()
        for seg in props.get("Chains", "").split(","):
            m = re.search(r"=(\d+)-(\d+)", seg)
            if m:
                cov |= set(range(int(m.group(1)), int(m.group(2)) + 1))
        ov = cov & ext
        out.append({"pdb": x["id"], "method": props.get("Method", ""), "cover": cov, "ecto_overlap": len(ov)})
    return out


def rcsb_post(chunk):
    q = "{entries(entry_ids:%s){rcsb_id polymer_entities{rcsb_polymer_entity{pdbx_description}}}}"
    body = json.dumps({"query": q % json.dumps(chunk)}).encode()
    hdr = dict(UA)
    hdr["Content-Type"] = "application/json"
    req = urllib.request.Request("https://data.rcsb.org/graphql", data=body, headers=hdr)
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def rcsb_descriptions(ids):
    out = {}
    todo = sorted(ids)
    for k in range(0, len(todo), 100):
        chunk = todo[k:k + 100]
        dest = os.path.join(RAW, "rcsb_%s_%d.json" % (chunk[0], len(chunk)))
        if not os.path.exists(dest):
            json.dump(rcsb_post(chunk), open(dest, "w"))
        d = json.load(open(dest))
        for e in (d.get("data") or {}).get("entries") or []:
            if e:
                ents = e.get("polymer_entities") or []
                out[e["rcsb_id"]] = [(x.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "" for x in ents]
    return out


def main():
    os.makedirs(RAW, exist_ok=True)
    gs = json.load(open("results/depmap_gene_sets.json"))
    groups = [("gated", gs["gated"]), ("background", gs["background"]), ("reference", gs["references"])]
    rows, ab_ids = [], set()
    for grp, genes in groups:
        for sym in genes:
            r = uniprot(sym)
            if not r:
                rows.append({"gene": sym, "group": grp, "found": 0})
                continue
            ext, topo, ntm = ectodomain(r)
            pl = alphafold_plddt(r["primaryAccession"]) if ext else None
            ents = pdb_entries(r, ext)
            for e in ents:
                if e["ecto_overlap"] >= 10:
                    ab_ids.add(e["pdb"])
            rows.append({"gene": sym, "group": grp, "found": 1, "acc": r["primaryAccession"],
                         "topology": topo, "n_tm": ntm, "ecto": sorted(ext), "plddt": pl,
                         "pdb": [{k: v for k, v in e.items() if k != "cover"} | {"cover": sorted(e["cover"] & ext)} for e in ents]})
            print(sym, grp, topo, len(ext), flush=True)
    desc = rcsb_descriptions(ab_ids)
    json.dump({"rows": rows, "rcsb_desc": desc}, open(os.path.join(RAW, "epitope_raw.json"), "w"))
    print("genes", len(rows), "rcsb entries", len(desc))


if __name__ == "__main__":
    main()

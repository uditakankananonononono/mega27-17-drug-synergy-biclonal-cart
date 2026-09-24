#!/usr/bin/env python3
"""Fetch CPTAC protein expression for the gated-antigen audit set via cBioPortal API.

Studies (curated CPTAC pan-cancer publications on cBioPortal):
  brca_cptac_2020 coad_cptac_2019 gbm_cptac_2021 luad_cptac_2020
  lusc_cptac_2021 paad_cptac_2021 ucec_cptac_2020
Profile per study: {study}_protein_quantification (LOG2-VALUE).

Output: results/cptac_protein_rows.csv  (study, sample_id, symbol, log2_protein)
        results/cptac_fetch_meta.json   (gene map, sample counts, missing symbols)
"""
import json, urllib.request, sys, csv, time

API = "https://www.cbioportal.org/api"
STUDIES = ["brca_cptac_2020","coad_cptac_2019","gbm_cptac_2021","luad_cptac_2020",
           "lusc_cptac_2021","paad_cptac_2021","ucec_cptac_2020"]

def post(path, payload):
    req = urllib.request.Request(API+path, data=json.dumps(payload).encode(),
        headers={"Content-Type":"application/json","accept":"application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)

def get(path):
    req = urllib.request.Request(API+path, headers={"accept":"application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)

def main():
    sets = json.load(open("results/depmap_gene_sets.json"))
    symbols = sorted(set(sets["gated"])|set(sets["background"])|set(sets["references"]))
    genes = post("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL&projection=SUMMARY",
                 symbols)
    sym2entrez = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in genes}
    missing = [s for s in symbols if s not in sym2entrez]
    entrez = sorted(sym2entrez.values())
    print(f"mapped {len(sym2entrez)}/{len(symbols)} symbols; missing: {missing}", flush=True)

    rows, sample_counts = [], {}
    for st in STUDIES:
        prof = f"{st}_protein_quantification"
        lists = get(f"/studies/{st}/sample-lists?projection=SUMMARY")
        sl = next((l for l in lists if l["sampleListId"] == f"{st}_protein_quantification"), None)
        if sl is None:
            sl = next((l for l in lists if "protein_quantification" in l["sampleListId"]), None)
        if sl is None:
            sl = next(l for l in lists if "protein" in l["sampleListId"]
                      and "phospho" not in l["sampleListId"] and "acetyl" not in l["sampleListId"])
        ids = get(f"/sample-lists/{sl['sampleListId']}/sample-ids")
        sample_counts[st] = len(ids)
        try:
            data = post(f"/molecular-profiles/{prof}/molecular-data/fetch?projection=SUMMARY",
                        {"entrezGeneIds": entrez, "sampleListId": sl["sampleListId"]})
        except Exception as e:
            print(f"{st}: FETCH FAILED {e}", flush=True); continue
        ent2sym = {v:k for k,v in sym2entrez.items()}
        n = 0
        for rec in data:
            if rec.get("value") in (None,"NA",""): continue
            rows.append((st, rec["sampleId"], ent2sym[rec["entrezGeneId"]], float(rec["value"])))
            n += 1
        print(f"{st}: {sample_counts[st]} samples, {n} quantified gene-sample values", flush=True)
        time.sleep(1)

    with open("results/cptac_protein_rows.csv","w",newline="") as f:
        w = csv.writer(f); w.writerow(["study","sample_id","symbol","log2_protein"]); w.writerows(rows)
    json.dump({"studies": STUDIES, "sample_counts": sample_counts,
               "symbols_queried": len(symbols), "symbols_mapped": len(sym2entrez),
               "missing_symbols": missing, "n_rows": len(rows),
               "source": "cBioPortal API /api/molecular-profiles/{study}_protein_quantification/molecular-data/fetch",
               "fetched": time.strftime("%Y-%m-%d")},
              open("results/cptac_fetch_meta.json","w"), indent=1)
    print("total rows:", len(rows))

if __name__ == "__main__":
    main()

"""Annotate NCI-ALMANAC drugs with ChEMBL assay-level bioactivity.

NSC -> name (PubChem PUG synonyms) -> ChEMBL molecule -> per-assay activities
(pchembl). Each ChEMBL assay (CHEMBL accession) used is one accession-level
dataset; the analysis uses them to annotate drug-target mechanisms for the
synergy model. Output: data/chembl/<nsc>.json + results/chembl_annotation.json
"""
import json, pathlib, time, urllib.request, urllib.parse

OUT = pathlib.Path("data/chembl"); OUT.mkdir(exist_ok=True)

def get(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)

import pandas as pd
df = pd.read_csv("data/almanac_synergy.tsv", sep="\t")
nscs = sorted(set(df.NSC1) | set(df.NSC2))
print(len(nscs), "unique NSC codes", flush=True)
manifest = {}
for nsc in nscs:
    f = OUT / f"{nsc}.json"
    if f.exists():
        manifest[nsc] = json.loads(f.read_text()).get("n_assays", 0)
        continue
    rec = {"nsc": int(nsc)}
    try:
        syn = get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/NSC{nsc}/synonyms/JSON")
        names = syn["InformationList"]["Information"][0].get("Synonym", [])
        rec["names"] = names[:5]
        rec["name"] = next((n for n in names if not n.replace("-", "").isdigit()), names[0] if names else None)
        if rec["name"]:
            q = urllib.parse.quote(rec["name"])
            mol = get(f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name__iexact={q}&limit=1")
            if not mol["molecules"]:
                mol = get(f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?q={q}&limit=1")
            if mol["molecules"]:
                m = mol["molecules"][0]
                rec["chembl_id"] = m["molecule_chembl_id"]
                rec["pref_name"] = m.get("pref_name")
                acts = get(f"https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id={m['molecule_chembl_id']}&pchembl_value__isnull=false&limit=1000")
                rec["n_assays"] = acts["page_meta"]["total_count"]
                seen = {}
                for a in acts["activities"]:
                    aid = a["assay_chembl_id"]
                    if aid not in seen:
                        seen[aid] = {"assay": aid, "target": a.get("target_pref_name"),
                                     "type": a.get("standard_type"),
                                     "pchembl_max": float(a["pchembl_value"])}
                    else:
                        seen[aid]["pchembl_max"] = max(seen[aid]["pchembl_max"], float(a["pchembl_value"]))
                rec["assays"] = list(seen.values())
        f.write_text(json.dumps(rec))
        manifest[nsc] = rec.get("n_assays", 0)
        print(nsc, rec.get("name"), rec.get("chembl_id"), rec.get("n_assays"), flush=True)
    except Exception as e:
        print("FAIL", nsc, repr(e)[:90], flush=True)
        time.sleep(3)
    time.sleep(0.4)
json.dump({"n_drugs_annotated": sum(1 for v in manifest.values() if v),
           "n_unique_assays_total": None,  # filled by analysis
           "per_drug_assays": manifest}, open("results/chembl_annotation.json", "w"), indent=1)
print("done")

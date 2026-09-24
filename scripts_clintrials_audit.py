"""ClinicalTrials.gov landscape audit (API v2) of the 205-gene item-17 set.
Question: how much verified, target-directed oncology development already exists
for the gated surfaceome and the seven AND-gate antigens?
Tiers per gene: strict (symbol token in intervention names), strict_onc (strict +
oncology condition filter; PRIMARY), broad (full-text + modality; sensitivity),
alias (curated aliases for 7 antigens + 4 references).
Acronym-collision curation: every gene with strict_onc>0 was title-inspected;
collisions/biomarker-only hits are zeroed via the committed results/clintrials_curation.csv.
Calibration gates BEFORE claims: CD19 verified >= 100 (positive control);
POLR2A + RPS3 verified = 0 (text-noise controls). PCNA/PSMA1 control collisions documented.
Writes results/clintrials_audit.json + results/clintrials_per_gene.csv.
"""
import csv, json, os
from scipy.stats import fisher_exact, mannwhitneyu

ACTIVE = {"RECRUITING", "ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION", "NOT_YET_RECRUITING"}
AND_GATE = ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6"]
CACHE = "data/clintrials"

def load(gene, tier):
    p = os.path.join(CACHE, f"{tier}_{gene}.json")
    if not os.path.exists(p):
        return {"total": 0, "records": [], "capped": False}
    return json.load(open(p))

def per_gene(gene):
    s = load(gene, "strict")
    so = load(gene, "strict_onc")
    b = load(gene, "broad")
    a = load(gene, "alias")
    act = lambda recs: sum(1 for r in recs if r.get("status") in ACTIVE)
    union = {r["nct"] for r in so["records"]} | {r["nct"] for r in b["records"]} | {r["nct"] for r in a["records"]}
    return {"gene": gene,
            "strict_total": s["total"], "strict_onc_total": so["total"],
            "strict_onc_active": act(so["records"]), "strict_onc_capped": so.get("capped", False),
            "broad_total": b["total"], "broad_capped": b.get("capped", False),
            "alias_total": a["total"], "alias_active": act(a["records"]),
            "alias_capped": a.get("capped", False),
            "n_nct_union": len(union), "union_ncts": sorted(union)}

def main():
    sets = json.load(open("results/depmap_gene_sets.json"))
    classes = {g: "gated" for g in sets["gated"]}
    classes.update({g: "background" for g in sets["background"]})
    classes.update({g: "reference" for g in sets["references"]})
    classes.update({g: "positive_control" for g in sets["positive_controls"]})
    cur = {}
    with open("results/clintrials_curation.csv") as f:
        for r in csv.DictReader(f):
            cur[r["gene"]] = r
    rows = []
    for g, cls in classes.items():
        r = per_gene(g)
        r["class"] = cls
        c = cur.get(g)
        r["curated"] = bool(c)
        r["verdict"] = c["verdict"] if c else ("silent" if r["strict_onc_total"] == 0 else "uninspected")
        r["verified"] = int(c["verified_count"]) if c else r["strict_onc_total"]
        if c and c["verdict"] in ("collision", "biomarker"):
            r["verified_active"] = 0
        elif c and c["verdict"] == "mixed":
            r["verified_active"] = 0  # conservative: only AOH1996 genuine; control gene, not in stats
        else:
            r["verified_active"] = r["strict_onc_active"]
        rows.append(r)
    assert all(r["verdict"] != "uninspected" for r in rows), "inspection gap"

    gated = [r for r in rows if r["class"] == "gated"]
    bg = [r for r in rows if r["class"] == "background"]
    refs = {r["gene"]: r for r in rows if r["class"] in ("reference", "positive_control")}

    calibration = {
        "cd19_verified": refs["CD19"]["verified"], "cd19_gate": refs["CD19"]["verified"] >= 100,
        "noise_controls_verified_sum": refs["POLR2A"]["verified"] + refs["RPS3"]["verified"],
        "noise_controls": {"POLR2A": refs["POLR2A"]["verified"], "RPS3": refs["RPS3"]["verified"],
                           "PCNA": refs["PCNA"]["verified"], "PSMA1": refs["PSMA1"]["verified"]},
    }
    calibration["noise_gate"] = calibration["noise_controls_verified_sum"] == 0
    calibration["passed"] = calibration["cd19_gate"] and calibration["noise_gate"]

    def med(v):
        s = sorted(v); n = len(s)
        return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2

    def mw(key):
        x = [r[key] for r in gated]; y = [r[key] for r in bg]
        u = mannwhitneyu(x, y, alternative="two-sided")
        return {"p": u.pvalue, "gated_median": med(x), "bg_median": med(y),
                "gated_mean": sum(x) / len(x), "bg_mean": sum(y) / len(y)}

    def fisher_any(key):
        ta = sum(1 for r in gated if r[key] > 0); tb = sum(1 for r in bg if r[key] > 0)
        _, p = fisher_exact([[ta, len(gated) - ta], [tb, len(bg) - tb]])
        return {"p": p, "gated_any": f"{ta}/{len(gated)}", "bg_any": f"{tb}/{len(bg)}"}

    stats = {
        "verified_mw": mw("verified"),
        "verified_fisher": fisher_any("verified"),
        "verified_active_mw": mw("verified_active"),
        "strict_onc_raw_mw": mw("strict_onc_total"),
        "strict_raw_mw": mw("strict_total"),
        "broad_mw": mw("broad_total"),
    }

    alias_rows = [r for r in rows if r["gene"] in
                  ["CA9", "CA12", "CLDN18", "CLDN6", "MSLN", "PSCA", "SLC39A6",
                   "CD19", "ERBB2", "FOLR1", "TNFRSF17"]]
    alias_table = [{"gene": r["gene"], "class": ("and_gate" if r["gene"] in AND_GATE else "reference"),
                    "total": r["alias_total"], "active": r["alias_active"], "capped": r["alias_capped"],
                    "verified_strict": r["verified"]}
                   for r in sorted(alias_rows, key=lambda r: -r["alias_total"])]

    all_ncts = set()
    for r in rows:
        all_ncts |= set(r["union_ncts"])

    n_curated = sum(1 for r in rows if r["curated"])
    n_collision = sum(1 for r in rows if r["verdict"] in ("collision", "biomarker"))
    audit = {
        "question": "how much verified target-directed oncology development exists for the gated surfaceome and AND-gate antigens?",
        "source": "ClinicalTrials.gov API v2; strict_onc (intervention-name token + oncology condition) primary, curated by title inspection; broad + alias tiers as sensitivity/depth",
        "n_genes": len(rows), "classes": {c: sum(1 for r in rows if r["class"] == c) for c in set(classes.values())},
        "calibration": calibration,
        "stats": stats,
        "alias_table": alias_table,
        "curation": {"n_inspected": n_curated, "n_collision_or_biomarker": n_collision,
                     "n_genuine": sum(1 for r in rows if r["verdict"] == "genuine"),
                     "n_mixed": sum(1 for r in rows if r["verdict"] == "mixed"),
                     "genuine_gated": sorted(r["gene"] for r in gated if r["verified"] > 0),
                     "genuine_background": sorted(r["gene"] for r in bg if r["verified"] > 0)},
        "and_gate_summary": {
            "n": 7,
            "verified_strict": {r["gene"]: r["verified"] for r in rows if r["gene"] in AND_GATE},
            "alias_total": {r["gene"]: r["alias_total"] for r in rows if r["gene"] in AND_GATE},
            "alias_active": {r["gene"]: r["alias_active"] for r in rows if r["gene"] in AND_GATE}},
        "datasets": {"unique_nct_globally": len(all_ncts),
                     "per_gene_union_sum": sum(r["n_nct_union"] for r in rows)},
    }
    with open("results/clintrials_audit.json", "w") as f:
        json.dump(audit, f, indent=1)
    with open("results/clintrials_per_gene.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "class", "strict_total", "strict_onc_total", "verified",
                    "strict_onc_active", "broad_total", "alias_total", "alias_active",
                    "verdict", "n_nct_union"])
        for r in sorted(rows, key=lambda r: (r["class"], r["gene"])):
            w.writerow([r["gene"], r["class"], r["strict_total"], r["strict_onc_total"],
                        r["verified"], r["strict_onc_active"], r["broad_total"],
                        r["alias_total"], r["alias_active"], r["verdict"], r["n_nct_union"]])
    print(json.dumps({"calibration": calibration, "stats": stats,
                      "unique_nct": len(all_ncts), "curation": audit["curation"]}, indent=1))

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""CAR-epitope structural audit: ectodomain size, AlphaFold confidence, experimental
coverage and antibody co-structures for gated antigens vs background vs CAR-T references.
Input: data/epitope/epitope_raw.json (scripts_epitope_fetch.py).
Writes results/epitope_audit.json + results/epitope_per_gene.csv."""
import json, re
import numpy as np, pandas as pd
from scipy import stats

AB = re.compile(r"antibod|\bfab\b|scfv|nanobod|heavy chain|light chain|\bvhh\b", re.I)
AND_GATE = ["CLDN18", "MSLN", "CA9", "CA12", "CLDN6", "PSCA", "SLC39A6"]


def longest_run(pos):
    best = cur = 0
    prev = None
    for p in pos:
        cur = cur + 1 if prev is not None and p == prev + 1 else 1
        best = max(best, cur)
        prev = p
    return best


def main():
    d = json.load(open("data/epitope/epitope_raw.json"))
    desc = d["rcsb_desc"]
    rows = []
    for r in d["rows"]:
        if not r.get("found"):
            rows.append({"gene": r["gene"], "group": r["group"], "found": 0})
            continue
        ext = sorted(r["ecto"])
        pl = r.get("plddt") or {}
        ext_pl = [pl[str(p)] for p in ext if str(p) in pl]
        cov = set()
        ab = 0
        nstruct = 0
        for e in r.get("pdb", []):
            if e["ecto_overlap"] >= 10:
                nstruct += 1
                cov |= set(e["cover"])
                if any(AB.search(x or "") for x in desc.get(e["pdb"], [])):
                    ab += 1
        rows.append({"gene": r["gene"], "group": r["group"], "found": 1, "topology": r["topology"],
                     "n_tm": r["n_tm"], "ecto_len": len(ext), "ecto_longest": longest_run(ext),
                     "mean_plddt_ecto": float(np.mean(ext_pl)) if ext_pl else None,
                     "exp_cov_ecto": len(cov & set(ext)) / len(ext) if ext else None,
                     "n_pdb_ecto": nstruct, "n_ab_structures": ab,
                     "has_structure": int(nstruct > 0), "has_ab_structure": int(ab > 0)})
    df = pd.DataFrame(rows)
    df.to_csv("results/epitope_per_gene.csv", index=False)
    out = {"n_genes": int(df.found.sum()), "n_unannotated": int((df.topology == "unannotated").sum())}
    ann = df[(df.found == 1) & (df.topology != "unannotated")]
    g = ann[ann.group == "gated"]; b = ann[ann.group == "background"]; ref = ann[ann.group == "reference"]
    out["groups"] = {}
    for name, s in (("gated", g), ("background", b), ("reference", ref)):
        out["groups"][name] = {"n": int(len(s)),
                               "ecto_len_median": float(s.ecto_len.median()),
                               "mean_plddt_median": float(s.mean_plddt_ecto.dropna().median()),
                               "exp_cov_median": float(s.exp_cov_ecto.dropna().median()),
                               "frac_has_structure": float(s.has_structure.mean()),
                               "frac_has_ab_structure": float(s.has_ab_structure.mean()),
                               "ab_structures_total": int(s.n_ab_structures.sum())}
    tests = {}
    for col in ("ecto_len", "mean_plddt_ecto", "exp_cov_ecto", "n_ab_structures"):
        x, y = g[col].dropna(), b[col].dropna()
        u = stats.mannwhitneyu(x, y, alternative="two-sided")
        tests[col] = {"gated_median": float(x.median()), "bg_median": float(y.median()),
                      "auroc_gated_vs_bg": float(u.statistic / (len(x) * len(y))), "p": float(u.pvalue)}
    for col in ("has_structure", "has_ab_structure"):
        t = pd.crosstab(ann.group.isin(["gated"]).map({True: "gated", False: "bg"}), ann[col])
        o, p = stats.fisher_exact(t.values)
        tests[col] = {"table": t.values.tolist(), "or": float(o), "p": float(p)}
    out["gated_vs_background"] = tests
    # positive controls: the four approved-CAR reference antigens must show antibody co-structures
    out["reference_controls"] = df[df.group == "reference"][["gene", "has_ab_structure", "n_ab_structures"]].to_dict("records")
    out["and_gate"] = df[df.gene.isin(AND_GATE)][["gene", "topology", "ecto_len", "mean_plddt_ecto",
                                                  "exp_cov_ecto", "n_pdb_ecto", "n_ab_structures"]].to_dict("records")
    json.dump(out, open("results/epitope_audit.json", "w"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str)[:4000])


if __name__ == "__main__":
    main()

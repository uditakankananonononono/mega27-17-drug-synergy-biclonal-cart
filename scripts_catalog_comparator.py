#!/usr/bin/env python3
"""Tier-4 #11 secondary endpoint: published-catalog comparator (prespec results/blinded_benchmark_prespec.json).

Catalog source: PMC12374352 ("Advancements and challenges in CAR-T cell therapy for
solid tumors", Cancer Cell International), full-text XML via Europe PMC:
https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12374352/fullTextXML

EXTRACTION RULE (pre-declared, deterministic):
  1. Antigen tokens from Table 1 ("Common target antigens of solid tumors") rows,
     citation brackets stripped.
  2. Antigen subsection headings under "Common antigen targets in solid tumors" and
     "Novel antigen targets in CAR-T therapy for solid tumors".
  3. Alias map (below) normalizes review names to corpus gene symbols. Non-gene
     antigens (glycolipid GD2, glycoforms), receptor-axis entries (NKG2D), and
     non-specific families (ERBB family) are excluded and logged.
  4. Catalog membership is evaluated only over the 205-gene corpus.

TEST: Fisher exact on catalog membership x clintrials verdict (genuine vs not),
over corpus genes present in both resources. One-sided enrichment direction
(catalog genes more often genuine), alpha per prespec; descriptive only, no
multiplicity family (single secondary test).
"""
import csv, io, json, re, urllib.request
from xml.etree import ElementTree as ET

URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12374352/fullTextXML"
ALIAS = {
    "HER2": "ERBB2", "CEA": "CEACAM5", "CAIX": "CA9", "FRalpha": "FOLR1",
    "PSMA": "FOLH1", "cMET": "MET", "VEGFR2": "KDR", "VEGFR": "KDR",
    "CLDN18.2": "CLDN18", "EGFR VIII": "EGFR", "CD44V6": "CD44",
    "CS1": "SLAMF7", "TEM8/ANTXR1": "ANTXR1", "IL13RA": "IL13RA2",
    "GD2": None, "NKG2D": None, "ERBB family": None, "Glycoforms of antigens": None,
}

def text_of(e):
    return "".join(e.itertext())

def norm(tok):
    tok = tok.strip()
    tok = tok.replace("FR\u03b1", "FRalpha")  # FR alpha
    return ALIAS.get(tok, tok)

xml = urllib.request.urlopen(URL, timeout=60).read()
root = ET.fromstring(xml)
raw, excluded = set(), set()

# 1. Table 1 rows
for tw in root.iter():
    if tw.tag.split("}")[-1] == "table-wrap":
        lab = tw.find("label")
        if lab is not None and text_of(lab).strip() == "Table 1":
            for tr in tw.iter():
                if tr.tag.split("}")[-1] == "tr":
                    cells = [text_of(c) for c in tr if c.tag.split("}")[-1] in ("td", "th")]
                    if len(cells) >= 2 and cells[0].strip() not in ("Tissues",):
                        # 3-cell rows: [tissue, antigen, stage]; 2-cell continuation rows: [antigen, stage]
                        cell = cells[1] if len(cells) == 3 else cells[0]
                        tok = re.sub(r"\[[^\]]*\]", "", cell).strip()
                        if tok and tok != "Targeted antigens":
                            g = norm(tok)
                            (raw if g else excluded).add(g or tok)

# 2. Antigen subsection headings
def walk(e):
    for ch in e:
        if ch.tag.split("}")[-1] == "sec":
            t = ch.find("title")
            title = text_of(t) if t is not None else ""
            if title in ("Common antigen targets in solid tumors",
                         "Novel antigen targets in CAR-T therapy for solid tumors"):
                for sub in ch:
                    if sub.tag.split("}")[-1] == "sec":
                        st = sub.find("title")
                        if st is not None:
                            h = text_of(st).strip()
                            m = re.search(r"\(([^)]+)\)\s*$", h)  # "Mesothelin (MSLN)"
                            tok = m.group(1) if m else h
                            g = norm(tok)
                            (raw if g else excluded).add(g or h)
        walk(ch)
walk(root)

# corpus labels
labels = {}
with open("results/clintrials_per_gene.csv") as f:
    for row in csv.DictReader(f):
        labels[row["gene"]] = row["verdict"]
catalog = {g for g in raw if g in labels}
missed = sorted(g for g in raw if g not in labels)

# Fisher exact (catalog membership vs genuine)
in_cat_gen = sum(1 for g in catalog if labels[g] == "genuine")
in_cat_not = len(catalog) - in_cat_gen
out_gen = sum(1 for g, v in labels.items() if g not in catalog and v == "genuine")
out_not = len(labels) - len(catalog) - out_gen
a, b, c, d = in_cat_gen, in_cat_not, out_gen, out_not
# one-sided Fisher (enrichment) via hypergeometric tail
from math import comb
n = a + b + c + d
K = a + c
N = a + b
def hyper(x):
    return comb(K, x) * comb(n - K, N - x) / comb(n, N)
p = sum(hyper(x) for x in range(a, min(K, N) + 1))
or_ = (a * d) / (b * c) if b * c else None

out = {
    "catalog_source": URL,
    "catalog_size_mapped": len(raw),
    "catalog_in_corpus": len(catalog),
    "catalog_not_in_corpus": missed,
    "excluded_non_gene": sorted(excluded),
    "table": {"catalog_genuine": a, "catalog_not_genuine": b,
              "noncatalog_genuine": c, "noncatalog_not_genuine": d},
    "fisher_one_sided_p": p,
    "odds_ratio": or_,
    "catalog_gene_verdicts": {g: labels[g] for g in sorted(catalog)},
}
json.dump(out, open("results/catalog_comparator.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "catalog_gene_verdicts"}, indent=1))

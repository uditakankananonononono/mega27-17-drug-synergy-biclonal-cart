"""Stream-aggregate 1GB HPA rna_cancer_sample -> per (gene, cancer) mean nTPM."""
import zipfile, csv
from collections import defaultdict

acc = defaultdict(lambda: [0.0, 0])
with zipfile.ZipFile("data/rna_cancer_sample.tsv.zip") as z:
    name = next(n for n in z.namelist() if n.endswith(".tsv"))
    with z.open(name) as fh:
        rdr = csv.DictReader((l.decode("utf-8", "replace") for l in fh), delimiter="\t")
        cols = rdr.fieldnames
        print("columns:", cols, flush=True)
        for row in rdr:
            try:
                tpm = float(row.get("FPKM") or 0)
            except ValueError:
                continue
            key = (row["Gene"], row["Cancer"])
            e = acc[key]; e[0] += tpm; e[1] += 1
with open("data/cancer_rna_mean.tsv", "w") as out:
    out.write("Gene\tCancer\tmean_nTPM\tn_samples\n")
    for (g, c), (s, n) in sorted(acc.items()):
        out.write(f"{g}\t{c}\t{s/n:.3f}\t{n}\n")
print("gene-cancer rows:", len(acc))

"""Stream-aggregate HPA v23 rna_cancer_sample (7 GB TSV inside zip) into
per (gene, cancer) mean FPKM. Fast split-based parser with progress logging."""
import io, zipfile
from collections import defaultdict

acc = defaultdict(lambda: [0.0, 0])
n = 0
with zipfile.ZipFile("data/rna_cancer_sample.tsv.zip") as z:
    name = next(x for x in z.namelist() if x.endswith(".tsv"))
    with z.open(name) as fh:
        f = io.TextIOWrapper(fh, encoding="utf-8", errors="replace")
        cols = f.readline().rstrip("\n").split("\t")
        print("columns:", cols, flush=True)
        gi, ci, vi = cols.index("Gene"), cols.index("Cancer"), cols.index("FPKM")
        for line in f:
            p = line.split("\t")
            try:
                v = float(p[vi])
            except (ValueError, IndexError):
                continue
            e = acc[(p[gi], p[ci].strip())]; e[0] += v; e[1] += 1
            n += 1
            if n % 10_000_000 == 0:
                print("rows", n, flush=True)
with open("data/cancer_rna_mean.tsv", "w") as out:
    out.write("Gene\tCancer\tmean_FPKM\tn_samples\n")
    for (g, c), (s, k) in sorted(acc.items()):
        out.write(f"{g}\t{c}\t{s/k:.3f}\t{k}\n")
print("done rows", n, "gene-cancer", len(acc), flush=True)

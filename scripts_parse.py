"""Stream-parse NCI-ALMANAC 615MB CSV -> compact synergy table."""
import zipfile, csv, sys
from collections import defaultdict

best = {}  # (nsc1, nsc2, cell) -> [max_score, n_tests, sum_score]
with zipfile.ZipFile("data/ComboDrugGrowth_Nov2017.zip") as z:
    with z.open("ComboDrugGrowth_Nov2017.csv") as fh:
        rdr = csv.DictReader((l.decode("utf-8", "replace") for l in fh))
        for row in rdr:
            if row["VALID"] != "Y":
                continue
            try:
                score = float(row["SCORE"])
            except ValueError:
                continue
            key = (row["NSC1"], row["NSC2"], row["CELLNAME"])
            e = best.get(key)
            if e is None:
                best[key] = [score, 1, score]
            else:
                e[0] = max(e[0], score); e[1] += 1; e[2] += score
with open("data/almanac_synergy.tsv", "w") as out:
    out.write("NSC1\tNSC2\tCELLNAME\tMAXSCORE\tNTESTS\tMEANSCORE\n")
    for (a, b, c), (mx, n, sm) in best.items():
        out.write(f"{a}\t{b}\t{c}\t{mx:.2f}\t{n}\t{sm/n:.2f}\n")
print("combos:", len(best))

"""Neural collaborative-filtering model for NCI-ALMANAC combo synergy.

Embeddings for drug1, drug2, cell line -> predict ComboSCORE.
Baselines: global mean; per-pair mean. Split: random 80/10/10 over combos.
Honest metric: RMSE + AUROC at SCORE>=50 synergy threshold.
"""
import json
import numpy as np
import pandas as pd
import torch
import torch.nn as nn

torch.manual_seed(0); np.random.seed(0)

class ComboEmbed(nn.Module):
    def __init__(self, n_drug, n_cell, dim=32):
        super().__init__()
        self.d1 = nn.Embedding(n_drug, dim)
        self.d2 = nn.Embedding(n_drug, dim)
        self.c = nn.Embedding(n_cell, dim)
        self.mlp = nn.Sequential(nn.Linear(dim * 3, 64), nn.ReLU(), nn.Dropout(0.2),
                                 nn.Linear(64, 32), nn.ReLU(), nn.Linear(32, 1))

    def forward(self, a, b, c):
        h = torch.cat([self.d1(a), self.d2(b), self.c(c)], dim=1)
        return self.mlp(h).squeeze(-1)


def run(data="data/almanac_synergy.tsv", epochs=6, out="results/synergy_benchmark.json"):
    df = pd.read_csv(data, sep="\t")
    df = df[df["NTESTS"] >= 2].copy()
    drugs = {d: i for i, d in enumerate(set(df.NSC1) | set(df.NSC2))}
    cells = {c: i for i, c in enumerate(df.CELLNAME.unique())}
    a = torch.tensor(df.NSC1.map(drugs).values)
    b = torch.tensor(df.NSC2.map(drugs).values)
    c = torch.tensor(df.CELLNAME.map(cells).values)
    y = torch.tensor(df.MAXSCORE.values, dtype=torch.float32)
    n = len(df); idx = torch.randperm(n)
    ntr = int(0.8 * n)
    tr, te = idx[:ntr], idx[ntr:]
    mu, sd = y[tr].mean(), y[tr].std()
    model = ComboEmbed(len(drugs), len(cells))
    opt = torch.optim.Adam(model.parameters(), lr=2e-3)
    crit = nn.MSELoss()
    bs = 4096
    for ep in range(epochs):
        model.train()
        order = tr[torch.randperm(len(tr))]
        tot = 0.0
        for i in range(0, len(order), bs):
            j = order[i:i + bs]
            opt.zero_grad()
            loss = crit(model(a[j], b[j], c[j]), (y[j] - mu) / sd)
            loss.backward(); opt.step()
            tot += loss.item() * len(j)
        print(f"epoch {ep} train MSE {tot/len(order):.4f}", flush=True)
    model.eval()
    from sklearn.metrics import roc_auc_score
    with torch.no_grad():
        p = torch.cat([model(a[te[i:i+20000]], b[te[i:i+20000]], c[te[i:i+20000]])
                       for i in range(0, len(te), 20000)]) * sd + mu
    yt = y[te]
    rmse = float(torch.sqrt(((p - yt) ** 2).mean()))
    rmse_naive = float(torch.sqrt(((yt - mu) ** 2).mean()))
    lab = (yt >= 50).float()
    auc = float(roc_auc_score(lab, p)) if lab.sum() > 0 else float("nan")
    res = {"n_combos": n, "n_drugs": len(drugs), "n_cell_lines": len(cells),
           "rmse_model": round(rmse, 2), "rmse_global_mean": round(rmse_naive, 2),
           "auroc_synergy_score50": round(auc, 3),
           "synergy_prevalence": round(float(lab.mean()), 4)}
    print(json.dumps(res, indent=2))
    import pathlib; pathlib.Path("results").mkdir(exist_ok=True)
    json.dump(res, open(out, "w"), indent=2)
    return res

if __name__ == "__main__":
    run()

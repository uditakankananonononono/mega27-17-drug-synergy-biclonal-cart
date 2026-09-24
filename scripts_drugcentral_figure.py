"""Figure for the DrugCentral engagement audit -> paper/figs/fig_drugcentral.pdf."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = json.load(open("results/drugcentral_engagement.json"))
T, S = D["tests"], D["summary"]
AG = D["and_gate"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.6))
metrics = [("any", "any activity"), ("tclin", "Tclin target"), ("moa", "MOA record")]
x = range(len(metrics)); w = 0.36
gv = [S["gated"][m] / S["gated"]["n"] for m, _ in metrics]
bv = [S["background"][m] / S["background"]["n"] for m, _ in metrics]
ax1.bar([i - w/2 for i in x], gv, w, label=f"gated (n={S['gated']['n']})", color="#c44e52")
ax1.bar([i + w/2 for i in x], bv, w, label=f"background (n={S['background']['n']})", color="#4c72b0")
for i, (m, lab) in enumerate(metrics):
    p = T[m]["p"]
    ax1.text(i, max(gv[i], bv[i]) + 0.012, f"p={p:.2f}", ha="center", fontsize=8)
ax1.set_xticks(list(x)); ax1.set_xticklabels([l for _, l in metrics])
ax1.set_ylabel("fraction of genes"); ax1.set_ylim(0, 0.28)
ax1.set_title("Gated vs random surfaceome: DrugCentral engagement (null)")
ax1.legend(fontsize=8, frameon=False)

genes = list(AG); vals = [AG[g]["n_drugs"] for g in genes]
colors = ["#55a868" if v > 0 else "#c44e52" for v in vals]
ax2.bar(genes, vals, color=colors)
for i, v in enumerate(vals):
    ax2.text(i, v + 0.8, str(v), ha="center", fontsize=8)
ax2.set_ylabel("unique DrugCentral drugs"); ax2.set_title("AND-gate antigens: quantitative engagement")
ax2.tick_params(axis="x", rotation=35)
fig.tight_layout(); fig.savefig("paper/figs/fig_drugcentral.pdf")
print("wrote paper/figs/fig_drugcentral.pdf")

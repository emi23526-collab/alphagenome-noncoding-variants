"""Analysis of AlphaGenome scores for non-coding BRCA1 / INS ClinVar variants. No API calls."""
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve

P = "data/processed/"
FIG = "results/figures/"
os.makedirs(FIG, exist_ok=True)
pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 30)
SCORE = "score_absmax"  # primary score

# ---------- load and merge ----------
meta = pd.read_csv(P + "clinvar_noncoding_hg38.csv")
cols = ["key", "status", "score_mean", "score_absmax"]
sc = pd.concat([pd.read_csv(P + "alphagenome_scores_hg38.csv")[cols],
                pd.read_csv(P + "alphagenome_scores_noncoding.csv")[cols]])
sc = sc[sc["status"] == "ok"].drop_duplicates("key")
df = meta.merge(sc, on="key", how="inner")
df = df[df["group"].isin(["P/LP", "B/LB", "VUS", "Conflicting"])].copy()
df[SCORE] = df[SCORE].astype(float)


def dist(name):
    m = re.search(r"c\.[-*]?\d+([+-])(\d+)", str(name))
    return int(m.group(2)) if m else np.nan


df["dist"] = df["clinvar_name"].apply(dist)  # distance from exon edge (intronic variants)
df.to_csv("results/scored_noncoding_merged.csv", index=False)

print("scored non-coding variants used:", len(df))
print("\n=== gene x group ===")
print(pd.crosstab(df["gene"], df["group"]).to_string())
print("\n=== class x group (BRCA1) ===")
print(pd.crosstab(df[df.gene == "BRCA1"]["reg_class"], df[df.gene == "BRCA1"]["group"]).to_string())
print("\n=== median", SCORE, "by gene and group ===")
print(df.groupby(["gene", "group"])[SCORE].median().unstack().to_string())


# ---------- AUROC / AUPRC ----------
def auc_ci(y, s, n=1000, seed=0):
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(set(y[i])) == 2:
            vals.append(roc_auc_score(y[i], s[i]))
    return np.percentile(vals, [2.5, 97.5])


def evaluate(name, d, score=SCORE):
    d = d[d["group"].isin(["P/LP", "B/LB"])]
    y = (d["group"] == "P/LP").astype(int).values
    n1, n0 = int(y.sum()), int((1 - y).sum())
    if n1 < 5 or n0 < 5:
        print(f"  {name}: skipped (P/LP={n1}, B/LB={n0})")
        return None
    s = d[score].values
    lo, hi = auc_ci(y, s)
    return dict(subset=name, n_PLP=n1, n_BLB=n0, auroc=roc_auc_score(y, s), ci_lo=lo, ci_hi=hi,
                auprc=average_precision_score(y, s), auprc_baseline=n1 / len(y),
                mannwhitney_p=mannwhitneyu(s[y == 1], s[y == 0]).pvalue)


brca = df[df.gene == "BRCA1"]
ins = df[df.gene == "INS"]
intr = brca[brca["dist"].notna()].copy()
intr["neg_dist"] = -intr["dist"]
splice = brca[brca["reg_class"] == "Splice region (<=50bp)"]
nonsplice = brca[brca["reg_class"] != "Splice region (<=50bp)"]

rows = [
    evaluate("BRCA1 all non-coding", brca),
    evaluate("BRCA1 intronic only (model)", intr),
    evaluate("BRCA1 intronic only (baseline: closeness to exon)", intr, "neg_dist"),
    evaluate("BRCA1 splice-region only", splice),
    evaluate("BRCA1 not splice-region (deep intronic/UTR/promoter)", nonsplice),
    evaluate("INS non-coding", ins),
]
res = pd.DataFrame([r for r in rows if r])
res.to_csv("results/performance.csv", index=False)
print("\n=== performance (P/LP vs B/LB) ===")
print(res.round(3).to_string(index=False))

# ---------- VUS candidates above the benign 95th percentile ----------
cand = []
for g, d in df.groupby("gene"):
    ben = d[d["group"] == "B/LB"][SCORE]
    if len(ben) < 5:
        continue
    thr = ben.quantile(0.95)
    c = d[d["group"].isin(["VUS", "Conflicting"]) & (d[SCORE] > thr)].copy()
    c["threshold_benign_p95"] = thr
    cand.append(c)
    print(f"\n{g}: benign 95th percentile = {thr:.4f}; VUS/Conflicting above it: {len(c)}")
if cand:
    cand = pd.concat(cand).sort_values(SCORE, ascending=False)
    cand[["gene", "clinvar_name", "clinvar_sig", "reg_class", "dist", SCORE,
          "threshold_benign_p95"]].to_csv("results/top_candidates.csv", index=False)
    print(cand[["gene", "clinvar_name", "group", "reg_class", SCORE]].head(15).to_string(index=False))

# ---------- figures ----------
ORDER = ["B/LB", "VUS", "Conflicting", "P/LP"]
COL = {"B/LB": "#4c9f70", "VUS": "#e0a030", "Conflicting": "#8a8a8a", "P/LP": "#c44e52"}
rng = np.random.default_rng(1)


def box(ax, d, title):
    data = [np.clip(d[d.group == g][SCORE].values, 1e-6, None) for g in ORDER]
    ax.boxplot(data, positions=range(len(ORDER)), widths=0.5, showfliers=False)
    for i, g in enumerate(ORDER):
        v = data[i]
        ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=6, alpha=0.4, color=COL[g])
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([f"{g}\n(n={len(x)})" for g, x in zip(ORDER, data)])
    ax.set_yscale("log")
    ax.set_ylabel("score_absmax (log)")
    ax.set_title(title)


fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
box(axs[0], brca, "BRCA1 (breast)")
box(axs[1], ins, "INS (beta cell)")
plt.tight_layout()
plt.savefig(FIG + "fig1_scores_by_clinvar_group.png", dpi=150)
plt.close()

fig, ax = plt.subplots(figsize=(10, 4.5))
classes = [c for c in brca["reg_class"].unique()]
for k, c in enumerate(classes):
    for j, g in enumerate(["B/LB", "P/LP"]):
        v = np.clip(brca[(brca.reg_class == c) & (brca.group == g)][SCORE].values, 1e-6, None)
        if len(v) == 0:
            continue
        x = k + (j - 0.5) * 0.35
        ax.boxplot([v], positions=[x], widths=0.3, showfliers=False)
        ax.scatter(x + rng.uniform(-0.08, 0.08, len(v)), v, s=6, alpha=0.4, color=COL[g])
        ax.text(x, ax.get_ylim()[0], f"n={len(v)}", ha="center", va="bottom", fontsize=7)
ax.set_xticks(range(len(classes)))
ax.set_xticklabels(classes, rotation=15)
ax.set_yscale("log")
ax.set_ylabel("score_absmax (log)")
ax.set_title("BRCA1 by class: B/LB (green, left) vs P/LP (red, right)")
plt.tight_layout()
plt.savefig(FIG + "fig2_scores_by_class.png", dpi=150)
plt.close()

fig, ax = plt.subplots(figsize=(5.5, 5))
for name, d, s in [("AlphaGenome score", intr, SCORE), ("Closeness to exon only", intr, "neg_dist")]:
    dd = d[d.group.isin(["P/LP", "B/LB"])]
    y = (dd.group == "P/LP").astype(int)
    if y.sum() >= 5 and (1 - y).sum() >= 5:
        fpr, tpr, _ = roc_curve(y, dd[s])
        ax.plot(fpr, tpr, label=f"{name} (AUROC {roc_auc_score(y, dd[s]):.2f})")
ax.plot([0, 1], [0, 1], "k--", lw=0.8)
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title("BRCA1 intronic: P/LP vs B/LB")
ax.legend(loc="lower right")
plt.tight_layout()
plt.savefig(FIG + "fig3_roc_brca1.png", dpi=150)
plt.close()

print("\nSaved figures to", FIG)
print("Saved results/performance.csv, results/top_candidates.csv, results/scored_noncoding_merged.csv")
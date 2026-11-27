"""Compare RNA-seq score, splice score and distance. No API calls."""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 30)

base = pd.read_csv("results/scored_noncoding_merged.csv")
base = base.rename(columns={"score_absmax": "rna"})
sp = pd.read_csv("data/processed/alphagenome_scores_splice.csv")
sp = sp[sp["status"] == "ok"].drop_duplicates("key")[["key", "score_absmax"]]
sp = sp.rename(columns={"score_absmax": "splice"})
df = base.merge(sp, on="key", how="inner")
df["rna"] = df["rna"].astype(float)
df["splice"] = df["splice"].astype(float)
df["neg_dist"] = -df["dist"]
print("splice-region variants with both scores:", len(df))
print(df.groupby(["gene", "group"]).size().unstack(fill_value=0).to_string())

d = df[(df["gene"] == "BRCA1") & df["dist"].notna() & df["group"].isin(["P/LP", "B/LB"])].copy()
d["y"] = (d["group"] == "P/LP").astype(int)


def auc_ci(y, s, n=1000, seed=0):
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(set(y[i])) == 2:
            vals.append(roc_auc_score(y[i], s[i]))
    return np.percentile(vals, [2.5, 97.5])


subsets = [
    ("3-10 bp", d[(d["dist"] >= 3) & (d["dist"] <= 10)]),
    ("11-50 bp", d[(d["dist"] >= 11) & (d["dist"] <= 50)]),
    ("3-50 bp", d[(d["dist"] >= 3) & (d["dist"] <= 50)]),
    ("all intronic splice-region", d),
]
predictors = [("distance", "neg_dist"), ("RNA-seq score", "rna"), ("splice score", "splice")]

rows = []
print("\n=== AUROC with bootstrap 95% CI (P/LP vs B/LB, BRCA1) ===")
for name, sub in subsets:
    y = sub["y"].values
    n1, n0 = int(y.sum()), int((1 - y).sum())
    if n1 < 5 or n0 < 5:
        print(f"  {name}: skipped (P/LP={n1}, B/LB={n0})")
        continue
    for pname, col in predictors:
        s = sub[col].values
        lo, hi = auc_ci(y, s)
        rows.append({"subset": name, "n_PLP": n1, "n_BLB": n0, "predictor": pname,
                     "auroc": roc_auc_score(y, s), "ci_lo": lo, "ci_hi": hi,
                     "mannwhitney_p": mannwhitneyu(s[y == 1], s[y == 0]).pvalue})
res = pd.DataFrame(rows)
if len(res):
    print(res.round(3).to_string(index=False))
    res.to_csv("results/splice_comparison.csv", index=False)


def cv_auc(X, y, repeats=10):
    aucs = []
    for seed in range(repeats):
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
        model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))
        p = cross_val_predict(model, X, y, cv=cv, method="predict_proba")[:, 1]
        aucs.append(roc_auc_score(y, p))
    return float(np.mean(aucs))


print("\n=== combined models, 5-fold cross-validation x 10 repeats ===")
comb = []
for name, sub in subsets:
    y = sub["y"].values
    if y.sum() < 10 or (1 - y).sum() < 10:
        print(f"  {name}: skipped")
        continue
    dist = np.log1p(sub["dist"].values)
    rna = np.log10(sub["rna"].values + 1e-6)
    spl = np.log10(sub["splice"].values + 1e-6)
    feats = {
        "distance": np.column_stack([dist]),
        "distance + RNA-seq": np.column_stack([dist, rna]),
        "distance + splice": np.column_stack([dist, spl]),
        "distance + RNA-seq + splice": np.column_stack([dist, rna, spl]),
    }
    for k, X in feats.items():
        comb.append({"subset": name, "model": k, "cv_auroc": cv_auc(X, y)})
if comb:
    comb = pd.DataFrame(comb)
    print(comb.round(3).to_string(index=False))
    comb.to_csv("results/splice_combined.csv", index=False)

print("\n=== VUS / Conflicting above the benign 95th percentile of the splice score ===")
cand = []
for g, dd in df.groupby("gene"):
    ben = dd[dd["group"] == "B/LB"]["splice"]
    if len(ben) < 5:
        print(f"  {g}: too few benign variants, skipped")
        continue
    thr = ben.quantile(0.95)
    c = dd[dd["group"].isin(["VUS", "Conflicting"]) & (dd["splice"] > thr)].copy()
    c["threshold_benign_p95"] = thr
    cand.append(c)
    print(f"  {g}: threshold {thr:.4f}; candidates above it: {len(c)}")
if cand:
    cand = pd.concat(cand).sort_values("splice", ascending=False)
    cand[["gene", "clinvar_name", "group", "dist", "rna", "splice",
          "threshold_benign_p95"]].to_csv("results/splice_candidates.csv", index=False)
    print(cand[["gene", "clinvar_name", "group", "dist", "rna", "splice"]].head(15).to_string(index=False))

print("\nSaved results/splice_comparison.csv, splice_combined.csv, splice_candidates.csv")
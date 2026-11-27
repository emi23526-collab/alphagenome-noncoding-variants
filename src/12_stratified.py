"""Stratified test: does AlphaGenome add anything beyond distance to the exon? No API calls."""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

SCORE = "score_absmax"
pd.set_option("display.width", 200)

df = pd.read_csv("results/scored_noncoding_merged.csv")
d = df[(df["gene"] == "BRCA1") & df["dist"].notna() & df["group"].isin(["P/LP", "B/LB"])].copy()
d["y"] = (d["group"] == "P/LP").astype(int)
d["neg_dist"] = -d["dist"]
print("BRCA1 intronic variants (P/LP + B/LB):", len(d),
      "| P/LP:", int(d["y"].sum()), "| B/LB:", int((1 - d["y"]).sum()))

# ---------- 1. Within distance bins ----------
bins = [(1, 2, "canonical (1-2 bp)"), (3, 10, "3-10 bp"), (11, 50, "11-50 bp"),
        (51, 10**9, "deep intronic (>50 bp)")]
rows = []
print("\n=== within distance bins (P/LP vs B/LB) ===")
for lo, hi, name in bins:
    b = d[(d["dist"] >= lo) & (d["dist"] <= hi)]
    n1, n0 = int(b["y"].sum()), int((1 - b["y"]).sum())
    rec = {"bin": name, "n_PLP": n1, "n_BLB": n0, "auroc_model": np.nan, "auroc_distance": np.nan}
    if n1 >= 5 and n0 >= 5:
        rec["auroc_model"] = roc_auc_score(b["y"], b[SCORE])
        rec["auroc_distance"] = roc_auc_score(b["y"], b["neg_dist"])
    rows.append(rec)
bin_res = pd.DataFrame(rows)
print(bin_res.round(3).to_string(index=False))
print("(empty AUROC = fewer than 5 P/LP or 5 B/LB in that bin, so it is not tested)")
bin_res.to_csv("results/stratified_bins.csv", index=False)


# ---------- 2. Combined model (cross-validated) ----------
def cv_auc(X, y, repeats=10):
    aucs = []
    for seed in range(repeats):
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
        model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))
        p = cross_val_predict(model, X, y, cv=cv, method="predict_proba")[:, 1]
        aucs.append(roc_auc_score(y, p))
    return float(np.mean(aucs)), float(np.min(aucs)), float(np.max(aucs))


def compare(name, sub):
    y = sub["y"].values
    if y.sum() < 10 or (1 - y).sum() < 10:
        print(f"  {name}: skipped (P/LP={int(y.sum())}, B/LB={int((1 - y).sum())})")
        return []
    feats = {
        "distance only": np.log1p(sub[["dist"]].values),
        "AlphaGenome only": np.log10(sub[[SCORE]].values + 1e-6),
        "distance + AlphaGenome": np.column_stack([np.log1p(sub["dist"].values),
                                                   np.log10(sub[SCORE].values + 1e-6)]),
    }
    out = []
    for k, X in feats.items():
        m, lo, hi = cv_auc(X, y)
        out.append({"subset": name, "model": k, "n_PLP": int(y.sum()),
                    "n_BLB": int((1 - y).sum()), "cv_auroc_mean": m,
                    "cv_auroc_min": lo, "cv_auroc_max": hi})
    return out


print("\n=== combined model, 5-fold cross-validation x 10 repeats ===")
res = compare("all intronic", d) + compare("distance 3-50 bp only", d[(d["dist"] >= 3) & (d["dist"] <= 50)])
if res:
    comb = pd.DataFrame(res)
    print(comb.round(3).to_string(index=False))
    comb.to_csv("results/stratified_combined.csv", index=False)

print("\nSaved results/stratified_bins.csv and results/stratified_combined.csv")
"""Overlap scored variants with ENCODE candidate cis-regulatory elements (cCREs). No API calls.

Needs the file data/external/GRCh38-cCREs.bed (see the download step).
"""
import os

import pandas as pd

BED = "data/external/GRCh38-cCREs.bed"
SCORE = "score_absmax"
if not os.path.exists(BED):
    raise SystemExit(f"{BED} not found - do the download step first.")

# ---------- load cCREs (BED: 0-based start, end exclusive; class is the last column) ----------
raw = pd.read_csv(BED, sep="\t", header=None, comment="#")
cc = raw.iloc[:, [0, 1, 2, raw.shape[1] - 1]].copy()
cc.columns = ["chrom", "start", "end", "ccre_class"]
cc["start"] = cc["start"].astype(int)
cc["end"] = cc["end"].astype(int)
print("cCREs loaded:", len(cc))

# ---------- variants ----------
df = pd.read_csv("results/scored_noncoding_merged.csv")
df["pos"] = df["key"].str.split(":").str[1].astype(int)  # key looks like chr17:43067696:C:T:BRCA1
df["chrom"] = df["key"].str.split(":").str[0]


def classes_at(chrom, pos):
    sub = cc[(cc["chrom"] == chrom) & (cc["start"] < pos) & (cc["end"] >= pos)]
    return ",".join(sorted(set(sub["ccre_class"]))) if len(sub) else ""


# only look at cCREs near the genes, to keep it fast
parts = []
for chrom, g in df.groupby("chrom"):
    lo, hi = g["pos"].min() - 5000, g["pos"].max() + 5000
    parts.append(cc[(cc["chrom"] == chrom) & (cc["end"] >= lo) & (cc["start"] <= hi)])
cc = pd.concat(parts)
print("cCREs near the scored variants:", len(cc))

df["ccre_class"] = [classes_at(c, p) for c, p in zip(df["chrom"], df["pos"])]
df["in_ccre"] = df["ccre_class"] != ""

print("\n=== share of variants overlapping a cCRE, by group ===")
print(df.groupby(["gene", "group"])["in_ccre"].agg(["mean", "sum", "count"]).round(3).to_string())

# ---------- candidates above the benign 95th percentile ----------
out = []
for g, d in df.groupby("gene"):
    ben = d[d["group"] == "B/LB"][SCORE].astype(float)
    if len(ben) < 5:
        continue
    thr = ben.quantile(0.95)
    c = d[d["group"].isin(["VUS", "Conflicting"]) & (d[SCORE].astype(float) > thr)].copy()
    c["threshold_benign_p95"] = thr
    out.append(c)
if out:
    out = pd.concat(out).sort_values(SCORE, ascending=False)
    cols = ["gene", "clinvar_name", "group", "reg_class", SCORE, "in_ccre", "ccre_class"]
    out[cols].to_csv("results/candidates_with_ccre.csv", index=False)
    print("\n=== candidates NOT in the splice region (UTR / promoter / deep intronic) ===")
    ns = out[out["reg_class"] != "Splice region (<=50bp)"]
    print(ns[cols].head(25).to_string(index=False))
    print("\nSaved results/candidates_with_ccre.csv")
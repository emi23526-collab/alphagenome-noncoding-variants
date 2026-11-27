"""Real-data check: what is in the 1,500 scored variants?"""
import re
import pandas as pd

pd.set_option("display.width", 200)
pd.set_option("display.max_colwidth", 80)

sc = pd.read_csv("data/processed/alphagenome_scores_hg38.csv")
meta = pd.read_csv("data/processed/clinvar_hg38.csv")[["key", "clinvar_name"]].drop_duplicates("key")
df = sc.merge(meta, on="key", how="left")
df["clinvar_name"] = df["clinvar_name"].fillna("")
df = df[df["status"] == "ok"].copy()
print("scored OK:", len(df), "\n")


def group(sig):
    s = str(sig).lower()
    if "conflicting" in s: return "Conflicting"
    if "uncertain" in s: return "VUS"
    if "pathogenic" in s and "benign" not in s: return "P/LP"
    if "benign" in s and "pathogenic" not in s: return "B/LB"
    return "Other"


def reg_class(name):
    n = str(name)
    if re.search(r"\(p\.", n): return "Coding (protein change)"
    m = re.search(r"c\.[-*]?\d+([+-])(\d+)", n)
    if m: return "Splice region (<=50bp)" if int(m.group(2)) <= 50 else "Deep intronic (>50bp)"
    if re.search(r"c\.-\d+", n): return "5' UTR / promoter"
    if re.search(r"c\.\*\d+", n): return "3' UTR"
    return "Unclassified"


df["group"] = df["clinvar_sig"].apply(group)
df["reg_class"] = df["clinvar_name"].apply(reg_class)
df["noncoding"] = df["reg_class"] != "Coding (protein change)"

print("=== raw clinvar_sig values ===")
print(df.groupby(["gene", "clinvar_sig"]).size().to_string(), "\n")
print("=== cohort ===")
print(df.groupby(["gene", "cohort"]).size().to_string(), "\n")
print("=== group by gene ===")
print(pd.crosstab(df["gene"], df["group"]).to_string(), "\n")
print("=== regulatory class by gene ===")
print(pd.crosstab(df["gene"], df["reg_class"]).to_string(), "\n")
print("=== NON-CODING only: group by gene ===")
nc = df[df["noncoding"]]
print(pd.crosstab(nc["gene"], nc["group"]).to_string(), "\n")

for gene in ["BRCA1", "INS", "pooled"]:
    s = nc if gene == "pooled" else nc[nc["gene"] == gene]
    n = int((s["group"] == "P/LP").sum())
    branch = "A" if n >= 20 else ("B" if n >= 10 else "C")
    print(f"{gene}: non-coding P/LP = {n} -> branch {branch}")

print("\n=== 5 example variant names ===")
print(df[["gene", "clinvar_sig", "clinvar_name", "reg_class"]].sample(5, random_state=1).to_string(index=False))
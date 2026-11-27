"""Pull NON-CODING BRCA1 / INS variants from ClinVar's bulk table (GRCh38 only). No API calls."""
import os
import re
import pandas as pd

RAW = "data/raw/variant_summary.txt.gz"
OUT = "data/processed/clinvar_noncoding_hg38.csv"
SCORED = "data/processed/alphagenome_scores_hg38.csv"
GENES = ("BRCA1", "INS")
MAX_ALLELE_LEN = 50

if not os.path.exists(RAW):
    raise SystemExit(f"Missing {RAW} - run the download command first")

WANT = ["Type", "Name", "GeneSymbol", "ClinicalSignificance", "ReviewStatus", "Assembly",
        "Chromosome", "VariationID", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"]

parts, total = [], 0
for chunk in pd.read_csv(RAW, sep="\t", compression="gzip", usecols=lambda c: c in WANT,
                         dtype=str, chunksize=500_000):
    missing = [c for c in WANT if c not in chunk.columns]
    if missing:
        raise SystemExit(f"Expected columns missing: {missing}\nColumns found: {list(chunk.columns)}")
    total += len(chunk)
    keep = chunk["GeneSymbol"].str.contains("BRCA1|INS", na=False) & (chunk["Assembly"] == "GRCh38")
    parts.append(chunk[keep])
print("rows read from ClinVar:", f"{total:,}")

df = pd.concat(parts, ignore_index=True)


def pick_gene(s):
    toks = {x.strip() for x in re.split(r"[;|,]", str(s))}
    for g in GENES:
        if g in toks:
            return g
    return None


df["gene"] = df["GeneSymbol"].apply(pick_gene)
df = df[df["gene"].notna()].copy()
print("GRCh38 rows for BRCA1 / INS:", len(df))

df["ref"] = df["ReferenceAlleleVCF"].fillna("").str.upper()
df["alt"] = df["AlternateAlleleVCF"].fillna("").str.upper()
ok = (df["ref"].str.fullmatch(r"[ACGT]+") & df["alt"].str.fullmatch(r"[ACGT]+") &
      (df["ref"].str.len() <= MAX_ALLELE_LEN) & (df["alt"].str.len() <= MAX_ALLELE_LEN) &
      df["PositionVCF"].str.fullmatch(r"\d+", na=False))
print("dropped (no usable VCF alleles or >50 bp):", int((~ok).sum()))
df = df[ok].copy()

df["pos"] = df["PositionVCF"].astype(int)
df["chrom"] = "chr" + df["Chromosome"].astype(str)
df["clinvar_name"] = df["Name"].fillna("")
df["clinvar_sig"] = df["ClinicalSignificance"].fillna("")


def group(sig):
    s = sig.lower()
    if "conflicting" in s: return "Conflicting"
    if "uncertain" in s: return "VUS"
    if "pathogenic" in s and "benign" not in s: return "P/LP"
    if "benign" in s and "pathogenic" not in s: return "B/LB"
    return "Other"


def reg_class(name):
    if re.search(r"\(p\.", name): return "Coding (protein change)"
    m = re.search(r"c\.[-*]?\d+([+-])(\d+)", name)
    if m: return "Splice region (<=50bp)" if int(m.group(2)) <= 50 else "Deep intronic (>50bp)"
    if re.search(r"c\.-\d+", name): return "5' UTR / promoter"
    if re.search(r"c\.\*\d+", name): return "3' UTR"
    return "Unclassified"


df["group"] = df["clinvar_sig"].apply(group)
df["reg_class"] = df["clinvar_name"].apply(reg_class)
df["key"] = (df["chrom"] + ":" + df["pos"].astype(str) + ":" + df["ref"] + ":" +
             df["alt"] + ":" + df["gene"])
df = df.drop_duplicates("key")

print("\nexcluded as coding:", int((df["reg_class"] == "Coding (protein change)").sum()))
print("excluded as unclassified:", int((df["reg_class"] == "Unclassified").sum()))
nc = df[~df["reg_class"].isin(["Coding (protein change)", "Unclassified"])].copy()

done = set()
if os.path.exists(SCORED):
    s = pd.read_csv(SCORED)
    done = set(s.loc[s["status"] == "ok", "key"])
nc["already_scored"] = nc["key"].isin(done)

nc = nc.rename(columns={"ReviewStatus": "review_status", "VariationID": "variation_id"})
cols = ["key", "chrom", "pos", "ref", "alt", "gene", "clinvar_sig", "clinvar_name", "reg_class",
        "group", "review_status", "variation_id", "already_scored"]
os.makedirs("data/processed", exist_ok=True)
nc[cols].to_csv(OUT, index=False)

print(f"\nNON-CODING variants saved: {len(nc)} -> {OUT}")
print("\n=== group by gene ===")
print(pd.crosstab(nc["gene"], nc["group"]).to_string())
print("\n=== regulatory class by gene ===")
print(pd.crosstab(nc["gene"], nc["reg_class"]).to_string())
print("\n=== class x group (all genes) ===")
print(pd.crosstab(nc["reg_class"], nc["group"]).to_string())

todo = int((~nc["already_scored"]).sum())
print(f"\nalready scored: {int(nc['already_scored'].sum())} | still to score: {todo}")
print(f"estimated scoring time: about {todo / 0.84 / 60:.0f} minutes")
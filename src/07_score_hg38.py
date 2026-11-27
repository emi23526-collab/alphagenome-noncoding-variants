"""Convert mixed hg19/hg38 variants to hg38, then score them with AlphaGenome.

Usage:
    python -u src\\07_score_hg38.py dry   # convert only, no API calls
    python -u src\\07_score_hg38.py 3     # score 3 variants
    python -u src\\07_score_hg38.py       # full run (safe to stop and restart)
"""
import csv
import os
import sys
import time

import pandas as pd
from pyliftover import LiftOver

IN_CSV = "data/processed/clinvar_rarity.csv"
CONV_CSV = "data/processed/clinvar_hg38.csv"
OUT_CSV = "data/processed/alphagenome_scores_hg38.csv"
HALF_WINDOW = 524_288
MAX_CONSECUTIVE_ERRORS = 10

TISSUE = {
    "BRCA1": "UBERON:0000310",  # breast
    "INS": "CL:0000169",        # pancreatic beta cell
}


def build_of(gene, pos):
    """Guess the genome build from the position (genes sit at different places in each build)."""
    if gene == "BRCA1":   # hg19 ~41.2M, hg38 ~43.1M
        return "hg19" if pos < 42_100_000 else "hg38"
    if gene == "INS":     # hg19 ~2.182M, hg38 ~2.160M
        return "hg19" if pos > 2_170_000 else "hg38"
    return "unknown"


# ---------- 1. Convert everything to hg38 ----------
lo = LiftOver("data/hg19ToHg38.over.chain.gz")
df = pd.read_csv(IN_CSV)
df = df[df["gene"].isin(TISSUE)].copy()
df["build"] = [build_of(g, int(p)) for g, p in zip(df["gene"], df["pos"])]
df["orig_pos"] = df["pos"].astype(int)

new_pos = []
for b, c, p in zip(df["build"], df["chrom"], df["orig_pos"]):
    if b == "hg38":
        new_pos.append(p)
    else:
        res = lo.convert_coordinate(c, int(p) - 1)  # liftover is 0-based
        new_pos.append(res[0][1] + 1 if res else None)
df["pos"] = new_pos

print("Build counts per gene:")
print(df.groupby(["gene", "build"]).size().to_string())
failed = int(df["pos"].isna().sum())
print("liftover failures (dropped):", failed)
df = df.dropna(subset=["pos"]).copy()
df["pos"] = df["pos"].astype(int)

df["key"] = (df["chrom"].astype(str) + ":" + df["pos"].astype(str) + ":" +
             df["ref"].astype(str) + ":" + df["alt"].astype(str) + ":" + df["gene"].astype(str))
before = len(df)
df = df.drop_duplicates("key")
print(f"duplicates removed after conversion: {before - len(df)} | unique variants: {len(df)}")

lifted = df[df["build"] == "hg19"].head(5)
print("\nExample conversions (hg19 -> hg38):")
print(lifted[["gene", "orig_pos", "pos"]].to_string(index=False))
print("\nPosition range per gene after conversion:")
print(df.groupby("gene")["pos"].agg(["min", "max"]).to_string())

df.to_csv(CONV_CSV, index=False)
print("\nSaved", CONV_CSV)

mode = sys.argv[1] if len(sys.argv) > 1 else None
if mode == "dry":
    raise SystemExit("Dry run finished - no API calls made.")

# ---------- 2. Score ----------
from alphagenome.data.genome import Interval, Variant
from alphagenome.models.v2 import dna_client, variant_scorers
from alphagenome.models.v2.dna_model import Organism

FIELDS = ["key", "chrom", "pos", "orig_pos", "build", "ref", "alt", "gene", "clinvar_sig",
          "gnomad_af", "cohort", "condition", "status", "score_mean", "score_absmax", "error"]

api_key = os.environ.get("ALPHA_GENOME_API_KEY", "")
if not api_key:
    raise SystemExit("API KEY MISSING - load the key first")
model = dna_client.create(api_key)
scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]

done = set()
if os.path.exists(OUT_CSV):
    prev = pd.read_csv(OUT_CSV)
    done = set(prev.loc[prev["status"] != "error", "key"])

todo = df[~df["key"].isin(done)]
limit = int(mode) if mode else None
if limit:
    todo = todo.head(limit)
print(f"\nvariants: {len(df)} | already done: {len(done)} | to do now: {len(todo)}")

new_file = not os.path.exists(OUT_CSV)
out = open(OUT_CSV, "a", newline="", encoding="utf-8")
writer = csv.DictWriter(out, fieldnames=FIELDS)
if new_file:
    writer.writeheader()

errors_in_row = 0
t0 = time.time()
for n, (_, r) in enumerate(todo.iterrows(), 1):
    cond = TISSUE[r["gene"]]
    pos = int(r["pos"])
    rec = {
        "key": r["key"], "chrom": r["chrom"], "pos": pos, "orig_pos": int(r["orig_pos"]),
        "build": r["build"], "ref": r["ref"], "alt": r["alt"], "gene": r["gene"],
        "clinvar_sig": r.get("clinvar_sig", ""), "gnomad_af": r.get("gnomad_af", ""),
        "cohort": r.get("cohort", ""), "condition": cond,
        "status": "", "score_mean": "", "score_absmax": "", "error": "",
    }
    try:
        variant = Variant(chromosome=r["chrom"], position=pos,
                          reference_bases=str(r["ref"]), alternate_bases=str(r["alt"]))
        interval = Interval(r["chrom"], max(0, pos - HALF_WINDOW), pos + HALF_WINDOW)
        res = model.score_variant(interval=interval, variant=variant,
                                  variant_scorers=[scorer], condition=cond,
                                  organism=Organism.HOMO_SAPIENS)
        ad = res[0]
        rows = [i for i, g in enumerate(ad.obs["gene_name"]) if g == r["gene"]]
        if not rows:
            rec["status"] = "gene_not_in_window"
        else:
            vals = ad.X[rows[0]]
            rec["score_mean"] = float(vals.mean())
            rec["score_absmax"] = float(abs(vals).max())
            rec["status"] = "ok"
        errors_in_row = 0
    except Exception as e:
        rec["status"] = "error"
        rec["error"] = f"{type(e).__name__}: {str(e)[:200]}".replace("\n", " ")
        errors_in_row += 1
        print(f"  error on {r['key']}: {rec['error']}")
        time.sleep(2)

    writer.writerow(rec)
    out.flush()

    if n % 10 == 0 or n == len(todo):
        rate = n / (time.time() - t0)
        print(f"[{n}/{len(todo)}] last={rec['status']}  {rate:.2f} variants/s")

    if errors_in_row >= MAX_CONSECUTIVE_ERRORS:
        print(f"\nStopping: {MAX_CONSECUTIVE_ERRORS} errors in a row (quota, key or network). "
              "Progress is saved. Fix the cause and run again.")
        break

out.close()
print("\nDone. Results in", OUT_CSV)

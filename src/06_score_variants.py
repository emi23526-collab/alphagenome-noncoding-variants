"""Score every variant in clinvar_rarity.csv with AlphaGenome (RNA-seq effect on its own gene).

Usage:
    python -u src\\06_score_variants.py 3      # test on the first 3 variants
    python -u src\\06_score_variants.py        # full run (safe to stop and restart)
"""
import csv
import os
import sys
import time

import pandas as pd
from alphagenome.data.genome import Interval, Variant
from alphagenome.models.v2 import dna_client, variant_scorers
from alphagenome.models.v2.dna_model import Organism

IN_CSV = "data/processed/clinvar_rarity.csv"
OUT_CSV = "data/processed/alphagenome_scores.csv"
HALF_WINDOW = 524_288  # interval width = 1,048,576 (largest supported)
MAX_CONSECUTIVE_ERRORS = 10

TISSUE = {
    "BRCA1": "UBERON:0000310",  # breast
    "INS": "CL:0000169",        # pancreatic beta cell
}

FIELDS = ["key", "chrom", "pos", "ref", "alt", "gene", "clinvar_sig", "gnomad_af",
          "cohort", "condition", "status", "score_mean", "score_absmax", "error"]

key_ = os.environ.get("ALPHA_GENOME_API_KEY", "")
if not key_:
    raise SystemExit("API KEY MISSING - load the key first")
model = dna_client.create(key_)
scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]

df = pd.read_csv(IN_CSV)
df["key"] = (df["chrom"].astype(str) + ":" + df["pos"].astype(str) + ":" +
             df["ref"].astype(str) + ":" + df["alt"].astype(str) + ":" + df["gene"].astype(str))
df = df.drop_duplicates("key")
df = df[df["gene"].isin(TISSUE)]

done = set()
if os.path.exists(OUT_CSV):
    prev = pd.read_csv(OUT_CSV)
    done = set(prev.loc[prev["status"] != "error", "key"])

todo = df[~df["key"].isin(done)]
limit = int(sys.argv[1]) if len(sys.argv) > 1 else None
if limit:
    todo = todo.head(limit)
print(f"variants in file: {len(df)} | already done: {len(done)} | to do now: {len(todo)}")

new_file = not os.path.exists(OUT_CSV)
out = open(OUT_CSV, "a", newline="", encoding="utf-8")
writer = csv.DictWriter(out, fieldnames=FIELDS)
if new_file:
    writer.writeheader()

errors_in_row = 0
t0 = time.time()
for n, (_, r) in enumerate(todo.iterrows(), 1):
    cond = TISSUE[r["gene"]]
    rec = {
        "key": r["key"], "chrom": r["chrom"], "pos": int(r["pos"]), "ref": r["ref"],
        "alt": r["alt"], "gene": r["gene"], "clinvar_sig": r.get("clinvar_sig", ""),
        "gnomad_af": r.get("gnomad_af", ""), "cohort": r.get("cohort", ""),
        "condition": cond, "status": "", "score_mean": "", "score_absmax": "", "error": "",
    }
    try:
        pos = int(r["pos"])
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
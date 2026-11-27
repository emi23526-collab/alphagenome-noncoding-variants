"""Score the non-coding ClinVar variants not yet scored. Safe to stop and restart.

Usage:
    python -u src\\10_score_noncoding.py 3    # test on 3 variants
    python -u src\\10_score_noncoding.py      # full run
"""
import csv
import os
import sys
import time

import pandas as pd
from alphagenome.data.genome import Interval, Variant
from alphagenome.models.v2 import dna_client, variant_scorers
from alphagenome.models.v2.dna_model import Organism

IN_CSV = "data/processed/clinvar_noncoding_hg38.csv"
OUT_CSV = "data/processed/alphagenome_scores_noncoding.csv"
HALF_WINDOW = 524_288
MAX_CONSECUTIVE_ERRORS = 10
KEEP_GROUPS = {"P/LP", "B/LB", "VUS", "Conflicting"}

TISSUE = {
    "BRCA1": "UBERON:0000310",  # breast
    "INS": "CL:0000169",        # pancreatic beta cell
}

FIELDS = ["key", "chrom", "pos", "ref", "alt", "gene", "clinvar_sig", "clinvar_name",
          "reg_class", "group", "condition", "status", "score_mean", "score_absmax", "error"]

api_key = os.environ.get("ALPHA_GENOME_API_KEY", "")
if not api_key:
    raise SystemExit("API KEY MISSING - load the key first")
model = dna_client.create(api_key)
scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]

df = pd.read_csv(IN_CSV)
df = df[(~df["already_scored"]) & (df["ref"] != df["alt"]) &
        df["group"].isin(KEEP_GROUPS) & df["gene"].isin(TISSUE)]
df = df.sample(frac=1, random_state=0)  # shuffle so a partial run still mixes all classes

done = set()
if os.path.exists(OUT_CSV):
    prev = pd.read_csv(OUT_CSV)
    done = set(prev.loc[prev["status"] != "error", "key"])

todo = df[~df["key"].isin(done)]
limit = int(sys.argv[1]) if len(sys.argv) > 1 else None
if limit:
    todo = todo.head(limit)
print(f"eligible: {len(df)} | already done: {len(done)} | to do now: {len(todo)}")
print(f"estimated time: about {len(todo) / 0.84 / 60:.0f} minutes")

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
        "key": r["key"], "chrom": r["chrom"], "pos": pos, "ref": r["ref"], "alt": r["alt"],
        "gene": r["gene"], "clinvar_sig": r["clinvar_sig"], "clinvar_name": r["clinvar_name"],
        "reg_class": r["reg_class"], "group": r["group"], "condition": cond,
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

    if n % 25 == 0 or n == len(todo):
        rate = n / (time.time() - t0)
        print(f"[{n}/{len(todo)}] last={rec['status']}  {rate:.2f} variants/s")

    if errors_in_row >= MAX_CONSECUTIVE_ERRORS:
        print(f"\nStopping: {MAX_CONSECUTIVE_ERRORS} errors in a row (quota, key or network). "
              "Progress is saved. Fix the cause and run again.")
        break

out.close()
print("\nDone. Results in", OUT_CSV)
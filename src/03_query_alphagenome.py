"""
Step 3: Query AlphaGenome for variant effect scores.
Input:  data/processed/clinvar_rarity.csv
Output: data/processed/alphagenome_scores.csv
Resumes from successful scores only.
"""

import os
import csv
import time

from alphagenome.data import genome
from alphagenome.models.v2 import dna_client
from alphagenome.models.v2 import variant_scorers

IN_CSV = "data/processed/clinvar_rarity.csv"
OUT_CSV = "data/processed/alphagenome_scores.csv"

CHECKPOINT_EVERY = 25
SLEEP_BETWEEN = 0.5

FIELDNAMES = [
    "chrom", "pos", "ref", "alt", "gene", "clinvar_sig",
    "cohort", "tissue", "score", "status"
]

TISSUES = {
    "BRCA1": ["UBERON:0001911", "UBERON:0000017"],
    "INS": ["UBERON:0001264"],
}

os.makedirs("data/processed", exist_ok=True)


def load_done():
    if not os.path.exists(OUT_CSV):
        return set(), []

    done = set()
    rows = []

    with open(OUT_CSV, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("status") == "ok" and r.get("score") not in ("", None):
                key = (r["chrom"], r["pos"], r["ref"], r["alt"], r["tissue"])
                done.add(key)
                rows.append(r)

    print(f"[OK] Resuming: {len(done):,} successful scores already done", flush=True)
    return done, rows


def save_rows(rows):
    tmp = OUT_CSV + ".tmp"

    with open(tmp, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    os.replace(tmp, OUT_CSV)


def score_variant(model, chrom, pos, ref, alt, tissue):
    try:
        variant = genome.Variant(
            chromosome=chrom,
            position=int(pos),
            reference_bases=ref,
            alternate_bases=alt,
        )

        interval = variant.reference_interval.resize(
            dna_client.SEQUENCE_LENGTH_1MB
        )

        scores = model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=[
                variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]
            ],
            condition=tissue,
        )

        vals = []

        for track_scores in scores.values():
            try:
                vals.extend(
                    abs(float(x))
                    for x in track_scores.flatten()
                    if x == x
                )
            except Exception:
                continue

        return max(vals) if vals else None

    except Exception as e:
        print(
            f"    [!] score failed for {chrom}:{pos} "
            f"{ref}>{alt} [{tissue}]: {e}",
            flush=True,
        )
        return None


def main():
    api_key = os.environ.get("ALPHA_GENOME_API_KEY")

    if not api_key:
        raise ValueError("ALPHA_GENOME_API_KEY not set in environment")

    model = dna_client.create(api_key)
    print("[OK] AlphaGenome client created", flush=True)

    print("[*] Loading candidates...", flush=True)

    candidates = []

    with open(IN_CSV, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["cohort"] != "rare_candidate":
                continue

            gene = r["gene"]

            if gene not in TISSUES:
                continue

            candidates.append({
                "chrom": r["chrom"],
                "pos": r["pos"],
                "ref": r["ref"],
                "alt": r["alt"],
                "gene": gene,
                "clinvar_sig": r["clinvar_sig"],
                "cohort": r["cohort"],
            })

    print(f"[OK] Candidates: {len(candidates):,}", flush=True)

    tasks = []

    for c in candidates:
        for tissue in TISSUES[c["gene"]]:
            tasks.append((c, tissue))

    print(f"[OK] Total variant x tissue tasks: {len(tasks):,}", flush=True)

    done_keys, done_rows = load_done()

    remaining = [
        (c, t)
        for c, t in tasks
        if (c["chrom"], c["pos"], c["ref"], c["alt"], t) not in done_keys
    ]

    print(f"[OK] Remaining tasks: {len(remaining):,}", flush=True)

    if not remaining:
        print("[OK] Nothing to do — already complete.", flush=True)
        return

    all_rows = list(done_rows)
    processed = 0

    for c, tissue in remaining:
        score = score_variant(
            model,
            c["chrom"],
            c["pos"],
            c["ref"],
            c["alt"],
            tissue,
        )

        all_rows.append({
            "chrom": c["chrom"],
            "pos": c["pos"],
            "ref": c["ref"],
            "alt": c["alt"],
            "gene": c["gene"],
            "clinvar_sig": c["clinvar_sig"],
            "cohort": c["cohort"],
            "tissue": tissue,
            "score": score if score is not None else "",
            "status": "ok" if score is not None else "failed",
        })

        processed += 1

        if processed % CHECKPOINT_EVERY == 0:
            save_rows(all_rows)
            print(
                f"    [checkpoint] saved {len(all_rows):,} rows",
                flush=True,
            )

        if processed % 10 == 0:
            print(
                f"    processed {processed:,}/{len(remaining):,}",
                flush=True,
            )

        time.sleep(SLEEP_BETWEEN)

    save_rows(all_rows)

    print(f"\n[OK] Saved: {OUT_CSV}", flush=True)

    ok = sum(1 for r in all_rows if r["status"] == "ok")
    failed = sum(1 for r in all_rows if r["status"] == "failed")

    print("\n--- Summary ---", flush=True)
    print(f"  ok: {ok:,}", flush=True)
    print(f"  failed: {failed:,}", flush=True)


if __name__ == "__main__":
    main()
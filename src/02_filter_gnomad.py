"""
Step 3: Query AlphaGenome for variant effect scores.
Input:  data/processed/clinvar_rarity.csv
Output: data/processed/alphagenome_scores.csv
Checkpoint every 100 variants. Resumes on restart.
"""

import os
import csv
import time

from alphagenome.data import genome
from alphagenome.models.v2 import dna_client

IN_CSV = "data/processed/clinvar_rarity.csv"
OUT_CSV = "data/processed/alphagenome_scores.csv"

CHECKPOINT_EVERY = 100
SLEEP_BETWEEN = 0.5

FIELDNAMES = [
    "chrom", "pos", "ref", "alt", "gene", "clinvar_sig",
    "cohort", "tissue", "score", "status"
]

# Tissue ontology per gene
TISSUES = {
    "BRCA1": ["UBERON:0001911", "UBERON:0000017"],  # mammary gland, ovary
    "INS":   ["UBERON:0001264"],                    # pancreatic islet
}

os.makedirs("data/processed", exist_ok=True)


def load_done():
    if not os.path.exists(OUT_CSV):
        return set(), []
    done = set()
    rows = []
    with open(OUT_CSV, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            key = (r["chrom"], r["pos"], r["ref"], r["alt"], r["tissue"])
            done.add(key)
            rows.append(r)
    print(f"[OK] Resuming: {len(done):,} scores already done")
    return done, rows


def save_rows(rows):
    tmp = OUT_CSV + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, OUT_CSV)


def score_variant(model, chrom, pos, ref, alt, tissue):
    """Query AlphaGenome for one variant in one tissue. Returns scalar or None."""
    try:
        variant = genome.Variant(
            chromosome=chrom,
            position=int(pos),
            reference_bases=ref,
            alternate_bases=alt,
        )
        interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

        scores = model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=[dna_client.RecommendedVariantScorers.RNA_SEQ],
            condition=tissue,
        )

        # Extract a single scalar: max absolute effect across tracks
        vals = []
        for track_scores in scores.values():
            try:
                vals.extend([abs(float(x)) for x in track_scores.flatten() if x == x])
            except Exception:
                continue
        return max(vals) if vals else None
    except Exception as e:
        print(f"    [!] score failed for {chrom}:{pos} {ref}>{alt} [{tissue}]: {e}")
        return None


def main():
    api_key = os.environ.get("ALPHA_GENOME_API_KEY")
    if not api_key:
        raise ValueError("ALPHA_GENOME_API_KEY not set in environment")

    model = dna_client.create(api_key)
    print("[OK] AlphaGenome client created")

    # Load candidates
    print("[*] Loading candidates...")
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

    print(f"[OK] Candidates: {len(candidates):,}")

    # Build task list: (variant, tissue)
    tasks = []
    for c in candidates:
        for tissue in TISSUES[c["gene"]]:
            tasks.append((c, tissue))
    print(f"[OK] Total (variant x tissue) tasks: {len(tasks):,}")

    # Resume
    done_keys, done_rows = load_done()
    remaining = [
        (c, t) for c, t in tasks
        if (c["chrom"], c["pos"], c["ref"], c["alt"], t) not in done_keys
    ]
    print(f"[OK] Remaining tasks: {len(remaining):,}")

    if not remaining:
        print("[OK] Nothing to do — already complete.")
        return

    all_rows = list(done_rows)
    processed = 0

    for c, tissue in remaining:
        score = score_variant(model, c["chrom"], c["pos"], c["ref"], c["alt"], tissue)
        all_rows.append({
            "chrom": c["chrom"], "pos": c["pos"], "ref": c["ref"], "alt": c["alt"],
            "gene": c["gene"], "clinvar_sig": c["clinvar_sig"], "cohort": c["cohort"],
            "tissue": tissue,
            "score": score if score is not None else "",
            "status": "ok" if score is not None else "failed",
        })
        processed += 1

        if processed % CHECKPOINT_EVERY == 0:
            save_rows(all_rows)
            print(f"    [checkpoint] saved {len(all_rows):,} total rows")

        if processed % 10 == 0:
            print(f"    processed {processed:,}/{len(remaining):,}")

        time.sleep(SLEEP_BETWEEN)

    save_rows(all_rows)
    print(f"\n[OK] Saved: {OUT_CSV}")

    ok = sum(1 for r in all_rows if r["status"] == "ok")
    failed = sum(1 for r in all_rows if r["status"] == "failed")
    print(f"\n--- Summary ---")
    print(f"  ok: {ok:,}")
    print(f"  failed: {failed:,}")


if __name__ == "__main__":
    main()
"""
Step 1: Pull non-coding BRCA1 and INS variants from ClinVar.
Pure-Python streaming (no pandas during load).
"""

import os
import csv
import gzip
import urllib.request

GENES = {"BRCA1", "INS"}
RAW_DIR = "data/raw"
OUT_CSV = os.path.join(RAW_DIR, "clinvar_extracted.csv")
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
LOCAL_GZ = os.path.join(RAW_DIR, "variant_summary.txt.gz")

CODING_TERMS = [
    "missense", "nonsense", "frameshift", "stop_gained", "stop_lost",
    "start_lost", "inframe_deletion", "inframe_insertion",
    "protein_altering", "coding_sequence_variant",
]

KEEP_COLS = [
    "AlleleID", "Type", "Name", "GeneSymbol", "ClinicalSignificance",
    "LastEvaluated", "RS# (dbSNP)", "Chromosome", "Start", "Stop",
    "ReferenceAlleleVCF", "AlternateAlleleVCF", "PositionVCF",
    "Assembly", "MolecularConsequence", "ReviewStatus", "PhenotypeList",
]


def download_clinvar():
    os.makedirs(RAW_DIR, exist_ok=True)
    if os.path.exists(LOCAL_GZ):
        print(f"[OK] ClinVar already downloaded: {LOCAL_GZ}")
        return
    print("[*] Downloading ClinVar...")
    urllib.request.urlretrieve(CLINVAR_URL, LOCAL_GZ)
    print(f"[OK] Downloaded: {LOCAL_GZ}")


def stream_filter():
    print("[*] Streaming ClinVar line-by-line (pure python)...")

    out_rows = []
    total = 0
    kept_genes = 0

    with gzip.open(LOCAL_GZ, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames
        missing = [c for c in KEEP_COLS if c not in fieldnames]
        if missing:
            print(f"[!] Missing columns in ClinVar: {missing}")

        for row in reader:
            total += 1
            if total % 500_000 == 0:
                print(f"    rows seen: {total:,}, kept: {kept_genes:,}")

            gene = row.get("GeneSymbol", "")
            if gene not in GENES:
                continue
            kept_genes += 1

            mc = (row.get("MolecularConsequence") or "").lower()
            if any(t in mc for t in CODING_TERMS):
                continue

            out_rows.append({k: row.get(k, "") for k in KEEP_COLS})

    print(f"[OK] Streamed. Total rows: {total:,}, BRCA1/INS rows: {kept_genes:,}")
    print(f"[OK] Non-coding BRCA1/INS rows: {len(out_rows):,}")

    os.makedirs(RAW_DIR, exist_ok=True)
    with open(OUT_CSV, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=KEEP_COLS)
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"[OK] Saved: {OUT_CSV}")

    from collections import Counter
    gene_counts = Counter(r["GeneSymbol"] for r in out_rows)
    sig_counts = Counter(r["ClinicalSignificance"] for r in out_rows)

    print("\n--- By gene ---")
    for g, c in gene_counts.most_common():
        print(f"  {g}: {c:,}")

    print("\n--- Top clinical significance ---")
    for s, c in sig_counts.most_common(10):
        print(f"  {s}: {c:,}")

    return out_rows


if __name__ == "__main__":
    download_clinvar()
    stream_filter()
    print("\n[OK] Step 1 complete.")
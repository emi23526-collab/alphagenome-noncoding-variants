"""Probe v3: find a supported condition, then dump the score_variant result structure."""
import inspect
import os

from alphagenome.data.genome import Interval, Variant
from alphagenome.models.v2 import dna_client, variant_scorers
from alphagenome.models.v2.dna_model import Organism

# --- 1. Client ---
key = os.environ.get("ALPHA_GENOME_API_KEY", "")
if not key:
    raise SystemExit("API KEY MISSING - reload the key first")
model = dna_client.create(key)
print("OK: client created")

# --- 2. Try to list supported conditions from model metadata (free, not a scoring call) ---
found = []
try:
    meta = model.output_metadata(Organism.HOMO_SAPIENS)
    print("OK: output_metadata loaded:", type(meta))
    for attr in [a for a in dir(meta) if not a.startswith("_")]:
        df = getattr(meta, attr, None)
        if hasattr(df, "columns") and "ontology_curie" in df.columns:
            name_col = "biosample_name" if "biosample_name" in df.columns else None
            print(f"  table '{attr}': {df.shape[0]} rows, columns={list(df.columns)[:8]}")
            if name_col:
                hits = df[df[name_col].astype(str).str.contains(
                    "breast|mammary|pancrea|islet|beta", case=False, na=False)]
                for _, row in hits.drop_duplicates("ontology_curie").head(8).iterrows():
                    print(f"    MATCH: {row['ontology_curie']}  {row[name_col]}")
                    found.append(str(row["ontology_curie"]))
            sample = df["ontology_curie"].dropna().unique()[:10]
            print("    sample curies:", list(sample))
except Exception as e:
    print("output_metadata not available:", type(e).__name__, str(e)[:150])
    print("public model methods:", [m for m in dir(model) if not m.startswith("_")])

# --- 3. Build the test query ---
BRCA1_POS = 43092028
variant = Variant(
    chromosome="chr17",
    position=BRCA1_POS,
    reference_bases="G",
    alternate_bases="A",
)
ref_interval = Interval("chr17", BRCA1_POS - 524_288, BRCA1_POS + 524_288)  # width 1,048,576
scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]

# --- 4. Try conditions until one is accepted ---
candidates = found[:4] + [
    "UBERON:0000310",  # breast
    "UBERON:0001264",  # pancreas
    "EFO:0001203",     # MCF-7 (breast cancer cell line)
    "EFO:0001187",     # HepG2
    "EFO:0002067",     # K562
    "EFO:0002784",     # GM12878
]

output = None
used = None
for curie in candidates:
    try:
        output = model.score_variant(
            interval=ref_interval,
            variant=variant,
            variant_scorers=[scorer],
            condition=curie,
            organism=Organism.HOMO_SAPIENS,
        )
        used = curie
        print(f"\nACCEPTED condition: {curie}")
        break
    except Exception as e:
        print(f"rejected {curie}: {type(e).__name__}: {str(e)[:110]}")

if output is None:
    raise SystemExit("\nNo candidate accepted. Paste everything above.")

# --- 5. Dump result structure ---
print("\n=== RESULT STRUCTURE ===")
print("condition used:", used)
print("type(output):", type(output), "len:", len(output))
for i, ad in enumerate(output):
    print(f"\n--- output[{i}] (AnnData) ---")
    print("  obs:", ad.obs.shape, list(ad.obs.columns)[:10])
    print("  var:", ad.var.shape, list(ad.var.columns)[:10])
    print("  X type:", type(ad.X), "shape:", getattr(ad.X, "shape", None))
    try:
        print("  X[:3]:", ad.X[:3])
    except Exception as e:
        print("  X preview failed:", e)
    try:
        print("  obs head:\n", ad.obs.head(3).to_string())
        print("  var head:\n", ad.var.head(3).to_string())
    except Exception as e:
        print("  head failed:", e)
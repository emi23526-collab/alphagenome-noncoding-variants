"""Check which pancreas conditions have RNA-seq and contain INS in the results."""
import os
import pandas as pd
from alphagenome.data.genome import Interval, Variant
from alphagenome.models.v2 import dna_client, variant_scorers
from alphagenome.models.v2.dna_model import Organism

pd.set_option("display.width", 250)

CANDIDATES = [
    "UBERON:0001264",
    "UBERON:0000016",
    "CL:0000169",
    "UBERON:0000006_CL:0000169",
]

# --- 1. What the table says about RNA-seq support (free) ---
df = pd.read_csv("data/processed/supported_conditions.csv")
sub = df[df["condition"].isin(CANDIDATES)]
print(sub[["condition", "biosample_type", "biosample_name", "rna_seq"]].to_string(index=False))

# --- 2. Test-score an INS variant with each candidate ---
model = dna_client.create(os.environ["ALPHA_GENOME_API_KEY"])
scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]

INS_POS = 2160500  # inside INS (chr11, hg38)
variant = Variant(chromosome="chr11", position=INS_POS,
                  reference_bases="G", alternate_bases="A")
interval = Interval("chr11", INS_POS - 524_288, INS_POS + 524_288)

for cond in CANDIDATES:
    print(f"\n=== {cond} ===")
    try:
        out = model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=[scorer],
            condition=cond,
            organism=Organism.HOMO_SAPIENS,
        )
        ad = out[0]
        print("  X shape:", ad.X.shape, "| tracks:", ad.var.shape[0])
        print("  track names:", list(ad.var["name"])[:3])
        rows = ad.obs.index[ad.obs["gene_name"] == "INS"]
        if len(rows) == 0:
            print("  INS NOT in results")
        else:
            i = ad.obs.index.get_loc(rows[0])
            print("  INS found. mean score over tracks:", float(ad.X[i].mean()))
    except Exception as e:
        print(f"  FAILED: {type(e).__name__}: {str(e)[:150]}")
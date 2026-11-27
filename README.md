# AlphaGenome scoring of non-coding BRCA1 and INS ClinVar variants

Computational study of how well AlphaGenome variant-effect scores separate pathogenic from benign non-coding ClinVar variants in two genes: BRCA1 (breast tissue) and INS (pancreatic beta cell). All results are computational flags only, not clinical evidence. This analysis does not follow the original project protocol in full: it was not pre-registered, and it used the RNA_SEQ and SPLICE_SITES scorers instead of DNase/ChIP features.

## Data

- Variants: non-coding BRCA1 and INS variants from ClinVar's bulk table (GRCh38 coordinates only).
- Groups used: P/LP (pathogenic or likely pathogenic), B/LB (benign or likely benign), VUS, Conflicting.
- Variants with identical reference and alternate base were removed (invalid entries).
- RNA-seq analysis: 2,461 scored non-coding variants (BRCA1: 1,603 B/LB, 279 P/LP, 352 VUS, 152 Conflicting. INS: 33 B/LB, 5 P/LP, 24 VUS, 13 Conflicting).
- Splice analysis: 1,457 splice-region variants (within 50 bp of an exon) scored with both scorers (BRCA1: 798 B/LB, 275 P/LP, 224 VUS, 132 Conflicting. INS: 15 B/LB, 1 P/LP, 8 VUS, 4 Conflicting).
- Distance to the exon was read from the HGVS name (for example c.671-2A>G gives 2 bp) and is only available for intronic variants.

## Methods

- Model: AlphaGenome, 524,288 bp on each side of the variant.
- Scorers: RNA_SEQ (all 2,461 variants) and SPLICE_SITES (splice-region variants only), both from the recommended scorers.
- Tissue context: BRCA1 uses breast (UBERON:0000310). INS uses pancreatic beta cell (CL:0000169).
- Score: score_absmax, the maximum absolute predicted effect across tracks for the target gene.
- Evaluation: AUROC and AUPRC for P/LP vs B/LB, bootstrap 95% confidence intervals (1,000 resamples), Mann-Whitney test.
- Baseline: closeness to the exon (negative distance) as the only predictor.
- Combined models: logistic regression on log distance and log scores, 5-fold cross-validation, 10 repeats.
- ENCODE: overlap of each variant's listed position with ENCODE candidate cis-regulatory elements (Registry V3, GRCh38).

## Results: RNA-seq scorer (2,461 variants)

Median score_absmax was higher for P/LP than B/LB in both genes (BRCA1 0.025 vs 0.0017; INS 0.339 vs 0.018).

BRCA1 intronic variants (277 P/LP vs 1,513 B/LB, including deep intronic benign variants):

| Predictor | AUROC | Notes |
|---|---|---|
| RNA-seq score | 0.931 (95% CI 0.919 to 0.943) | AUPRC 0.655 (chance 0.155) |
| Closeness to exon only | 0.972 (95% CI 0.959 to 0.983) | AUPRC 0.931 |
| Distance + RNA-seq score (cross-validated) | 0.985 | |

Distance wins overall because 202 of the 277 intronic P/LP variants are at the canonical +/-1 or +/-2 splice sites, where only 2 benign variants exist.

INS (5 P/LP vs 33 B/LB): AUROC 0.82 (95% CI 0.54 to 1.00). Descriptive only.

## Results: splice scorer (BRCA1 splice-region variants)

AUROC for P/LP vs B/LB (95% bootstrap CI):

| Subset | n P/LP | n B/LB | Distance only | RNA-seq score | Splice score |
|---|---|---|---|---|---|
| All splice-region | 275 | 798 | 0.960 (0.942 to 0.975) | 0.918 (0.901 to 0.933) | 0.994 (0.987 to 0.998) |
| 3 to 50 bp | 73 | 796 | 0.852 (0.789 to 0.910) | 0.914 (0.887 to 0.937) | 0.986 (0.963 to 0.997) |
| 3 to 10 bp | 60 | 215 | 0.836 (0.773 to 0.890) | 0.926 (0.891 to 0.957) | 0.995 (0.988 to 1.000) |
| 11 to 50 bp | 13 | 581 | 0.510 (0.318 to 0.689) | 0.884 (0.834 to 0.931) | 0.946 (0.832 to 0.999) |

Cross-validated combined models (AUROC):

| Subset | Distance | Distance + RNA-seq | Distance + splice | All three |
|---|---|---|---|---|
| All splice-region | 0.959 | 0.981 | 0.994 | 0.995 |
| 3 to 50 bp | 0.849 | 0.941 | 0.980 | 0.986 |

The splice scorer separates P/LP from B/LB better than distance alone, also outside the canonical +/-1/+/-2 sites (3 to 50 bp). The cross-validated distance-only AUROC in the 11 to 50 bp bin (0.339, below chance) comes from having only 13 P/LP variants and should be read as noise.

Candidate variants: VUS or Conflicting variants scoring above the 95th percentile of the benign scores. RNA-seq score: BRCA1 81, INS 4. Splice score: BRCA1 168, INS 5. Top splice-score VUS include variants at c.548-1, c.548-2 and c.4358-2. See results/top_candidates.csv and results/splice_candidates.csv. These are computational flags only and need manual review of the ClinVar evidence.

## Results: ENCODE cCRE overlap (non-splice candidates)

Variants were overlapped with ENCODE candidate cis-regulatory elements (Registry V3, GRCh38, 1,063,878 elements) using the variant's listed position (src\16_ccre_overlap.py).

Share of variants inside a cCRE: BRCA1 B/LB 50 of 1,603, Conflicting 4 of 152, P/LP 1 of 279, VUS 18 of 352. INS B/LB 14 of 33, Conflicting 6 of 13, P/LP 2 of 5, VUS 8 of 24. Overlap does not separate pathogenic from benign variants in this data; it is context for the UTR and promoter candidates.

Candidates (VUS or Conflicting, RNA-seq score above the benign 95th percentile):

- BRCA1 5' UTR / promoter variants inside a promoter-like, CTCF-bound cCRE: c.-106C>G, c.-106C>T, c.-166A>T, c.-70A>G, c.-125C>T, c.-62G>A, c.-61C>T.
- Highest-scoring INS candidates, not inside a cCRE: c.-187_-164del (score 1.28), c.-88A>G (0.56), c.-153C>G (0.26). Only the listed position was checked, so a deletion could still partly overlap an element.
- BRCA1 3' UTR c.*1363A>T and c.*1358A>G, not inside a cCRE.

GTEx and MPRA data were not checked, so all of these remain computational-only flags. See results/candidates_with_ccre.csv.

## Limitations

- ClinVar P/LP and B/LB labels for splice variants can themselves rely on splice-prediction tools, so the AUROCs may be optimistic (circularity). Overlap of AlphaGenome's training data with ClinVar was not assessed.
- INS has only 5 P/LP variants (1 in the splice region), so it is a descriptive case study, not a performance estimate.
- Almost all BRCA1 P/LP variants (275 of 279) are in the splice region. Only 4 are outside it, so deep intronic, UTR and promoter variants cannot be evaluated.
- The 11 to 50 bp bin has only 13 P/LP variants.
- The splice scorer was only run on splice-region variants.
- The benign 95th-percentile cutoff for candidates is lenient and gives many candidates; they are flags, not classifications.
- score_absmax (maximum absolute effect across tracks) is one of several possible summaries of the scores.
- cCRE overlap uses the listed variant position only; GTEx eQTL and MPRA validation were not performed.
- The analysis was not pre-registered.

## Reproducing

Run from the project folder in an environment with the alphagenome package and an API key loaded in ALPHA_GENOME_API_KEY. The ENCODE step needs data\external\GRCh38-cCREs.bed (downloaded from ENCODE SCREEN, Registry V3).

    python -u src\10_score_noncoding.py
    python -u src\11_analysis.py
    python -u src\12_stratified.py
    python -u src\13_score_splice.py
    python -u src\14_splice_analysis.py
    python -u src\16_ccre_overlap.py

Scripts 11, 12, 14 and 16 make no API calls. Outputs go to results\ (figures, performance.csv, top_candidates.csv, stratified_bins.csv, stratified_combined.csv, splice_comparison.csv, splice_combined.csv, splice_candidates.csv, candidates_with_ccre.csv, scored_noncoding_merged.csv).

## Pipeline: all scripts in src\

    01_pull_clinvar.py        pull BRCA1 and INS variants from ClinVar
    02_filter_gnomad.py       filter variants using gnomAD
    03_probe.py               probe the AlphaGenome API (tissues, settings)
    03_query_alphagenome.py   first AlphaGenome queries
    04_conditions.py          choose tissue / cell-type conditions
    05_ins_test.py            test scoring on INS
    06_score_variants.py      first scoring run
    07_score_hg38.py          scoring after the hg19/hg38 coordinate fix (1,500 variants)
    08_check_data.py          data checks (most scored variants turned out to be coding)
    09_pull_noncoding.py      non-coding variants from ClinVar's bulk table (GRCh38)
    10_score_noncoding.py     score the non-coding variants with RNA_SEQ (API calls)
    11_analysis.py            merge scores, AUROC/AUPRC, figures, candidates
    12_stratified.py          distance-bin test and combined model
    13_score_splice.py        score splice-region variants with SPLICE_SITES (API calls)
    14_splice_analysis.py     compare distance, RNA-seq and splice scores
    15_check_tasks.py         count scored variants per file
    16_ccre_overlap.py        overlap variants with ENCODE cCREs

Early scores used a hg19/hg38 coordinate mix; this was fixed and all reported results use GRCh38.
# app-pip

A scalable Snakemake pipeline for multi-layer genetic characterisation of the
**antigen presentation pathway (APP)**. It turns phased, imputed genome-wide
genotypes and HLA allele calls into a single, analysis-ready per-sample feature
matrix spanning the pathway, together with a feature dictionary (provenance)
and a QC summary.

## What it produces

| Layer | Source | Example columns |
| --- | --- | --- |
| Variant features | catalogue variants, additive ALT-allele dosage (0/1/2) | `rs2549782`, `rs1057141` |
| Allotypes | `allotype.R`, genes with ≥2 defining missense residues | `TAP1.GPVSVVCVLIDQRK` (one residue per defining site) |
| HLA functional features | [`hlafactor`](https://github.com/QiJingS/hlafactor) | `TD.HLA-B`, `Exp.HLA-C`, `Het.A`, `Supertype.B07`, `NK.C1C2` |

Outputs (under `results/`):

- `app_feature_matrix.tsv`: one row per sample, one column per feature
- `feature_dictionary.tsv`: definition and source of every column
- `qc_summary.txt`: variants extracted and retained, allotype sites used, feature counts
- `allotype_sites.tsv`: status of every allotype-defining site (used / flipped / absent after QC)

## Variant catalogue

`resources/app_variants_v2.tsv` lists 281 variants in 26 APP genes:

| Category | Genes |
| --- | --- |
| HLA class I processing and loading | PSMB5–10, PSME1, PSME2, TAP1, TAP2, ERAP1, ERAP2, TAPBP, TAPBPL, CALR, CANX, PDIA3 |
| HLA class II processing and loading | CD74, CTSS, CTSL, HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB |
| Transcriptional regulators | NLRC5 (class I), CIITA (class II) |

Coding and splice-region variants on the canonical transcript with gnomAD v4.1
exome allele frequency ≥ 0.001, with GRCh38 and GRCh37 positions, plus ten
literature-curated non-coding variants (see the `source` column). B2M has no
coding variant above this threshold.

`resources/allotype_definitions.tsv` lists, for every gene with ≥2 missense
residues in the catalogue, the defining sites with their protein positions;
residues are read N→C in that order. Which genes are reconstructed is set by
`allotype_genes` in the config.

Both tables are plain TSVs: edit them to add variants or change
allotype-defining sites, no code change needed.

## Requirements

- `snakemake`, `bcftools` (with `tabix`/`bgzip`)
- `R` (base installation only; no extra R packages are needed)
- [`hlafactor`](https://github.com/QiJingS/hlafactor): download the prebuilt
  binary for Linux x86_64 or macOS arm64 from its releases page and put
  `bin/hlafactor` on `PATH`, or set `hlafactor.binary` in the config to its path

```bash
conda env create -f environment.yml && conda activate app-pip
```

## HLA functional features

`hlafactor` is run once per module listed in `hlafactor.modules`; its outputs
are converted to numeric columns:

| Module | Columns | Encoding |
| --- | --- | --- |
| `-het` | `Het.HLA-A`, `Het.HLA-DRB1`, … | 1 = heterozygous, 0 = homozygous |
| `-hetf` | `FH.HLA-A/B/C` | functional heterozygosity (continuous) |
| `-exp` | `Exp.HLA-A`, `Exp.HLA-C` | expression proxy (continuous) |
| `-tapasin` | `TD.HLA-A/B/C`, `TD.global` | tapasin dependence (continuous) |
| `-sup` | `Supertype.A01`, `Supertype.B07`, … | supertype copies (0/1/2) |
| `-lignk` | `NK.C1`, `NK.C2`, `NK.Bw4`, `NK.Bw6` | ligand-category copies (0/1/2) |

Missing values (`NA`, `nan`, or a locus with only one called allele) become `NA`.
HLA QC (region, MAF, imputation R²) is done with `bcftools` so that any INFO
field (`R2`, `DR2`, …) can be used.

## Inputs

1. **Phased, imputed genome-wide genotypes**: VCF/VCF.GZ, per chromosome or a
   single file, GRCh38 or GRCh37, contigs named `6` or `chr6`.
   **Genotypes must be phased** (`GT` separator `|`); if any unphased call is
   found the pipeline stops and asks you to phase the input first.
2. **HLA allele calls**: imputed (e.g. BEAGLE/Minimac, `HLA_A*01:01` style IDs)
   or molecularly typed, as VCF/VCF.GZ on the same build. Both 2-field and
   4-field calls are retained, as `hlafactor` uses both.

## Usage

```bash
# 1. edit paths and thresholds
$EDITOR config/config.yaml

# 2. run
snakemake --configfile config/config.yaml --cores 8
```

## Testing

An end-to-end test on simulated data (with a stand-in for `hlafactor`) is
included:

```bash
bash test/run_test.sh
```

## Layout

```
app-pip/
├── config/
│   ├── config.yaml            # user choices
│   └── config.schema.yaml     # validated at start-up; unknown keys are rejected
├── workflow/
│   ├── Snakefile
│   └── scripts/
│       ├── snp_features.R      # additive SNP dosages
│       ├── allotype.R          # multi-variant allotype reconstruction
│       └── combine_matrix.R    # merge layers + dictionary + QC
├── resources/
│   ├── app_variants_v2.tsv           # versioned variant catalogue
│   └── allotype_definitions.tsv      # defining missense sites, ordered N→C
└── test/                             # simulated-data smoke test
```

## Processing notes

- Variants are extracted at the exact catalogue positions; multi-allelic
  records are split, and MAF is recomputed from the data before filtering.
- Variants lacking the imputation-quality INFO field (e.g. directly typed) are kept.
- Allotype sites are matched on CHROM, POS, REF and ALT; a site with REF/ALT
  swapped in the input is matched and flipped. A missing call at any defining
  site sets that sample's allotypes for the gene to NA.
- A requested `hlafactor` module whose output is missing stops the run, so no
  feature is dropped silently.
- Sample IDs are preserved verbatim.

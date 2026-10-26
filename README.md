# app-pip

A scalable Snakemake pipeline for multi-layer genetic characterisation of the
**antigen presentation pathway (APP)**. It turns phased, imputed genome-wide
genotypes and HLA allele calls into a single, analysis-ready per-sample feature
matrix spanning the pathway, together with a feature dictionary (provenance)
and a QC summary.

## What it produces

| Layer | Source | Example columns |
| --- | --- | --- |
| SNP features | APP genes, dosage under a genetic model | `rs2549794`, `TAP1:…` |
| Allotypes | `allotype` (R), genes with ≥2 missense sites | `TAP1.IDAVRVG`, `ERAP2.TN` |
| HLA functional features | [`hlafactor`](https://github.com/QiJingS/hlafactor) | `TD.HLA-B`, `Exp.HLA-C`, `Het.A`, `Supertype.B07`, `NK.C1C2` |

Outputs (under `results/`):

- `app_feature_matrix.tsv` — one row per sample, one column per feature
- `feature_dictionary.tsv` — definition and source of every column
- `qc_summary.txt` — feature counts and run summary

## Requirements

- `snakemake`, `bcftools` (with `tabix`/`bgzip`)
- `R` with the `vcfR` package
- the `hlafactor` binary on `PATH` (or set its path in the config)

## Inputs

1. **Phased, imputed WGS genotypes** — VCF/VCF.GZ, per chromosome or a single
   file. **Genotypes must be phased** (`GT` separator `|`). If any unphased
   call is found the pipeline stops and asks you to phase the input first.
2. **HLA allele calls** — imputed (e.g. BEAGLE/Minimac, `HLA_A*01:01` style IDs)
   or molecularly typed, as VCF/VCF.GZ.

## Usage

```bash
# 1. edit paths and thresholds
$EDITOR config/config.yaml

# 2. run
snakemake --configfile config/config.yaml --cores 8
```

## Layout

```
app-pip/
├── config/
│   ├── config.yaml            # user choices
│   └── config.schema.yaml     # validated at start-up
├── workflow/
│   ├── Snakefile
│   └── scripts/
│       ├── snp_features.R      # SNP dosage under the genetic model
│       ├── allotype.R          # multi-variant allotype reconstruction
│       └── combine_matrix.R    # merge layers + dictionary + QC
└── resources/
    ├── app_variants.tsv              # versioned variant reference table
    └── allotype_definitions.tsv      # defining missense sites, ordered N→C
```

## Reference tables

`resources/app_variants.tsv` and `resources/allotype_definitions.tsv` are plain,
versioned TSVs derived from the curated APP variant catalogue. Edit them to add
variants, change allotype-defining sites, or update annotations — no code change
or recompilation needed. `allotype_definitions.tsv` lists the defining missense
sites per gene with their protein positions; sites are read N→C in that order.

## Notes / changes from the prototype

- HLA functional features are computed by `hlafactor` (not re-implemented here),
  keeping this pipeline and `hlafactor` consistent.
- The genetic model in the config is now honoured.
- SNP filtering uses minor allele frequency (not raw ALT frequency); the
  imputation-quality INFO field is configurable (`R2`, `DR2`, …).
- Extraction is region-based (one pass per file) rather than per-SNP.
- Sample IDs are preserved verbatim (no automatic `X` prefix).
- Allotype sites are ordered by protein position and missing calls are handled
  explicitly rather than mis-assigned.

#!/usr/bin/env bash
# End-to-end smoke test on simulated data (no real genotypes needed).
# Usage (from the repository root):  bash test/run_test.sh
set -euo pipefail
python3 test/make_test.py
bgzip -f test/wgs.vcf && tabix -f -p vcf test/wgs.vcf.gz
bgzip -f test/hla.vcf && tabix -f -p vcf test/hla.vcf.gz
chmod +x test/fake_hlafactor
sed -e 's|/path/to/wgs/chr{chrom}.phased.vcf.gz|test/wgs.vcf.gz|' \
    -e 's|/path/to/hla.phased.vcf.gz|test/hla.vcf.gz|' \
    -e 's|binary: "hlafactor"|binary: "test/fake_hlafactor"|' \
    -e 's|output_dir: "results"|output_dir: "test/results"|' \
    config/config.yaml > test/config.test.yaml
snakemake -s workflow/Snakefile --configfile test/config.test.yaml --cores 2 -q --forceall
cat test/results/qc_summary.txt

# Reconstruct multi-variant protein allotypes from PHASED genotypes.
#
# For each gene with >= 2 defining missense sites present after QC, the amino
# acid encoded on each haplotype is read at every defining site, ordered from
# the N- to the C-terminus (by protein position), and concatenated into an
# allotype string. Each sample is assigned a count (0/1/2) per observed allotype.
#
# Robustness fixes vs the original script:
#   * sites are ordered by protein position, not spreadsheet row order
#   * missing / unphased calls at a defining site -> that sample's allotypes
#     for the gene are set to NA (not silently mis-assigned)
#   * matching uses CHROM + POS + REF/ALT, not POS alone
#   * sample IDs are preserved verbatim
suppressMessages(library(vcfR))

vcf_file <- snakemake@input[["vcf"]]
defs_file <- snakemake@input[["defs"]]
build    <- as.character(snakemake@params[["build"]])
min_freq <- as.numeric(snakemake@params[["min_freq"]])
out_file <- snakemake@output[["tsv"]]

defs <- read.table(defs_file, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)
poscol <- if (build == "38") "pos38" else "pos37"

vcf <- read.vcfR(vcf_file, verbose = FALSE)
fix <- as.data.frame(getFIX(vcf), stringsAsFactors = FALSE)
fix$POS <- as.character(fix$POS)
gt_all <- extract.gt(vcf, element = "GT")
samples <- colnames(gt_all)

allo_tab <- data.frame(SampleID = samples, stringsAsFactors = FALSE,
                       check.names = FALSE)

for (gene in unique(defs$gene)) {
  sites <- defs[defs$gene == gene, ]
  sites <- sites[order(suppressWarnings(as.integer(sites$aa_position))), ]  # N -> C
  sites$pos <- as.character(sites[[poscol]])

  # keep defining sites that are present in the QC'd VCF (match on chrom+pos)
  key_vcf  <- paste(fix$CHROM, fix$POS)
  key_site <- paste(sites$chrom, sites$pos)
  present  <- key_site %in% key_vcf
  sites <- sites[present, ]
  if (nrow(sites) < 2) next           # not a multi-variant allotype gene here

  idx <- match(paste(sites$chrom, sites$pos), key_vcf)
  gt  <- gt_all[idx, , drop = FALSE]

  # build the two haplotype amino-acid strings per sample
  hap1 <- character(length(samples)); hap2 <- character(length(samples))
  valid <- rep(TRUE, length(samples))
  for (s in seq_len(nrow(sites))) {
    g <- gt[s, ]
    a <- sub("[|/].*$", "", g); b <- sub("^.*[|/]", "", g)
    aa_a <- ifelse(a == "0", sites$ref_aa[s],
             ifelse(a == "1", sites$alt_aa[s], NA))
    aa_b <- ifelse(b == "0", sites$ref_aa[s],
             ifelse(b == "1", sites$alt_aa[s], NA))
    valid <- valid & !is.na(aa_a) & !is.na(aa_b)
    hap1 <- paste0(hap1, ifelse(is.na(aa_a), "?", aa_a))
    hap2 <- paste0(hap2, ifelse(is.na(aa_b), "?", aa_b))
  }

  # per-sample allotype counts
  all_types <- sort(unique(c(hap1[valid], hap2[valid])))
  grid <- matrix(0, nrow = length(samples), ncol = length(all_types),
                 dimnames = list(samples, all_types))
  for (i in seq_along(samples)) {
    if (!valid[i]) { grid[i, ] <- NA; next }
    grid[i, hap1[i]] <- grid[i, hap1[i]] + 1
    grid[i, hap2[i]] <- grid[i, hap2[i]] + 1
  }

  # frequency filter (ignoring NA samples)
  denom <- sum(valid) * 2
  if (denom > 0) {
    keep <- colSums(grid, na.rm = TRUE) / denom >= min_freq
    grid <- grid[, keep, drop = FALSE]
  }
  if (ncol(grid) == 0) next
  colnames(grid) <- paste0(gene, ".", colnames(grid))
  allo_tab <- cbind(allo_tab, as.data.frame(grid, check.names = FALSE))
}

write.table(allo_tab, out_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

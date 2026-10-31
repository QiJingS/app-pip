# Recode phased genotypes into additive ALT-allele dosages (0/1/2).
# Multi-allelic sites have been split upstream, so each record has one ALT.
# Sample IDs are preserved verbatim (no automatic "X" prefixing).
suppressMessages(library(vcfR))

vcf_file <- snakemake@input[["vcf"]]
out_file <- snakemake@output[["tsv"]]

vcf <- read.vcfR(vcf_file, verbose = FALSE)
gt  <- extract.gt(vcf, element = "GT", IDtoRowNames = FALSE)          # rows = variants, cols = samples
fix <- getFIX(vcf)

# Variant identifier: rsID (or CHROM:POS if absent); records sharing an
# identifier (split multi-allelic sites) get their ALT allele appended.
vid   <- fix[, "ID"]
nomid <- is.na(vid) | vid == "." | vid == ""
vid[nomid] <- paste0(sub("^chr", "", fix[nomid, "CHROM"]), ":", fix[nomid, "POS"])
dup <- vid %in% vid[duplicated(vid)]
vid[dup] <- paste0(vid[dup], "_", fix[dup, "ALT"])
rownames(gt) <- make.unique(vid)

recode <- function(g) {
  a <- suppressWarnings(as.integer(sub("[|/].*$", "", g)))
  b <- suppressWarnings(as.integer(sub("^.*[|/]", "", g)))
  a + b                                          # NA if either allele missing
}

mat <- apply(gt, 1, recode)                      # rows = samples, cols = variants
if (is.null(dim(mat))) mat <- matrix(mat, ncol = nrow(gt),
                                     dimnames = list(colnames(gt), rownames(gt)))
rownames(mat) <- colnames(gt)
mat <- as.data.frame(mat, check.names = FALSE)
out <- cbind(SampleID = rownames(mat), mat)

write.table(out, out_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

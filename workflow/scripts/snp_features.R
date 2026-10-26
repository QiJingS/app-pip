# Recode phased SNP genotypes into dosages under the chosen genetic model.
# Sample IDs are preserved verbatim (no automatic "X" prefixing).
suppressMessages(library(vcfR))

vcf_file <- snakemake@input[["vcf"]]
model    <- snakemake@params[["model"]]
out_file <- snakemake@output[["tsv"]]

vcf <- read.vcfR(vcf_file, verbose = FALSE)
gt  <- extract.gt(vcf, element = "GT")          # rows = variants, cols = samples

# Variant identifier: prefer rsID, fall back to CHROM:POS
fix   <- getFIX(vcf)
vid   <- fix[, "ID"]
nomid <- is.na(vid) | vid == "." | vid == ""
vid[nomid] <- paste0(fix[nomid, "CHROM"], ":", fix[nomid, "POS"])
rownames(gt) <- make.unique(vid)

recode <- function(g, model) {
  # accept phased ("a|b") or unphased ("a/b"); count ALT alleles
  a <- suppressWarnings(as.integer(sub("[|/].*$", "", g)))
  b <- suppressWarnings(as.integer(sub("^.*[|/]", "", g)))
  dose <- a + b
  if (model == "additive") {
    out <- dose
  } else if (model == "dominant") {
    out <- ifelse(dose >= 1, 1, 0)
  } else if (model == "recessive") {
    out <- ifelse(dose == 2, 1, 0)
  } else {
    stop(sprintf("unknown genetic_model: %s", model))
  }
  out[is.na(dose)] <- NA
  out
}

mat <- apply(gt, 1, recode, model = model)       # rows = samples, cols = variants
mat <- as.data.frame(mat, check.names = FALSE)
out <- cbind(SampleID = rownames(mat), mat)

write.table(out, out_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

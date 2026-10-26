# Merge every feature layer into one per-sample matrix, and emit a feature
# dictionary (provenance) and a QC summary. Sample IDs preserved verbatim.

snps_file <- snakemake@input[["snps"]]
allo_file <- snakemake@input[["allo"]]
hladir    <- snakemake@params[["hladir"]]
modules   <- snakemake@params[["modules"]]

out_matrix <- snakemake@output[["matrix"]]
out_dict   <- snakemake@output[["dict"]]
out_qc     <- snakemake@output[["qc"]]

rd <- function(f) {
  d <- read.table(f, header = TRUE, sep = "	", stringsAsFactors = FALSE,
                  check.names = FALSE, quote = "", comment.char = "", fill = TRUE)
  nm <- names(d)
  d[, !(is.na(nm) | trimws(nm) == ""), drop = FALSE]
}

snps <- rd(snps_file)
allo <- rd(allo_file)

dict <- data.frame(feature = character(), source = character(),
                   description = character(), stringsAsFactors = FALSE)
add_dict <- function(cols, source, desc) {
  if (length(cols) == 0) return(invisible())
  dict <<- rbind(dict, data.frame(feature = cols, source = source,
                                  description = desc, stringsAsFactors = FALSE))
}

matrix <- snps
add_dict(setdiff(colnames(snps), "SampleID"), "app-pip:snp_features",
         "SNP dosage under the configured genetic model")

if (ncol(allo) > 1) {
  matrix <- merge(matrix, allo, by = "SampleID", all = TRUE)
  add_dict(setdiff(colnames(allo), "SampleID"), "app-pip:allotype",
           "Multi-variant protein allotype count (0/1/2)")
}

# ---- helper: read a hlafactor output, rename the sample column to SampleID,
#      keep selected value columns renamed with a prefix
take <- function(path, value_cols, newnames, source, desc) {
  if (!file.exists(path)) return(invisible())
  d <- rd(path)
  scol <- intersect(c("SampleID", "sampleid"), colnames(d))[1]
  colnames(d)[colnames(d) == scol] <- "SampleID"
  value_cols <- intersect(value_cols, colnames(d))
  if (length(value_cols) == 0) return(invisible())
  sub <- d[, c("SampleID", value_cols), drop = FALSE]
  colnames(sub) <- c("SampleID", newnames[seq_along(value_cols)])
  matrix <<- merge(matrix, sub, by = "SampleID", all = TRUE)
  add_dict(newnames[seq_along(value_cols)], source, desc)
}

p <- function(m, suffix) file.path(hladir, paste0(m, suffix))

for (m in modules) {
  if (m == "het") {
    f <- p(m, "_heterozygosity.txt")
    if (file.exists(f)) {
      d <- rd(f); colnames(d)[1] <- "SampleID"
      loci <- setdiff(colnames(d), "SampleID")
      colnames(d) <- c("SampleID", paste0("Het.", loci))
      matrix <- merge(matrix, d, by = "SampleID", all = TRUE)
      add_dict(paste0("Het.", loci), "hlafactor:-het",
               "HLA gene heterozygosity (Het/Hom)")
    }
  } else if (m == "hetf") {
    take(p(m, "_functional_zygosity.txt"),
         c("HLA_A_FH", "HLA_B_FH", "HLA_C_FH"),
         c("FH.HLA-A", "FH.HLA-B", "FH.HLA-C"),
         "hlafactor:-hetf", "HLA functional heterozygosity (allotype divergence)")
  } else if (m == "exp") {
    take(p(m, "_exp_A.txt"), "expression_value", "Exp.HLA-A",
         "hlafactor:-exp", "HLA-A allele-specific expression proxy")
    take(p(m, "_exp_C.txt"), "expression_value", "Exp.HLA-C",
         "hlafactor:-exp", "HLA-C allele-specific expression proxy")
  } else if (m == "tapasin") {
    take(p(m, "_taps_A.txt"), "tapsA_value", "TD.HLA-A", "hlafactor:-tapasin", "Tapasin dependence HLA-A")
    take(p(m, "_taps_B.txt"), "tapsB_value", "TD.HLA-B", "hlafactor:-tapasin", "Tapasin dependence HLA-B")
    take(p(m, "_taps_C.txt"), "tapsC_value", "TD.HLA-C", "hlafactor:-tapasin", "Tapasin dependence HLA-C")
    take(p(m, "_taps_global.txt"), "taps_global_value", "TD.global", "hlafactor:-tapasin", "Tapasin dependence global")
  } else if (m == "sup") {
    f <- p(m, "_supertypes.txt")
    if (file.exists(f)) {
      d <- rd(f); colnames(d)[1] <- "SampleID"
      st <- setdiff(colnames(d), "SampleID")
      colnames(d) <- c("SampleID", paste0("Supertype.", st))
      matrix <- merge(matrix, d, by = "SampleID", all = TRUE)
      add_dict(paste0("Supertype.", st), "hlafactor:-sup", "HLA supertype count")
    }
  } else if (m == "lignk") {
    take(p(m, "_HLA_variant_NK.txt"),
         c("HLA_C_ligand", "HLA_B_ligand"),
         c("NK.C1C2", "NK.Bw4Bw6"),
         "hlafactor:-lignk", "NK cell ligand group / HLA-B -21 M/T")
  } else if (m == "amac") {
    # residue-level table is large & transposed; kept as a separate file
    message("[combine] -amac output kept separate (not merged into main matrix)")
  }
}

write.table(matrix, out_matrix, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
write.table(dict, out_dict, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

qc <- c(
  sprintf("app-pip feature matrix"),
  sprintf("samples:            %d", nrow(matrix)),
  sprintf("total features:     %d", ncol(matrix) - 1),
  sprintf("  SNP features:     %d", sum(dict$source == "app-pip:snp_features")),
  sprintf("  allotypes:        %d", sum(dict$source == "app-pip:allotype")),
  sprintf("  HLA features:     %d", sum(grepl("^hlafactor", dict$source))),
  sprintf("hlafactor modules:  %s", paste(modules, collapse = ", "))
)
writeLines(qc, out_qc)

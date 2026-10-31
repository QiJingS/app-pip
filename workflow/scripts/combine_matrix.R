# Merge every feature layer into one per-sample matrix, and emit a feature
# dictionary (provenance) and a QC summary. Sample IDs preserved verbatim.

snps_file   <- snakemake@input[["snps"]]
allo_file   <- snakemake@input[["allo"]]
sites_file  <- snakemake@input[["sites"]]
counts_file <- snakemake@input[["counts"]]
hladir      <- snakemake@params[["hladir"]]
modules     <- unlist(snakemake@params[["modules"]])
build       <- snakemake@params[["build"]]
variants    <- snakemake@params[["variants"]]

out_matrix <- snakemake@output[["matrix"]]
out_dict   <- snakemake@output[["dict"]]
out_qc     <- snakemake@output[["qc"]]

rd <- function(f) {
  d <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
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
         "ALT-allele dosage (additive, 0/1/2)")

if (ncol(allo) > 1) {
  matrix <- merge(matrix, allo, by = "SampleID", all = TRUE)
  add_dict(setdiff(colnames(allo), "SampleID"), "app-pip:allotype",
           "Copies of the protein allotype (0/1/2); residues ordered N- to C-terminus")
}

# ---- hlafactor outputs. A requested module whose output is missing is an
#      error (features are never dropped silently).
need <- function(path, m) {
  if (!file.exists(path))
    stop(sprintf("hlafactor module '-%s' did not produce expected file: %s", m, path))
  path
}
# hlafactor writes missing values as "NA", "nan" or "."; every feature is
# returned as numeric so the matrix can go straight into association tools.
num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.nan(x)] <- NA
  x
}
merge_in <- function(d, newnames, source, desc, numeric = TRUE) {
  colnames(d) <- c("SampleID", newnames)
  if (numeric) for (k in newnames) d[[k]] <- num(d[[k]])
  matrix <<- merge(matrix, d, by = "SampleID", all = TRUE)
  add_dict(newnames, source, desc)
}
# "HLA_A" -> "HLA-A"
locus <- function(x) sub("^HLA_", "HLA-", x)
# Genotype strings such as "C1/C2" or "BW6/BW6" -> per-category copy counts
# (0/1/2). A sample with a missing allele ("." or absent) is NA for the locus.
count_categories <- function(g) {
  parts <- strsplit(ifelse(is.na(g), "", g), "/", fixed = TRUE)
  ok <- vapply(parts, function(p) length(p) == 2 && !any(p %in% c("", ".", "NA")),
               logical(1))
  cats <- sort(unique(unlist(parts[ok])))
  out <- sapply(cats, function(cc) vapply(parts, function(p) sum(p == cc), numeric(1)))
  out <- matrix(out, nrow = length(g), dimnames = list(NULL, cats))
  out[!ok, ] <- NA
  out
}
take <- function(path, m, value_cols, newnames, source, desc) {
  d <- rd(need(path, m))
  scol <- intersect(c("SampleID", "sampleid"), colnames(d))[1]
  if (is.na(scol)) stop(sprintf("no sample ID column in %s", path))
  miss <- setdiff(value_cols, colnames(d))
  if (length(miss)) stop(sprintf("columns %s missing in %s",
                                 paste(miss, collapse = ", "), path))
  merge_in(d[, c(scol, value_cols), drop = FALSE], newnames, source, desc)
}

p <- function(m, suffix) file.path(hladir, paste0(m, suffix))

for (m in modules) {
  if (m == "het") {
    d <- rd(need(p(m, "_heterozygosity.txt"), m))
    for (k in colnames(d)[-1])
      d[[k]] <- ifelse(d[[k]] == "Het", 1, ifelse(d[[k]] == "Hom", 0, NA))
    merge_in(d, paste0("Het.", locus(colnames(d)[-1])), "hlafactor:-het",
             "HLA locus heterozygosity (1 = heterozygous, 0 = homozygous)")
  } else if (m == "hetf") {
    take(p(m, "_functional_zygosity.txt"), m,
         c("HLA_A_FH", "HLA_B_FH", "HLA_C_FH"),
         c("FH.HLA-A", "FH.HLA-B", "FH.HLA-C"),
         "hlafactor:-hetf", "HLA functional heterozygosity (allotype divergence)")
  } else if (m == "exp") {
    take(p(m, "_exp_A.txt"), m, "expression_value", "Exp.HLA-A",
         "hlafactor:-exp", "HLA-A allele-specific expression proxy")
    take(p(m, "_exp_C.txt"), m, "expression_value", "Exp.HLA-C",
         "hlafactor:-exp", "HLA-C allele-specific expression proxy")
  } else if (m == "tapasin") {
    take(p(m, "_taps_A.txt"), m, "tapsA_value", "TD.HLA-A", "hlafactor:-tapasin", "Tapasin dependence HLA-A")
    take(p(m, "_taps_B.txt"), m, "tapsB_value", "TD.HLA-B", "hlafactor:-tapasin", "Tapasin dependence HLA-B")
    take(p(m, "_taps_C.txt"), m, "tapsC_value", "TD.HLA-C", "hlafactor:-tapasin", "Tapasin dependence HLA-C")
    take(p(m, "_taps_global.txt"), m, "taps_global_value", "TD.global", "hlafactor:-tapasin", "Tapasin dependence global")
  } else if (m == "sup") {
    d <- rd(need(p(m, "_supertypes.txt"), m))
    merge_in(d, paste0("Supertype.", colnames(d)[-1]), "hlafactor:-sup",
             "HLA supertype count")
  } else if (m == "lignk") {
    d <- rd(need(p(m, "_HLA_variant_NK.txt"), m))
    scol <- intersect(c("SampleID", "sampleid"), colnames(d))[1]
    for (k in setdiff(colnames(d), scol)) {
      cnt <- count_categories(d[[k]])
      if (ncol(cnt) == 0) next
      cats <- sub("^BW", "Bw", colnames(cnt))
      nm <- paste0("NK.", cats)
      if (any(nm %in% dict$feature)) nm <- paste0("NK.", locus(k), ".", cats)
      merge_in(data.frame(d[[scol]], cnt, check.names = FALSE), nm,
               "hlafactor:-lignk",
               sprintf("Copies of NK-cell ligand category (%s; 0/1/2)", locus(k)))
    }
  } else if (m == "amac") {
    message("[combine] -amac output kept separate (not merged into main matrix)")
  } else {
    stop(sprintf("unknown hlafactor module: %s", m))
  }
}

write.table(matrix, out_matrix, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
write.table(dict, out_dict, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

counts <- rd(counts_file)
sites  <- rd(sites_file)
cnt <- function(step) counts$n[counts$step == step]
allo_genes <- unique(sites$gene[grepl("^used", sites$status)])
qc <- c(
  "app-pip feature matrix",
  sprintf("reference build:    GRCh%s", build),
  sprintf("variant catalogue:  %s", variants),
  sprintf("samples:            %d", nrow(matrix)),
  "",
  "variant extraction",
  sprintf("  catalogue sites:          %s", cnt("catalogue_sites")),
  sprintf("  records in input:         %s", cnt("extracted_records")),
  sprintf("  passed MAF / R2 filters:  %s", cnt("passed_maf_r2")),
  "",
  "allotype-defining sites",
  sprintf("  used:                     %d", sum(grepl("^used", sites$status))),
  sprintf("  absent after QC:          %d", sum(sites$status == "absent_after_qc")),
  sprintf("  genes with allotypes:     %s",
          if (length(allo_genes)) paste(allo_genes, collapse = ", ") else "none"),
  "",
  "features",
  sprintf("  total:                    %d", ncol(matrix) - 1),
  sprintf("  SNP features:             %d", sum(dict$source == "app-pip:snp_features")),
  sprintf("  allotypes:                %d", sum(dict$source == "app-pip:allotype")),
  sprintf("  HLA features:             %d", sum(grepl("^hlafactor", dict$source))),
  sprintf("hlafactor modules:  %s", paste(modules, collapse = ", "))
)
writeLines(qc, out_qc)

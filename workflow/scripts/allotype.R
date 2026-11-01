# Reconstruct multi-variant protein allotypes from PHASED genotypes.
#
# For each gene with >= 2 defining missense residues present after QC, the
# amino acid encoded on each haplotype is read at every defining residue,
# ordered from the N- to the C-terminus (by protein position), and the residues
# are concatenated into an allotype string. Each sample is assigned the number
# of copies (0/1/2) of every observed allotype.
#
#   * sites are matched on CHROM + POS + REF + ALT; a site whose REF/ALT are
#     swapped in the input VCF is matched and its alleles flipped
#   * several ALT alleles at the same residue are combined into one position
#   * a missing call at any defining site sets that sample's allotypes for the
#     gene to NA (never silently mis-assigned)
#   * every defining site is reported as used / absent in `allotype_sites.tsv`
#   * sample IDs are preserved verbatim
vcf_file  <- snakemake@input[["gt"]]
defs_file <- snakemake@input[["defs"]]
build     <- as.character(snakemake@params[["build"]])
genes_sel <- unlist(snakemake@params[["genes"]])
target    <- unlist(snakemake@params[["target"]])
min_freq  <- as.numeric(snakemake@params[["min_freq"]])
out_file  <- snakemake@output[["tsv"]]
site_file <- snakemake@output[["sites"]]

defs <- read.table(defs_file, header = TRUE, sep = "\t", quote = "",
                   stringsAsFactors = FALSE, check.names = FALSE,
                   colClasses = "character")
poscol <- if (build == "38") "pos38" else "pos37"
if (length(genes_sel)) defs <- defs[defs$gene %in% genes_sel, ]
if (length(target))    defs <- defs[defs$gene %in% target, ]

# Read the genotype table written by bcftools query (rule genotype_table):
# columns CHROM POS ID REF ALT, then one GT column per sample.
read_gt_table <- function(f) {
  d <- read.delim(f, header = TRUE, colClasses = "character", check.names = FALSE,
                  quote = "", comment.char = "", na.strings = character(0))
  list(fix = d[, 1:5, drop = FALSE],
       gt  = as.matrix(d[, -(1:5), drop = FALSE]))
}

tab <- read_gt_table(vcf_file)
fix <- tab$fix
fix$CHROM <- sub("^chr", "", fix$CHROM)
gt_all  <- tab$gt
samples <- colnames(gt_all)
key_fwd <- paste(fix$CHROM, fix$POS, fix$REF, fix$ALT, sep = ":")

allo_tab <- data.frame(SampleID = samples, stringsAsFactors = FALSE,
                       check.names = FALSE)
site_log <- list()

allele <- function(g, side) {
  x <- if (side == 1) sub("[|/].*$", "", g) else sub("^.*[|/]", "", g)
  suppressWarnings(as.integer(x))                # NA for "." / missing
}

for (gene in unique(defs$gene)) {
  sites <- defs[defs$gene == gene, ]
  sites$aa_position <- as.integer(sites$aa_position)
  sites <- sites[order(sites$aa_position), ]
  k_fwd <- paste(sites$chrom, sites[[poscol]], sites$ref, sites$alt, sep = ":")
  k_rev <- paste(sites$chrom, sites[[poscol]], sites$alt, sites$ref, sep = ":")
  idx   <- match(k_fwd, key_fwd)
  flip  <- is.na(idx) & !is.na(match(k_rev, key_fwd))
  idx[flip] <- match(k_rev[flip], key_fwd)
  sites$vcf_idx <- idx; sites$flipped <- flip

  residues <- sort(unique(sites$aa_position[!is.na(idx)]))
  used <- !is.na(idx) & length(residues) >= 2
  site_log[[gene]] <- data.frame(
    gene = gene, rsid = sites$rsid, aa_position = sites$aa_position,
    change = paste0(sites$ref_aa, sites$aa_position, sites$alt_aa),
    status = ifelse(is.na(idx), "absent_after_qc",
                    ifelse(used, ifelse(flip, "used_flipped", "used"),
                           "unused_lt2_residues")),
    stringsAsFactors = FALSE)
  if (length(residues) < 2) next                 # not a multi-variant gene here
  sites <- sites[!is.na(idx), ]

  n <- length(samples)
  hap <- list(character(n), character(n))
  valid <- rep(TRUE, n)
  for (r in residues) {
    s_r <- sites[sites$aa_position == r, ]
    for (h in 1:2) {
      aa <- rep(s_r$ref_aa[1], n)
      n_alt <- integer(n)
      for (j in seq_len(nrow(s_r))) {
        a <- allele(gt_all[s_r$vcf_idx[j], ], h)
        is_alt <- if (s_r$flipped[j]) a == 0 else a == 1
        valid <- valid & !is.na(is_alt)
        is_alt[is.na(is_alt)] <- FALSE
        aa[is_alt] <- s_r$alt_aa[j]
        n_alt <- n_alt + is_alt
      }
      valid <- valid & n_alt <= 1                # two ALTs on one haplotype: invalid
      hap[[h]] <- paste0(hap[[h]], aa)
    }
  }

  all_types <- sort(unique(c(hap[[1]][valid], hap[[2]][valid])))
  if (length(all_types) == 0) next
  grid <- matrix(0, nrow = n, ncol = length(all_types),
                 dimnames = list(samples, all_types))
  for (h in 1:2) {
    ix <- cbind(which(valid), match(hap[[h]][valid], all_types))
    grid[ix] <- grid[ix] + 1
  }
  grid[!valid, ] <- NA

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
site_df <- if (length(site_log)) do.call(rbind, site_log) else
  data.frame(gene = character(), rsid = character(), aa_position = integer(),
             change = character(), status = character())
write.table(site_df, site_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

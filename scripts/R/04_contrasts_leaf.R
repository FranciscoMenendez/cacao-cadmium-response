# scripts/R/04_contrasts_leaf.R
# Define leaf contrasts + the exact library indices used on each side (A,B).
# Output: data/interim/contrasts_leaf.rda
# Also writes human-readable tables under results/tables/.

source("scripts/R/00_setup.R")
load(file.path(PATHS$interim, "dd_objects.rda"))      # dd_leaf.le (from step 03)
load(file.path(PATHS$interim, "expsum_built.rda"))    # expsum_leaf (from step 01)

# model matrix (ensures names match current design in the object)
mm_leaf <- model.matrix(object = design(dd_leaf.le), data = as.data.frame(colData(dd_leaf.le)))
rownames(mm_leaf) <- as.character(dd_leaf.le$library)

# tiny helpers
mmean <- function(i) colMeans(mm_leaf[i, , drop = FALSE])    # mean of rows (libraries)
idx   <- function(g=NULL, t=NULL) {
  with(as.data.frame(colData(dd_leaf.le)),
       which((is.null(g) | genotype %in% g) &
               (is.null(t) | time     %in% t)))
}

# group indices (by genotype x time)
i_PA0  <- idx("PA121","0");   i_PA24 <- idx("PA121","24");  i_PA48 <- idx("PA121","48") # PA0 = libraries of PA121 at 0 hours, etc.
i_TS0  <- idx("TSH660","0");  i_TS24 <- idx("TSH660","24"); i_TS48 <- idx("TSH660","48") # TS0 = libraries of TSH660 at 0 hours, etc.
i_0    <- idx(t="0");         i_24   <- idx(t="24");        i_48   <- idx(t="48")
i_cd   <- idx(t=c("24","48")) # any Cd exposure (24 or 48)

# numeric contrast vectors (12)
c_PA_24_vs_0   <- mmean(i_PA24) - mmean(i_PA0)
c_PA_48_vs_0   <- mmean(i_PA48) - mmean(i_PA0)
c_TS_24_vs_0   <- mmean(i_TS24) - mmean(i_TS0)
c_TS_48_vs_0   <- mmean(i_TS48) - mmean(i_TS0)

c_PA_vs_TS_0   <- mmean(i_PA0)  - mmean(i_TS0)
c_PA_vs_TS_24  <- mmean(i_PA24) - mmean(i_TS24)
c_PA_vs_TS_48  <- mmean(i_PA48) - mmean(i_TS48)

c_PA_cd_vs_no  <- mmean(intersect(i_cd,  idx("PA121", c("24","48")))) - mmean(i_PA0)
c_TS_cd_vs_no  <- mmean(intersect(i_cd,  idx("TSH660", c("24","48")))) - mmean(i_TS0)
c_all_24_vs_0  <- mmean(i_24) - mmean(i_0)
c_all_48_vs_0  <- mmean(i_48) - mmean(i_0)

# interaction (ΔΔ) terms taken as named coefficients from the fit (no numeric vector needed)
coef_24  <- "genotypeTSH660.time24"
coef_48  <- "genotypeTSH660.time48"

# naming
leaf_names <- c(
  "leaf_pa121_24h_0h","leaf_pa121_48h_0h",
  "leaf_tsh660_24h_0h","leaf_tsh660_48h_0h",
  "leaf_pa121_tsh660_0h_0h","leaf_pa121_tsh660_24h_24h","leaf_pa121_tsh660_48h_48h",
  "leaf_pa121_cd_no_cd","leaf_tsh660_cd_no_cd",
  "leaf_24h_0h","leaf_48h_0h",
  # interaction terms (by name)
  "leaf_tsh660.24h","leaf_tsh660.48h"
)

# assemble
leaf_contrasts <- list(
  c_PA_24_vs_0, c_PA_48_vs_0, c_TS_24_vs_0, c_TS_48_vs_0,
  c_PA_vs_TS_0, c_PA_vs_TS_24, c_PA_vs_TS_48,
  c_PA_cd_vs_no, c_TS_cd_vs_no,
  c_all_24_vs_0, c_all_48_vs_0,
  coef_24, coef_48
)
names(leaf_contrasts) <- leaf_names

# indices per comparison (A then B, concatenated), for base mean calcs later
index_leaf <- list(
  leaf_pa121_24h_0h        = c(i_PA24, i_PA0),
  leaf_pa121_48h_0h        = c(i_PA48, i_PA0),
  leaf_tsh660_24h_0h       = c(i_TS24, i_TS0),
  leaf_tsh660_48h_0h       = c(i_TS48, i_TS0),
  leaf_pa121_tsh660_0h_0h  = c(i_PA0,  i_TS0),
  leaf_pa121_tsh660_24h_24h= c(i_PA24, i_TS24),
  leaf_pa121_tsh660_48h_48h= c(i_PA48, i_TS48),
  leaf_pa121_cd_no_cd      = c(intersect(i_cd, idx("PA121", c("24","48"))), i_PA0),
  leaf_tsh660_cd_no_cd     = c(intersect(i_cd, idx("TSH660", c("24","48"))), i_TS0),
  leaf_24h_0h              = c(i_24, i_0),
  leaf_48h_0h              = c(i_48, i_0)
  # note: the two interaction-name contrasts have no A/B index by construction
)

# write compact artifacts (numeric contrasts only)
dir.create(PATHS$tables, recursive = TRUE, showWarnings = FALSE)
num_mat <- do.call(rbind, leaf_contrasts[1:11])  # first 11 are numeric here
rownames(num_mat) <- names(leaf_contrasts)[1:11]
write.table(num_mat,
            file = file.path(PATHS$tables, "leaf_contrasts_matrix.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# index table (A|B library names for readability)
libN <- as.character(dd_leaf.le$library)
leaf_idx_tbl <- data.table::data.table(
  comparison = names(index_leaf),
  A_libs = vapply(index_leaf, function(v) paste(libN[v][seq_len(length(v)/2)], collapse=", "), ""),
  B_libs = vapply(index_leaf, function(v) paste(libN[v][-seq_len(length(v)/2)], collapse=", "), "")
)
data.table::fwrite(leaf_idx_tbl, file.path(PATHS$tables, "leaf_contrasts_index.txt"), sep = "\t")

# save for downstream runners
save(leaf_contrasts, index_leaf, mm_leaf,
     file = file.path(PATHS$interim, "contrasts_leaf.rda"))

cat("[04_contrasts_leaf] wrote contrasts_leaf.rda and tables/leaf_contrasts_*.txt\n")

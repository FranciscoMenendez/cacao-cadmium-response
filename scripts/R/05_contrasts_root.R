# scripts/R/05_contrasts_root.R
# Define root contrasts + library indices.
# Output: data/interim/contrasts_root.rda
# Also writes human-readable tables under results/tables/.

source("scripts/R/00_setup.R")
load(file.path(PATHS$interim, "dd_objects.rda"))      # dd_root.le
load(file.path(PATHS$interim, "expsum_built.rda"))    # expsum_root

mm_root <- model.matrix(object = design(dd_root.le), data = as.data.frame(colData(dd_root.le)))
rownames(mm_root) <- as.character(dd_root.le$library)

mmean <- function(i) colMeans(mm_root[i, , drop = FALSE])
idx   <- function(g=NULL, t=NULL) {
  with(as.data.frame(colData(dd_root.le)),
       which((is.null(g) | genotype %in% g) &
               (is.null(t) | time     %in% t)))
}

# groups (root has 0 & 48)
i_PA0 <- idx("PA121","0");  i_PA48 <- idx("PA121","48")
i_TS0 <- idx("TSH660","0"); i_TS48 <- idx("TSH660","48")
i_0   <- idx(t="0");        i_48   <- idx(t="48")

# contrasts (all numeric for roots)
c_PA_vs_TS_0   <- mmean(i_PA0)  - mmean(i_TS0)
c_PA_vs_TS_48  <- mmean(i_PA48) - mmean(i_TS48)
c_PA_48_vs_0   <- mmean(i_PA48) - mmean(i_PA0)
c_TS_48_vs_0   <- mmean(i_TS48) - mmean(i_TS0)
c_all_48_vs_0  <- mmean(i_48)   - mmean(i_0)

root_names <- c(
  "root_pa121_tsh660_0h_0h",
  "root_pa121_tsh660_48h_48h",
  "root_pa121_48h_0h",
  "root_tsh660_48h_0h",
  "root_48h_0h"
)

root_contrasts <- list(
  c_PA_vs_TS_0, c_PA_vs_TS_48, c_PA_48_vs_0, c_TS_48_vs_0, c_all_48_vs_0
)
names(root_contrasts) <- root_names

index_root <- list(
  root_pa121_tsh660_0h_0h  = c(i_PA0,  i_TS0),
  root_pa121_tsh660_48h_48h= c(i_PA48, i_TS48),
  root_pa121_48h_0h        = c(i_PA48, i_PA0),
  root_tsh660_48h_0h       = c(i_TS48, i_TS0),
  root_48h_0h              = c(i_48,  i_0)
)

# write compact artifacts
dir.create(PATHS$tables, recursive = TRUE, showWarnings = FALSE)
num_mat <- do.call(rbind, root_contrasts)
rownames(num_mat) <- names(root_contrasts)
write.table(num_mat,
            file = file.path(PATHS$tables, "root_contrasts_matrix.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

libN <- as.character(dd_root.le$library)
root_idx_tbl <- data.table::data.table(
  comparison = names(index_root),
  A_libs = vapply(index_root, function(v) paste(libN[v][seq_len(length(v)/2)], collapse=", "), ""),
  B_libs = vapply(index_root, function(v) paste(libN[v][-seq_len(length(v)/2)], collapse=", "), "")
)
data.table::fwrite(root_idx_tbl, file.path(PATHS$tables, "root_contrasts_index.txt"), sep = "\t")

save(root_contrasts, index_root, mm_root,
     file = file.path(PATHS$interim, "contrasts_root.rda"))

cat("[05_contrasts_root] wrote contrasts_root.rda and tables/root_contrasts_*.txt\n")

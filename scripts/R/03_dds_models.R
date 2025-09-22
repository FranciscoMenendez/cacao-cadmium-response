# scripts/R/03_dds_models.R
# Build DESeq2 datasets for leaf and root 
# No normalization or DE is done here; those happen in later steps.

source("scripts/R/00_setup.R")  # loads data.table, DESeq2, PATHS, etc.

# ---- Load inputs from step 01 (expsum & counts) ----
load(file.path(PATHS$interim, "expsum_built.rda"))      # expsum, expsum_leaf, expsum_root
if (file.exists(file.path(PATHS$interim, "count_leaf_built.rda")))
  load(file.path(PATHS$interim, "count_leaf_built.rda"))  # count_leaf
if (file.exists(file.path(PATHS$interim, "count_root_built.rda")))
  load(file.path(PATHS$interim, "count_root_built.rda")) # count_root

# ---- Factor levels (explicit/ref levels) ----
expsum$genotype <- factor(expsum$genotype, levels = c("PA121","TSH660"))
expsum$time     <- factor(expsum$time,     levels = c("0","24","48"))
expsum$block    <- factor(expsum$block,    levels = c("B1","B2","B3"))
expsum$tissue   <- factor(expsum$tissue,   levels = c("leaf","root"))

expsum_leaf$genotype <- factor(expsum_leaf$genotype, levels = c("PA121","TSH660"))
expsum_leaf$time     <- factor(expsum_leaf$time,     levels = c("0","24","48"))
expsum_leaf$block    <- factor(expsum_leaf$block,    levels = c("B1","B2","B3"))

expsum_root$genotype <- factor(expsum_root$genotype, levels = c("PA121","TSH660"))
expsum_root$time     <- factor(expsum_root$time,     levels = c("0","48"))
expsum_root$block    <- factor(expsum_root$block,    levels = c("B1","B2","B3"))

# ---- Minimal consistency checks (no file outputs) ----
stopifnot(identical(colnames(count_leaf),  as.character(expsum_leaf$library)))
stopifnot(identical(colnames(count_root), as.character(expsum_root$library)))
stopifnot(identical(rownames(count_leaf),  rownames(count_root)))  # same gene universe

# ---- Build DESeq2 datasets ----
# Leaf model: ~ genotype + time + block + genotype:time + block:genotype
dd_leaf.le <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_leaf,
  colData   = expsum_leaf,
  design    = ~ genotype + time + block + genotype:time + block:genotype
)

# Root-only model (no leaf libraries): ~ genotype + time + block + genotype:time + genotype:block
dd_root.le <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_root,
  colData   = expsum_root,
  design    = ~ genotype + time + block + genotype:time + genotype:block
)

# Optional: combined (leaf+root) object if you want cross-tissue designs downstream.
# It’s created but we won’t use it unless needed later.
make_dd_all <- TRUE
if (make_dd_all) {
  # The above is verbose; simpler way:
  counts_all <- cbind(count_leaf, count_root)[, libs_all]
  stopifnot(identical(colnames(counts_all), libs_all))
  
  dd_all.le <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_all,
    colData   = expsum,
    design    = ~ tissue + genotype + time + block + genotype:tissue + genotype:time + genotype:block
  )
} else {
  dd_all.le <- NULL
}

# ---- Filter low-count genes (same simple filter you used) ----
keep_leaf     <- rowSums(counts(dd_leaf.le))     >= 10L
keep_root  <- rowSums(counts(dd_root.le))  >= 10L
dd_leaf.le    <- dd_leaf.le[keep_leaf, ]
dd_root.le <- dd_root.le[keep_root, ]

if (!is.null(dd_all.le)) {
  keep_all  <- rowSums(counts(dd_all.le)) >= 10L
  dd_all.le <- dd_all.le[keep_all, ]
}

# ---- Save objects for downstream steps (no size factors estimated here) ----
save(dd_leaf.le, dd_root.le, dd_all.le,
     file = file.path(PATHS$interim, "dd_objects.rda"))

# ---- Export model matrices (human-readable, versionable) ----
dir.create(file.path(PATHS$tables), recursive = TRUE, showWarnings = FALSE)

mm_leaf <- as.data.frame(model.matrix(object = design(dd_leaf.le),    data = as.data.frame(colData(dd_leaf.le))))
mm_root <- as.data.frame(model.matrix(object = design(dd_root.le), data = as.data.frame(colData(dd_root.le))))

write.table(mm_leaf,   file = file.path(PATHS$tables, "model_matrix_leaf.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(mm_root, file = file.path(PATHS$tables, "model_matrix_root.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

cat("[03_dds_models] Built dd_leaf.le and dd_root.le; saved to data/interim/dd_objects.rda\n")

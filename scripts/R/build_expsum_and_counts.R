# Build experimental summary tables (expsum_*) and counts matrices (count_leaf, count_root1)
# Robust to naming quirks; avoids attach(); aligns counts columns to expsum row order
# Author: Francisco + ChatGPT
# Date: 2025-09-19

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# -----------------------------
# 1) Read raw counts
# -----------------------------
# Expect a tab-delimited file where the first column is gene IDs and the remaining are sample libraries
# Prefer fread() for speed; it will keep Geneid as character and counts as integer/numeric
counts_dt <- fread("data/raw/counts_genes_last_exon_new.txt")
stopifnot(names(counts_dt)[1] %in% c("Geneid","geneid","GeneID"))
setnames(counts_dt, names(counts_dt)[1], "geneid")

# Keep only the true library columns (numeric). This is safer than fixed column indices.
lib_cols <- grep("\\.bam$", names(counts_dt), value = TRUE, ignore.case = TRUE)

if (length(lib_cols) == 0) stop("No numeric library columns found. Check input file.")

# Optional: drop any T8 libraries if present (we'll also filter later by time)
# lib_cols <- lib_cols[!grepl("T8", lib_cols, ignore.case = TRUE)]

# -----------------------------
# 2) Melt to long and parse metadata from library names
# -----------------------------
ctmel <- melt(counts_dt[, c("geneid", lib_cols), with = FALSE],
              id.vars = "geneid", variable.name = "library", value.name = "count")

# Parse genotype, time, block, tissue from the library string
ctmel[, genotype := fifelse(grepl("TSH", library, ignore.case = TRUE), "TSH660",
                            fifelse(grepl("PA", library,  ignore.case = TRUE), "PA121", NA_character_))]
# Time: 
ctmel[, time := {
  m <- str_match(library, "T(0|24|48|8)(?:[^0-9]|$)")
  # m[,2] grabs the captured group; keep as character ('0','24','48','8')
  m[, 2]
}]

# Block: vector output; case-insensitive; returns e.g. "B1","B2","B3"
ctmel[, block := toupper(str_extract(as.character(library), "(?i)B[123]"))]
ctmel[, block := factor(block, levels = c("B1","B2","B3"))]
stopifnot(nrow(ctmel) == length(ctmel$block))
table(ctmel$block, useNA = "ifany") # all blocks have the same number of genes (qc step)

# Tissue
ctmel[, tissue := fifelse(grepl("raiz", library, ignore.case = TRUE), "root",
                          fifelse(grepl("hoja", library, ignore.case = TRUE), "leaf", NA_character_))] # since we're splitting later, no need to convert to factor

# Filter to the timepoints of interest
ctmel <- ctmel[time %in% c("0","24","48")]
table(ctmel$block, useNA = "ifany") # 24h libraries are fewer

# -----------------------------
# 3) Build experimental summary table (expsum) (one row per library)
# -----------------------------
expsum <- unique(ctmel[, .(library, genotype, time, block, tissue)])
setorder(expsum, tissue, genotype, block, time, library)

# Sanity: each library should have exactly one row
stopifnot(!any(duplicated(expsum$library)))

# Coerce to factors with explicit levels; set rownames = library for DESeq2 compatibility
expsum[, `:=`(
  genotype = factor(genotype, levels = c("PA121","TSH660")),
  time     = factor(time, levels = c("0","24","48")),
  block    = factor(block, levels = c("B1","B2","B3")),
  tissue   = factor(tissue, levels = c("leaf","root"))
)]

# Add plant ID (cell) = block × genotype (== actual plant here: 1 plant per cell)
expsum[, plant := interaction(block, genotype, drop = TRUE)]

# -----------------------------
# 4) QC summaries (helpful to eyeball design correctness)
# -----------------------------
# Counts of libraries per cell and time
qc_cell_time <- expsum[, .N, by = .(tissue, block, genotype, time)][order(tissue, block, genotype, time)]
# Expectation with your design: leaf has 3 timepoints per cell; root has 2
qc_cell_tot  <- expsum[, .N, by = .(tissue, block, genotype)][order(tissue, block, genotype)]

# Print quick summaries
print(qc_cell_time)
print(qc_cell_tot)

# Check for any NA in parsed fields
if (anyNA(expsum)) {
  warning("Some metadata fields are NA after parsing. Inspect rows:")
  print(expsum[!complete.cases(expsum)])
}

# -----------------------------
# 5) Split expsum by tissue with correct time ranges
# -----------------------------
expsum_leaf  <- expsum[tissue == "leaf"] # 0,24,48 by design
expsum_root <- expsum[tissue == "root" ] # roots only at 0 & 48

# -----------------------------
# 6) Build counts matrices aligned to expsum order
# -----------------------------
# Ensure column order of counts matches row order of expsum_* exactly

# Rebuild a wide counts matrix holding only the libraries we actually keep
counts_mat <- as.matrix(counts_dt[, ..lib_cols])
rownames(counts_mat) <- counts_dt$geneid

# Leaf
if (nrow(expsum_leaf)) {
  missing_leaf <- setdiff(expsum_leaf$library, colnames(counts_mat))
  if (length(missing_leaf)) stop("Leaf libraries not found in counts: ", paste(missing_leaf, collapse=","))
  count_leaf <- counts_mat[, expsum_leaf$library, drop = FALSE]
}


# Root (root-only libraries at 0/48)
if (nrow(expsum_root)) {
  missing_root <- setdiff(expsum_root$library, colnames(counts_mat))
  if (length(missing_root)) stop("Root libraries not found in counts: ", paste(missing_root, collapse=","))
  count_root <- counts_mat[, expsum_root$library, drop = FALSE]
}

# Count matrix check
stopifnot(!anyDuplicated(expsum_leaf$library))    # no dup library IDs
stopifnot(!anyDuplicated(expsum_root$library))    # no dup library IDs
stopifnot(all(rownames(counts_mat) != ""))        # gene IDs set

# -----------------------------
# 7) Additional QC: library size and duplicates
# -----------------------------
libsize_leaf <- if (exists("count_leaf")) data.table(library = colnames(count_leaf), libsize = colSums(count_leaf), detected_genes= colSums(count_leaf > 0)) else NULL
libsize_root <- if (exists("count_root")) data.table(library = colnames(count_root), libsize = colSums(count_root), detected_genes= colSums(count_root > 0)) else NULL

# Attach metadata for readability (tissue, genotype, time, block)
exmeta <- as.data.table(expsum)[, .(library, tissue, genotype, time, block)]
libsize_leaf <- exmeta[libsize_leaf, on = "library"]  # left join by library
libsize_root <- exmeta[libsize_root, on = "library"]  # left join by library
rm(exmeta)


# Check uniqueness of plant per (block×genotype) cell
plant_check <- expsum[, .N, by = .(plant, tissue)][order(plant, tissue)]
print(plant_check)

# -----------------------------
# 8) Save artifacts
# -----------------------------
# RDA (binary artifacts; ignored by git per your .gitignore)
save(expsum, expsum_leaf, expsum_root1, file = "data/interim/expsum_built.rda")
if (exists("count_leaf"))  save(count_leaf,  file = "data/interim/count_leaf_built.rda")
if (exists("count_root")) save(count_root1, file = "data/interim/count_root_built.rda")
save(libsize_leaf, libsize_root, file = 'data/interim/libsize_without_sizefactor.rda') # (add size factor column later, save as result)

# sample-sheet CSVs (version these)
data.table::fwrite(expsum,       "metadata/expsum_all.csv")
data.table::fwrite(expsum_leaf,  "metadata/expsum_leaf.csv")
data.table::fwrite(expsum_root1, "metadata/expsum_root.csv")


message("expsum_* and count_* built successfully. Columns in counts are aligned to expsum row order.")

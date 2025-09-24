# scripts/R/06_deseq_leaf_run.R
# Run DE for leaf contrasts.
# - Numeric contrasts use per-contrast normalization (subset on the fly).
# - Interaction coefficients (optional) come from a single global fit.
# Outputs go under results/res_leaf/* and data/interim/df_leaf_geno.rda

source("scripts/R/00_setup.R")

# ---- dependencies (your helpers) ----
source("scripts/R/lib/filtered_results.R")  # filtered_results(dds, name_contrast, lib_index, contr)
source("scripts/R/lib/dds2df.R")            # dds2df(res, name)
source("scripts/R/lib/baseMean_add.R")      # baseMean_add(res_df, des, index)
source("scripts/R/lib/make_vol.R")          # make.vol(res, res_name, out, file_path, file_name)
source("scripts/R/lib/heatmaps.R")         # write_go_en(), write.gsealist(), make_idx() if needed

# ---- inputs from earlier steps ----
load(file.path(PATHS$interim, "dd_objects.rda"))         # dd_leaf.le
load(file.path(PATHS$interim, "contrasts_leaf.rda"))     # leaf_contrasts, index_leaf, mm_leaf
gffb <- NULL
if (file.exists(file.path(PATHS$interim, "gff_b97.rda"))) {
  load(file.path(PATHS$interim, "gff_b97.rda"))           # goterm (id, annot)
}
id_map <- NULL
if (file.exists(file.path(PATHS$interim, "id_map_ncbi.rds"))) {
  id_map <- readRDS(file.path(PATHS$interim, "id_map_ncbi.rds"))  # named char vec: geneid -> ncbi
}

# ---- outputs (dirs) ----
res_dir  <- file.path(PATHS$results, "res_leaf")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(res_dir, "res_leaf_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(res_dir, "res_leaf_volcano"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(res_dir, "GOTermE"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(res_dir, "GSEA_lists"), recursive = TRUE, showWarnings = FALSE)

# ---- split contrasts: numeric vs. coefficient-name (interaction) ----
comp_names <- names(leaf_contrasts)
is_numeric <- vapply(leaf_contrasts, is.numeric, logical(1))
num_names  <- comp_names[is_numeric]
coef_names <- comp_names[!is_numeric]  # typically c("leaf_tsh660.24h","leaf_tsh660.48h")

# ---- run numeric contrasts with per-contrast normalization ----
say("Running %d numeric leaf contrasts with per-contrast normalization…", length(num_names))

res_num <- mapply(
  FUN = function(cn, v) {
    filtered_results(
      dds          = dd_leaf.le,
      name_contrast= leaf_contrasts[[cn]],  # numeric vector over model matrix columns
      lib_index    = index_leaf[[cn]],      # concatenated indices: c(A,B)
      contr        = TRUE                   # numeric contrast, not a coef name
    )
  },
  cn = num_names,
  v  = num_names,
  SIMPLIFY = FALSE
)

# ---- optional: run interaction coefficients from a single global fit ----
RUN_INTERACTIONS <- length(coef_names) > 0
res_coef <- list()
if (RUN_INTERACTIONS) {
  say("Fitting single global model on all leaf libraries for interaction coefficients…")
  dds_full <- DESeq2::DESeq(dd_leaf.le)  # one fit; uses all leaf libs
  for (nm in coef_names) {
    coef_id <- leaf_contrasts[[nm]]      # e.g., "genotypeTSH660.time24"
    res_coef[[nm]] <- DESeq2::results(dds_full, name = coef_id)
  }
}

# ---- volcano plots (numeric + coef); tables; GO; GSEA ----
# format -> df via your dds2df(); annotate; baseMeanA/B for numeric
df_leaf_geno <- list()

# helper to add annotations consistently
.annotate_df <- function(dfx) {
  if (!is.null(gffb) && all(c("id","annot") %in% names(gffb))) {
    dfx$annot2 <- gffb$annot[ match(dfx$geneid, goterm$id) ]
  }
  if (!is.null(id_map)) {
    dfx$ncbi <- unname(id_map[ dfx$geneid ])
  }
  dfx
}

# numeric contrasts
for (nm in num_names) {
  res_obj <- res_num[[nm]]
  # table/plot names
  base_name <- nm
  
  # volcano
  make.vol(res_obj, base_name, out = TRUE,
           file_path = file.path(res_dir, "res_leaf_volcano"),
           file_name = paste0(base_name, "_volcano"))
  
  # tidy df
  dfx <- dds2df(res_obj, base_name)
  dfx <- .annotate_df(dfx)
  
  # add BaseMeanA/B using your index (A first half, B second half)
  dfx <- baseMean_add(dfx, des = dd_leaf.le, index = index_leaf[[nm]])
  
  # write table
  data.table::fwrite(dfx,
                     file = file.path(res_dir, "res_leaf_tables", paste0(base_name, "_results_table.txt")),
                     sep = "\t", quote = FALSE
  )
  
  # GO up/down lists + GSEA ranking
  write_go_en(dfx, base_name, directory = file.path(res_dir, "GOTermE"), obj = dd_leaf.le)
  write.gsealist(dfx[order(dfx$log2FoldChange, decreasing = TRUE), ],
                 file_name = file.path(res_dir, "GSEA_lists", paste0(base_name, "_gsea_list.txt")))
  
  df_leaf_geno[[nm]] <- dfx
}

# interaction coefficients (if any) — no per-contrast baseMeans
if (RUN_INTERACTIONS) {
  for (nm in coef_names) {
    res_obj <- res_coef[[nm]]
    base_name <- nm
    
    make.vol(res_obj, base_name, out = TRUE,
             file_path = file.path(res_dir, "res_leaf_volcano"),
             file_name = paste0(base_name, "_volcano"))
    
    dfx <- dds2df(res_obj, base_name)
    dfx <- .annotate_df(dfx)
    
    data.table::fwrite(dfx,
                       file = file.path(res_dir, "res_leaf_tables", paste0(base_name, "_results_table.txt")),
                       sep = "\t", quote = FALSE
    )
    
    write_go_en(dfx, base_name, directory = file.path(res_dir, "GOTermE"), obj = dd_leaf.le)
    write.gsealist(dfx[order(dfx$log2FoldChange, decreasing = TRUE), ],
                   file_name = file.path(res_dir, "GSEA_lists", paste0(base_name, "_gsea_list.txt")))
    
    df_leaf_geno[[nm]] <- dfx
  }
}

# ---- save combined object for downstream merges/figures ----
save(df_leaf_geno, file = file.path(PATHS$interim, "df_leaf_geno.rda"))

say("Leaf DE complete. Wrote %d numeric contrasts%s.",
    length(num_names),
    if (RUN_INTERACTIONS) sprintf(" and %d interaction coefficients", length(coef_names)) else "")

# ---- diagnostics: one row per contrast with major stats + size factors ----

# helper: size factors for a subset of libraries (per-contrast normalization)
.sf_for <- function(dds, libs) {
  ds <- dds[, libs]
  ds <- DESeq2::estimateSizeFactors(ds)
  data.frame(library = colnames(ds),
             sizeFactor = as.numeric(sizeFactors(ds)),
             stringsAsFactors = FALSE)
}

diag_rows <- list()
libN <- as.character(dd_leaf.le$library)

# numeric contrasts (per-contrast size factors)
for (nm in num_names) {
  dfx <- df_leaf_geno[[nm]]             # already created above
  idx_vec <- index_leaf[[nm]]
  half <- length(idx_vec) / 2L
  libsA <- libN[ idx_vec[seq_len(half)] ]
  libsB <- libN[ idx_vec[(half + 1L):length(idx_vec)] ]
  sft   <- .sf_for(dd_leaf.le, c(libsA, libsB))
  
  sfa <- subset(sft, library %in% libsA)
  sfb <- subset(sft, library %in% libsB)
  
  diag_rows[[nm]] <- data.table::data.table(
    contrast         = nm,
    n_tested         = sum(!is.na(dfx$padj)),
    n_DE_padj_lt_0.05= sum(dfx$padj < 0.05, na.rm = TRUE),
    n_up             = sum(dfx$padj < 0.05 & dfx$log2FoldChange > 0, na.rm = TRUE),
    n_down           = sum(dfx$padj < 0.05 & dfx$log2FoldChange < 0, na.rm = TRUE),
    median_LFC       = stats::median(dfx$log2FoldChange, na.rm = TRUE),
    mean_LFC         = mean(dfx$log2FoldChange, na.rm = TRUE),
    min_padj         = suppressWarnings(min(dfx$padj, na.rm = TRUE)),
    median_baseMean  = stats::median(dfx$baseMean, na.rm = TRUE),
    libs_A           = paste(libsA, collapse = ","),
    libs_B           = paste(libsB, collapse = ","),
    sfA_median       = stats::median(sfa$sizeFactor),
    sfA_min          = min(sfa$sizeFactor),
    sfA_max          = max(sfa$sizeFactor),
    sfB_median       = stats::median(sfb$sizeFactor),
    sfB_min          = min(sfb$sizeFactor),
    sfB_max          = max(sfb$sizeFactor),
    sfA_detail       = paste(sprintf("%s:%.4f", sfa$library, sfa$sizeFactor), collapse = ";"),
    sfB_detail       = paste(sprintf("%s:%.4f", sfb$library, sfb$sizeFactor), collapse = ";")
  )
}

# interaction coefficients (global fit size factors)
if (RUN_INTERACTIONS) {
  # size factors from the single global leaf fit
  sf_full <- data.frame(library = colnames(dds_full),
                        sizeFactor = as.numeric(sizeFactors(dds_full)),
                        stringsAsFactors = FALSE)
  sf_full_str <- paste(sprintf("%s:%.4f", sf_full$library, sf_full$sizeFactor), collapse = ";")
  
  for (nm in coef_names) {
    dfx <- df_leaf_geno[[nm]]
    diag_rows[[nm]] <- data.table::data.table(
      contrast         = nm,
      n_tested         = sum(!is.na(dfx$padj)),
      n_DE_padj_lt_0.05= sum(dfx$padj < 0.05, na.rm = TRUE),
      n_up             = sum(dfx$padj < 0.05 & dfx$log2FoldChange > 0, na.rm = TRUE),
      n_down           = sum(dfx$padj < 0.05 & dfx$log2FoldChange < 0, na.rm = TRUE),
      median_LFC       = stats::median(dfx$log2FoldChange, na.rm = TRUE),
      mean_LFC         = mean(dfx$log2FoldChange, na.rm = TRUE),
      min_padj         = suppressWarnings(min(dfx$padj, na.rm = TRUE)),
      median_baseMean  = stats::median(dfx$baseMean, na.rm = TRUE),
      libs_A           = NA_character_,   # not applicable
      libs_B           = NA_character_,   # not applicable
      sfA_median       = stats::median(sf_full$sizeFactor),
      sfA_min          = min(sf_full$sizeFactor),
      sfA_max          = max(sf_full$sizeFactor),
      sfB_median       = NA_real_,        # not split A/B for interactions
      sfB_min          = NA_real_,
      sfB_max          = NA_real_,
      sfA_detail       = sf_full_str,
      sfB_detail       = NA_character_
    )
  }
}

# bind, write, and save
diag_tbl <- data.table::rbindlist(diag_rows, use.names = TRUE, fill = TRUE)
dir.create(PATHS$tables, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(diag_tbl,
                   file = file.path(PATHS$tables, "leaf_contrasts_diagnostics.tsv"),
                   sep = "\t", quote = FALSE
)
saveRDS(diag_tbl, file = file.path(PATHS$interim, "leaf_contrasts_diagnostics.rds"))
say("Wrote diagnostics: results/tables/leaf_contrasts_diagnostics.tsv")

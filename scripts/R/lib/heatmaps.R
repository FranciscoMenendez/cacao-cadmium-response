# scripts/R/lib/heatmaps.R

# ---- Sample distance heatmap (pheatmap) ----
# Filters genes detected in >= min_detected libraries (within the selected libs),
# applies VST, computes sample-sample distances, and draws a distance heatmap.
make.heatmap <- function(dds,
                         lib_index = NULL,
                         title = "",
                         min_detected = 2,
                         cluster_cols = FALSE,
                         n_colors = 100) {
  stopifnot(methods::is(dds, "DESeqDataSet"))
  
  if (!is.null(lib_index)) dds <- dds[, lib_index]
  
  dds <- DESeq2::estimateSizeFactors(dds)
  keep <- rowSums(DESeq2::counts(dds, normalized = TRUE) > 0) >= min_detected
  if (!any(keep)) stop("No genes pass filter for heatmap.")
  dds <- dds[keep, ]
  
  vsd <- DESeq2::vst(dds, blind = TRUE)
  smat <- as.matrix(stats::dist(t(SummarizedExperiment::assay(vsd))))
  
  # Default labels: genotype_time_block_tissue
  cd <- as.data.frame(SummarizedExperiment::colData(vsd))
  lbl_cols <- intersect(c("genotype","time","block","tissue"), colnames(cd))
  labs <- do.call(paste, c(cd[, lbl_cols, drop = FALSE], sep = "_"))
  rownames(smat) <- labs
  colnames(smat) <- labs
  
  pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"))(n_colors)
  
  pheatmap::pheatmap(
    smat,
    main = title,
    cluster_cols = cluster_cols,
    clustering_distance_rows = stats::as.dist(smat),
    clustering_distance_cols = stats::as.dist(smat),
    color = pal
  )
  
  invisible(smat)
}

# ---- Z-score heatmap of counts matrix (gplots::heatmap.2) ----
# Accepts a numeric matrix (rows=genes, cols=libraries). Rows are z-scored,
# top ngene by variance are plotted (or 'named_genes' subset). Optional
# RowSideColors from a data.frame mapping gene_id -> Category.
make_zscr_heatmap <- function(data,
                              file_name,
                              ngene = 150,
                              dendo = "row",          # "row", "column", or "both"
                              row_ordered = TRUE,     # if FALSE, rows not clustered
                              colsidecolors = NULL,   # kept for API compatibility (unused here)
                              colsp = NULL,
                              rowsp = NULL,
                              named_genes = NULL,
                              brks = NULL,
                              named_cat = NULL,       # data.frame with cols: gene_id, Category
                              ncat = 3,               # min palette size
                              HEIGHT = 26,
                              WIDTH = 12,
                              MARGINS = c(14, 12)) {
  
  # Ensure matrix
  mat <- as.matrix(data)
  storage.mode(mat) <- "numeric"
  
  # Tiny jitter to avoid zero-variance hiccups
  mat <- jitter(mat, factor = 1, amount = 1e-5)
  
  # Row-wise z-score
  z <- t(scale(t(mat)))
  
  # Order rows by variance (DESC) and select
  vars <- matrixStats::rowVars(z)
  ord  <- order(vars, decreasing = TRUE)
  
  if (is.null(named_genes)) {
    z <- z[ord, , drop = FALSE]
    if (nrow(z) > ngene) z <- z[seq_len(ngene), , drop = FALSE]
  } else {
    z <- z[rownames(z) %in% named_genes, , drop = FALSE]
    # keep current order of 'z'; if you want 'named_genes' order, reindex explicitly
    # z <- z[intersect(named_genes, rownames(z)), , drop = FALSE]
  }
  
  # Optional RowSideColors from categories matched by rownames
  row_side <- NULL
  if (!is.null(named_cat)) {
    stopifnot(all(c("gene_id","Category") %in% colnames(named_cat)))
    cat_map <- setNames(as.character(named_cat$Category), named_cat$gene_id)
    f <- factor(cat_map[rownames(z)])
    nlev <- max(ncat, nlevels(f))
    pal  <- RColorBrewer::brewer.pal(nlev, "Dark2")
    row_side <- pal[as.integer(f)]
  }
  
  # Draw PDF
  grDevices::pdf(file_name, width = WIDTH, height = HEIGHT)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  gplots::heatmap.2(
    z,
    col = gplots::greenred,               # palette function OK here
    density.info = "none",
    trace = "none",
    dendrogram = dendo,
    Rowv = if (isTRUE(row_ordered)) NULL else FALSE,
    Colv = NULL,                          # no column clustering by default (as in your original)
    cexRow = 0.6, cexCol = 0.6,
    margins = MARGINS,
    lhei = c(.45, 5),
    colsep = colsp,
    rowsep = rowsp,
    breaks = brks,
    RowSideColors = row_side,
    srtCol = 45
  )
  
  if (!is.null(row_side)) {
    # Simple legend (top-right)
    fshow <- factor(cat_map[rownames(z)])
    levs  <- levels(fshow)
    nlev  <- length(levs)
    pal   <- RColorBrewer::brewer.pal(max(ncat, nlev), "Dark2")[seq_len(nlev)]
    legend("topright", legend = levs, col = pal, pch = 15, bty = "n", cex = 0.8)
  }
  
  invisible(z)
}

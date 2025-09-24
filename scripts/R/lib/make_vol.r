# scripts/R/lib/make_vol.R

make.vol <- function(res,
                     res_name = "noname",
                     out = FALSE,
                     file_name = NULL,
                     file_path = NULL,
                     file_width = 7,
                     file_height = 7,
                     file_units = "in",
                     file_dpi = 300,
                     fdr = 0.05,
                     label_top = 10,
                     label_by = c("padj", "absLFC"),
                     label_sig_only = TRUE,
                     point_size = 0.4,
                     alpha = 0.5,
                     colors = c(sig = "#18A558", ns = "#050A30"),
                     xlim_sym = TRUE,
                     xlim_quantile = 0.995,
                     ymax = NULL,
                     label_col = NULL,   # e.g., "ncbi" to label with NCBI IDs
                     verbose = TRUE) {
  
  # Coerce DESeqResults / DataFrame -> data.frame
  df <- as.data.frame(res, stringsAsFactors = FALSE)
  if (!all(c("log2FoldChange", "padj") %in% names(df))) {
    stop("make.vol: 'res' must have columns log2FoldChange and padj")
  }
  
  # Safe padj for plotting: replace NA with 1, clamp lower bound to avoid Inf
  padj_safe <- df$padj
  padj_safe[is.na(padj_safe)] <- 1
  padj_safe <- pmax(padj_safe, .Machine$double.xmin)
  
  # Significance flag
  sig_flag <- ifelse(df$padj < fdr, "FDR<0.05", "Not Sig")
  
  # Choose label metric
  label_by <- match.arg(label_by)
  metric <- if (label_by == "padj") padj_safe else abs(df$log2FoldChange)
  
  # Who to label?
  lab_pool <- seq_len(nrow(df))
  if (isTRUE(label_sig_only)) lab_pool <- lab_pool[which(df$padj < fdr)]
  lab_ord <- order(metric[lab_pool], decreasing = (label_by == "absLFC"))
  lab_idx <- head(lab_pool[lab_ord], label_top)
  
  # Label text: prefer a provided column, else rownames
  txt <- if (!is.null(label_col) && label_col %in% names(df)) df[[label_col]] else rownames(df)
  lab_df <- df[lab_idx, , drop = FALSE]
  lab_txt <- txt[lab_idx]
  
  # Symmetric x-limits (quantile-based) if requested
  xlims <- NULL
  if (isTRUE(xlim_sym) && any(is.finite(df$log2FoldChange))) {
    q <- stats::quantile(abs(df$log2FoldChange), probs = xlim_quantile, na.rm = TRUE)
    if (is.finite(q) && q > 0) xlims <- c(-q, q)
  }
  
  # y-limits
  yvals <- -log10(padj_safe)
  if (!is.null(ymax)) ylims <- c(0, ymax) else ylims <- NULL
  
  # Build plot
  plt <- ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange, y = yvals)) +
    ggplot2::geom_point(ggplot2::aes(colour = sig_flag),
                        size = point_size, alpha = alpha) +
    ggplot2::scale_color_manual(values = c("FDR<0.05" = colors["sig"],
                                           "Not Sig"  = colors["ns"]),
                                name = NULL) +
    ggplot2::labs(title = res_name, x = "log2 fold-change", y = expression(-log[10](padj))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(face = "bold"))
  
  if (length(lab_idx) > 0) {
    plt <- plt +
      ggrepel::geom_text_repel(
        data = lab_df,
        ggplot2::aes(x = lab_df$log2FoldChange,
                     y = -log10(padj_safe[lab_idx]),
                     label = lab_txt),
        size = 3, max.overlaps = 50, box.padding = 0.35, point.padding = 0.2
      )
  }
  if (!is.null(xlims)) plt <- plt + ggplot2::coord_cartesian(xlim = xlims)
  if (!is.null(ylims)) plt <- plt + ggplot2::coord_cartesian(ylim = ylims)
  
  if (isTRUE(verbose)) {
    n_sig <- sum(df$padj < fdr, na.rm = TRUE)
    message(sprintf(
      "[make.vol] %s • n=%d • FDR<%.2g: %d • labels: %d by %s%s",
      res_name, nrow(df), fdr, n_sig, length(lab_idx), label_by,
      if (label_sig_only) " (sig only)" else ""
    ))
  }
  
  if (!isTRUE(out)) return(plt)
  
  # Save to disk (PNG by default)
  if (is.null(file_name)) file_name <- paste0(res_name, "_volcano.png")
  outfile <- if (!is.null(file_path)) file.path(file_path, file_name) else file_name
  ggplot2::ggsave(filename = outfile, plot = plt, device = "png",
                  width = file_width, height = file_height,
                  units = file_units, dpi = file_dpi)
  invisible(plt)
}

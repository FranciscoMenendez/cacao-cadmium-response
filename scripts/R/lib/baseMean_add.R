# scripts/R/lib/baseMean_add.R

baseMean_add <- function(res,
                         des,
                         index = NULL,          # concatenated: c(A, B)
                         index_a = NULL,        # if provided, use these (unequal A/B ok)
                         index_b = NULL,
                         gene_col = "geneid",
                         use_normalized = TRUE,
                         verbose = FALSE) {
  stopifnot(gene_col %in% names(res))
  
  # Resolve A/B library indices
  if (!is.null(index_a) && !is.null(index_b)) {
    idxA <- as.integer(index_a)
    idxB <- as.integer(index_b)
  } else {
    stopifnot(!is.null(index), length(index) >= 2L, length(index) %% 2L == 0L)
    half <- length(index) / 2L
    idxA <- as.integer(index[seq_len(half)])
    idxB <- as.integer(index[(half + 1L):length(index)])
  }
  libs <- c(idxA, idxB)
  
  # Subset DESeqDataSet to the libraries used in this contrast
  dsub <- des[, libs]
  dsub <- DESeq2::estimateSizeFactors(dsub)
  
  # Pull counts (normalized or raw)
  cm <- DESeq2::counts(dsub, normalized = isTRUE(use_normalized))
  # cm columns are now in order: c(idxA, idxB)
  
  # Align rows of counts to 'res' rows
  gid <- res[[gene_col]]
  m <- match(gid, rownames(cm))           # can be NA if gene missing
  # Build matrices for A and B (preserve row order, fill NAs where missing)
  A <- cm[m, seq_along(idxA), drop = FALSE]
  B <- cm[m, length(idxA) + seq_along(idxB), drop = FALSE]
  
  # Row means (NA if gene not found)
  baseMeanA <- rowMeans(A, na.rm = FALSE)
  baseMeanB <- rowMeans(B, na.rm = FALSE)
  
  # Attach to result frame (don’t reorder)
  res$baseMeanA <- baseMeanA
  res$baseMeanB <- baseMeanB
  # l2(B/A), NA if either side NA or zero
  res$l2B_A <- {
    num <- baseMeanB
    den <- baseMeanA
    out <- rep(NA_real_, length(num))
    ok <- is.finite(num) & is.finite(den) & den > 0
    out[ok] <- log2(num[ok] / den[ok])
    out
  }
  
  if (isTRUE(verbose)) {
    message(sprintf(
      "[baseMean_add] n=%d • libs: A=%d, B=%d • normalized=%s • NA rows=%d",
      nrow(res), ncol(A), ncol(B), use_normalized, sum(is.na(m))
    ))
  }
  
  res
}

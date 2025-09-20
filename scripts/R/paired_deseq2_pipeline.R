library(DESeq2); library(data.table)

# expsum_leaf must have: library, block, genotype, time
expsum_leaf <- as.data.table(expsum_leaf)
expsum_leaf[, block    := factor(block)]
expsum_leaf[, genotype := relevel(factor(genotype), "PA121")]
expsum_leaf[, time     := factor(time, levels = c("0","24","48"))]
# paired blocker = the cell factor (identical to plant here)
expsum_leaf[, cell := interaction(block, genotype, drop = TRUE)]

dds_leaf <- DESeqDataSetFromMatrix(count_leaf, expsum_leaf,
                                   design = ~ cell + time + genotype:time)
dds_leaf <- DESeq(dds_leaf)  # ONE fit

# Sanity: coefficients you’ll use
resultsNames(dds_leaf)
# "time_24_vs_0" "time_48_vs_0" "genotypeTSH660.time24" "genotypeTSH660.time48"

# Within-genotype time effects
res_leaf_pa121_24h_0h  <- results(dds_leaf, name = "time_24_vs_0")
res_leaf_pa121_48h_0h  <- results(dds_leaf, name = "time_48_vs_0")
res_leaf_tsh660_24h_0h <- results(dds_leaf,
                                  contrast = list(c("time_24_vs_0","genotypeTSH660.time24")), listValues = c(1,1))
res_leaf_tsh660_48h_0h <- results(dds_leaf,
                                  contrast = list(c("time_48_vs_0","genotypeTSH660.time48")), listValues = c(1,1))

# Interaction (ΔΔ) terms
res_leaf_tsh660_24h <- results(dds_leaf, name = "genotypeTSH660.time24")
res_leaf_tsh660_48h <- results(dds_leaf, name = "genotypeTSH660.time48")

# Pooled Cd vs no-Cd across genotypes (weighted by sample counts)
cd <- as.data.table(colData(dds_leaf))
ns24 <- cd[time=="24", .N, by=genotype][order(genotype)]$N
ns48 <- cd[time=="48", .N, by=genotype][order(genotype)]$N

res_leaf_cd_no_cd <- results(dds_leaf,
                             contrast = list(
                               c("time_24_vs_0", "time_24_vs_0", "genotypeTSH660.time24",
                                 "time_48_vs_0", "time_48_vs_0", "genotypeTSH660.time48")),
                             listValues = c(ns24[1], ns24[2], ns24[2], ns48[1], ns48[2], ns48[2]) / (sum(ns24)+sum(ns48))
)

# Consistent baseMeans from the same normalized counts
norm <- counts(dds_leaf, normalized = TRUE)

pick <- function(g=NULL,t=NULL){
  idx <- rep(TRUE, nrow(cd))
  if(!is.null(g)) idx <- idx & cd$genotype %in% g
  if(!is.null(t)) idx <- idx & cd$time %in% t
  which(idx)
}

bm <- function(A,B){
  Amean <- rowMeans(norm[,A,drop=FALSE]); Bmean <- rowMeans(norm[,B,drop=FALSE])
  data.table(baseMeanA=Amean, baseMeanB=Bmean, l2B_A = log2((Amean+1e-8)/(Bmean+1e-8)))
}

# Example: attach baseMeans to PA121 24h vs 0h
bm_pa121_24 <- bm(pick("PA121","24"), pick("PA121","0"))
res_pa121_24_dt <- cbind(as.data.table(res_leaf_pa121_24h_0h, keep.rownames="geneid"), bm_pa121_24)

# Build a 'plant' dds and a 'cell' dds; then compare key results:
dds_plant <- DESeqDataSetFromMatrix(count_leaf, within(expsum_leaf, {
  plant <- interaction(block, genotype, drop=TRUE); rm(cell)
}), design = ~ plant + time + genotype:time) |> DESeq()

cmp <- function(r1, r2){
  m1 <- as.data.table(r1, keep.rownames="gene")[,.(gene, lfc1=log2FoldChange, p1=padj)]
  m2 <- as.data.table(r2, keep.rownames="gene")[,.(gene, lfc2=log2FoldChange, p2=padj)]
  merge(m1,m2,by="gene")[, .(
    cor_lfc = cor(lfc1, lfc2, use="complete.obs"),
    max_abs_diff_lfc = max(abs(lfc1-lfc2), na.rm=TRUE)
  )]
}

cmp(results(dds_plant, name="time_24_vs_0"),
    results(dds_leaf,  name="time_24_vs_0"))
# Expect correlation ~1 and tiny max_abs_diff_lfc (floating-point noise)

# scripts/R/02_annotations_build.R
# Build unified annotation map (cr), chip file, and goterm index.

source("scripts/R/00_setup.R")
suppressPackageStartupMessages(library(ape))  # for read.FASTA

# ---- Inputs (paths) ----
FASTA_CDS   <- file.path(PATHS$raw, "Theobroma_cacaoV2_annot_cds.fna")
RBH_cr3     <- file.path(PATHS$raw, "RBBH/rbbh_criollov3")          # CriolloV2 <- CriolloV3
RBH_cr3_rev <- file.path(PATHS$raw, "RBBH/rbbh_CrV2onCrV3")         # CriolloV2 -> CriolloV3
RBH_ar1     <- file.path(PATHS$raw, "RBBH/rbbh_araport11_crv2")
RBH_ar2     <- file.path(PATHS$raw, "RBBH/rbbh_tc_araport11_e4")
RBH_ar3     <- file.path(PATHS$raw, "RBBH/rbbh_araport11_tc_5e2")
RBH_ar4     <- file.path(PATHS$raw, "RBBH/rbbh_tc_araport11_5e2")
RBH_ncbi    <- file.path(PATHS$raw, "RBBH/rbbh_NCBI_TcB97v3_cds")
PHY_crTC    <- file.path(PATHS$raw, "RBBH/phytozome_TC_criolloV2.txt")
PHY_crTC2   <- file.path(PATHS$raw, "RBBH/phytozome_criolloV2_TCphyto.txt")
PHY_TCB97   <- file.path(PATHS$raw, "RBBH/phytozome_TCB97.txt")
GFF3_PATH   <- file.path(PATHS$raw, "TC_B97_consensusV2_rev2021_annotated.gff3")

# ---- Small helpers ----
first_token <- function(x) sub(" .*", "", x)
to_gene_id <- function(x) {
  x <- as.character(x)
  # case 1: transcript form → keep digits after _t, drop optional .version
  x <- sub("(_)[tT](\\d+)(?:\\.\\d+)?$", "\\1g\\2", x, perl = TRUE)
  # case 2: if already in gene form but has a .version → strip version
  x <- sub("(_g\\d+)\\..*$", "\\1", x, perl = TRUE)
  x
}

read_rbh <- function(path) {
  # First 2 columns are query, target
  dt <- data.table::fread(file=path, fill = TRUE, header = FALSE, sep = "\t", quote = "")
  if (ncol(dt) < 2) stop("RBH file malformed: ", path)
  data.table::setnames(dt, 1:2, c("query","target"))
  dt[, .(query, target)]
}

# ---- Build 'cr' mapping ----
say("Reading CDS FASTA…")
cds <- read.FASTA(FASTA_CDS)
lab <- labels(cds); rm(cds)
cr  <- data.table::data.table(
  cr      = lab,
  lookup  = first_token(lab)          # e.g., CriolloV2 CDS ID (with isoform)
)

say("Reading RBH/orthology tables…")
ideq      <- read_rbh(RBH_cr3)        # CriolloV2 <- CriolloV3
ideq2     <- read_rbh(RBH_cr3_rev)    # CriolloV2 -> CriolloV3
araport   <- read_rbh(RBH_ar1)
araport2  <- read_rbh(RBH_ar2)
araport3  <- read_rbh(RBH_ar3)
araport4  <- read_rbh(RBH_ar4)
ncbi <- read_rbh(RBH_ncbi)
data.table::setnames(ncbi, c("query","target"), c("query","target"))  # keep as is
phy1      <- data.table::fread(file = PHY_crTC)   # TC -> CrV2
phy2      <- data.table::fread(file = PHY_crTC2)  # CrV2 -> TCphyto
phyB97    <- data.table::fread(file = PHY_TCB97)

# Normalize column names for phytozome tables (first two are query/target)
for (dt in list(phy1, phy2, phyB97)) {
  if (ncol(dt) >= 2) data.table::setnames(dt, 1:2, c("query","target"))
}

# Map CriolloV2 -> CriolloV3 (two directions available)
cr[, crv3   := ideq [match(lookup, ideq$target),   query]]
cr[, crv3_2 := ideq2[match(lookup, ideq2$query),   target]]
# Gene-level normalization

cr[, `:=`(crv3   = to_gene_id(crv3),
          crv3_2 = to_gene_id(crv3_2))]

# Araport mappings (multiple runs with different thresholds/directions)
cr[, araport  := araport [match(lookup, araport$target),  query]]
cr[, araport2 := araport2[match(lookup, araport2$query),  target]]
cr[, araport3 := araport3[match(lookup, araport3$target), query]]
cr[, araport4 := araport4[match(lookup, araport4$query),  target]]

# NCBI mapping (normalize TC gene ids)
ncbi_dt <- data.table::as.data.table(ncbi)
ncbi_dt[, tc_gen := to_gene_id(target)]
ncbi_dt[, tc_gen := substr(tc_gen, 1, 16)]
cr[, ncbi := ncbi_dt[match(crv3_2, tc_gen), query]]

# Phytozome mappings (two directions around CrV2)
cr[, phytozome  := phy1 [match(lookup, phy1$target),  query]]
cr[, phytozome2 := phy2 [match(lookup, phy2$query),   target]]

# Save 'cr'
data.table::setDT(cr)
save(cr, file = file.path(PATHS$interim, "cr.rda"))
say("Saved: data/interim/cr.rda  (rows: %d)", nrow(cr))

# ---- Write chip file (Probe Set ID, Gene Symbol, Gene Title) ----
chip <- cr[!is.na(crv3_2), .(ProbeSetID = crv3_2, GeneSymbol = lookup, GeneTitle = cr)]
data.table::fwrite(chip, file = file.path(PATHS$tables, "b97.chip"),
                   sep = "\t", quote = FALSE, col.names = TRUE)
say("Wrote chip: results/tables/b97.chip  (rows: %d)", nrow(chip))

# ---- Build annotation for b97 (from GFF3) ----
say("Parsing GFF3 for gene annotations…")
gff <- data.table::fread(GFF3_PATH, fill = TRUE, sep = "\t", header = FALSE, quote = "")
# Keep only gene rows (column 3 == "gene")
gff <- gff[V3 == "gene"]
# Parse ID from attributes (V9); keep the rest collapsed as 'annot'
gff[, id    := sub("^ID=([^;]+).*", "\\1", V9)]
data.table::setnames(gff, "V9", "annot")
gff[, annot := apply(.SD, 1, function(x) paste(x[9:length(x)], collapse = " ")), .SDcols = names(gff)]
gffb <- gff[, .(id, annot)] #gff b97 
save(gffb, file = file.path(PATHS$interim, "gff_b97.rda"))

# optional tiny CSV for browsing
data.table::fwrite(gffb, file = file.path(PATHS$tables, "gff_b97_index.csv"))
say("Saved: data/interim/gff_b97.rda; results/tables/ggff_b97_index.csv")

say("Annotations build complete.")

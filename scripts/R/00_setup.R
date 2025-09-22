# scripts/R/00_setup.R

# 1) Root sanity
if (!file.exists("cacao-cadmium-response.Rproj") && !file.exists("README.md")) {
  stop("Run scripts from the project root (where the .Rproj / README.md live).")
}

# 2) Core packages (now includes DESeq2)
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(DESeq2)  # <- added
})

# 3) here() helper (optional)
if (requireNamespace("here", quietly = TRUE)) {
  HERE <- here::here
  try(here::i_am("README.md"), silent = TRUE)
} else {
  HERE <- function(...) file.path(getwd(), ...)
}

# 4) Canonical paths
PATHS <- list(
  raw       = HERE("data","raw"),
  interim   = HERE("data","interim"),
  processed = HERE("data","processed"),
  meta      = HERE("metadata"),
  results   = HERE("results"),
  tables    = HERE("results","tables"),
  figures   = HERE("results","figures")
)

invisible(lapply(PATHS, dir.create, recursive = TRUE, showWarnings = FALSE))

# 5) Mild options
options(datatable.print.nrows = 50L, width = 120)

# 6) Tiny logger/seed (optional)
say  <- function(...) cat("[cacao-cd] ", sprintf(...), "\n")
set.seed(1337)
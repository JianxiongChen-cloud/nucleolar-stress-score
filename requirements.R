# Main R package dependencies for the nucleolar-stress-score analysis scripts.
# Run this file in R to install missing CRAN and Bioconductor packages.

cran_packages <- c(
  "broom",
  "babelgene",
  "data.table",
  "dplyr",
  "forcats",
  "ggplot2",
  "ggpubr",
  "ggrepel",
  "grid",
  "matrixStats",
  "msigdbr",
  "patchwork",
  "pbapply",
  "purrr",
  "RColorBrewer",
  "readr",
  "readxl",
  "rstatix",
  "scales",
  "Seurat",
  "stringr",
  "survival",
  "survminer",
  "tibble",
  "tidyr",
  "UpSetR",
  "viridis"
)

bioc_packages <- c(
  "AnnotationDbi",
  "biomaRt",
  "clusterProfiler",
  "cmapR",
  "enrichplot",
  "GEOquery",
  "GSVA",
  "GSEABase",
  "hgu133a2.db",
  "hgu133plus2.db",
  "hugene10sttranscriptcluster.db",
  "limma",
  "org.Hs.eg.db",
  "rhdf5",
  "UCSCXenaTools"
)

install_missing_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing)
  }
}

install_missing_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }
}

install_missing_cran(cran_packages)
install_missing_bioc(bioc_packages)

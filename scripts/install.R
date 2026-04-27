#!/usr/bin/env Rscript
## Install all ProteoForge dependencies.
## Usage: Rscript scripts/install.R

message("Installing ProteoForge dependencies...")

# CRAN packages
cran_pkgs <- c(
  "data.table", "arrow", "dplyr", "tidyr", "ggplot2", "ggrepel",
  "patchwork", "openxlsx2", "officer", "flextable", "yaml", "digest",
  "fs", "logger", "jsonlite", "optparse", "pheatmap", "igraph",
  "ggraph", "ComplexUpset", "ggVennDiagram", "ragg", "svglite",
  "shiny", "golem", "promises", "future", "waiter", "shinytest2",
  "DT", "plotly", "covr", "testthat", "knitr", "rmarkdown",
  "pkgdown", "roxygen2", "lintr", "styler", "pak", "renv", "targets"
)

installed <- rownames(installed.packages())
to_install <- setdiff(cran_pkgs, installed)

if (length(to_install) > 0) {
  message("Installing ", length(to_install), " CRAN packages...")
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  message("All CRAN packages already installed.")
}

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
bioc_pkgs <- c(
  "limma", "DEqMS", "QFeatures", "SummarizedExperiment",
  "msqrob2", "MSstats", "sva", "clusterProfiler", "enrichplot",
  "fgsea", "STRINGdb", "OmnipathR", "decoupleR"
)
to_install_bioc <- setdiff(bioc_pkgs, installed)
if (length(to_install_bioc) > 0) {
  message("Installing ", length(to_install_bioc), " Bioconductor packages...")
  BiocManager::install(to_install_bioc, ask = FALSE)
} else {
  message("All Bioconductor packages already installed.")
}

message("\nInstallation complete. Run R CMD check to verify.")

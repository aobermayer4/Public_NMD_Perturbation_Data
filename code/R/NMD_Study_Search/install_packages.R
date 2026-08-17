# Packages needed by the cleaned NMD metadata workflow.
# Run inside the project's renv environment.

cran_packages <- c(
  "data.table", "dplyr", "httr", "jsonlite", "purrr",
  "stringr", "tibble", "tidyr", "writexl", "xml2"
)

bioc_packages <- c(
  "GEOquery",
  "ArrayExpress"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_cran) > 0L) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_bioc) > 0L) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

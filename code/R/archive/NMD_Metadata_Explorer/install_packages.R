cran_packages <- c("shiny", "bslib", "DT", "data.table", "zip")
missing <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All required packages are already installed.")
}

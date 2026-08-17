

install.packages("renv")

renv::init(bare=TRUE)

renv::install(c("BiocManager","data.table","dplyr","jsonlite","purrr","remotes","stringr","tibble","tidyr","matrixStats","abind"))

renv::install(c("Bioc::GEOquery","Bioc::oligo","Bioc::ArrayExpress"))

renv::install(c("Bioc::SparseArray","Bioc::DelayedArray","Bioc::Biobase","Bioc::GenomeInfoDbData",
"Bioc::UCSC.utils","Bioc::GenomeInfoDb","Bioc::SummarizedExperiment","Bioc::limma","Bioc::GEOquery","Bioc::affxparser"
"Bioc::zlibbioc","Bioc::affyio","Bioc::Biostrings","Bioc::oligoClasses","Bioc::preprocessCore","Bioc::oligo","Bioc::ArrayExpress",,))

install.packages(c("BiocManager","data.table","dplyr","jsonlite","purrr","remotes","stringr","tibble","tidyr"))


install.packages("~/Projects/Active/Public_NMD_Perturbation_Data/local/BiocVersion_3.19.1.tar.gz", repos = NULL, type = "source")
install.packages("~/Projects/Active/Public_NMD_Perturbation_Data/local/BiocGenerics_0.58.1.tar.gz", repos = NULL, type = "source")
install.packages("matrixStats")
install.packages("~/Projects/Active/Public_NMD_Perturbation_Data/local/MatrixGenerics_1.24.0.tar.gz", repos = NULL, type = "source")
install.packages("abind")
install.packages("~/Projects/Active/Public_NMD_Perturbation_Data/local/IRanges_2.46.0.tar.gz", repos = NULL, type = "source")
install.packages("~/Projects/Active/Public_NMD_Perturbation_Data/local/S4Arrays_1.12.0.tar.gz", repos = NULL, type = "source")
install.packages("~/Projects/Active/Public_NMD_Perturbation_Data/local/Seqinfo_1.2.0.tar.gz", repos = NULL, type = "source")
install.packages("~/Projects/Active/Public_NMD_Perturbation_Data/local/GenomicRanges_1.64.0.tar.gz", repos = NULL, type = "source")
install.packages("~/Projects/Active/Public_NMD_Perturbation_Data/local/XVector_0.52.0.tar.gz", repos = NULL, type = "source")
renv::install("Bioc::SparseArray")
renv::install("Bioc::DelayedArray")
install.packages("~/Projects/Active/Public_NMD_Perturbation_Data/local/S4Vectors_0.50.1.tar.gz", repos = NULL, type = "source")
renv::install("Bioc::Biobase")
renv::install("Bioc::GenomeInfoDbData", prompt = F)
renv::install("Bioc::UCSC.utils", prompt = F)
renv::install("Bioc::IRanges", prompt = F)
renv::install("Bioc::Seqinfo", prompt = F)
renv::install("Bioc::GenomeInfoDb", prompt = F)
renv::install("Bioc::SummarizedExperiment", prompt = F)
renv::install("Bioc::limma", prompt = F)
renv::install("Bioc::GEOquery", prompt = F)
renv::install("Bioc::affxparser", prompt = F)
renv::install("Bioc::zlibbioc", prompt = F)
renv::install("Bioc::affyio", prompt = F)
renv::install("Bioc::Biostrings", prompt = F)
renv::install("Bioc::oligoClasses", prompt = F)
renv::install("Bioc::preprocessCore", prompt = F)
renv::install("Bioc::oligo", prompt = F)
renv::install("Bioc::ArrayExpress", prompt = F)

renv::install(c("Bioc::BiocVersion","Bioc::BiocGenerics","Bioc::MatrixGenerics","Bioc::IRanges","Bioc::S4Arrays","Bioc::Seqinfo","Bioc::GenomicRanges","Bioc::XVector"))
renv::install(c("bioc::S4Vectors"))
renv::install(c("Bioc::GEOquery","Bioc::oligo","Bioc::ArrayExpress"))

renv::install(c("bioc::SparseArray","bioc::DelayedArray","bioc::Biobase","bioc::GenomeInfoDbData",
"bioc::UCSC.utils","bioc::GenomeInfoDb","bioc::SummarizedExperiment","bioc::limma","bioc::affxparser",
"bioc::zlibbioc","bioc::affyio","bioc::Biostrings","bioc::oligoClasses","bioc::preprocessCore","bioc::oligo"))

BiocManager::install("zlibbioc")

renv::install("zlibbioc", prompt = F)

renv::install(c("Bioc::GEOquery","Bioc::oligo","Bioc::ArrayExpress"))

renv::install("Bioc::SparseArray")

The following package(s) were not installed successfully:
- [affxparser]: package 'affxparser' is not available
- [affyio]: package 'affyio' is not available
- [Biobase]: package 'Biobase' is not available
- [BiocGenerics]: package 'BiocGenerics' is not available
- [Biostrings]: package 'Biostrings' is not available
- [limma]: package 'limma' is not available
- [oligoClasses]: package 'oligoClasses' is not available
- [preprocessCore]: package 'preprocessCore' is not available
- [S4Vectors]: package 'S4Vectors' is not available
- [SummarizedExperiment]: package 'SummarizedExperiment' is not available
- [oligo]: install failed
- [GEOquery]: install failed
- [ArrayExpress]: dependency failed (Biobase, oligo, limma)
You may need to manually download and install these packages.

# oligo ----------------------------------------------------------------------
ERROR: dependencies ‘BiocGenerics’, ‘oligoClasses’, ‘Biobase’, ‘Biostrings’, ‘affyio’, ‘affxparser’, ‘preprocessCore’ are not available for package ‘oligo’
Perhaps try a variation of:
install.packages(c('BiocGenerics', 'oligoClasses', 'Biobase', 'Biostrings', 'affyio', 'affxparser', 'preprocessCore'))
* removing ‘/home/alyssa/Projects/Public_NMD_Perturbation_Data/renv/staging/1/oligo’

# GEOquery -------------------------------------------------------------------
ERROR: dependencies ‘Biobase’, ‘limma’, ‘SummarizedExperiment’, ‘S4Vectors’ are not available for package ‘GEOquery’
Perhaps try a variation of:
install.packages(c('Biobase', 'limma', 'SummarizedExperiment', 'S4Vectors'))
* removing ‘/home/alyssa/Projects/Public_NMD_Perturbation_Data/renv/staging/1/GEOquery’



cat("R version:\n")
print(R.version.string)

cat("\nR home:\n")
print(R.home())

cat("\nLibrary paths:\n")
print(.libPaths())

cat("\nRepositories:\n")
print(getOption("repos"))

cat("\nBiocManager installed:\n")
print(requireNamespace("BiocManager", quietly = TRUE))

if (requireNamespace("BiocManager", quietly = TRUE)) {
  cat("\nBioconductor version:\n")
  print(BiocManager::version())
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}
BiocManager::install(version = "3.23", ask = FALSE)
renv::project()
.libPaths()
renv::settings$bioconductor.version("3.23")
renv::settings$bioconductor.version()
renv::install(c(
  "bioc::GEOquery",
  "bioc::ArrayExpress"
))

package_check <- c(
  "GEOquery",
  "ArrayExpress",
  "Biobase",
  "limma",
  "SummarizedExperiment",
  "S4Vectors",
  "oligo"
)

data.frame(
  package = package_check,
  installed = vapply(
    package_check,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  ),
  version = vapply(
    package_check,
    function(x) {
      if (requireNamespace(x, quietly = TRUE)) {
        as.character(packageVersion(x))
      } else {
        NA_character_
      }
    },
    FUN.VALUE = character(1)
  )
)
library(GEOquery)
library(ArrayExpress)

renv::install("fs")
renv::install("bslib")
renv::install("DT")
renv::install("shiny")

renv::install("writexl")
renv::install("zip")


suppressPackageStartupMessages({
  library(ArrayExpress)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(jsonlite)
})

setwd("C:\\Users\\4474108\\Projects\\Active\\Public_NMD_Perturbation_Data")


source("code/R/arrayexpress_metadata_parser.R")

source("code/R/arrayexpress_one_row_samples_patch.R")

ae_meta_res <- readRDS("C:\\Users\\4474108\\Projects\\Active\\Public_NMD_Perturbation_Data\\data\\ArrayExpress_Metadata_v2\\ArrayExpress_metadata_results.rds")

ae_tables <- ae_combine_results(ae_meta_res)


sdrf <- ae_tables$samples

le

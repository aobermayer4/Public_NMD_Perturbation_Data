
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

#setwd("C:\\Users\\4474108\\Projects\\Active\\Public_NMD_Perturbation_Data")


source("code/R/arrayexpress_metadata_parser_with_one_row_samples.R")

#source("code/R/arrayexpress_one_row_samples_patch.R")

ae_meta_res <- readRDS("data/ArrayExpress_Metadata_v4/ArrayExpress_metadata_results.rds")

ae_tables <- ae_combine_results(ae_meta_res)

writexl::write_xlsx(ae_tables,"data/ArrayExpress_Parsed_Output_20260718.xlsx", format_headers = F)

sdrf <- ae_tables$samples_one_row



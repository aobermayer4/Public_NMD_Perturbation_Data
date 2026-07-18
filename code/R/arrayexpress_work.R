









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

setwd("~/Projects/Active/Public_NMD_Perturbation_Data/")


source("code/R/arrayexpress_metadata_parser.R")
source("code/R/arrayexpress_sdrf_character_patch.R")

study_ids <- c(
  "E-MTAB-9330","E-MTAB-10716","E-MTAB-8461","E-MTAB-13839","E-MTAB-13788",
  "E-MTAB-13789","E-MTAB-13787","E-MTAB-13829","E-MTAB-13829","E-MTAB-13836",
  "E-MTAB-16399","E-MTAB-13949","E-MTAB-14755","E-MTAB-14755","E-MTAB-14725",
  "E-MTAB-14725","E-MTAB-13837"
)
study_ids <- unique(study_ids)



ae_results <- ae_fetch_many(
  accessions = study_ids,
  base_dir = "data/ArrayExpress_Metadata_v2",
  overwrite = FALSE,
  use_api_fallback = TRUE,
  verbose = TRUE
)

ae_tables <- ae_combine_results(ae_results)

ae_write_results(
  ae_tables,
  out_dir = "data/ArrayExpress_Metadata_v2/ArrayExpress_Metadata_Tables"
)

ae_search_metadata(
  ae_tables,
  terms = c(
    "\\bUPF1\\b",
    "\\bSMG1\\b",
    "\\bSMG6\\b",
    "knockdown",
    "siRNA",
    "shRNA",
    "CRISPR",
    "inhibitor",
    "cycloheximide",
    "emetine"
  ),
  scope = "samples"
)


ae_search_metadata(
  ae_tables,
  terms = c(
    "nonsense-mediated decay",
    "\\bNMD\\b",
    "premature termination codon",
    "\\bUPF1\\b"
  ),
  scope = "idf"
)

ae_search_metadata(
  ae_tables,
  terms = c(
    "siRNA",
    "shRNA",
    "knockdown",
    "depletion",
    "CRISPR",
    "transfection",
    "transduction",
    "inhibitor",
    "cycloheximide",
    "emetine"
  ),
  scope = "protocols"
)


ae_search_metadata(
  ae_tables,
  terms = c("\\bUPF1\\b", "knockdown"),
  scope = "samples",
  require_all = TRUE
)

one_study <- ae_results[["E-MTAB-9330"]]


one_study$sdrf_long |>
  count(field_raw, sort = TRUE)


one_study$sdrf_long |>
  filter(
    str_detect(
      field_raw,
      regex(
        "Factor Value|Characteristics|Protocol REF|Parameter Value",
        ignore_case = TRUE
      )
    )
  ) |>
  arrange(Sample_ID, field_raw) |>
  View()

one_study$protocol_refs |>
  select(
    Study_ID,
    Sample_ID,
    field_raw,
    Protocol_REF
  ) |>
  View()

one_study$protocols |>
  View()








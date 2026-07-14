



library(data.table)
library(stringr)
library(GEOquery)
library(dplyr)


setwd("D:/R/NMD_Perturbation_Data")

gse_list <- fread("NMD_Perturbation_GSE_Studies.txt",data.table = F)

source("annotate_gse_ids.R")

GSE_IDs <- unique(toupper(gse_list$GSE_IDs))

gse_metas <- lapply(GSE_IDs, extract_gse_meta, verbose = F)
gse_metas <- setNames(gse_metas,GSE_IDs)
saveRDS(gse_metas,file = "NMD_Perturbation_GSE_Metas_20260712.rds")
results <- lapply(gse_metas, lookup_gse, verbose = F)
res_df <- bind_rows(results, .id = "GSE_ID")

gse_metas_int_cols <- Reduce(intersect,lapply(gse_metas,colnames))

gse_metas_cols <- as.data.frame(table(unname(unlist(lapply(gse_metas,colnames)))))
gse_metas_cols <- gse_metas_cols[order(gse_metas_cols$Freq, decreasing = T),]


source("nmd_geo_metadata_parser.R")

parsed_nmd <- parse_nmd_geo_metadata(gse_metas)

GSE_TopOpts <- c("GSE232185","GSE232333","GSE12928","GSE30499","GSE305669","GSE16856")


parsed_nmd_top <- lapply(parsed_nmd,function(df) {
  return(df[which(df$GSE_ID %in% GSE_TopOpts),])
})

parsed_nmd_samp_top <- parsed_nmd$samples[which(parsed_nmd$samples$GSE_ID %in% GSE_TopOpts),]

write.table(parsed_nmd_samp_top,"NMD_Perturbation_TopStudySamples_20260712.txt", sep = '\t', row.names = F)


# Study-level ranked summary
View(parsed_nmd$studies)

# Sample-level classifications
View(parsed_nmd$samples)

# Exact metadata fields that caused each match
View(parsed_nmd$evidence)

# Inventory of all metadata fields
View(parsed_nmd$fields)


write_nmd_geo_results(
  parsed_nmd,
  out_dir = "NMD_GEO_Metadata_Review_20260712"
)


parsed_nmd$studies |>
  select(
    GSE_ID,
    study_title,
    n_samples,
    n_candidate_perturbation,
    n_candidate_controls,
    nmd_targets,
    perturbation_mechanisms,
    agents,
    organism,
    disease_field_values,
    tissue_field_values,
    study_priority_score,
    review_category
  ) |>
  View()



parsed_nmd$studies |>
  filter(
    review_category ==
      "highest priority: candidate perturbation and control arms"
  )



parsed_nmd$evidence |>
  filter(GSE_ID == "GSE12345") |>
  View()


search_geo_metadata(
  parsed_nmd$long_metadata,
  terms = c("UPF1", "knockdown")
)



search_geo_metadata(
  parsed_nmd$long_metadata,
  terms = c("UPF1", "knockdown"),
  require_all_terms = TRUE
)

search_geo_metadata(
  parsed_nmd$long_metadata,
  terms = c("siRNA", "shRNA", "CRISPR", "inhibitor"),
  fields = c("treatment", "protocol", "description", "characteristics")
)


my_agents <- NMD_AGENT_PATTERNS

my_agents$custom_NMD_compound <- c(
  "\\bExactCompoundName\\b",
  "\\bAlternativeCompoundName\\b"
)

parsed_nmd <- parse_nmd_geo_metadata(
  gse_metas,
  agent_patterns = my_agents
)

gse_metas_v2 <- lapply(
  GSE_IDs,
  extract_gse_meta_v2,
  verbose = TRUE
)

gse_metas_v2 <- setNames(gse_metas_v2, GSE_IDs)

saveRDS(
  gse_metas_v2,
  file = "NMD_Perturbation_GSE_Metas_v2_20260712.rds"
)

parsed_nmd_v2 <- parse_nmd_geo_metadata(gse_metas_v2)











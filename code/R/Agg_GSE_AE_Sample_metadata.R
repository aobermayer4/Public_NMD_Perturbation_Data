# Try to combine GSE and Array Express sample data


source("code/R/arrayexpress_metadata_parser_with_one_row_samples.R")
source("code/R/annotate_gse_ids.R")
source("code/R/nmd_geo_metadata_parser.R")

ae_meta_res <- readRDS("data/ArrayExpress_Metadata_v4/ArrayExpress_metadata_results.rds")

ae_tables <- ae_combine_results(ae_meta_res)


gse_meta_res <- readRDS("code\\R\\NMD_Perturbation_GSE_Metas_20260716.rds")

GSE_IDs <- c("GSE232185","GSE232333","GSE12928","GSE30499",
"GSE305669","GSE16856","GSE289050","GSE24205","GSE37210","GSE5486",
"GSE109143","GSE152436","GSE152435","GSE176197","GSE134059","GSE162699",
"GSE59884","GSE204985","GSE185655","GSE162199","GSE60045","GSE61398")

gse_meta_parsed <- parse_nmd_geo_metadata(gse_meta_res)


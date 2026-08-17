# 08/16/2026
# Test NMD study search cleaned workflow

source("code/R/NMD_Study_Search/run_nmd_workflow.R")


result <- run_nmd_workflow(
  discover = FALSE,
  geo_ids = c("GSE305669"),
  arrayexpress_ids = c("E-MTAB-9330")
)

result <- run_nmd_workflow(
  discover = FALSE,
  geo_ids = c(
    "GSE305669"
  ),
  arrayexpress_ids = c(
    "E-MTAB-9330",
    "E-MTAB-10716"
  )
)

# =============================================================================
# NMD perturbation public-data screening workflow
# =============================================================================
# This is the main script to source/run.
#
# Typical use:
#   source("run_nmd_workflow.R")
#
#   result <- run_nmd_workflow(
#     repositories = c("GEO", "ArrayExpress"),
#     discover = TRUE,
#     geo_ids = c("GSE305669"),          # optional extra/manual IDs
#     arrayexpress_ids = c("E-MTAB-9330") # optional extra/manual IDs
#   )
#
# The workflow searches repositories, retrieves metadata, screens studies and
# samples for NMD perturbation evidence, and writes one Excel workbook.
# =============================================================================

nmd_workflow_root <- local({
  candidate <- file.path(
    getwd(),
    "code",
    "R",
    "NMD_Study_Search"
  )

  expected_file <- file.path(
    candidate,
    "R",
    "nmd_context.R"
  )

  if (!file.exists(expected_file)) {
    stop(
      "Could not locate the NMD workflow folder.\n",
      "Expected to find:\n  ",
      expected_file,
      "\n\nCurrent working directory:\n  ",
      getwd()
    )
  }

  normalizePath(
    candidate,
    winslash = "/",
    mustWork = TRUE
  )
})

source(file.path(nmd_workflow_root, "R", "nmd_context.R"))
source(file.path(nmd_workflow_root, "R", "geo_workflow.R"))
source(file.path(nmd_workflow_root, "R", "arrayexpress_workflow.R"))
source(file.path(nmd_workflow_root, "R", "export_workflow.R"))


run_nmd_workflow <- function(
  repositories = c("GEO", "ArrayExpress"),
  discover = TRUE,
  geo_ids = NULL,
  arrayexpress_ids = NULL,
  context = nmd_default_context(),
  include_broad_queries = TRUE,
  include_indirect_inhibitors = TRUE,
  include_stress_modulators = FALSE,
  output_dir = file.path("data", "NMD_metadata_screen"),
  workbook_name = NULL,
  geo_email = NULL,
  geo_api_key = NULL,
  overwrite_arrayexpress_metadata = FALSE,
  include_long_metadata_in_excel = TRUE,
  verbose = TRUE
) {
  repositories <- unique(tolower(trimws(repositories)))
  run_geo <- "geo" %in% repositories
  run_ae <- "arrayexpress" %in% repositories

  if (!run_geo && !run_ae) {
    stop("repositories must include 'GEO' and/or 'ArrayExpress'.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(
    file.path(output_dir, "cache"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  dir.create(
    file.path(output_dir, "ArrayExpress_metadata"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  query_manifest <- nmd_query_manifest(
    context = context,
    include_broad = include_broad_queries,
    include_indirect = include_indirect_inhibitors,
    include_stress = include_stress_modulators
  )
  query_manifest <- query_manifest[
    tolower(query_manifest$repository) %in% repositories,
    ,
    drop = FALSE
  ]

  stamp <- format(Sys.Date(), "%Y%m%d")
  if (is.null(workbook_name)) {
    workbook_name <- paste0("NMD_Perturbation_Metadata_Screen_", stamp, ".xlsx")
  }

  # --------------------------------------------------------------------------
  # GEO
  # --------------------------------------------------------------------------

  geo_search <- NULL
  geo_results <- NULL

  if (run_geo) {
    if (isTRUE(discover)) {
      if (verbose) {
        message("\n===== GEO discovery =====")
      }

      geo_search <- geo_search_nmd(
        context = context,
        include_broad = include_broad_queries,
        include_indirect = include_indirect_inhibitors,
        include_stress = include_stress_modulators,
        api_key = geo_api_key,
        email = geo_email,
        verbose = verbose
      )
    }

    discovered_geo_ids <- if (is.null(geo_search)) {
      character()
    } else {
      geo_search$accessions
    }
    all_geo_ids <- unique(c(discovered_geo_ids, geo_ids))
    all_geo_ids <- toupper(trimws(as.character(all_geo_ids)))
    all_geo_ids <- all_geo_ids[grepl("^GSE[0-9]+$", all_geo_ids)]

    if (length(all_geo_ids) > 0L) {
      if (verbose) {
        message(
          "\n===== GEO metadata + NMD screening: ",
          length(all_geo_ids),
          " studies ====="
        )
      }

      geo_results <- geo_run_nmd_screen(
        accessions = all_geo_ids,
        verbose = verbose,
        cache_file = file.path(output_dir, "cache", "GEO_metadata_latest.rds")
      )
    } else if (verbose) {
      message("No GEO accessions to screen.")
    }
  }

  # --------------------------------------------------------------------------
  # ArrayExpress / BioStudies
  # --------------------------------------------------------------------------

  ae_search <- NULL
  ae_results <- NULL

  if (run_ae) {
    if (isTRUE(discover)) {
      if (verbose) {
        message("\n===== ArrayExpress discovery =====")
      }

      ae_search <- ae_search_nmd(
        context = context,
        include_broad = include_broad_queries,
        include_indirect = include_indirect_inhibitors,
        include_stress = include_stress_modulators,
        verbose = verbose
      )
    }

    discovered_ae_ids <- if (is.null(ae_search)) {
      character()
    } else {
      ae_search$accessions
    }
    all_ae_ids <- unique(c(discovered_ae_ids, arrayexpress_ids))
    all_ae_ids <- toupper(trimws(as.character(all_ae_ids)))
    all_ae_ids <- all_ae_ids[grepl("^E-[A-Z0-9]+-[0-9]+$", all_ae_ids)]

    if (length(all_ae_ids) > 0L) {
      if (verbose) {
        message(
          "\n===== ArrayExpress metadata + NMD screening: ",
          length(all_ae_ids),
          " studies ====="
        )
      }

      ae_results <- ae_screen_nmd_accessions(
        accessions = all_ae_ids,
        base_dir = file.path(output_dir, "ArrayExpress_metadata"),
        overwrite = overwrite_arrayexpress_metadata,
        use_api_fallback = TRUE,
        query_biostudies_inventory = TRUE,
        verbose = verbose,
        prescreen = TRUE,
        prescreen_require_raw_evidence = TRUE,
        prescreen_require_perturbation_evidence = TRUE,
        prescreen_require_human_evidence = FALSE,
        prescreen_exclude_clear_nonhuman = TRUE,
        prescreen_require_cell_model_evidence = FALSE,
        prescreen_include_stress_perturbations = include_stress_modulators,
        prescreen_on_error = "continue"
      )
    } else if (verbose) {
      message("No ArrayExpress accessions to screen.")
    }
  }

  # --------------------------------------------------------------------------
  # Unified output
  # --------------------------------------------------------------------------

  workbook_path <- file.path(output_dir, workbook_name)

  exported <- export_nmd_workbook(
    path = workbook_path,
    context = context,
    geo_search = geo_search,
    geo_results = geo_results,
    ae_search = ae_search,
    ae_results = ae_results,
    query_manifest = query_manifest,
    include_long_metadata = include_long_metadata_in_excel
  )

  result <- list(
    context = context,
    query_manifest = query_manifest,
    geo_search = geo_search,
    geo = geo_results,
    arrayexpress_search = ae_search,
    arrayexpress = ae_results,
    unified_ranking = exported$ranking,
    workbook = workbook_path
  )

  saveRDS(
    result,
    file = file.path(output_dir, paste0("NMD_workflow_result_", stamp, ".rds"))
  )

  if (verbose) {
    message("\n===== Complete =====")
    message(
      "Excel workbook: ",
      normalizePath(workbook_path, winslash = "/", mustWork = FALSE)
    )
    message("Studies in unified ranking: ", nrow(result$unified_ranking))
  }

  # 'invisible()' prevents R from dumping output to console
  invisible(result)
}

# =============================================================================
# Minimal examples
# =============================================================================
#
# 1) Discover and screen both repositories:
# result <- run_nmd_workflow()
#
# 2) Only inspect specific new IDs, without doing broad database discovery:
# result <- run_nmd_workflow(
#   discover = FALSE,
#   geo_ids = c("GSE305669"),
#   arrayexpress_ids = c("E-MTAB-9330", "E-MTAB-10716")
# )
#
# 3) Search ArrayExpress only:
# result <- run_nmd_workflow(
#   repositories = "ArrayExpress"
# )
#
# 4) Search GEO only:
# result <- run_nmd_workflow(
#   repositories = "GEO"
# )

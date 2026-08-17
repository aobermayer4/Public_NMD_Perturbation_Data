# =============================================================================
# Unified study ranking and Excel export
# =============================================================================

nmd_empty_df <- function() data.frame(stringsAsFactors = FALSE)


nmd_col <- function(df, name, default = NA) {
  if (is.null(df) || !is.data.frame(df) || !name %in% names(df)) {
    return(rep(default, if (is.data.frame(df)) nrow(df) else 0L))
  }
  df[[name]]
}


nmd_nonempty_flag <- function(x) {
  x <- trimws(as.character(x))
  !is.na(x) & nzchar(x)
}


nmd_contains_any_literal <- function(text, terms) {
  terms <- unique(trimws(as.character(terms)))
  terms <- terms[!is.na(terms) & nzchar(terms)]

  vapply(
    as.character(text),
    function(x) {
      if (is.na(x) || !nzchar(x) || length(terms) == 0L) return(FALSE)
      any(vapply(
        terms,
        function(term) grepl(term, x, ignore.case = TRUE, fixed = TRUE),
        logical(1)
      ))
    },
    logical(1)
  )
}


nmd_fit_category <- function(
    perturbation,
    control,
    raw_rnaseq,
    rnaseq,
    nmd_target,
    explicit_nmd,
    score) {

  dplyr::case_when(
    perturbation & control & raw_rnaseq & rnaseq ~
      "Very strong candidate",
    perturbation & raw_rnaseq & rnaseq ~
      "Strong candidate; control review needed",
    perturbation & rnaseq ~
      "Candidate; raw-data availability needs review",
    perturbation ~
      "NMD perturbation candidate; assay/raw data needs review",
    nmd_target | explicit_nmd ~
      "NMD-related study requiring manual review",
    score >= 20 ~
      "Weak metadata match",
    TRUE ~
      "Low fit"
  )
}


build_geo_fit_table <- function(
    geo_results,
    context = nmd_default_context()) {

  studies <- geo_results$studies
  samples <- geo_results$samples

  if (is.null(studies) || !is.data.frame(studies) || nrow(studies) == 0L) {
    return(nmd_empty_df())
  }

  sample_flags <- if (!is.null(samples) && is.data.frame(samples) && nrow(samples) > 0L) {
    samples %>%
      dplyr::mutate(
        .raw_rnaseq = nmd_contains_any_literal(
          paste(
            dplyr::coalesce(as.character(.data$all_metadata), ""),
            dplyr::coalesce(as.character(.data$protocol_metadata), "")
          ),
          context$raw_rnaseq
        ),
        .cell_model = nmd_contains_any_literal(
          paste(
            dplyr::coalesce(as.character(.data$all_metadata), ""),
            dplyr::coalesce(as.character(.data$tissue_field_values), "")
          ),
          context$cell_models
        )
      ) %>%
      dplyr::group_by(.data$GSE_ID) %>%
      dplyr::summarise(
        raw_rnaseq = any(.data$.raw_rnaseq, na.rm = TRUE),
        cell_model = any(.data$.cell_model, na.rm = TRUE),
        explicit_nmd = any(.data$explicit_nmd_term, na.rm = TRUE),
        nmd_target = any(!is.na(.data$nmd_targets), na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    data.frame(
      GSE_ID = studies$GSE_ID,
      raw_rnaseq = FALSE,
      cell_model = FALSE,
      explicit_nmd = FALSE,
      nmd_target = nmd_nonempty_flag(studies$nmd_targets),
      stringsAsFactors = FALSE
    )
  }

  out <- studies %>%
    dplyr::left_join(sample_flags, by = "GSE_ID") %>%
    dplyr::mutate(
      repository = "GEO",
      Study_ID = .data$GSE_ID,
      title = .data$study_title,
      repository_tier = .data$review_category,
      perturbation = .data$n_candidate_perturbation > 0L,
      control = .data$n_candidate_controls > 0L,
      rnaseq = .data$has_rna_seq_evidence,
      raw_rnaseq = dplyr::coalesce(.data$raw_rnaseq, FALSE),
      nmd_target = dplyr::coalesce(.data$nmd_target, FALSE),
      explicit_nmd = dplyr::coalesce(.data$explicit_nmd, FALSE),
      cell_model = dplyr::coalesce(.data$cell_model, FALSE),
      standardized_fit_score =
        30L * as.integer(.data$perturbation) +
        20L * as.integer(.data$control) +
        20L * as.integer(.data$raw_rnaseq) +
        10L * as.integer(.data$rnaseq) +
        10L * as.integer(.data$nmd_target) +
        5L * as.integer(.data$explicit_nmd) +
        5L * as.integer(.data$cell_model),
      standardized_fit_category = nmd_fit_category(
        .data$perturbation,
        .data$control,
        .data$raw_rnaseq,
        .data$rnaseq,
        .data$nmd_target,
        .data$explicit_nmd,
        .data$standardized_fit_score
      ),
      n_perturbation_samples = .data$n_candidate_perturbation,
      n_control_samples = .data$n_candidate_controls,
      targets = .data$nmd_targets,
      perturbation_methods = .data$perturbation_mechanisms,
      model_or_tissue = .data$tissue_field_values,
      technical_summary = ifelse(
        .data$has_rna_seq_evidence,
        "RNA-seq evidence present in GEO metadata",
        NA_character_
      )
    ) %>%
    dplyr::select(
      .data$repository,
      .data$Study_ID,
      .data$title,
      .data$standardized_fit_score,
      .data$standardized_fit_category,
      .data$repository_tier,
      .data$perturbation,
      .data$control,
      .data$rnaseq,
      .data$raw_rnaseq,
      .data$nmd_target,
      .data$explicit_nmd,
      .data$cell_model,
      .data$n_samples,
      .data$n_perturbation_samples,
      .data$n_control_samples,
      .data$targets,
      .data$perturbation_methods,
      .data$model_or_tissue,
      .data$organism,
      .data$technical_summary,
      repository_priority_score = .data$study_priority_score
    )

  out
}


build_ae_fit_table <- function(
    ae_results,
    context = nmd_default_context()) {

  studies <- ae_results$study_summary

  if (is.null(studies) || !is.data.frame(studies) || nrow(studies) == 0L) {
    return(nmd_empty_df())
  }

  studies %>%
    dplyr::mutate(
      repository = "ArrayExpress",
      title = .data$project_title,
      repository_tier = .data$screening_tier,
      perturbation = .data$n_nmd_perturbation_samples > 0L,
      control = .data$n_control_samples > 0L,
      rnaseq =
        .data$raw_rnaseq_available |
        nmd_nonempty_flag(.data$library_strategy) |
        nmd_contains_any_literal(.data$project_title, context$rnaseq),
      raw_rnaseq = .data$raw_rnaseq_available,
      nmd_target = nmd_nonempty_flag(.data$nmd_targets),
      explicit_nmd = nmd_nonempty_flag(.data$direct_nmd_terms),
      cell_model = nmd_nonempty_flag(.data$cell_models),
      standardized_fit_score =
        30L * as.integer(.data$perturbation) +
        20L * as.integer(.data$control) +
        20L * as.integer(.data$raw_rnaseq) +
        10L * as.integer(.data$rnaseq) +
        10L * as.integer(.data$nmd_target) +
        5L * as.integer(.data$explicit_nmd) +
        5L * as.integer(.data$cell_model),
      standardized_fit_category = nmd_fit_category(
        .data$perturbation,
        .data$control,
        .data$raw_rnaseq,
        .data$rnaseq,
        .data$nmd_target,
        .data$explicit_nmd,
        .data$standardized_fit_score
      ),
      n_perturbation_samples = .data$n_nmd_perturbation_samples,
      targets = .data$nmd_targets,
      model_or_tissue = .data$cell_models,
      technical_summary = paste(
        "strategy:", dplyr::coalesce(as.character(.data$library_strategy), "unknown"),
        "| selection:", dplyr::coalesce(as.character(.data$library_selection), "unknown"),
        "| layout:", dplyr::coalesce(as.character(.data$library_layout), "unknown"),
        "| instrument:", dplyr::coalesce(as.character(.data$sequencing_instrument), "unknown")
      )
    ) %>%
    dplyr::select(
      .data$repository,
      .data$Study_ID,
      .data$title,
      .data$standardized_fit_score,
      .data$standardized_fit_category,
      .data$repository_tier,
      .data$perturbation,
      .data$control,
      .data$rnaseq,
      .data$raw_rnaseq,
      .data$nmd_target,
      .data$explicit_nmd,
      .data$cell_model,
      .data$n_samples,
      .data$n_perturbation_samples,
      .data$n_control_samples,
      .data$targets,
      perturbation_methods = .data$perturbation_methods,
      .data$model_or_tissue,
      organism = NA_character_,
      .data$technical_summary,
      repository_priority_score = NA_real_
    )
}


build_unified_study_ranking <- function(
    geo_results = NULL,
    ae_results = NULL,
    context = nmd_default_context()) {

  geo <- if (is.null(geo_results)) nmd_empty_df() else {
    build_geo_fit_table(geo_results, context = context)
  }

  ae <- if (is.null(ae_results)) nmd_empty_df() else {
    build_ae_fit_table(ae_results, context = context)
  }

  out <- dplyr::bind_rows(geo, ae)

  if (nrow(out) == 0L) return(out)

  tier_order <- c(
    "Very strong candidate",
    "Strong candidate; control review needed",
    "Candidate; raw-data availability needs review",
    "NMD perturbation candidate; assay/raw data needs review",
    "NMD-related study requiring manual review",
    "Weak metadata match",
    "Low fit"
  )

  out %>%
    dplyr::arrange(
      factor(.data$standardized_fit_category, levels = tier_order),
      dplyr::desc(.data$standardized_fit_score),
      .data$repository,
      .data$Study_ID
    )
}


# ----------------------------------------------------------------------------
# Excel helpers
# ----------------------------------------------------------------------------

nmd_excel_safe_df <- function(df, max_chars = 30000L) {
  if (is.null(df) || !is.data.frame(df)) return(nmd_empty_df())

  out <- as.data.frame(df, stringsAsFactors = FALSE)

  out[] <- lapply(out, function(x) {
    if (is.list(x)) {
      x <- vapply(
        x,
        function(z) paste(as.character(unlist(z, recursive = TRUE)), collapse = "; "),
        character(1)
      )
    }

    if (is.factor(x)) x <- as.character(x)

    if (is.character(x)) {
      x <- enc2utf8(x)
      too_long <- !is.na(x) & nchar(x, type = "chars") > max_chars
      x[too_long] <- paste0(substr(x[too_long], 1L, max_chars - 20L), " [truncated]")
    }

    x
  })

  out
}


nmd_workbook_readme <- function() {
  data.frame(
    item = c(
      "Purpose",
      "Unified_Ranking",
      "Repository-specific tiers",
      "Standardized score",
      "Raw RNA-seq caution",
      "Evidence tables",
      "Manual review"
    ),
    description = c(
      "Screen public GEO and ArrayExpress studies for NMD perturbation RNA-seq candidates.",
      "One row per study across repositories, sorted by a transparent 0-100 metadata fit score.",
      "GEO review_category and ArrayExpress screening_tier are preserved; they remain the repository-native classifications.",
      "30 perturbation + 20 control + 20 raw RNA-seq + 10 RNA-seq + 10 NMD target + 5 explicit NMD language + 5 cell-model evidence.",
      "ArrayExpress raw-file evidence is explicit. GEO raw-file evidence is inferred from terms/accessions retained in GEO metadata and should be confirmed before download.",
      "Use GEO_Evidence and AE_Evidence to see the metadata that produced matches rather than trusting the summary blindly.",
      "A high score is a screening result, not proof of the experimental design. Controls, indirect inhibitors, rescue experiments, and selection drugs still need scientific review."
    ),
    stringsAsFactors = FALSE
  )
}


export_nmd_workbook <- function(
    path,
    context = nmd_default_context(),
    geo_search = NULL,
    geo_results = NULL,
    ae_search = NULL,
    ae_results = NULL,
    query_manifest = NULL,
    include_long_metadata = TRUE) {

  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop("Package 'writexl' is required for Excel export: install.packages('writexl').")
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  ranking <- build_unified_study_ranking(
    geo_results = geo_results,
    ae_results = ae_results,
    context = context
  )

  sheets <- list(
    README = nmd_workbook_readme(),
    Search_Context = nmd_context_table(context),
    Search_Queries = if (is.null(query_manifest)) nmd_query_manifest(context) else query_manifest,
    Unified_Ranking = ranking
  )

  if (!is.null(geo_search)) {
    sheets$GEO_Search <- geo_search$unique_studies
    sheets$GEO_All_Query_Hits <- geo_search$all_query_hits
  }

  if (!is.null(geo_results)) {
    sheets$GEO_Studies <- geo_results$studies
    sheets$GEO_Samples <- geo_results$samples
    sheets$GEO_Evidence <- geo_results$evidence
    sheets$GEO_Fields <- geo_results$fields
    sheets$GEO_Errors <- geo_results$errors

    if (isTRUE(include_long_metadata)) {
      sheets$GEO_Long <- geo_results$long_metadata
    }
  }

  if (!is.null(ae_search)) {
    sheets$AE_Search <- ae_search$unique_studies
    sheets$AE_All_Query_Hits <- ae_search$all_query_hits
  }

  if (!is.null(ae_results)) {
    ae_map <- c(
      AE_Prescreen = "prescreen_summary",
      AE_PreEvidence = "prescreen_evidence",
      AE_Skipped = "skipped_prescreen",
      AE_Studies = "study_summary",
      AE_Samples = "sample_summary",
      AE_Eligible = "eligible_sample_summary",
      AE_Perturb = "perturbation_sample_summary",
      AE_Controls = "control_sample_summary",
      AE_Evidence = "evidence",
      AE_Files = "file_inventory",
      AE_Errors = "errors"
    )

    for (sheet_name in names(ae_map)) {
      table_name <- unname(ae_map[[sheet_name]])
      sheets[[sheet_name]] <- ae_results[[table_name]]
    }

    if (isTRUE(include_long_metadata)) {
      sheets$AE_SDRF_Long <- ae_results$sdrf_long
      sheets$AE_IDF_Long <- ae_results$idf_long
    }
  }

  sheets <- lapply(sheets, nmd_excel_safe_df)

  # Remove entirely NULL entries while preserving useful zero-row tables.
  sheets <- sheets[!vapply(sheets, is.null, logical(1))]

  writexl::write_xlsx(
    x = sheets,
    path = path,
    format_headers = TRUE
  )

  invisible(list(path = path, ranking = ranking, sheets = names(sheets)))
}

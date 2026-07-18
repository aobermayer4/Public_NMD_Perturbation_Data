build_geo_index <- function(geo) {
  samples <- data.table::as.data.table(data.table::copy(geo$samples %||% empty_dt()))
  studies <- data.table::as.data.table(data.table::copy(geo$studies %||% empty_dt()))
  metadata <- data.table::as.data.table(data.table::copy(geo$metadata %||% empty_dt()))
  evidence <- data.table::as.data.table(data.table::copy(geo$evidence %||% empty_dt()))

  if (nrow(samples) == 0L || !all(c("GSE_ID", "Sample_ID") %in% names(samples))) {
    return(list(samples = empty_dt(), studies = empty_dt()))
  }

  samples[, annotation_text := combine_row_text(.SD), .SDcols = names(samples)]

  if (nrow(metadata) > 0L && all(c("GSE_ID", "Sample_ID", "field_raw", "value") %in% names(metadata))) {
    metadata[, meta_piece := paste0(field_raw, ": ", value)]
    meta_agg <- metadata[, .(
      metadata_text = collapse_unique_text(meta_piece, max_values = Inf)
    ), by = .(GSE_ID, Sample_ID)]
  } else {
    meta_agg <- data.table::data.table(GSE_ID = character(), Sample_ID = character(), metadata_text = character())
  }

  if (nrow(evidence) > 0L && all(c("GSE_ID", "Sample_ID", "evidence_type", "field_raw", "value") %in% names(evidence))) {
    evidence[, evidence_piece := paste0(evidence_type, " | ", field_raw, ": ", value)]
    evidence_agg <- evidence[, .(
      evidence_text = collapse_unique_text(evidence_piece, max_values = Inf)
    ), by = .(GSE_ID, Sample_ID)]
  } else {
    evidence_agg <- data.table::data.table(GSE_ID = character(), Sample_ID = character(), evidence_text = character())
  }

  sample_index <- merge(samples, meta_agg, by = c("GSE_ID", "Sample_ID"), all.x = TRUE, sort = FALSE)
  sample_index <- merge(sample_index, evidence_agg, by = c("GSE_ID", "Sample_ID"), all.x = TRUE, sort = FALSE)
  sample_index[, sample_text := paste(annotation_text, metadata_text, evidence_text, sep = " | ")]
  if ("study_title" %in% names(sample_index)) {
    sample_index[, sample_text := mapply(
      function(text, title) {
        if (is.na(title) || !nzchar(title)) text else gsub(title, " ", text, fixed = TRUE)
      },
      sample_text,
      study_title,
      USE.NAMES = FALSE
    )]
  }

  if (nrow(studies) > 0L && "GSE_ID" %in% names(studies)) {
    studies[, study_text := combine_row_text(.SD), .SDcols = names(studies)]
  }

  list(samples = sample_index, studies = studies)
}

run_geo_search <- function(geo, index, filters) {
  samples_raw <- data.table::as.data.table(data.table::copy(geo$samples %||% empty_dt()))
  studies_raw <- data.table::as.data.table(data.table::copy(geo$studies %||% empty_dt()))
  metadata <- data.table::as.data.table(data.table::copy(geo$metadata %||% empty_dt()))
  evidence <- data.table::as.data.table(data.table::copy(geo$evidence %||% empty_dt()))

  sample_index <- data.table::as.data.table(data.table::copy(index$samples))
  study_index <- data.table::as.data.table(data.table::copy(index$studies))

  if (nrow(sample_index) == 0L) stop("No GEO sample annotation table is loaded.")

  allowed_studies <- unique(sample_index$GSE_ID)
  if (nrow(study_index) > 0L) allowed_studies <- unique(study_index$GSE_ID)

  if (length(filters$study_ids) > 0L) {
    allowed_studies <- intersect(allowed_studies, filters$study_ids)
  }

  if (length(filters$organisms) > 0L && nrow(study_index) > 0L && "organism" %in% names(study_index)) {
    allowed_studies <- intersect(
      allowed_studies,
      study_index[grepl(or_pattern(filters$organisms), organism, ignore.case = TRUE, perl = TRUE), GSE_ID]
    )
  }

  if (length(filters$review_categories) > 0L && nrow(study_index) > 0L && "review_category" %in% names(study_index)) {
    allowed_studies <- intersect(
      allowed_studies,
      study_index[review_category %in% filters$review_categories, GSE_ID]
    )
  }

  sample_index <- sample_index[GSE_ID %in% allowed_studies]

  if (length(filters$targets) > 0L) {
    sample_index <- sample_index[match_selected_concepts(sample_text, filters$targets, NMD_TARGETS)]
  }
  if (length(filters$mechanisms) > 0L) {
    sample_index <- sample_index[match_selected_concepts(sample_text, filters$mechanisms, NMD_MECHANISMS)]
  }
  if (length(filters$agents) > 0L) {
    sample_index <- sample_index[match_selected_concepts(sample_text, filters$agents, NMD_AGENTS)]
  }
  if (length(filters$controls) > 0L) {
    sample_index <- sample_index[match_selected_concepts(sample_text, filters$controls, CONTROL_CONCEPTS)]
  }
  if (length(filters$model_presets) > 0L) {
    sample_index <- sample_index[match_selected_concepts(sample_text, filters$model_presets, CELL_MODEL_CONCEPTS)]
  }
  if (nzchar(trimws(filters$model_text %||% ""))) {
    sample_index <- sample_index[grepl(regex_escape(filters$model_text), sample_text, ignore.case = TRUE, perl = TRUE)]
  }
  if (length(filters$sample_roles) > 0L && "sample_role" %in% names(sample_index)) {
    sample_index <- sample_index[sample_role %in% filters$sample_roles]
  }
  if (length(filters$confidence) > 0L && "evidence_confidence" %in% names(sample_index)) {
    sample_index <- sample_index[evidence_confidence %in% filters$confidence]
  }

  query_info <- smart_query_groups(
    filters$query,
    use_synonyms = isTRUE(filters$use_synonyms),
    mode = filters$match_mode
  )

  query_present <- nzchar(trimws(filters$query %||% ""))
  sample_query_hit <- if (query_present) {
    match_search_groups(sample_index$sample_text, query_info$groups, filters$match_mode)
  } else {
    rep(TRUE, nrow(sample_index))
  }

  study_query_ids <- character()
  if (query_present && nrow(study_index) > 0L && "study_text" %in% names(study_index)) {
    study_query_ids <- study_index[
      GSE_ID %in% allowed_studies &
        match_search_groups(study_text, query_info$groups, filters$match_mode),
      GSE_ID
    ]
  }

  if (query_present) {
    if (isTRUE(filters$include_study_context)) {
      sample_index[, query_match_source := ifelse(
        sample_query_hit & GSE_ID %in% study_query_ids,
        "sample and study text",
        ifelse(sample_query_hit, "sample text", ifelse(GSE_ID %in% study_query_ids, "study context", NA_character_))
      )]
      sample_index <- sample_index[sample_query_hit | GSE_ID %in% study_query_ids]
    } else {
      sample_index[, query_match_source := ifelse(sample_query_hit, "sample text", NA_character_)]
      sample_index <- sample_index[sample_query_hit]
    }
  } else {
    sample_index[, query_match_source := "structured filters"]
  }

  sample_keys <- unique(sample_index[, .(GSE_ID, Sample_ID)])
  result_study_ids <- unique(c(sample_keys$GSE_ID, study_query_ids))
  result_study_ids <- intersect(result_study_ids, allowed_studies)

  studies_out <- if (nrow(studies_raw) > 0L) studies_raw[GSE_ID %in% result_study_ids] else empty_dt()

  samples_out <- if (nrow(sample_keys) > 0L) {
    merge(
      samples_raw,
      unique(sample_index[, .(GSE_ID, Sample_ID, query_match_source)]),
      by = c("GSE_ID", "Sample_ID"),
      all = FALSE,
      sort = FALSE
    )
  } else {
    samples_raw[0]
  }

  evidence_out <- if (nrow(sample_keys) > 0L && nrow(evidence) > 0L) {
    merge(evidence, sample_keys, by = c("GSE_ID", "Sample_ID"), all = FALSE, sort = FALSE)
  } else {
    evidence[0]
  }

  metadata_out <- if (nrow(sample_keys) > 0L && nrow(metadata) > 0L) {
    merge(metadata, sample_keys, by = c("GSE_ID", "Sample_ID"), all = FALSE, sort = FALSE)
  } else {
    metadata[0]
  }

  fields_out <- make_field_summary(metadata_out)

  list(
    tables = list(
      GEO_studies = studies_out,
      GEO_samples = samples_out,
      GEO_evidence = evidence_out,
      GEO_metadata_long = metadata_out,
      GEO_field_summary = fields_out
    ),
    interpretation = query_info$display,
    counts = c(
      Studies = nrow(studies_out),
      Samples = nrow(samples_out),
      Evidence = nrow(evidence_out),
      `Metadata rows` = nrow(metadata_out)
    )
  )
}

geo_explorer_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::sidebarLayout(
    shiny::sidebarPanel(
      width = 3,
      shiny::h4("Plain-language search"),
      shiny::textAreaInput(
        ns("query"),
        "Describe what you are looking for",
        placeholder = "Example: UPF1 knockdown in HEK293 cells",
        rows = 3
      ),
      shiny::selectInput(
        ns("match_mode"),
        "How should free-text concepts combine?",
        choices = c(
          "Match all concepts" = "all",
          "Match any concept" = "any",
          "Match the exact phrase" = "exact"
        ),
        selected = "all"
      ),
      shiny::checkboxInput(ns("use_synonyms"), "Expand common aliases and synonyms", TRUE),
      shiny::checkboxInput(ns("include_study_context"), "Include all samples when study text matches", TRUE),
      shiny::hr(),
      shiny::h4("Preset filters"),
      shiny::selectizeInput(ns("targets"), "NMD factor", choices = names(NMD_TARGETS), multiple = TRUE),
      shiny::selectizeInput(ns("mechanisms"), "Perturbation mechanism", choices = names(NMD_MECHANISMS), multiple = TRUE),
      shiny::selectizeInput(ns("agents"), "Drug or treatment", choices = names(NMD_AGENTS), multiple = TRUE),
      shiny::selectizeInput(ns("controls"), "Control type", choices = names(CONTROL_CONCEPTS), multiple = TRUE),
      shiny::selectizeInput(ns("model_presets"), "Common biological model", choices = names(CELL_MODEL_CONCEPTS), multiple = TRUE),
      shiny::textInput(ns("model_text"), "Other cell line, tissue, or model", placeholder = "HeLa, nasopharyngeal carcinoma, T cells..."),
      shiny::selectizeInput(ns("organisms"), "Organism", choices = NULL, multiple = TRUE),
      shiny::selectizeInput(ns("sample_roles"), "Sample classification", choices = NULL, multiple = TRUE),
      shiny::selectizeInput(ns("confidence"), "Evidence confidence", choices = NULL, multiple = TRUE),
      shiny::selectizeInput(ns("review_categories"), "Study review category", choices = NULL, multiple = TRUE),
      shiny::selectizeInput(ns("study_ids"), "Limit to GSE IDs", choices = NULL, multiple = TRUE),
      shiny::div(
        class = "action-row",
        shiny::actionButton(ns("run"), "Run search", class = "btn-primary"),
        shiny::actionButton(ns("reset"), "Reset")
      ),
      shiny::br(),
      shiny::downloadButton(ns("download_zip"), "Download filtered ZIP", class = "btn-success")
    ),
    shiny::mainPanel(
      width = 9,
      shiny::uiOutput(ns("summary")),
      shiny::tabsetPanel(
        id = ns("results_tabs"),
        shiny::tabPanel("Studies", DT::DTOutput(ns("studies"))),
        shiny::tabPanel("Samples", DT::DTOutput(ns("samples"))),
        shiny::tabPanel("Evidence", DT::DTOutput(ns("evidence"))),
        shiny::tabPanel("Metadata", DT::DTOutput(ns("metadata"))),
        shiny::tabPanel("Field summary", DT::DTOutput(ns("fields"))),
        shiny::tabPanel(
          "Search interpretation",
          shiny::div(class = "interpretation-box", shiny::verbatimTextOutput(ns("interpretation")))
        )
      )
    )
  )
}

geo_explorer_server <- function(id, geo_data) {
  shiny::moduleServer(id, function(input, output, session) {
    index <- shiny::reactive({
      geo <- geo_data()
      shiny::req(geo)
      build_geo_index(geo)
    })

    shiny::observeEvent(geo_data(), {
      geo <- geo_data()
      studies <- data.table::as.data.table(geo$studies %||% empty_dt())
      samples <- data.table::as.data.table(geo$samples %||% empty_dt())

      shiny::updateSelectizeInput(session, "study_ids", choices = safe_unique(studies$GSE_ID), server = TRUE)
      shiny::updateSelectizeInput(session, "organisms", choices = split_multivalue(c(studies$organism, samples$organism)), server = TRUE)
      shiny::updateSelectizeInput(session, "sample_roles", choices = safe_unique(samples$sample_role), server = TRUE)
      shiny::updateSelectizeInput(session, "confidence", choices = safe_unique(samples$evidence_confidence), server = TRUE)
      shiny::updateSelectizeInput(session, "review_categories", choices = safe_unique(studies$review_category), server = TRUE)
    }, ignoreInit = FALSE)

    shiny::observeEvent(input$reset, {
      shiny::updateTextAreaInput(session, "query", value = "")
      shiny::updateSelectInput(session, "match_mode", selected = "all")
      shiny::updateCheckboxInput(session, "use_synonyms", value = TRUE)
      shiny::updateCheckboxInput(session, "include_study_context", value = TRUE)
      for (id2 in c("targets", "mechanisms", "agents", "controls", "model_presets", "organisms", "sample_roles", "confidence", "review_categories", "study_ids")) {
        shiny::updateSelectizeInput(session, id2, selected = character())
      }
      shiny::updateTextInput(session, "model_text", value = "")
    })

    filters <- shiny::eventReactive(input$run, {
      list(
        query = input$query %||% "",
        match_mode = input$match_mode %||% "all",
        use_synonyms = isTRUE(input$use_synonyms),
        include_study_context = isTRUE(input$include_study_context),
        targets = input$targets %||% character(),
        mechanisms = input$mechanisms %||% character(),
        agents = input$agents %||% character(),
        controls = input$controls %||% character(),
        model_presets = input$model_presets %||% character(),
        model_text = input$model_text %||% "",
        organisms = input$organisms %||% character(),
        sample_roles = input$sample_roles %||% character(),
        confidence = input$confidence %||% character(),
        review_categories = input$review_categories %||% character(),
        study_ids = input$study_ids %||% character()
      )
    }, ignoreInit = TRUE)

    result <- shiny::eventReactive(input$run, {
      geo <- geo_data()
      validation <- validate_geo_bundle(geo)
      if (!is.null(validation)) stop(validation)
      run_geo_search(geo, index(), filters())
    }, ignoreInit = TRUE)

    output$summary <- shiny::renderUI({
      if (input$run == 0L) {
        return(shiny::div(class = "empty-state", "Choose filters and click Run search."))
      }
      res <- result()
      summary_cards(as.list(res$counts))
    })

    output$studies <- DT::renderDT({ shiny::req(result()); result_table(result()$tables$GEO_studies) }, server = TRUE)
    output$samples <- DT::renderDT({ shiny::req(result()); result_table(result()$tables$GEO_samples) }, server = TRUE)
    output$evidence <- DT::renderDT({ shiny::req(result()); result_table(result()$tables$GEO_evidence) }, server = TRUE)
    output$metadata <- DT::renderDT({ shiny::req(result()); result_table(result()$tables$GEO_metadata_long) }, server = TRUE)
    output$fields <- DT::renderDT({ shiny::req(result()); result_table(result()$tables$GEO_field_summary) }, server = TRUE)

    output$interpretation <- shiny::renderText({
      shiny::req(result())
      paste0(
        "Free-text interpretation:\n", result()$interpretation,
        "\n\nSelected structured filters:\n",
        "Targets: ", compact_filter_summary(filters()$targets), "\n",
        "Mechanisms: ", compact_filter_summary(filters()$mechanisms), "\n",
        "Agents: ", compact_filter_summary(filters()$agents), "\n",
        "Controls: ", compact_filter_summary(filters()$controls), "\n",
        "Models: ", compact_filter_summary(c(filters()$model_presets, filters()$model_text)), "\n",
        "Organisms: ", compact_filter_summary(filters()$organisms), "\n",
        "Sample roles: ", compact_filter_summary(filters()$sample_roles), "\n",
        "Confidence: ", compact_filter_summary(filters()$confidence), "\n",
        "Review category: ", compact_filter_summary(filters()$review_categories), "\n",
        "GSE IDs: ", compact_filter_summary(filters()$study_ids)
      )
    })

    output$download_zip <- shiny::downloadHandler(
      filename = function() paste0("GEO_NMD_subset_", format(Sys.Date(), "%Y%m%d"), ".zip"),
      content = function(file) {
        shiny::req(result(), filters())
        write_result_zip(
          result(),
          platform = "GEO",
          filters = filters(),
          interpretation = result()$interpretation,
          output_file = file
        )
      }
    )
  })
}

build_ae_index <- function(ae) {
  samples <- data.table::as.data.table(data.table::copy(ae$samples %||% empty_dt()))
  sdrf_long <- data.table::as.data.table(data.table::copy(ae$sdrf_long %||% empty_dt()))
  studies <- data.table::as.data.table(data.table::copy(ae$studies %||% empty_dt()))
  idf_long <- data.table::as.data.table(data.table::copy(ae$idf_long %||% empty_dt()))
  protocols <- data.table::as.data.table(data.table::copy(ae$protocols %||% empty_dt()))
  refs <- data.table::as.data.table(data.table::copy(ae$protocol_refs %||% empty_dt()))

  key_cols <- c("Study_ID", "Sample_ID", "SDRF_Row")
  if (nrow(samples) == 0L || !all(key_cols %in% names(samples))) {
    return(list(samples = empty_dt(), studies = empty_dt(), protocol_links = empty_dt()))
  }

  samples[, wide_text := combine_row_text(.SD), .SDcols = names(samples)]

  if (nrow(sdrf_long) > 0L && all(c(key_cols, "field_raw", "value") %in% names(sdrf_long))) {
    sdrf_long[, metadata_piece := paste0(field_raw, ": ", value)]
    sample_meta <- sdrf_long[, .(
      sdrf_text = collapse_unique_text(metadata_piece, max_values = Inf)
    ), by = key_cols]
  } else {
    sample_meta <- data.table::data.table(Study_ID = character(), Sample_ID = character(), SDRF_Row = character(), sdrf_text = character())
  }

  protocol_links <- empty_dt()
  protocol_agg <- data.table::data.table(Study_ID = character(), Sample_ID = character(), SDRF_Row = character(), protocol_text = character())

  if (nrow(refs) > 0L && nrow(protocols) > 0L &&
      all(c(key_cols, "Protocol_REF") %in% names(refs)) &&
      all(c("Study_ID", "protocol_name") %in% names(protocols))) {
    protocol_links <- merge(
      refs,
      protocols,
      by.x = c("Study_ID", "Protocol_REF"),
      by.y = c("Study_ID", "protocol_name"),
      all.x = TRUE,
      sort = FALSE
    )
    protocol_links[, protocol_piece := combine_row_text(.SD), .SDcols = names(protocol_links)]
    protocol_agg <- protocol_links[, .(
      protocol_text = collapse_unique_text(protocol_piece, max_values = Inf)
    ), by = key_cols]
  }

  sample_index <- merge(samples, sample_meta, by = key_cols, all.x = TRUE, sort = FALSE)
  sample_index <- merge(sample_index, protocol_agg, by = key_cols, all.x = TRUE, sort = FALSE)

  if (nrow(idf_long) > 0L && all(c("Study_ID", "field_raw", "value") %in% names(idf_long))) {
    idf_long[, idf_piece := paste0(field_raw, ": ", value)]
    idf_agg <- idf_long[, .(
      idf_text = collapse_unique_text(idf_piece, max_values = Inf)
    ), by = Study_ID]
  } else {
    idf_agg <- data.table::data.table(Study_ID = character(), idf_text = character())
  }

  if (nrow(studies) > 0L && "Study_ID" %in% names(studies)) {
    studies[, summary_text := combine_row_text(.SD), .SDcols = names(studies)]
    study_index <- merge(studies, idf_agg, by = "Study_ID", all.x = TRUE, sort = FALSE)
    study_index[, study_text := paste(summary_text, idf_text, sep = " | ")]
    sample_index <- merge(sample_index, study_index[, .(Study_ID, study_text)], by = "Study_ID", all.x = TRUE, sort = FALSE)
  } else {
    study_index <- empty_dt()
    sample_index[, study_text := NA_character_]
  }

  sample_index[, sample_text := paste(wide_text, sdrf_text, protocol_text, sep = " | ")]
  sample_index[, structured_text := paste(sample_text, study_text, sep = " | ")]

  list(samples = sample_index, studies = study_index, protocol_links = protocol_links)
}

run_ae_search <- function(ae, index, filters) {
  studies_raw <- data.table::as.data.table(data.table::copy(ae$studies %||% empty_dt()))
  samples_raw <- data.table::as.data.table(data.table::copy(ae$samples %||% empty_dt()))
  sdrf_long <- data.table::as.data.table(data.table::copy(ae$sdrf_long %||% empty_dt()))
  idf_long <- data.table::as.data.table(data.table::copy(ae$idf_long %||% empty_dt()))
  idf_summary <- data.table::as.data.table(data.table::copy(ae$idf_summary %||% empty_dt()))
  protocols <- data.table::as.data.table(data.table::copy(ae$protocols %||% empty_dt()))
  refs <- data.table::as.data.table(data.table::copy(ae$protocol_refs %||% empty_dt()))
  errors <- data.table::as.data.table(data.table::copy(ae$errors %||% empty_dt()))

  sample_index <- data.table::as.data.table(data.table::copy(index$samples))
  study_index <- data.table::as.data.table(data.table::copy(index$studies))
  protocol_links <- data.table::as.data.table(data.table::copy(index$protocol_links))
  key_cols <- c("Study_ID", "Sample_ID", "SDRF_Row")

  if (nrow(sample_index) == 0L) stop("No ArrayExpress SDRF sample table is loaded.")

  allowed_studies <- unique(sample_index$Study_ID)
  if (nrow(study_index) > 0L) allowed_studies <- unique(study_index$Study_ID)

  if (length(filters$study_ids) > 0L) {
    allowed_studies <- intersect(allowed_studies, filters$study_ids)
  }
  if (length(filters$organisms) > 0L && nrow(study_index) > 0L && "organism" %in% names(study_index)) {
    allowed_studies <- intersect(
      allowed_studies,
      study_index[grepl(or_pattern(filters$organisms), organism, ignore.case = TRUE, perl = TRUE), Study_ID]
    )
  }

  sample_index <- sample_index[Study_ID %in% allowed_studies]

  if (length(filters$targets) > 0L) {
    sample_index <- sample_index[match_selected_concepts(structured_text, filters$targets, NMD_TARGETS)]
  }
  if (length(filters$mechanisms) > 0L) {
    sample_index <- sample_index[match_selected_concepts(structured_text, filters$mechanisms, NMD_MECHANISMS)]
  }
  if (length(filters$agents) > 0L) {
    sample_index <- sample_index[match_selected_concepts(structured_text, filters$agents, NMD_AGENTS)]
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

  if (length(filters$protocol_types) > 0L) {
    if (nrow(protocol_links) == 0L || !"protocol_type" %in% names(protocol_links)) {
      sample_index <- sample_index[0]
    } else {
      protocol_keys <- unique(protocol_links[protocol_type %in% filters$protocol_types, ..key_cols])
      sample_index <- merge(sample_index, protocol_keys, by = key_cols, all = FALSE, sort = FALSE)
    }
  }

  if (nzchar(filters$characteristic_field %||% "")) {
    characteristic_rows <- sdrf_long[field_raw == filters$characteristic_field]
    if (length(filters$characteristic_values) > 0L) {
      characteristic_rows <- characteristic_rows[value %in% filters$characteristic_values]
    }
    characteristic_keys <- unique(characteristic_rows[, ..key_cols])
    sample_index <- merge(sample_index, characteristic_keys, by = key_cols, all = FALSE, sort = FALSE)
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
      Study_ID %in% allowed_studies &
        match_search_groups(study_text, query_info$groups, filters$match_mode),
      Study_ID
    ]
  }

  if (query_present) {
    if (isTRUE(filters$include_study_context)) {
      sample_index[, query_match_source := ifelse(
        sample_query_hit & Study_ID %in% study_query_ids,
        "sample and study text",
        ifelse(sample_query_hit, "sample or protocol text", ifelse(Study_ID %in% study_query_ids, "study context", NA_character_))
      )]
      sample_index <- sample_index[sample_query_hit | Study_ID %in% study_query_ids]
    } else {
      sample_index[, query_match_source := ifelse(sample_query_hit, "sample or protocol text", NA_character_)]
      sample_index <- sample_index[sample_query_hit]
    }
  } else {
    sample_index[, query_match_source := "structured filters"]
  }

  sample_keys <- unique(sample_index[, ..key_cols])
  result_study_ids <- unique(c(sample_keys$Study_ID, study_query_ids))
  result_study_ids <- intersect(result_study_ids, allowed_studies)

  studies_out <- if (nrow(studies_raw) > 0L) studies_raw[Study_ID %in% result_study_ids] else empty_dt()

  samples_out <- if (nrow(sample_keys) > 0L) {
    merge(
      samples_raw,
      unique(sample_index[, c(key_cols, "query_match_source"), with = FALSE]),
      by = key_cols,
      all = FALSE,
      sort = FALSE
    )
  } else {
    samples_raw[0]
  }

  sdrf_out <- if (nrow(sample_keys) > 0L && nrow(sdrf_long) > 0L) {
    merge(sdrf_long, sample_keys, by = key_cols, all = FALSE, sort = FALSE)
  } else {
    sdrf_long[0]
  }

  refs_out <- if (nrow(sample_keys) > 0L && nrow(refs) > 0L) {
    merge(refs, sample_keys, by = key_cols, all = FALSE, sort = FALSE)
  } else {
    refs[0]
  }

  referenced_protocols <- if (nrow(refs_out) > 0L && "Protocol_REF" %in% names(refs_out)) unique(refs_out[, .(Study_ID, Protocol_REF)]) else empty_dt()

  if (nrow(protocols) > 0L && length(result_study_ids) > 0L) {
    protocols_out <- protocols[Study_ID %in% result_study_ids]
    if (nrow(referenced_protocols) > 0L && "protocol_name" %in% names(protocols_out)) {
      protocols_out <- merge(
        protocols_out,
        referenced_protocols,
        by.x = c("Study_ID", "protocol_name"),
        by.y = c("Study_ID", "Protocol_REF"),
        all = FALSE,
        sort = FALSE
      )
    }
  } else {
    protocols_out <- protocols[0]
  }

  idf_long_out <- if (nrow(idf_long) > 0L) idf_long[Study_ID %in% result_study_ids] else idf_long[0]
  idf_summary_out <- if (nrow(idf_summary) > 0L) idf_summary[Study_ID %in% result_study_ids] else idf_summary[0]
  errors_out <- if (nrow(errors) > 0L && "Study_ID" %in% names(errors)) errors[Study_ID %in% unique(c(result_study_ids, filters$study_ids))] else errors
  fields_out <- make_field_summary(sdrf_out)

  list(
    tables = list(
      ArrayExpress_studies = studies_out,
      ArrayExpress_samples = samples_out,
      ArrayExpress_SDRF_long = sdrf_out,
      ArrayExpress_protocols = protocols_out,
      ArrayExpress_sample_protocol_refs = refs_out,
      ArrayExpress_IDF_long = idf_long_out,
      ArrayExpress_IDF_summary = idf_summary_out,
      ArrayExpress_SDRF_field_summary = fields_out,
      ArrayExpress_download_errors = errors_out
    ),
    interpretation = query_info$display,
    counts = c(
      Studies = nrow(studies_out),
      `SDRF rows` = nrow(samples_out),
      `Sample metadata rows` = nrow(sdrf_out),
      Protocols = nrow(protocols_out)
    )
  )
}

ae_explorer_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::sidebarLayout(
    shiny::sidebarPanel(
      width = 3,
      shiny::h4("Plain-language search"),
      shiny::textAreaInput(
        ns("query"),
        "Describe what you are looking for",
        placeholder = "Example: SMG7 knockout with SMG6 knockdown",
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
      shiny::textInput(ns("model_text"), "Other cell line, tissue, or model", placeholder = "Flp-In T-REx 293, HeLa, patient tissue..."),
      shiny::selectizeInput(ns("organisms"), "Organism", choices = NULL, multiple = TRUE),
      shiny::selectizeInput(ns("protocol_types"), "Protocol type", choices = NULL, multiple = TRUE),
      shiny::selectizeInput(ns("characteristic_field"), "Sample characteristic or factor", choices = NULL, multiple = FALSE),
      shiny::selectizeInput(ns("characteristic_values"), "Characteristic value", choices = NULL, multiple = TRUE),
      shiny::selectizeInput(ns("study_ids"), "Limit to E-MTAB IDs", choices = NULL, multiple = TRUE),
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
        shiny::tabPanel("Sample metadata", DT::DTOutput(ns("sdrf"))),
        shiny::tabPanel("Protocols", DT::DTOutput(ns("protocols"))),
        shiny::tabPanel("Protocol links", DT::DTOutput(ns("refs"))),
        shiny::tabPanel("IDF metadata", DT::DTOutput(ns("idf"))),
        shiny::tabPanel("IDF summary", DT::DTOutput(ns("idf_summary"))),
        shiny::tabPanel("Field summary", DT::DTOutput(ns("fields"))),
        shiny::tabPanel("Download errors", DT::DTOutput(ns("errors"))),
        shiny::tabPanel(
          "Search interpretation",
          shiny::div(class = "interpretation-box", shiny::verbatimTextOutput(ns("interpretation")))
        )
      )
    )
  )
}

ae_explorer_server <- function(id, ae_data) {
  shiny::moduleServer(id, function(input, output, session) {
    index <- shiny::reactive({
      ae <- ae_data()
      shiny::req(ae)
      build_ae_index(ae)
    })

    shiny::observeEvent(ae_data(), {
      ae <- ae_data()
      studies <- data.table::as.data.table(ae$studies %||% empty_dt())
      sdrf <- data.table::as.data.table(ae$sdrf_long %||% empty_dt())
      protocols <- data.table::as.data.table(ae$protocols %||% empty_dt())

      shiny::updateSelectizeInput(session, "study_ids", choices = safe_unique(studies$Study_ID), server = TRUE)
      shiny::updateSelectizeInput(session, "organisms", choices = split_multivalue(studies$organism), server = TRUE)
      shiny::updateSelectizeInput(session, "protocol_types", choices = safe_unique(protocols$protocol_type), server = TRUE)

      characteristic_fields <- if (nrow(sdrf) > 0L && "field_raw" %in% names(sdrf)) {
        safe_unique(sdrf[grepl("^(Characteristics|Factor Value|Parameter Value)\\[", field_raw, ignore.case = TRUE), field_raw])
      } else character()
      shiny::updateSelectizeInput(session, "characteristic_field", choices = c("Any characteristic" = "", stats::setNames(characteristic_fields, characteristic_fields)), selected = "", server = TRUE)
    }, ignoreInit = FALSE)

    shiny::observeEvent(input$characteristic_field, {
      ae <- ae_data()
      sdrf <- data.table::as.data.table(ae$sdrf_long %||% empty_dt())
      values <- character()
      if (nzchar(input$characteristic_field %||% "") && nrow(sdrf) > 0L) {
        values <- safe_unique(sdrf[field_raw == input$characteristic_field, value])
      }
      shiny::updateSelectizeInput(session, "characteristic_values", choices = values, selected = character(), server = TRUE)
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$reset, {
      shiny::updateTextAreaInput(session, "query", value = "")
      shiny::updateSelectInput(session, "match_mode", selected = "all")
      shiny::updateCheckboxInput(session, "use_synonyms", value = TRUE)
      shiny::updateCheckboxInput(session, "include_study_context", value = TRUE)
      for (id2 in c("targets", "mechanisms", "agents", "controls", "model_presets", "organisms", "protocol_types", "characteristic_values", "study_ids")) {
        shiny::updateSelectizeInput(session, id2, selected = character())
      }
      shiny::updateSelectizeInput(session, "characteristic_field", selected = "")
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
        protocol_types = input$protocol_types %||% character(),
        characteristic_field = input$characteristic_field %||% "",
        characteristic_values = input$characteristic_values %||% character(),
        study_ids = input$study_ids %||% character()
      )
    }, ignoreInit = TRUE)

    result <- shiny::eventReactive(input$run, {
      ae <- ae_data()
      validation <- validate_ae_bundle(ae)
      if (!is.null(validation)) stop(validation)
      run_ae_search(ae, index(), filters())
    }, ignoreInit = TRUE)

    output$summary <- shiny::renderUI({
      if (input$run == 0L) {
        return(shiny::div(class = "empty-state", "Choose filters and click Run search."))
      }
      res <- result()
      summary_cards(as.list(res$counts))
    })

    output$studies <- DT::renderDT({ shiny::req(result()); result_table(result()$tables$ArrayExpress_studies) }, server = TRUE)
    output$samples <- DT::renderDT({ shiny::req(result()); result_table(result()$tables$ArrayExpress_samples) }, server = TRUE)
    output$sdrf <- DT::renderDT({ shiny::req(result()); result_table(result()$tables$ArrayExpress_SDRF_long) }, server = TRUE)
    output$protocols <- DT::renderDT({ shiny::req(result()); result_table(result()$tables$ArrayExpress_protocols) }, server = TRUE)
    output$refs <- DT::renderDT({ shiny::req(result()); result_table(result()$tables$ArrayExpress_sample_protocol_refs) }, server = TRUE)
    output$idf <- DT::renderDT({ shiny::req(result()); result_table(result()$tables$ArrayExpress_IDF_long) }, server = TRUE)
    output$idf_summary <- DT::renderDT({ shiny::req(result()); result_table(result()$tables$ArrayExpress_IDF_summary) }, server = TRUE)
    output$fields <- DT::renderDT({ shiny::req(result()); result_table(result()$tables$ArrayExpress_SDRF_field_summary) }, server = TRUE)
    output$errors <- DT::renderDT({ shiny::req(result()); result_table(result()$tables$ArrayExpress_download_errors) }, server = TRUE)

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
        "Protocol types: ", compact_filter_summary(filters()$protocol_types), "\n",
        "Characteristic field: ", compact_filter_summary(filters()$characteristic_field), "\n",
        "Characteristic values: ", compact_filter_summary(filters()$characteristic_values), "\n",
        "Study IDs: ", compact_filter_summary(filters()$study_ids)
      )
    })

    output$download_zip <- shiny::downloadHandler(
      filename = function() paste0("ArrayExpress_NMD_subset_", format(Sys.Date(), "%Y%m%d"), ".zip"),
      content = function(file) {
        shiny::req(result(), filters())
        write_result_zip(
          result(),
          platform = "ArrayExpress",
          filters = filters(),
          interpretation = result()$interpretation,
          output_file = file
        )
      }
    )
  })
}

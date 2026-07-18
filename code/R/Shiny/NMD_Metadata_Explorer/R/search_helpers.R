regex_escape <- function(x) {
  x <- as.character(x)
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x, perl = TRUE)
}

or_pattern <- function(values, word_boundaries = FALSE) {
  values <- unique(trimws(as.character(values)))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0L) return("(?!)")

  escaped <- regex_escape(values)
  if (isTRUE(word_boundaries)) {
    escaped <- paste0("\\b", escaped, "\\b")
  }
  paste0("(?:", paste(escaped, collapse = "|"), ")")
}

concept_pattern <- function(selected, dictionary) {
  selected <- intersect(as.character(selected), names(dictionary))
  if (length(selected) == 0L) return(NULL)
  or_pattern(unlist(dictionary[selected], use.names = FALSE))
}

split_multivalue <- function(x, separators = ";") {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (length(x) == 0L) return(character())
  values <- unlist(strsplit(x, separators, perl = TRUE), use.names = FALSE)
  sort(unique(trimws(values[nzchar(trimws(values))])))
}

combine_row_text <- function(dt, columns = NULL) {
  if (is.null(dt) || nrow(dt) == 0L) return(character())
  if (is.null(columns)) columns <- names(dt)
  columns <- intersect(columns, names(dt))
  if (length(columns) == 0L) return(rep("", nrow(dt)))

  pieces <- lapply(dt[, ..columns], function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  })
  do.call(paste, c(pieces, sep = " | "))
}

collapse_unique_text <- function(x, max_values = Inf, sep = " | ") {
  x <- trimws(as.character(x))
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x) == 0L) return(NA_character_)
  if (is.finite(max_values) && length(x) > max_values) {
    extra <- length(x) - max_values
    x <- c(x[seq_len(max_values)], paste0("[+", extra, " additional value(s)]"))
  }
  paste(x, collapse = sep)
}

smart_query_groups <- function(query, use_synonyms = TRUE, mode = "all") {
  query <- trimws(as.character(query %||% ""))
  if (!nzchar(query)) {
    return(list(groups = character(), labels = character(), display = "No free-text query"))
  }

  if (identical(mode, "exact")) {
    return(list(
      groups = regex_escape(query),
      labels = paste0("Exact phrase: ", query),
      display = paste0("Exact phrase: ", query)
    ))
  }

  remaining <- query
  groups <- character()
  labels <- character()

  if (isTRUE(use_synonyms)) {
    query_lower <- tolower(query)
    for (label in names(SMART_QUERY_CONCEPTS)) {
      terms <- SMART_QUERY_CONCEPTS[[label]]
      matched_terms <- terms[vapply(
        terms,
        function(term) grepl(tolower(term), query_lower, fixed = TRUE),
        logical(1)
      )]

      if (length(matched_terms) > 0L) {
        groups <- c(groups, or_pattern(terms))
        labels <- c(labels, paste0(label, " [", paste(terms, collapse = ", "), "]"))
        for (term in matched_terms) {
          remaining <- gsub(term, " ", remaining, ignore.case = TRUE, fixed = TRUE)
        }
      }
    }
  }

  remaining <- gsub("[\"']", " ", remaining)
  tokens <- unlist(strsplit(remaining, "[,;[:space:]]+", perl = TRUE), use.names = FALSE)
  tokens <- trimws(tokens)
  stop_words <- c(
    "and", "or", "with", "in", "of", "the", "a", "an", "to", "for",
    "from", "study", "studies", "sample", "samples", "cell", "cells",
    "dataset", "datasets", "experiment", "experiments", "show", "find"
  )
  tokens <- tokens[nzchar(tokens) & !tolower(tokens) %in% stop_words]

  if (length(tokens) > 0L) {
    groups <- c(groups, vapply(tokens, regex_escape, character(1)))
    labels <- c(labels, paste0("Literal term: ", tokens))
  }

  if (length(groups) == 0L) {
    groups <- regex_escape(query)
    labels <- paste0("Literal query: ", query)
  }

  list(
    groups = unique(groups),
    labels = unique(labels),
    display = paste(unique(labels), collapse = if (identical(mode, "all")) " AND " else " OR ")
  )
}

match_search_groups <- function(text, groups, mode = "all") {
  text <- as.character(text)
  text[is.na(text)] <- ""
  if (length(text) == 0L) return(logical())
  if (length(groups) == 0L) return(rep(TRUE, length(text)))

  hit_matrix <- vapply(
    groups,
    function(pattern) grepl(pattern, text, ignore.case = TRUE, perl = TRUE),
    logical(length(text))
  )

  if (is.null(dim(hit_matrix))) hit_matrix <- matrix(hit_matrix, ncol = 1L)
  if (identical(mode, "any")) {
    rowSums(hit_matrix) > 0L
  } else {
    rowSums(hit_matrix) == ncol(hit_matrix)
  }
}

match_selected_concepts <- function(text, selected, dictionary) {
  if (length(selected) == 0L) return(rep(TRUE, length(text)))
  pattern <- concept_pattern(selected, dictionary)
  grepl(pattern, as.character(text), ignore.case = TRUE, perl = TRUE)
}

safe_unique <- function(x) {
  x <- trimws(as.character(x))
  sort(unique(x[!is.na(x) & nzchar(x)]))
}

make_field_summary <- function(metadata) {
  metadata <- data.table::as.data.table(data.table::copy(metadata))
  required <- c("field", "field_raw", "value")
  if (nrow(metadata) == 0L || !all(required %in% names(metadata))) {
    return(data.table::data.table())
  }

  study_col <- if ("GSE_ID" %in% names(metadata)) "GSE_ID" else if ("Study_ID" %in% names(metadata)) "Study_ID" else NULL
  sample_col <- if ("Sample_ID" %in% names(metadata)) "Sample_ID" else NULL

  metadata[, value := as.character(value)]
  out <- metadata[, .(
    n_studies = if (!is.null(study_col)) data.table::uniqueN(get(study_col)) else NA_integer_,
    n_samples = if (!is.null(sample_col) && !is.null(study_col)) {
      data.table::uniqueN(paste(get(study_col), get(sample_col), sep = "::"))
    } else if (!is.null(sample_col)) {
      data.table::uniqueN(get(sample_col))
    } else {
      NA_integer_
    },
    n_nonempty_values = .N,
    n_unique_values = data.table::uniqueN(value),
    example_values = collapse_unique_text(value, max_values = 5L, sep = "; ")
  ), by = .(field, field_raw)]

  data.table::setorder(out, -n_studies, -n_samples, field)
  out
}

compact_filter_summary <- function(values) {
  values <- as.character(values)
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0L) "Any" else paste(values, collapse = ", ")
}

format_filter_log <- function(platform, filters, interpretation, counts) {
  lines <- c(
    paste0("NMD Metadata Explorer - ", platform, " subset"),
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "Search interpretation:",
    interpretation %||% "No free-text query",
    "",
    "Filters:"
  )

  for (nm in names(filters)) {
    value <- filters[[nm]]
    if (length(value) == 0L || all(is.na(value)) || identical(value, "")) value <- "Any"
    lines <- c(lines, paste0("- ", nm, ": ", paste(value, collapse = ", ")))
  }

  lines <- c(lines, "", "Exported table row counts:")
  for (nm in names(counts)) {
    lines <- c(lines, paste0("- ", nm, ": ", counts[[nm]]))
  }
  paste(lines, collapse = "\n")
}

write_result_zip <- function(result, platform, filters, interpretation, output_file) {
  export_dir <- tempfile(paste0("nmd_", tolower(platform), "_export_"))
  dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(export_dir, recursive = TRUE, force = TRUE), add = TRUE)

  tables <- result$tables
  counts <- integer()
  for (nm in names(tables)) {
    tab <- tables[[nm]]
    if (is.null(tab)) next
    tab <- data.table::as.data.table(tab)
    counts[[nm]] <- nrow(tab)
    if (ncol(tab) == 0L) next
    filename <- paste0(gsub("[^A-Za-z0-9_]+", "_", nm), ".tsv")
    data.table::fwrite(tab, file.path(export_dir, filename), sep = "\t", quote = FALSE, na = "")
  }

  writeLines(
    format_filter_log(platform, filters, interpretation, counts),
    file.path(export_dir, "SEARCH_FILTERS_AND_COUNTS.txt"),
    useBytes = TRUE
  )

  zip::zipr(
    zipfile = output_file,
    files = list.files(export_dir, full.names = FALSE),
    root = export_dir
  )
}

result_table <- function(data, page_length = 15L) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (ncol(data) == 0L) {
    data <- data.frame(Status = "No matching rows or this optional table was not loaded.", stringsAsFactors = FALSE)
  }
  DT::datatable(
    data,
    rownames = FALSE,
    filter = "top",
    extensions = "Buttons",
    escape = TRUE,
    options = list(
      dom = "Bfrtip",
      buttons = c("copy", "csv", "colvis"),
      pageLength = page_length,
      lengthMenu = c(10, 15, 25, 50, 100),
      scrollX = TRUE,
      autoWidth = TRUE,
      deferRender = TRUE,
      searchHighlight = TRUE
    ),
    class = "compact stripe hover"
  )
}

summary_cards <- function(items) {
  shiny::fluidRow(lapply(names(items), function(label) {
    shiny::column(
      width = max(2L, floor(12L / length(items))),
      shiny::div(
        class = "summary-card",
        shiny::div(class = "summary-value", format(items[[label]], big.mark = ",")),
        shiny::div(class = "summary-label", label)
      )
    )
  }))
}

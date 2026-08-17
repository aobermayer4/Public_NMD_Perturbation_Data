# =============================================================================
# ArrayExpress / BioStudies discovery + generic MAGE-TAB parser + NMD screen
# =============================================================================
# This file consolidates three older layers:
#   1. BioStudies ArrayExpress search/pagination
#   2. Generic MAGE-TAB download/read helpers (formerly buried in AML code)
#   3. NMD-specific pre-screening and sample/study classification
#
# Public functions used by the driver:
#   biostudies_search_many()
#   ae_search_nmd()
#   ae_screen_nmd_accessions()
#   ae_write_nmd_screening_results()
# =============================================================================

required_ae_packages <- c(
  "ArrayExpress", "httr", "jsonlite", "data.table",
  "dplyr", "tidyr", "purrr", "stringr", "tibble", "xml2"
)

missing_ae_packages <- required_ae_packages[
  !vapply(required_ae_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_ae_packages) > 0L) {
  stop(
    "Install the following packages before running the ArrayExpress workflow: ",
    paste(missing_ae_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(ArrayExpress)
  library(httr)
  library(jsonlite)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
})

# =============================================================================
# A. BioStudies ArrayExpress discovery
# =============================================================================
# =============================================================================
# Search the BioStudies ArrayExpress collection through its REST API
#
# This script retrieves every page of matching studies and returns a clean
# accession table. It does NOT download study data or MAGE-TAB files.
#
# API endpoint:
#   https://www.ebi.ac.uk/biostudies/api/v1/search
#
# Collection filter:
#   collection=arrayexpress
#
# Required packages:
#   httr, jsonlite, data.table
#
# Version 1.0.1: corrected two invalid split double-bracket expressions.
# Written without the native |> pipe for compatibility with R 4.0.5.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------

required_packages <- c(
  "httr",
  "jsonlite",
  "data.table"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(data.table)
})


# -----------------------------------------------------------------------------
# 1. General helpers
# -----------------------------------------------------------------------------

biostudies_api_base <- function() {
  "https://www.ebi.ac.uk/biostudies/api/v1"
}


biostudies_first_nonempty <- function(x) {
  x <- as.character(x)
  x <- trimws(x)

  keep <- !is.na(x) & nzchar(x)

  if (!any(keep)) {
    return(NA_character_)
  }

  x[which(keep)[1]]
}


biostudies_safe_collapse <- function(x) {
  x <- unique(
    trimws(
      as.character(x)
    )
  )

  x <- x[
    !is.na(x) &
      nzchar(x)
  ]

  if (length(x) == 0L) {
    return(NA_character_)
  }

  paste(x, collapse = "; ")
}


biostudies_empty_hits <- function() {
  data.table::data.table(
    accession = character(),
    type = character(),
    title = character(),
    author = character(),
    links = integer(),
    files = integer(),
    release_date = character(),
    views = integer(),
    isPublic = logical()
  )
}


# -----------------------------------------------------------------------------
# 2. Retrieve one search-results page
# -----------------------------------------------------------------------------

biostudies_search_page <- function(
    query,
    page = 1L,
    page_size = 100L,
    collection = "arrayexpress",
    sort_by = "relevance",
    sort_order = "descending",
    timeout_seconds = 90,
    retries = 4L,
    verbose = FALSE) {

  page <- as.integer(page)
  page_size <- as.integer(page_size)

  if (page < 1L) {
    stop("page must be at least 1.")
  }

  if (page_size < 1L) {
    stop("page_size must be at least 1.")
  }

  endpoint <- paste0(
    biostudies_api_base(),
    "/search"
  )

  request_query <- list(
    collection = collection,
    page = page,
    pageSize = page_size,
    sortBy = sort_by,
    sortOrder = sort_order
  )

  if (!is.null(query) &&
      length(query) == 1L &&
      !is.na(query) &&
      nzchar(trimws(query))) {

    request_query$query <- query
  }

  if (verbose) {
    message(
      "BioStudies search page ",
      page,
      ": ",
      ifelse(
        is.null(query),
        "<all studies>",
        query
      )
    )
  }

  response <- httr::RETRY(
    verb = "GET",
    url = endpoint,
    query = request_query,
    httr::accept_json(),
    httr::user_agent(
      "R BioStudies ArrayExpress pagination utility"
    ),
    httr::timeout(timeout_seconds),
    times = retries,
    pause_min = 1,
    pause_cap = 10,
    terminate_on = c(
      400L,
      401L,
      403L,
      404L
    ),
    quiet = !verbose
  )

  httr::stop_for_status(
    response,
    task = "search BioStudies"
  )

  response_text <- httr::content(
    response,
    as = "text",
    encoding = "UTF-8"
  )

  payload <- jsonlite::fromJSON(
    response_text,
    simplifyVector = TRUE,
    simplifyDataFrame = TRUE,
    flatten = TRUE
  )

  hits <- payload$hits

  if (is.null(hits) ||
      length(hits) == 0L) {

    hits <- biostudies_empty_hits()
  } else {
    hits <- data.table::as.data.table(
      hits
    )
  }

  list(
    page = as.integer(
      payload$page
    ),
    page_size = as.integer(
      payload$pageSize
    ),
    total_hits = as.numeric(
      payload$totalHits
    ),
    is_total_hits_exact =
      as.logical(
        payload$isTotalHitsExact
      ),
    sort_by =
      as.character(
        payload$sortBy
      ),
    sort_order =
      as.character(
        payload$sortOrder
      ),
    query_returned =
      as.character(
        payload$query
      ),
    facets =
      payload$facets,
    hits = hits,
    request_url =
      httr::build_url(
        httr::parse_url(
          xml2::url_absolute(
            endpoint,
            endpoint
          )
        )
      )
  )
}


# -----------------------------------------------------------------------------
# 3. Retrieve every page for one query
# -----------------------------------------------------------------------------

biostudies_search_all <- function(
    query,
    query_name = NULL,
    page_size = 100L,
    collection = "arrayexpress",
    sort_by = "relevance",
    sort_order = "descending",
    max_pages = Inf,
    max_hits = Inf,
    pause_seconds = 0.15,
    verbose = TRUE) {

  query_is_empty <- (
    is.null(query) ||
      length(query) == 0L ||
      is.na(query[1]) ||
      !nzchar(trimws(as.character(query[1])))
  )

  if (is.null(query_name) ||
      length(query_name) == 0L ||
      is.na(query_name[1]) ||
      !nzchar(trimws(as.character(query_name[1])))) {

    query_name <- if (query_is_empty) {
      "All ArrayExpress studies"
    } else {
      as.character(query[1])
    }
  }

  query_text <- if (query_is_empty) {
    NA_character_
  } else {
    as.character(query[1])
  }

  first_page <- biostudies_search_page(
    query = query,
    page = 1L,
    page_size = page_size,
    collection = collection,
    sort_by = sort_by,
    sort_order = sort_order,
    verbose = verbose
  )

  total_hits <- first_page$total_hits

  returned_page_size <- first_page$page_size

  if (is.na(returned_page_size) ||
      returned_page_size < 1L) {

    returned_page_size <- page_size
  }

  total_pages <- if (
    is.na(total_hits) ||
      !is.finite(total_hits)
  ) {
    Inf
  } else {
    ceiling(
      total_hits /
        returned_page_size
    )
  }

  pages_to_request <- min(
    total_pages,
    max_pages
  )

  page_tables <- list(
    first_page$hits
  )

  collected <- nrow(
    first_page$hits
  )

  if (verbose) {
    message(
      "Total matching studies reported by API: ",
      format(
        total_hits,
        big.mark = ",",
        scientific = FALSE
      )
    )
  }

  if (is.finite(max_hits) &&
      collected >= max_hits) {

    page_tables[[1]] <-
      page_tables[[1]][
        seq_len(max_hits)
      ]
  } else if (
    is.infinite(pages_to_request) ||
      pages_to_request >= 2L
  ) {

    page <- 2L

    repeat {
      if (is.finite(pages_to_request) &&
          page > pages_to_request) {
        break
      }

      if (is.finite(max_hits) &&
          collected >= max_hits) {
        break
      }

      if (pause_seconds > 0) {
        Sys.sleep(pause_seconds)
      }

      current <- biostudies_search_page(
        query = query,
        page = page,
        page_size = page_size,
        collection = collection,
        sort_by = sort_by,
        sort_order = sort_order,
        verbose = FALSE
      )

      if (nrow(current$hits) == 0L) {
        break
      }

      page_tables[[length(page_tables) + 1L]] <-
        current$hits

      collected <- collected +
        nrow(current$hits)

      if (verbose &&
          (
            page %% 10L == 0L ||
              (
                is.finite(pages_to_request) &&
                  page == pages_to_request
              )
          )) {

        message(
          "Retrieved page ",
          page,
          if (
            is.finite(pages_to_request)
          ) {
            paste0(
              " of ",
              pages_to_request
            )
          } else {
            ""
          },
          "; ",
          format(
            collected,
            big.mark = ","
          ),
          " rows collected."
        )
      }

      page <- page + 1L
    }
  }

  hits <- data.table::rbindlist(
    page_tables,
    fill = TRUE,
    use.names = TRUE
  )

  if (is.finite(max_hits) &&
      nrow(hits) > max_hits) {

    hits <- hits[
      seq_len(max_hits)
    ]
  }

  hits[
    ,
    `:=`(
      matched_query_name =
        as.character(query_name[1]),
      matched_query =
        query_text
    )
  ]

  preferred_order <- c(
    "accession",
    "type",
    "title",
    "author",
    "release_date",
    "files",
    "links",
    "views",
    "isPublic",
    "matched_query_name",
    "matched_query"
  )

  data.table::setcolorder(
    hits,
    c(
      intersect(
        preferred_order,
        names(hits)
      ),
      setdiff(
        names(hits),
        preferred_order
      )
    )
  )

  list(
    query = query_text,
    query_name =
      as.character(query_name[1]),
    collection = collection,
    total_hits_reported =
      total_hits,
    is_total_hits_exact =
      first_page$is_total_hits_exact,
    pages_retrieved =
      length(page_tables),
    hits = hits,
    facets = first_page$facets
  )
}


# -----------------------------------------------------------------------------
# 4. Run multiple searches and make a deduplicated accession table
# -----------------------------------------------------------------------------

biostudies_search_many <- function(
    queries,
    page_size = 100L,
    collection = "arrayexpress",
    sort_by = "relevance",
    sort_order = "descending",
    max_pages = Inf,
    max_hits_per_query = Inf,
    pause_seconds = 0.15,
    accession_pattern =
      "^E-(MTAB|GEOD)-[0-9]+$",
    verbose = TRUE) {

  queries <- as.character(
    queries
  )

  if (is.null(names(queries))) {
    names(queries) <- paste0(
      "query_",
      seq_along(queries)
    )
  }

  query_results <- vector(
    "list",
    length(queries)
  )

  names(query_results) <-
    names(queries)

  for (i in seq_along(queries)) {
    if (verbose) {
      message(
        "\n========== ",
        names(queries)[i],
        " =========="
      )
    }

    query_results[[i]] <-
      biostudies_search_all(
        query = queries[[i]],
        query_name =
          names(queries)[i],
        page_size = page_size,
        collection = collection,
        sort_by = sort_by,
        sort_order = sort_order,
        max_pages = max_pages,
        max_hits =
          max_hits_per_query,
        pause_seconds =
          pause_seconds,
        verbose = verbose
      )
  }

  all_hits <- data.table::rbindlist(
    lapply(
      query_results,
      function(x) x$hits
    ),
    fill = TRUE,
    use.names = TRUE
  )

  all_hits[
    ,
    accession :=
      toupper(
        trimws(
          as.character(
            accession
          )
        )
      )
  ]

  all_hits <- all_hits[
    !is.na(accession) &
      grepl(
        accession_pattern,
        accession
      )
  ]

  unique_hits <- all_hits[
    ,
    .(
      type =
        biostudies_first_nonempty(
          type
        ),
      title =
        biostudies_first_nonempty(
          title
        ),
      author =
        biostudies_first_nonempty(
          author
        ),
      release_date =
        biostudies_first_nonempty(
          release_date
        ),
      files = suppressWarnings(
        max(
          as.numeric(files),
          na.rm = TRUE
        )
      ),
      links = suppressWarnings(
        max(
          as.numeric(links),
          na.rm = TRUE
        )
      ),
      views = suppressWarnings(
        max(
          as.numeric(views),
          na.rm = TRUE
        )
      ),
      isPublic = any(
        isPublic %in% TRUE,
        na.rm = TRUE
      ),
      matched_query_names =
        biostudies_safe_collapse(
          matched_query_name
        ),
      matched_queries =
        biostudies_safe_collapse(
          matched_query
        ),
      n_query_matches =
        data.table::uniqueN(
          matched_query_name
        )
    ),
    by = accession
  ]

  # max(..., na.rm = TRUE) returns -Inf when every value is missing.
  numeric_columns <- c(
    "files",
    "links",
    "views"
  )

  for (column_name in numeric_columns) {
    column_values <- unique_hits[[column_name]]
    infinite_rows <- which(is.infinite(column_values))

    if (length(infinite_rows) > 0L) {
      data.table::set(
        unique_hits,
        i = infinite_rows,
        j = column_name,
        value = NA_real_
      )
    }
  }

  data.table::setorderv(
    unique_hits,
    cols = c(
      "n_query_matches",
      "release_date",
      "accession"
    ),
    order = c(-1L, -1L, 1L),
    na.last = TRUE
  )

  list(
    queries = queries,
    per_query = query_results,
    all_query_hits = all_hits,
    unique_studies = unique_hits,
    accessions =
      unique_hits$accession
  )
}


# -----------------------------------------------------------------------------
# 5. Retrieve detailed JSON for one accession
# -----------------------------------------------------------------------------

biostudies_get_study_json <- function(
    accession,
    info = FALSE,
    timeout_seconds = 90,
    retries = 4L) {

  accession <- toupper(
    trimws(
      as.character(accession)
    )
  )

  path <- paste0(
    biostudies_api_base(),
    "/studies/",
    utils::URLencode(
      accession,
      reserved = TRUE
    ),
    if (isTRUE(info)) {
      "/info"
    } else {
      ""
    }
  )

  response <- httr::RETRY(
    verb = "GET",
    url = path,
    httr::accept_json(),
    httr::user_agent(
      "R BioStudies study metadata utility"
    ),
    httr::timeout(
      timeout_seconds
    ),
    times = retries,
    pause_min = 1,
    pause_cap = 10,
    terminate_on = c(
      400L,
      401L,
      403L,
      404L
    ),
    quiet = TRUE
  )

  httr::stop_for_status(
    response,
    task = paste(
      "retrieve",
      accession
    )
  )

  jsonlite::fromJSON(
    httr::content(
      response,
      as = "text",
      encoding = "UTF-8"
    ),
    simplifyVector = TRUE,
    flatten = FALSE
  )
}


# -----------------------------------------------------------------------------
# 6. Write search results
# -----------------------------------------------------------------------------

write_biostudies_search_results <- function(
    search_result,
    output_prefix =
      "ArrayExpress_BioStudies_search") {

  data.table::fwrite(
    search_result$unique_studies,
    file = paste0(
      output_prefix,
      "_unique_studies.tsv"
    ),
    sep = "\t",
    quote = TRUE,
    na = ""
  )

  data.table::fwrite(
    search_result$all_query_hits,
    file = paste0(
      output_prefix,
      "_all_query_hits.tsv"
    ),
    sep = "\t",
    quote = TRUE,
    na = ""
  )

  writeLines(
    search_result$accessions,
    con = paste0(
      output_prefix,
      "_accessions.txt"
    ),
    useBytes = TRUE
  )

  saveRDS(
    search_result,
    file = paste0(
      output_prefix,
      "_complete_result.rds"
    )
  )

  invisible(
    c(
      unique_studies = paste0(
        output_prefix,
        "_unique_studies.tsv"
      ),
      all_query_hits = paste0(
        output_prefix,
        "_all_query_hits.tsv"
      ),
      accessions = paste0(
        output_prefix,
        "_accessions.txt"
      ),
      complete_result = paste0(
        output_prefix,
        "_complete_result.rds"
      )
    )
  )
}


# =============================================================================
# B. Generic ArrayExpress / MAGE-TAB helpers
# =============================================================================
# These functions were extracted from the older AML screening script because
# they are repository plumbing, not AML biology. Their ae_core_* names make
# that distinction explicit.
ae_core_safe_utf8 <- function(x) {
  x <- as.character(x)
  original_na <- is.na(x)
  out <- rep(NA_character_, length(x))

  keep <- !original_na

  if (any(keep)) {
    out[keep] <- suppressWarnings(
      iconv(
        x[keep],
        from = "",
        to = "UTF-8",
        sub = NA
      )
    )

    bad <- keep & is.na(out)

    if (any(bad)) {
      out[bad] <- suppressWarnings(
        iconv(
          x[bad],
          from = "CP1252",
          to = "UTF-8",
          sub = NA
        )
      )
    }

    bad <- keep & is.na(out)

    if (any(bad)) {
      out[bad] <- suppressWarnings(
        iconv(
          x[bad],
          from = "latin1",
          to = "UTF-8",
          sub = NA
        )
      )
    }

    bad <- keep & is.na(out)

    if (any(bad)) {
      out[bad] <- suppressWarnings(
        iconv(
          x[bad],
          from = "",
          to = "UTF-8",
          sub = "?"
        )
      )
    }
  }

  out[original_na] <- NA_character_
  out
}


ae_core_clean_text <- function(x) {
  x <- ae_core_safe_utf8(x)
  x <- stringr::str_squish(x)

  x[
    !is.na(x) &
      tolower(x) %in% c(
        "", "na", "n/a", "null", "none",
        "not applicable", "unknown"
      )
  ] <- NA_character_

  x
}


ae_core_normalize_field <- function(x) {
  x <- ae_core_clean_text(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}


ae_core_collapse_unique <- function(
    x,
    sep = "; ",
    max_values = Inf) {

  x <- ae_core_clean_text(x)
  x <- unique(x[!is.na(x) & nzchar(x)])

  if (length(x) == 0L) {
    return(NA_character_)
  }

  if (is.finite(max_values) && length(x) > max_values) {
    n_extra <- length(x) - max_values
    x <- c(
      x[seq_len(max_values)],
      paste0("[+", n_extra, " additional value(s)]")
    )
  }

  paste(x, collapse = sep)
}


ae_core_first_nonempty <- function(x) {
  x <- ae_core_clean_text(x)
  x <- x[!is.na(x) & nzchar(x)]

  if (length(x) == 0L) {
    return(NA_character_)
  }

  x[1]
}


ae_core_dictionary_hits <- function(text, dictionary) {
  text <- ae_core_clean_text(text)
  text <- paste(text[!is.na(text)], collapse = " ; ")

  if (!nzchar(text)) {
    return(character())
  }

  hit <- vapply(
    dictionary,
    function(pattern) {
      stringr::str_detect(
        text,
        stringr::regex(pattern, ignore_case = TRUE)
      )
    },
    logical(1)
  )

  names(dictionary)[hit]
}


ae_core_dictionary_hit_text <- function(text, dictionary) {
  hits <- ae_core_dictionary_hits(text, dictionary)

  if (length(hits) == 0L) {
    return(NA_character_)
  }

  paste(hits, collapse = "; ")
}


ae_core_has_dictionary_hit <- function(text, dictionary) {
  length(ae_core_dictionary_hits(text, dictionary)) > 0L
}


ae_core_retry <- function(
    expression_function,
    n = 3L,
    wait_seconds = 4,
    verbose = TRUE) {

  last_error <- NULL

  for (attempt in seq_len(n)) {
    answer <- tryCatch(
      expression_function(),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )

    if (!is.null(answer)) {
      return(answer)
    }

    if (verbose && !is.null(last_error)) {
      message(
        "Attempt ",
        attempt,
        "/",
        n,
        " failed: ",
        conditionMessage(last_error)
      )
    }

    if (attempt < n) {
      Sys.sleep(wait_seconds)
    }
  }

  if (!is.null(last_error)) {
    stop(last_error)
  }

  stop("Operation failed without an error object.")
}


# -----------------------------------------------------------------------------
# 3. MAGE-TAB metadata download
# -----------------------------------------------------------------------------

ae_core_download_file <- function(
    url,
    destination,
    overwrite = FALSE) {

  if (file.exists(destination) && !overwrite) {
    return(
      normalizePath(
        destination,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  dir.create(
    dirname(destination),
    recursive = TRUE,
    showWarnings = FALSE
  )

  ae_core_retry(
    function() {
      status <- suppressWarnings(
        utils::download.file(
          url = URLencode(url),
          destfile = destination,
          mode = "wb",
          quiet = TRUE,
          method = "libcurl"
        )
      )

      if (!identical(status, 0L) ||
          !file.exists(destination) ||
          file.info(destination)$size == 0) {
        stop("Download failed or returned an empty file: ", url)
      }

      normalizePath(
        destination,
        winslash = "/",
        mustWork = TRUE
      )
    }
  )
}


ae_core_download_with_arrayexpress <- function(
    accession,
    study_dir,
    overwrite = FALSE,
    verbose = TRUE) {

  getae_formals <- names(formals(ArrayExpress::getAE))

  args <- list(
    accession,
    path = study_dir,
    type = "mage",
    extract = TRUE
  )

  # The first argument has been called "input" in some older releases and
  # "accession" in newer releases. Passing it positionally avoids that drama.
  if ("overwrite" %in% getae_formals) {
    args$overwrite <- overwrite
  }

  if (verbose) {
    message(
      "Trying ArrayExpress::getAE(type = 'mage') for ",
      accession
    )
  }

  result <- do.call(
    ArrayExpress::getAE,
    args
  )

  if (is.null(result)) {
    stop("ArrayExpress::getAE() returned NULL.")
  }

  result
}


ae_core_biostudies_json <- function(accession) {
  study_url <- paste0(
    "https://www.ebi.ac.uk/biostudies/api/v1/studies/",
    accession
  )

  list(
    study = ae_core_retry(
      function() {
        jsonlite::fromJSON(
          study_url,
          simplifyVector = TRUE
        )
      }
    ),
    info = ae_core_retry(
      function() {
        jsonlite::fromJSON(
          paste0(study_url, "/info"),
          simplifyVector = TRUE
        )
      }
    )
  )
}


# Convert the PageTab study JSON and /info JSON into a small field/value table.
# This is much cheaper than downloading and expanding the MAGE-TAB archive.
ae_core_flatten_biostudies_json <- function(
    biostudies_json,
    accession) {

  flatten_one <- function(x, source_name) {
    flat <- unlist(
      x,
      recursive = TRUE,
      use.names = TRUE
    )

    if (length(flat) == 0L) {
      return(
        tibble::tibble(
          Study_ID = character(),
          json_source = character(),
          field = character(),
          value = character()
        )
      )
    }

    field_names <- names(flat)

    if (is.null(field_names)) {
      field_names <- paste0(
        source_name,
        "_value_",
        seq_along(flat)
      )
    }

    tibble::tibble(
      Study_ID = accession,
      json_source = source_name,
      field = ae_core_safe_utf8(field_names),
      value = ae_core_clean_text(
        as.character(flat)
      )
    ) %>%
      dplyr::filter(
        !is.na(.data$value),
        nzchar(.data$value)
      )
  }

  dplyr::bind_rows(
    flatten_one(
      biostudies_json$study,
      "study"
    ),
    flatten_one(
      biostudies_json$info,
      "info"
    )
  ) %>%
    dplyr::distinct()
}


ae_core_extract_prescreen_title <- function(flat_json) {
  if (nrow(flat_json) == 0L) {
    return(NA_character_)
  }

  title_rows <- stringr::str_detect(
    flat_json$field,
    stringr::regex(
      "(^|\\.)title($|\\.)|section.*title|attributes.*title",
      ignore_case = TRUE
    )
  )

  title <- ae_core_first_nonempty(
    flat_json$value[title_rows]
  )

  if (is.na(title)) {
    title <- ae_core_first_nonempty(
      flat_json$value
    )
  }

  title
}


# Fast gate used before any MAGE-TAB download.
#
# The gate intentionally does not require PBMC/blood/bone-marrow terms because
# these are often only present in the SDRF. It can, however, reject studies that
# are clearly non-human when no human evidence is present.
ae_core_download_from_biostudies <- function(
    accession,
    study_dir,
    overwrite = FALSE,
    verbose = TRUE,
    biostudies_json = NULL) {

  if (verbose) {
    message(
      "Using BioStudies metadata fallback for ",
      accession
    )
  }

  json <- biostudies_json

  if (is.null(json)) {
    json <- ae_core_biostudies_json(accession)
  }

  all_strings <- unique(
    ae_core_safe_utf8(
      as.character(
        unlist(
          c(json$study, json$info),
          recursive = TRUE,
          use.names = FALSE
        )
      )
    )
  )

  magetab_paths <- all_strings[
    !is.na(all_strings) &
      stringr::str_detect(
        basename(all_strings),
        stringr::regex(
          "(idf|sdrf)(\\.txt)?$",
          ignore_case = TRUE
        )
      )
  ]

  magetab_paths <- unique(magetab_paths)

  if (length(magetab_paths) == 0L) {
    stop(
      "No IDF or SDRF paths were found for ",
      accession,
      "."
    )
  }

  ftp_link <- NA_character_

  if (!is.null(json$info$ftpLink)) {
    ftp_link <- as.character(json$info$ftpLink)[1]
  }

  if (is.na(ftp_link) || !nzchar(ftp_link)) {
    stop(
      "BioStudies did not provide an ftpLink for ",
      accession,
      "."
    )
  }

  ftp_link <- sub("/+$", "", ftp_link)

  downloaded <- vapply(
    magetab_paths,
    function(relative_path) {
      relative_path <- sub(
        "^/?Files/",
        "",
        relative_path,
        ignore.case = TRUE
      )

      file_url <- paste0(
        ftp_link,
        "/Files/",
        relative_path
      )

      destination <- file.path(
        study_dir,
        basename(relative_path)
      )

      if (verbose) {
        message(
          "Downloading metadata file: ",
          basename(relative_path)
        )
      }

      ae_core_download_file(
        url = file_url,
        destination = destination,
        overwrite = overwrite
      )
    },
    character(1)
  )

  list(
    path = study_dir,
    mageTabFiles = downloaded,
    download_method = "BioStudies API fallback"
  )
}


ae_core_download_magetab <- function(
    accession,
    base_dir = "ArrayExpress_NMD_Metadata",
    overwrite = FALSE,
    use_api_fallback = TRUE,
    verbose = TRUE,
    biostudies_json = NULL) {

  accession <- toupper(trimws(accession))

  if (!stringr::str_detect(
    accession,
    "^E-[A-Z0-9]+-[0-9]+$"
  )) {
    stop(
      "Invalid-looking ArrayExpress accession: ",
      accession
    )
  }

  study_dir <- file.path(base_dir, accession)

  dir.create(
    study_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  package_result <- tryCatch(
    ae_core_download_with_arrayexpress(
      accession = accession,
      study_dir = study_dir,
      overwrite = overwrite,
      verbose = verbose
    ),
    error = function(e) e
  )

  if (!inherits(package_result, "error")) {
    package_result$path <- study_dir
    package_result$download_method <-
      "ArrayExpress::getAE(type = 'mage')"
    return(package_result)
  }

  if (!use_api_fallback) {
    stop(
      "ArrayExpress package download failed: ",
      conditionMessage(package_result)
    )
  }

  if (verbose) {
    message(
      "Package download failed: ",
      conditionMessage(package_result)
    )
  }

  ae_core_download_from_biostudies(
    accession = accession,
    study_dir = study_dir,
    overwrite = overwrite,
    verbose = verbose,
    biostudies_json = biostudies_json
  )
}


# -----------------------------------------------------------------------------
# 4. IDF and SDRF parsing
# -----------------------------------------------------------------------------

ae_core_make_sample_id <- function(
    data,
    original_names) {

  preferred_fields <- c(
    "Sample Name",
    "Source Name",
    "Assay Name",
    "Extract Name",
    "Scan Name",
    "Comment[ENA_SAMPLE]",
    "Comment[BIOSAMPLE]"
  )

  sample_id <- rep(
    NA_character_,
    nrow(data)
  )

  for (field in preferred_fields) {
    matching_columns <- names(data)[
      tolower(original_names) == tolower(field)
    ]

    for (column in matching_columns) {
      candidate <- ae_core_clean_text(
        data[[column]]
      )

      replace <- (
        is.na(sample_id) |
          !nzchar(sample_id)
      ) &
        !is.na(candidate) &
        nzchar(candidate)

      sample_id[replace] <- candidate[replace]
    }
  }

  missing <- is.na(sample_id) | !nzchar(sample_id)

  sample_id[missing] <- paste0(
    "SDRF_row_",
    which(missing)
  )

  sample_id
}


ae_core_read_sdrf <- function(
    file,
    accession) {

  sdrf <- data.table::fread(
    file,
    sep = "\t",
    header = TRUE,
    fill = TRUE,
    quote = "",
    na.strings = c(
      "", "NA", "N/A", "null", "NULL"
    ),
    data.table = FALSE,
    check.names = FALSE,
    colClasses = "character",
    encoding = "UTF-8",
    showProgress = FALSE
  )

  sdrf[] <- lapply(
    sdrf,
    ae_core_clean_text
  )

  original_names <- names(sdrf)

  unique_names <- make.unique(
    original_names,
    sep = "__duplicate_"
  )

  names(sdrf) <- unique_names

  sample_id <- ae_core_make_sample_id(
    data = sdrf,
    original_names = original_names
  )

  wide <- tibble::as_tibble(sdrf) %>%
    dplyr::mutate(
      Study_ID = accession,
      Sample_ID = sample_id,
      SDRF_Row = dplyr::row_number(),
      SDRF_File = basename(file),
      .before = 1
    )

  name_map <- stats::setNames(
    original_names,
    unique_names
  )

  long <- wide %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(unique_names),
      names_to = "field_unique",
      values_to = "value",
      values_transform = list(
        value = as.character
      )
    ) %>%
    dplyr::mutate(
      field_raw = unname(
        name_map[.data$field_unique]
      ),
      field = ae_core_normalize_field(
        .data$field_raw
      ),
      value = ae_core_clean_text(
        .data$value
      )
    ) %>%
    dplyr::filter(
      !is.na(.data$value),
      nzchar(.data$value)
    ) %>%
    dplyr::select(
      .data$Study_ID,
      .data$Sample_ID,
      .data$SDRF_Row,
      .data$SDRF_File,
      .data$field_raw,
      .data$field,
      .data$value
    )

  list(
    wide = wide,
    long = long
  )
}


ae_core_read_idf <- function(
    file,
    accession) {

  idf <- data.table::fread(
    file,
    sep = "\t",
    header = FALSE,
    fill = TRUE,
    quote = "",
    blank.lines.skip = FALSE,
    na.strings = NULL,
    data.table = FALSE,
    colClasses = "character",
    encoding = "UTF-8",
    showProgress = FALSE
  )

  idf[] <- lapply(
    idf,
    ae_core_clean_text
  )

  keep_column <- vapply(
    idf,
    function(x) {
      any(!is.na(x) & nzchar(x))
    },
    logical(1)
  )

  idf <- idf[, keep_column, drop = FALSE]

  if (ncol(idf) < 2L) {
    stop(
      "IDF file has fewer than two populated columns: ",
      file
    )
  }

  names(idf) <- paste0(
    "V",
    seq_len(ncol(idf))
  )

  field_values <- idf[[1]]
  field_values[
    is.na(field_values) |
      !nzchar(field_values)
  ] <- "Unlabeled IDF field"

  idf_long <- purrr::map_dfr(
    seq_len(nrow(idf)),
    function(i) {
      values <- as.character(
        unlist(
          idf[i, -1, drop = FALSE],
          use.names = FALSE
        )
      )

      tibble::tibble(
        Study_ID = accession,
        IDF_File = basename(file),
        IDF_Row = i,
        field_raw = field_values[i],
        field = ae_core_normalize_field(
          field_values[i]
        ),
        value_index = seq_along(values),
        value = ae_core_clean_text(values)
      )
    }
  ) %>%
    dplyr::filter(
      !is.na(.data$value),
      nzchar(.data$value)
    )

  idf_long
}


# -----------------------------------------------------------------------------
# 5. File inventory
# -----------------------------------------------------------------------------


# Generic file-type dictionaries used only for metadata/file inventory.
AE_CORE_IDAT_TERMS <- c(
  "IDAT" = "\\.idat(\\.gz)?([?#].*)?$|(^|[^A-Za-z0-9])IDAT([^A-Za-z0-9]|$)"
)

AE_CORE_RAW_SEQUENCE_TERMS <- c(
  "FASTQ" = "\\.(fastq|fq)(\\.gz)?([?#].*)?$",
  "SRA file" = "\\.sra([?#].*)?$",
  "ENA/SRA run accession" =
    "(^|[^A-Za-z0-9])(ERR|SRR|DRR)[0-9]+([^A-Za-z0-9]|$)"
)

AE_CORE_ALIGNED_TERMS <- c(
  "BAM" = "\\.bam([?#].*)?$",
  "CRAM" = "\\.cram([?#].*)?$"
)

AE_CORE_PROCESSED_FILE_TERMS <- c(
  "Expression or count matrix" =
    "count[s]?[ _-]*matrix|expression[ _-]*matrix|normalized[ _-]*matrix",
  "Tabular file" = "\\.(csv|tsv|txt)(\\.gz)?([?#].*)?$",
  "Spreadsheet" = "\\.(xls|xlsx)([?#].*)?$",
  "Genomic interval or signal" =
    "\\.(bed|bedgraph|bigwig|bw)(\\.gz)?([?#].*)?$"
)
ae_core_classify_file_value <- function(value) {
  value <- ae_core_clean_text(value)

  if (is.na(value) || !nzchar(value)) {
    return("Other metadata")
  }

  if (ae_core_has_dictionary_hit(
    value,
    AE_CORE_IDAT_TERMS
  )) {
    return("Raw methylation array file")
  }

  if (ae_core_has_dictionary_hit(
    value,
    AE_CORE_RAW_SEQUENCE_TERMS
  )) {
    return("Raw sequencing file or run")
  }

  if (ae_core_has_dictionary_hit(
    value,
    AE_CORE_ALIGNED_TERMS
  )) {
    return("Aligned sequencing file")
  }

  if (ae_core_has_dictionary_hit(
    value,
    AE_CORE_PROCESSED_FILE_TERMS
  )) {
    return("Processed or derived file")
  }

  "Other file or accession"
}


ae_core_sdrf_file_inventory <- function(sdrf_long) {
  if (nrow(sdrf_long) == 0L) {
    return(
      tibble::tibble(
        Study_ID = character(),
        Sample_ID = character(),
        inventory_source = character(),
        field_raw = character(),
        value = character(),
        file_class = character()
      )
    )
  }

  possible_file_row <- stringr::str_detect(
    sdrf_long$field,
    stringr::regex(
      "file|uri|url|scan_name|ena_run|sra_run|fastq|raw_data|array_data|submitted_file|bam|cram|idat",
      ignore_case = TRUE
    )
  ) |
    stringr::str_detect(
      sdrf_long$value,
      stringr::regex(
        "\\.(fastq|fq|sra|bam|cram|idat|csv|tsv|txt|xls|xlsx|bed|bigwig|bw)(\\.gz)?([?#].*)?$|(^|[^A-Za-z0-9])(ERR|SRR|DRR)[0-9]+([^A-Za-z0-9]|$)",
        ignore_case = TRUE
      )
    )

  sdrf_long[possible_file_row, , drop = FALSE] %>%
    dplyr::transmute(
      Study_ID = .data$Study_ID,
      Sample_ID = .data$Sample_ID,
      inventory_source = "SDRF",
      field_raw = .data$field_raw,
      value = .data$value,
      file_class = vapply(
        .data$value,
        ae_core_classify_file_value,
        character(1)
      )
    ) %>%
    dplyr::distinct()
}


ae_core_biostudies_file_inventory <- function(
    accession,
    verbose = TRUE,
    biostudies_json = NULL) {

  answer <- biostudies_json

  if (is.null(answer)) {
    answer <- tryCatch(
      ae_core_biostudies_json(accession),
      error = function(e) e
    )
  }

  if (inherits(answer, "error")) {
    if (verbose) {
      message(
        "Could not retrieve BioStudies inventory for ",
        accession,
        ": ",
        conditionMessage(answer)
      )
    }

    return(
      tibble::tibble(
        Study_ID = character(),
        Sample_ID = character(),
        inventory_source = character(),
        field_raw = character(),
        value = character(),
        file_class = character()
      )
    )
  }

  strings <- unique(
    ae_core_clean_text(
      as.character(
        unlist(
          c(answer$study, answer$info),
          recursive = TRUE,
          use.names = FALSE
        )
      )
    )
  )

  strings <- strings[
    !is.na(strings) &
      (
        stringr::str_detect(
          strings,
          stringr::regex(
            "\\.(fastq|fq|sra|bam|cram|idat|csv|tsv|txt|xls|xlsx|bed|bedgraph|bigwig|bw|zip|tar|gz)([?#].*)?$",
            ignore_case = TRUE
          )
        ) |
          stringr::str_detect(
            strings,
            stringr::regex(
              "(^|[^A-Za-z0-9])(ERR|SRR|DRR)[0-9]+([^A-Za-z0-9]|$)",
              ignore_case = TRUE
            )
          )
      )
  ]

  tibble::tibble(
    Study_ID = accession,
    Sample_ID = NA_character_,
    inventory_source = "BioStudies API",
    field_raw = "BioStudies file or accession inventory",
    value = strings,
    file_class = vapply(
      strings,
      ae_core_classify_file_value,
      character(1)
    )
  ) %>%
    dplyr::distinct()
}


# -----------------------------------------------------------------------------
# 6. Sample-level screening
# -----------------------------------------------------------------------------

ae_core_extract_project_title <- function(
    idf_long = NULL,
    prescreen_summary = NULL) {

  if (
    !is.null(prescreen_summary) &&
      nrow(prescreen_summary) > 0L &&
      "title_or_first_value" %in%
        names(prescreen_summary)
  ) {
    title <- ae_core_first_nonempty(
      prescreen_summary$title_or_first_value
    )

    if (!is.na(title) && nzchar(title)) {
      return(title)
    }
  }

  if (
    is.null(idf_long) ||
      nrow(idf_long) == 0L
  ) {
    return(NA_character_)
  }

  field_normalized <- ae_core_normalize_field(
    idf_long$field_raw
  )

  preferred_fields <- c(
    "investigation_title",
    "study_title",
    "experiment_title",
    "project_title",
    "title"
  )

  for (preferred_field in preferred_fields) {
    values <- idf_long$value[
      field_normalized == preferred_field
    ]

    title <- ae_core_first_nonempty(
      values
    )

    if (!is.na(title) && nzchar(title)) {
      return(title)
    }
  }

  title_rows <- stringr::str_detect(
    field_normalized,
    "(^|_)(investigation|study|experiment|project)?_?title($|_)"
  )

  ae_core_first_nonempty(
    idf_long$value[title_rows]
  )
}


# Retrieve study titles for the project IDs in a data frame.
#
# Parameters
# ----------
# data:
#   A data.frame or tibble containing ArrayExpress/BioStudies accessions.
#
# id_col:
#   Character string naming the accession column, for example "Study_ID".
#
# search_result:
#   Optional object returned by biostudies_search_many(). Titles found in
#   search_result$unique_studies are used before making any API requests.
#
# query_api:
#   Query the BioStudies detailed-study API for accessions not found in
#   search_result.
#
# Returns
# -------
# A character vector aligned to the rows of data. The vector is named with the
# accession values.

# =============================================================================
# C. NMD-specific ArrayExpress screening
# =============================================================================
# =============================================================================
# ArrayExpress / BioStudies NMD perturbation screening extension
#
# PURPOSE
# -------
# Screen ArrayExpress/BioStudies studies for human cell-line RNA-seq experiments
# that perturb nonsense-mediated mRNA decay (NMD), then parse study- and
# sample-level experimental and technical metadata.
#
# This script intentionally REUSES the generic BioStudies/MAGE-TAB functions
# from:
#
#   arrayexpress_workflow.R
#
# The underlying download / IDF / SDRF functions are generic despite their
# historical ae_core_* names. The biological screening below is NMD-specific.
#
# Recommended source order
# ------------------------
# source("arrayexpress_workflow.R")
# source("arrayexpress_nmd_perturbation_screening.R")
#
# Main functions
# --------------
# ae_nmd_quick_prescreen()
# ae_screen_one_nmd_study()
# ae_screen_nmd_accessions()
# ae_write_nmd_screening_results()
#
# Example
# -------
# nmd_screen <- ae_screen_nmd_accessions(
#   accessions = nmd_search$accessions,
#   prescreen = TRUE,
#   prescreen_require_raw_evidence = TRUE,
#   prescreen_require_perturbation_evidence = TRUE,
#   prescreen_require_cell_model_evidence = FALSE,
#   prescreen_on_error = "continue"
# )
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Confirm that the generic ArrayExpress helpers have already been sourced
# -----------------------------------------------------------------------------

ae_nmd_required_base_functions <- c(
  "ae_core_biostudies_json",
  "ae_core_flatten_biostudies_json",
  "ae_core_extract_prescreen_title",
  "ae_core_dictionary_hit_text",
  "ae_core_has_dictionary_hit",
  "ae_core_collapse_unique",
  "ae_core_first_nonempty",
  "ae_core_clean_text",
  "ae_core_normalize_field",
  "ae_core_download_magetab",
  "ae_core_read_sdrf",
  "ae_core_read_idf",
  "ae_core_sdrf_file_inventory",
  "ae_core_biostudies_file_inventory",
  "ae_core_extract_project_title"
)

ae_nmd_missing_base_functions <- ae_nmd_required_base_functions[
  !vapply(
    ae_nmd_required_base_functions,
    exists,
    logical(1),
    mode = "function",
    inherits = TRUE
  )
]

if (length(ae_nmd_missing_base_functions) > 0L) {
  stop(
    paste0(
      "The generic ArrayExpress core helpers must be available before this NMD ",
      "extension.\n\nSource:\n",
      "  arrayexpress_workflow.R",
      "\n\nMissing functions:\n  ",
      paste(
        ae_nmd_missing_base_functions,
        collapse = "\n  "
      )
    )
  )
}


# -----------------------------------------------------------------------------
# 1. NMD dictionaries
# -----------------------------------------------------------------------------

AE_NMD_DIRECT_TERMS <- c(
  "Nonsense-mediated decay" =
    "nonsense[ -]+mediated[ -]+(?:mRNA[ -]+)?decay",
  "NMD" =
    "(^|[^A-Za-z0-9])NMD([^A-Za-z0-9]|$)"
)


AE_NMD_TARGET_TERMS <- c(
  "UPF1" =
    "(^|[^A-Za-z0-9])(UPF1|hUPF1|RENT1)([^A-Za-z0-9]|$)",
  "UPF2" =
    "(^|[^A-Za-z0-9])(UPF2|hUPF2|RENT2)([^A-Za-z0-9]|$)",
  "UPF3A" =
    "(^|[^A-Za-z0-9])UPF3A([^A-Za-z0-9]|$)",
  "UPF3B" =
    "(^|[^A-Za-z0-9])(UPF3B|UPF3X)([^A-Za-z0-9]|$)",
  "SMG1" =
    "(^|[^A-Za-z0-9])SMG[ -]?1([^A-Za-z0-9]|$)",
  "SMG5" =
    "(^|[^A-Za-z0-9])SMG[ -]?5([^A-Za-z0-9]|$)",
  "SMG6" =
    "(^|[^A-Za-z0-9])SMG[ -]?6([^A-Za-z0-9]|$)",
  "SMG7" =
    "(^|[^A-Za-z0-9])SMG[ -]?7([^A-Za-z0-9]|$)",
  "SMG8" =
    "(^|[^A-Za-z0-9])SMG[ -]?8([^A-Za-z0-9]|$)",
  "SMG9" =
    "(^|[^A-Za-z0-9])SMG[ -]?9([^A-Za-z0-9]|$)",
  "DHX34" =
    "(^|[^A-Za-z0-9])DHX34([^A-Za-z0-9]|$)",
  "NBAS" =
    "(^|[^A-Za-z0-9])NBAS([^A-Za-z0-9]|$)"
)


AE_NMD_GENETIC_PERTURBATION_TERMS <- c(
  "siRNA" =
    "(^|[^A-Za-z0-9])siRNA([^A-Za-z0-9]|$)|small[ -]+interfering[ -]+RNA",
  "shRNA" =
    "(^|[^A-Za-z0-9])shRNA([^A-Za-z0-9]|$)|short[ -]+hairpin[ -]+RNA",
  "RNAi" =
    "(^|[^A-Za-z0-9])RNAi([^A-Za-z0-9]|$)",
  "Knockdown" =
    "knock[ -]?down|deplet(?:e|ed|ion)|silenc(?:e|ed|ing)",
  "Knockout" =
    "knock[ -]?out|(^|[^A-Za-z0-9])KO([^A-Za-z0-9]|$)",
  "CRISPR" =
    "CRISPR|Cas9",
  "Degron" =
    "degron|auxin[ -]+inducible|(^|[^A-Za-z0-9])AID(?:[0-9]+)?([^A-Za-z0-9]|$)",
  "Overexpression" =
    "over[ -]?express(?:ion|ed|ing)?|(^|[^A-Za-z0-9])OE([^A-Za-z0-9]|$)",
  "Inducible expression" =
    "inducible|doxycycline|(^|[^A-Za-z0-9])dox([^A-Za-z0-9]|$)"
)


AE_NMD_DIRECT_INHIBITOR_TERMS <- c(
  "SMG1 inhibitor" =
    "SMG[ -]?1[ -]+inhibitor|SMG1i",
  "SMG1 inhibitor 11j" =
    "(^|[^A-Za-z0-9])11j([^A-Za-z0-9]|$)",
  "NMDI-1" =
    "(^|[^A-Za-z0-9])NMDI[ -]?1([^A-Za-z0-9]|$)",
  "VG1" =
    "(^|[^A-Za-z0-9])VG1([^A-Za-z0-9]|$)"
)


AE_NMD_TRANSLATION_INHIBITOR_TERMS <- c(
  "Emetine" =
    "(^|[^A-Za-z0-9])emetine([^A-Za-z0-9]|$)",
  "Cycloheximide" =
    "cycloheximide|(^|[^A-Za-z0-9])CHX([^A-Za-z0-9]|$)",
  "Puromycin" =
    "(^|[^A-Za-z0-9])puromycin([^A-Za-z0-9]|$)",
  "Caffeine" =
    "(^|[^A-Za-z0-9])caffeine([^A-Za-z0-9]|$)"
)


AE_NMD_STRESS_TERMS <- c(
  "Thapsigargin" =
    "(^|[^A-Za-z0-9])thapsigargin([^A-Za-z0-9]|$)",
  "Tunicamycin" =
    "(^|[^A-Za-z0-9])tunicamycin([^A-Za-z0-9]|$)",
  "Arsenite" =
    "sodium[ -]+arsenite|(^|[^A-Za-z0-9])arsenite([^A-Za-z0-9]|$)",
  "ER stress" =
    "ER[ -]+stress|endoplasmic[ -]+reticulum[ -]+stress",
  "Integrated stress response" =
    "integrated[ -]+stress[ -]+response|(^|[^A-Za-z0-9])ISR([^A-Za-z0-9]|$)"
)


AE_NMD_CONTROL_TERMS <- c(
  "Untreated control" =
    "untreated|no[ -]+treatment|mock[ -]+treated",
  "Control" =
    "(^|[^A-Za-z0-9])controls?([^A-Za-z0-9]|$)",
  "Negative control" =
    "negative[ -]+control",
  "Scramble control" =
    "scrambl(?:e|ed)|non[ -]?targeting|nontargeting",
  "siRNA control" =
    "si(?:RNA)?[ -]?(?:control|ctrl)|siCTRL|siNC",
  "shRNA control" =
    "sh(?:RNA)?[ -]?(?:control|ctrl)",
  "DMSO vehicle" =
    "(^|[^A-Za-z0-9])DMSO([^A-Za-z0-9]|$)|vehicle",
  "Wild type" =
    "wild[ -]?type|(^|[^A-Za-z0-9])WT([^A-Za-z0-9]|$)",
  "Parental" =
    "(^|[^A-Za-z0-9])parental([^A-Za-z0-9]|$)",
  "Empty vector" =
    "empty[ -]+vector|vector[ -]+control",
  "Non-induced" =
    "non[ -]?induced|uninduced|without[ -]+dox(?:ycycline)?"
)


AE_NMD_HUMAN_TERMS <- c(
  "Homo sapiens" =
    "homo[ _-]*sapiens",
  "Human" =
    "(^|[^A-Za-z])human([^A-Za-z]|$)"
)


AE_NMD_NONHUMAN_TERMS <- c(
  "Mus musculus" =
    "mus[ _-]*musculus",
  "Mouse or murine" =
    "(^|[^A-Za-z])(mouse|mice|murine)([^A-Za-z]|$)",
  "Rattus norvegicus" =
    "rattus[ _-]*norvegicus",
  "Rat" =
    "(^|[^A-Za-z])rats?([^A-Za-z]|$)"
)


AE_NMD_CELL_MODEL_TERMS <- c(
  "Cell line" =
    "cell[ -]+line|cellline",
  "Cultured cells" =
    "cultured[ -]+cells?",
  "HEK293" =
    "HEK[ -]?293(?:T|TO)?|293T|293TO|Flp[ -]?In[ -]?T[ -]?REx[ -]?293",
  "HeLa" =
    "(^|[^A-Za-z0-9])HeLa([^A-Za-z0-9]|$)",
  "K562" =
    "(^|[^A-Za-z0-9])K562([^A-Za-z0-9]|$)",
  "Huh7" =
    "(^|[^A-Za-z0-9])Huh[ -]?7([^A-Za-z0-9]|$)",
  "HT1080" =
    "(^|[^A-Za-z0-9])HT[ -]?1080([^A-Za-z0-9]|$)",
  "A498" =
    "(^|[^A-Za-z0-9])A[ -]?498([^A-Za-z0-9]|$)",
  "A704" =
    "(^|[^A-Za-z0-9])A[ -]?704([^A-Za-z0-9]|$)",
  "HK1" =
    "(^|[^A-Za-z0-9])HK[ -]?1([^A-Za-z0-9]|$)",
  "SUNE1" =
    "(^|[^A-Za-z0-9])SUNE[ -]?1([^A-Za-z0-9]|$)"
)


AE_NMD_RNASEQ_TERMS <- c(
  "RNA-seq" =
    "RNA[ -]?seq|RNASeq",
  "RNA sequencing" =
    "RNA[ -]+sequencing",
  "Transcriptome sequencing" =
    "transcriptome[ -]+sequenc",
  "Expression profiling by high-throughput sequencing" =
    "expression[ -]+profiling[ -]+by[ -]+high[ -]+throughput[ -]+sequenc",
  "mRNA sequencing" =
    "mRNA[ -]+sequenc",
  "Total RNA sequencing" =
    "total[ -]+RNA[ -]+sequenc"
)


AE_NMD_RAW_SEQUENCE_TERMS <- c(
  "FASTQ" =
    "\\.(fastq|fq)(\\.gz)?([?#].*)?$",
  "SRA file" =
    "\\.sra([?#].*)?$",
  "ENA/SRA run accession" =
    "(^|[^A-Za-z0-9])(ERR|SRR|DRR)[0-9]+([^A-Za-z0-9]|$)",
  "SRA experiment accession" =
    "(^|[^A-Za-z0-9])(ERX|SRX|DRX)[0-9]+([^A-Za-z0-9]|$)",
  "BAM" =
    "\\.bam([?#].*)?$",
  "CRAM" =
    "\\.cram([?#].*)?$"
)


AE_NMD_GENERIC_RAW_TERMS <- c(
  "Raw data wording" =
    "(^|[^A-Za-z])raw[ -]+(data|file|files|read|reads)([^A-Za-z]|$)",
  "ENA" =
    "(^|[^A-Za-z0-9])ENA([^A-Za-z0-9]|$)|European[ -]+Nucleotide[ -]+Archive",
  "SRA" =
    "(^|[^A-Za-z0-9])SRA([^A-Za-z0-9]|$)|Sequence[ -]+Read[ -]+Archive",
  "BioProject accession" =
    "(^|[^A-Za-z0-9])PRJ(E|N|D)[A-Z]?[0-9]+([^A-Za-z0-9]|$)",
  "Sequencing project accession" =
    "(^|[^A-Za-z0-9])(ERP|SRP|DRP)[0-9]+([^A-Za-z0-9]|$)"
)


# -----------------------------------------------------------------------------
# 2. Classification helpers
# -----------------------------------------------------------------------------

ae_nmd_hits <- function(text, dictionary) {
  ae_core_dictionary_hit_text(
    text,
    dictionary
  )
}


ae_nmd_has <- function(text, dictionary) {
  ae_core_has_dictionary_hit(
    text,
    dictionary
  )
}


ae_nmd_classify_perturbation <- function(text) {

  text <- ae_core_clean_text(text)

  if (is.na(text) || !nzchar(text)) {
    return("Unclear")
  }

  if (ae_nmd_has(text, AE_NMD_CONTROL_TERMS)) {
    return("Control")
  }

  if (
    grepl(
      "degron|auxin[ -]+inducible|(^|[^A-Za-z0-9])AID(?:[0-9]+)?([^A-Za-z0-9]|$)",
      text,
      ignore.case = TRUE,
      perl = TRUE
    )
  ) {
    return("Inducible degradation")
  }

  if (
    grepl(
      "knock[ -]?out|CRISPR|Cas9|(^|[^A-Za-z0-9])KO([^A-Za-z0-9]|$)",
      text,
      ignore.case = TRUE,
      perl = TRUE
    )
  ) {
    return("Knockout / CRISPR")
  }

  if (
    grepl(
      "(^|[^A-Za-z0-9])siRNA([^A-Za-z0-9]|$)|small[ -]+interfering[ -]+RNA",
      text,
      ignore.case = TRUE,
      perl = TRUE
    )
  ) {
    return("siRNA knockdown")
  }

  if (
    grepl(
      "(^|[^A-Za-z0-9])shRNA([^A-Za-z0-9]|$)|short[ -]+hairpin[ -]+RNA",
      text,
      ignore.case = TRUE,
      perl = TRUE
    )
  ) {
    return("shRNA knockdown")
  }

  if (
    grepl(
      "knock[ -]?down|deplet(?:e|ed|ion)|silenc(?:e|ed|ing)|(^|[^A-Za-z0-9])RNAi([^A-Za-z0-9]|$)",
      text,
      ignore.case = TRUE,
      perl = TRUE
    )
  ) {
    return("Other knockdown / depletion")
  }

  if (ae_nmd_has(text, AE_NMD_DIRECT_INHIBITOR_TERMS)) {
    return("Direct NMD / SMG1 inhibitor")
  }

  if (ae_nmd_has(text, AE_NMD_TRANSLATION_INHIBITOR_TERMS)) {
    return("Translation inhibitor")
  }

  if (ae_nmd_has(text, AE_NMD_STRESS_TERMS)) {
    return("Stress-mediated NMD modulation")
  }

  if (
    grepl(
      "over[ -]?express(?:ion|ed|ing)?|(^|[^A-Za-z0-9])OE([^A-Za-z0-9]|$)",
      text,
      ignore.case = TRUE,
      perl = TRUE
    )
  ) {
    return("NMD-factor overexpression")
  }

  if (
    ae_nmd_has(text, AE_NMD_TARGET_TERMS) ||
      ae_nmd_has(text, AE_NMD_DIRECT_TERMS)
  ) {
    return("NMD perturbation; method unclear")
  }

  "Other / unclear"
}


ae_nmd_expected_effect <- function(
    perturbation_type,
    text = NA_character_) {

  if (is.na(perturbation_type)) {
    return("Unclear")
  }

  if (perturbation_type == "Control") {
    return("Control / baseline NMD")
  }

  if (
    perturbation_type %in% c(
      "Inducible degradation",
      "Knockout / CRISPR",
      "siRNA knockdown",
      "shRNA knockdown",
      "Other knockdown / depletion",
      "Direct NMD / SMG1 inhibitor"
    )
  ) {
    return("Likely NMD inhibition")
  }

  if (
    perturbation_type ==
      "Translation inhibitor"
  ) {
    return("Indirect NMD inhibition via translation inhibition")
  }

  if (
    perturbation_type ==
      "Stress-mediated NMD modulation"
  ) {
    return("Potential stress-mediated NMD inhibition")
  }

  if (
    perturbation_type ==
      "NMD-factor overexpression"
  ) {
    return("NMD-factor overexpression / altered NMD activity")
  }

  "NMD effect unclear"
}


ae_nmd_control_type <- function(text) {

  hits <- ae_nmd_hits(
    text,
    AE_NMD_CONTROL_TERMS
  )

  if (is.na(hits)) {
    return(NA_character_)
  }

  hits
}


# -----------------------------------------------------------------------------
# 3. Fast BioStudies pre-screen
# -----------------------------------------------------------------------------

ae_nmd_quick_prescreen <- function(
    accession,
    biostudies_json = NULL,
    require_raw_evidence = TRUE,
    require_perturbation_evidence = TRUE,
    require_human_evidence = FALSE,
    exclude_clear_nonhuman = TRUE,
    require_cell_model_evidence = FALSE,
    include_stress_perturbations = FALSE,
    verbose = TRUE) {

  accession <- toupper(
    trimws(accession)
  )

  if (is.null(biostudies_json)) {
    if (verbose) {
      message(
        "Quick NMD BioStudies pre-screen: ",
        accession
      )
    }

    biostudies_json <-
      ae_core_biostudies_json(
        accession
      )
  }

  flat_json <-
    ae_core_flatten_biostudies_json(
      biostudies_json =
        biostudies_json,
      accession = accession
    )

  # Avoid using publication/citation titles as biological evidence.
  screening_rows <-
    !stringr::str_detect(
      flat_json$field,
      stringr::regex(
        "reference|publication|pubmed|citation|bibliograph",
        ignore_case = TRUE
      )
    )

  screening_text <-
    ae_core_collapse_unique(
      paste(
        flat_json$field[screening_rows],
        flat_json$value[screening_rows],
        sep = ": "
      ),
      max_values = 2500L
    )

  all_text <-
    ae_core_collapse_unique(
      paste(
        flat_json$field,
        flat_json$value,
        sep = ": "
      ),
      max_values = 3500L
    )

  direct_nmd_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_DIRECT_TERMS
    )

  target_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_TARGET_TERMS
    )

  genetic_perturbation_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_GENETIC_PERTURBATION_TERMS
    )

  direct_inhibitor_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_DIRECT_INHIBITOR_TERMS
    )

  translation_inhibitor_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_TRANSLATION_INHIBITOR_TERMS
    )

  stress_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_STRESS_TERMS
    )

  cell_model_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_CELL_MODEL_TERMS
    )

  human_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_HUMAN_TERMS
    )

  nonhuman_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_NONHUMAN_TERMS
    )

  rnaseq_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_RNASEQ_TERMS
    )

  raw_sequence_terms <-
    ae_nmd_hits(
      all_text,
      AE_NMD_RAW_SEQUENCE_TERMS
    )

  generic_raw_terms <-
    ae_nmd_hits(
      all_text,
      AE_NMD_GENERIC_RAW_TERMS
    )

  nmd_flag <-
    !is.na(direct_nmd_terms) ||
    !is.na(target_terms) ||
    !is.na(direct_inhibitor_terms)

  perturbation_flag <-
    !is.na(genetic_perturbation_terms) ||
    !is.na(direct_inhibitor_terms) ||
    !is.na(translation_inhibitor_terms) ||
    (
      isTRUE(include_stress_perturbations) &&
        !is.na(stress_terms)
    )

  human_flag <- !is.na(human_terms)
  nonhuman_flag <- !is.na(nonhuman_terms)
  cell_model_flag <- !is.na(cell_model_terms)
  rnaseq_flag <- !is.na(rnaseq_terms)

  raw_sequence_flag <-
    !is.na(raw_sequence_terms) ||
    !is.na(generic_raw_terms)

  human_ok <- TRUE

  if (isTRUE(require_human_evidence)) {
    human_ok <- human_flag
  } else if (
    isTRUE(exclude_clear_nonhuman) &&
      nonhuman_flag &&
      !human_flag
  ) {
    human_ok <- FALSE
  }

  cell_model_ok <-
    !isTRUE(require_cell_model_evidence) ||
    cell_model_flag

  perturbation_ok <-
    !isTRUE(require_perturbation_evidence) ||
    perturbation_flag

  raw_ok <-
    !isTRUE(require_raw_evidence) ||
    raw_sequence_flag

  continue_full_parse <-
    nmd_flag &&
    perturbation_ok &&
    human_ok &&
    cell_model_ok &&
    rnaseq_flag &&
    raw_ok

  prescreen_status <- dplyr::case_when(
    !nmd_flag ~
      "Skip: NMD pathway or NMD factor not identified",

    !perturbation_ok ~
      "Skip: NMD perturbation evidence not identified",

    !human_ok ~
      "Skip: clearly non-human or human evidence required",

    !cell_model_ok ~
      "Skip: cell-line/model evidence not identified",

    !rnaseq_flag ~
      "Skip: RNA-seq not identified",

    !raw_ok ~
      "Skip: raw RNA-seq evidence not identified",

    TRUE ~
      "Continue: NMD perturbation RNA-seq candidate"
  )

  summary <- tibble::tibble(
    Study_ID = accession,
    project_title =
      ae_core_extract_prescreen_title(
        flat_json
      ),
    prescreen_status =
      prescreen_status,
    continue_full_parse =
      continue_full_parse,

    nmd_flag = nmd_flag,
    perturbation_flag =
      perturbation_flag,
    human_flag = human_flag,
    nonhuman_flag =
      nonhuman_flag,
    cell_model_flag =
      cell_model_flag,
    rnaseq_flag =
      rnaseq_flag,
    raw_sequence_flag =
      raw_sequence_flag,

    direct_nmd_terms =
      direct_nmd_terms,
    nmd_target_terms =
      target_terms,
    genetic_perturbation_terms =
      genetic_perturbation_terms,
    direct_inhibitor_terms =
      direct_inhibitor_terms,
    translation_inhibitor_terms =
      translation_inhibitor_terms,
    stress_terms =
      stress_terms,
    cell_model_terms =
      cell_model_terms,
    human_terms =
      human_terms,
    nonhuman_terms =
      nonhuman_terms,
    rnaseq_terms =
      rnaseq_terms,
    raw_sequence_terms =
      raw_sequence_terms,
    generic_raw_terms =
      generic_raw_terms
  )

  evidence <- flat_json %>%
    dplyr::mutate(
      row_text = paste(
        .data$field,
        .data$value,
        sep = ": "
      ),
      evidence_class =
        dplyr::case_when(
          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_DIRECT_TERMS
          ) ~ "NMD terminology",

          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_TARGET_TERMS
          ) ~ "NMD factor",

          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_GENETIC_PERTURBATION_TERMS
          ) ~ "Genetic perturbation",

          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_DIRECT_INHIBITOR_TERMS
          ) ~ "Direct NMD inhibitor",

          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_TRANSLATION_INHIBITOR_TERMS
          ) ~ "Translation inhibitor",

          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_STRESS_TERMS
          ) ~ "Stress perturbation",

          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_CELL_MODEL_TERMS
          ) ~ "Cell model",

          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_RNASEQ_TERMS
          ) ~ "RNA-seq",

          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_RAW_SEQUENCE_TERMS
          ) ~ "Raw sequencing evidence",

          TRUE ~ NA_character_
        )
    ) %>%
    dplyr::filter(
      !is.na(.data$evidence_class)
    ) %>%
    dplyr::select(
      .data$Study_ID,
      .data$json_source,
      .data$evidence_class,
      .data$field,
      .data$value
    ) %>%
    dplyr::distinct()

  list(
    summary = summary,
    evidence = evidence,
    json = biostudies_json
  )
}


# -----------------------------------------------------------------------------
# 4. Sample-level NMD metadata parser
# -----------------------------------------------------------------------------

ae_nmd_collapse_where <- function(
    values,
    keep,
    max_values = 30L) {

  ae_core_collapse_unique(
    values[keep],
    max_values = max_values
  )
}


ae_nmd_make_sample_summary <- function(
    accession,
    sdrf_long,
    idf_long,
    prescreen_summary = NULL) {

  if (nrow(sdrf_long) == 0L) {
    return(tibble::tibble())
  }

  study_context <-
    ae_core_collapse_unique(
      paste(
        idf_long$field_raw,
        idf_long$value,
        sep = ": "
      ),
      max_values = 300L
    )

  study_nmd_terms <-
    ae_nmd_hits(
      study_context,
      AE_NMD_DIRECT_TERMS
    )

  study_target_terms <-
    ae_nmd_hits(
      study_context,
      AE_NMD_TARGET_TERMS
    )

  study_rnaseq_terms <-
    ae_nmd_hits(
      study_context,
      AE_NMD_RNASEQ_TERMS
    )

  tagged <- sdrf_long %>%
    dplyr::mutate(
      field =
        ae_core_normalize_field(
          .data$field_raw
        ),

      is_organism_field =
        stringr::str_detect(
          .data$field,
          "organism|species"
        ),

      is_source_field =
        stringr::str_detect(
          .data$field,
          "source_name|cell_line|cell_type|model|biosource|material|organism_part|tissue"
        ),

      is_characteristic_field =
        stringr::str_detect(
          .data$field,
          "characteristic|factor_value|genotype|phenotype|disease|cell"
        ),

      is_perturbation_field =
        stringr::str_detect(
          .data$field,
          "factor_value|treatment|compound|drug|perturb|transfect|transduc|RNAi|siRNA|shRNA|CRISPR|knock|deplet|genotype|construct|vector|induc|dox|auxin|protocol"
        ),

      is_control_field =
        stringr::str_detect(
          .data$field,
          "factor_value|treatment|control|group|condition|genotype|construct|vector"
        ),

      is_assay_field =
        stringr::str_detect(
          .data$field,
          "assay|technology|sequenc|library|rna"
        ),

      is_library_strategy_field =
        stringr::str_detect(
          .data$field,
          "library_strategy|sequencing_strategy"
        ),

      is_library_source_field =
        stringr::str_detect(
          .data$field,
          "library_source"
        ),

      is_library_selection_field =
        stringr::str_detect(
          .data$field,
          "library_selection"
        ),

      is_library_layout_field =
        stringr::str_detect(
          .data$field,
          "library_layout|paired|single_end|singleend"
        ),

      is_library_prep_field =
        stringr::str_detect(
          .data$field,
          "library_prep|library_preparation|library_construction|protocol_ref"
        ),

      is_instrument_field =
        stringr::str_detect(
          .data$field,
          "instrument|sequencer|sequencing_platform|platform"
        ),

      is_stranded_field =
        stringr::str_detect(
          .data$field,
          "stranded|strand_specific|strand_specificity"
        ),

      is_read_length_field =
        stringr::str_detect(
          .data$field,
          "read_length|readlength|cycle"
        ),

      is_file_field =
        stringr::str_detect(
          .data$field,
          "file|uri|url|fastq|fq|bam|cram|ena|sra|run_accession|raw"
        )
    )

  sample_summary <- tagged %>%
    dplyr::group_by(
      .data$Study_ID,
      .data$Sample_ID
    ) %>%
    dplyr::summarise(
      all_sample_metadata =
        ae_core_collapse_unique(
          paste(
            .data$field_raw,
            .data$value,
            sep = ": "
          ),
          max_values = 250L
        ),

      organism_values =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_organism_field
        ),

      source_model_values =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_source_field
        ),

      characteristic_values =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_characteristic_field,
          max_values = 60L
        ),

      perturbation_values =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_perturbation_field,
          max_values = 80L
        ),

      control_condition_values =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_control_field,
          max_values = 60L
        ),

      assay_values =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_assay_field
        ),

      library_strategy =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_library_strategy_field
        ),

      library_source =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_library_source_field
        ),

      library_selection =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_library_selection_field
        ),

      library_layout =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_library_layout_field
        ),

      library_preparation =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_library_prep_field,
          max_values = 40L
        ),

      sequencing_instrument =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_instrument_field
        ),

      strandedness =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_stranded_field
        ),

      read_length =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_read_length_field
        ),

      file_values =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_file_field,
          max_values = 80L
        ),

      .groups = "drop"
    ) %>%
    dplyr::mutate(
      study_context =
        study_context,

      sample_search_text =
        paste(
          dplyr::coalesce(
            .data$all_sample_metadata,
            ""
          ),
          dplyr::coalesce(
            .data$source_model_values,
            ""
          ),
          dplyr::coalesce(
            .data$characteristic_values,
            ""
          ),
          dplyr::coalesce(
            .data$perturbation_values,
            ""
          )
        ),

      perturbation_search_text =
        paste(
          dplyr::coalesce(
            .data$perturbation_values,
            ""
          ),
          dplyr::coalesce(
            .data$control_condition_values,
            ""
          ),
          dplyr::coalesce(
            .data$characteristic_values,
            ""
          )
        ),

      assay_search_text =
        paste(
          dplyr::coalesce(
            .data$assay_values,
            ""
          ),
          dplyr::coalesce(
            .data$library_strategy,
            ""
          ),
          dplyr::coalesce(
            .data$file_values,
            ""
          ),
          dplyr::coalesce(
            .data$study_context,
            ""
          )
        ),

      human_terms =
        vapply(
          .data$sample_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_HUMAN_TERMS
        ),

      nonhuman_terms =
        vapply(
          .data$sample_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_NONHUMAN_TERMS
        ),

      cell_model_terms =
        vapply(
          .data$sample_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_CELL_MODEL_TERMS
        ),

      nmd_direct_terms =
        vapply(
          .data$sample_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_DIRECT_TERMS
        ),

      nmd_target_terms =
        vapply(
          .data$sample_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_TARGET_TERMS
        ),

      genetic_perturbation_terms =
        vapply(
          .data$perturbation_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_GENETIC_PERTURBATION_TERMS
        ),

      direct_inhibitor_terms =
        vapply(
          .data$perturbation_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_DIRECT_INHIBITOR_TERMS
        ),

      translation_inhibitor_terms =
        vapply(
          .data$perturbation_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_TRANSLATION_INHIBITOR_TERMS
        ),

      stress_terms =
        vapply(
          .data$perturbation_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_STRESS_TERMS
        ),

      control_terms =
        vapply(
          .data$perturbation_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_CONTROL_TERMS
        ),

      rnaseq_sample_terms =
        vapply(
          .data$assay_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_RNASEQ_TERMS
        ),

      raw_sequence_terms =
        vapply(
          dplyr::coalesce(
            .data$file_values,
            ""
          ),
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_RAW_SEQUENCE_TERMS
        ),

      human_flag =
        !is.na(.data$human_terms),

      nonhuman_flag =
        !is.na(.data$nonhuman_terms),

      cell_model_flag =
        !is.na(.data$cell_model_terms) |
        (
          !is.na(.data$source_model_values) &
            nzchar(.data$source_model_values)
        ),

      rnaseq_flag =
        !is.na(.data$rnaseq_sample_terms) |
        !is.na(study_rnaseq_terms),

      raw_sequence_flag =
        !is.na(.data$raw_sequence_terms),

      perturbation_type =
        vapply(
          .data$perturbation_search_text,
          ae_nmd_classify_perturbation,
          character(1)
        ),

      control_type =
        vapply(
          .data$perturbation_search_text,
          ae_nmd_control_type,
          character(1)
        ),

      control_flag =
        .data$perturbation_type ==
          "Control" |
        !is.na(.data$control_terms),

      direct_nmd_evidence_flag =
        !is.na(.data$nmd_direct_terms) |
        !is.na(.data$nmd_target_terms) |
        !is.na(.data$direct_inhibitor_terms),

      nmd_perturbation_flag =
        !.data$control_flag &
        (
          !is.na(.data$nmd_target_terms) |
          !is.na(.data$direct_inhibitor_terms) |
          (
            !is.na(.data$nmd_direct_terms) &
              .data$perturbation_type !=
                "Other / unclear"
          )
        ),

      expected_nmd_effect =
        mapply(
          ae_nmd_expected_effect,
          perturbation_type =
            .data$perturbation_type,
          text =
            .data$perturbation_search_text,
          USE.NAMES = FALSE
        )
    )

  # Controls in an NMD study are useful even though they do not contain the
  # NMD-factor name themselves.
  study_has_perturbation <- any(
    sample_summary$nmd_perturbation_flag,
    na.rm = TRUE
  )

  study_human_from_prescreen <- FALSE

  if (
    !is.null(prescreen_summary) &&
      nrow(prescreen_summary) > 0L &&
      "human_flag" %in%
        names(prescreen_summary)
  ) {
    study_human_from_prescreen <-
      isTRUE(
        prescreen_summary$human_flag[1]
      )
  }

  sample_summary <- sample_summary %>%
    dplyr::mutate(
      human_flag =
        .data$human_flag |
        study_human_from_prescreen,

      study_has_nmd_perturbation =
        study_has_perturbation,

      nmd_control_sample_flag =
        .data$control_flag &
        study_has_perturbation,

      eligible_nmd_sample_flag =
        .data$human_flag &
        .data$cell_model_flag &
        .data$rnaseq_flag &
        (
          .data$nmd_perturbation_flag |
          .data$nmd_control_sample_flag
        ),

      sample_role =
        dplyr::case_when(
          .data$nmd_perturbation_flag ~
            "NMD perturbation",

          .data$nmd_control_sample_flag ~
            "Control",

          TRUE ~
            "Other / unclear"
        ),

      nmd_target_or_study_target =
        dplyr::coalesce(
          .data$nmd_target_terms,
          study_target_terms
        )
    ) %>%
    dplyr::select(
      .data$Study_ID,
      .data$Sample_ID,
      .data$sample_role,

      .data$nmd_target_or_study_target,
      .data$nmd_target_terms,
      .data$perturbation_type,
      .data$expected_nmd_effect,
      .data$control_type,

      .data$perturbation_values,
      .data$control_condition_values,

      .data$organism_values,
      .data$source_model_values,
      .data$characteristic_values,

      .data$human_flag,
      .data$nonhuman_flag,
      .data$cell_model_flag,

      .data$nmd_perturbation_flag,
      .data$nmd_control_sample_flag,
      .data$eligible_nmd_sample_flag,

      .data$rnaseq_flag,
      .data$raw_sequence_flag,

      .data$assay_values,
      .data$library_strategy,
      .data$library_source,
      .data$library_selection,
      .data$library_layout,
      .data$library_preparation,
      .data$sequencing_instrument,
      .data$strandedness,
      .data$read_length,
      .data$file_values,

      .data$nmd_direct_terms,
      .data$genetic_perturbation_terms,
      .data$direct_inhibitor_terms,
      .data$translation_inhibitor_terms,
      .data$stress_terms,
      .data$control_terms
    )

  sample_summary
}


# -----------------------------------------------------------------------------
# 5. NMD evidence table
# -----------------------------------------------------------------------------

ae_nmd_make_evidence <- function(
    idf_long,
    sdrf_long,
    file_inventory) {

  tag_rows <- function(
      table,
      source,
      dictionary,
      evidence_class,
      sample_column = NULL) {

    if (nrow(table) == 0L) {
      return(tibble::tibble())
    }

    row_text <- paste(
      table$field_raw,
      table$value,
      sep = ": "
    )

    matched <- vapply(
      row_text,
      ae_nmd_hits,
      character(1),
      dictionary = dictionary
    )

    keep <- !is.na(matched)

    if (!any(keep)) {
      return(tibble::tibble())
    }

    sample_id <- rep(
      NA_character_,
      sum(keep)
    )

    if (
      !is.null(sample_column) &&
        sample_column %in% names(table)
    ) {
      sample_id <-
        as.character(
          table[[sample_column]][keep]
        )
    }

    tibble::tibble(
      Study_ID =
        table$Study_ID[keep],
      Sample_ID =
        sample_id,
      source = source,
      evidence_class =
        evidence_class,
      field =
        table$field_raw[keep],
      value =
        table$value[keep],
      matched_terms =
        matched[keep]
    )
  }

  evidence <- dplyr::bind_rows(
    tag_rows(
      idf_long,
      "IDF",
      AE_NMD_DIRECT_TERMS,
      "NMD terminology"
    ),
    tag_rows(
      sdrf_long,
      "SDRF",
      AE_NMD_DIRECT_TERMS,
      "NMD terminology",
      "Sample_ID"
    ),
    tag_rows(
      idf_long,
      "IDF",
      AE_NMD_TARGET_TERMS,
      "NMD factor"
    ),
    tag_rows(
      sdrf_long,
      "SDRF",
      AE_NMD_TARGET_TERMS,
      "NMD factor",
      "Sample_ID"
    ),
    tag_rows(
      idf_long,
      "IDF",
      AE_NMD_GENETIC_PERTURBATION_TERMS,
      "Genetic perturbation"
    ),
    tag_rows(
      sdrf_long,
      "SDRF",
      AE_NMD_GENETIC_PERTURBATION_TERMS,
      "Genetic perturbation",
      "Sample_ID"
    ),
    tag_rows(
      idf_long,
      "IDF",
      AE_NMD_DIRECT_INHIBITOR_TERMS,
      "Direct NMD inhibitor"
    ),
    tag_rows(
      sdrf_long,
      "SDRF",
      AE_NMD_DIRECT_INHIBITOR_TERMS,
      "Direct NMD inhibitor",
      "Sample_ID"
    ),
    tag_rows(
      sdrf_long,
      "SDRF",
      AE_NMD_CONTROL_TERMS,
      "Control condition",
      "Sample_ID"
    ),
    tag_rows(
      sdrf_long,
      "SDRF",
      AE_NMD_CELL_MODEL_TERMS,
      "Cell model",
      "Sample_ID"
    ),
    tag_rows(
      idf_long,
      "IDF",
      AE_NMD_RNASEQ_TERMS,
      "RNA-seq"
    ),
    tag_rows(
      sdrf_long,
      "SDRF",
      AE_NMD_RNASEQ_TERMS,
      "RNA-seq",
      "Sample_ID"
    )
  ) %>%
    dplyr::distinct()

  if (
    !is.null(file_inventory) &&
      nrow(file_inventory) > 0L
  ) {
    file_columns <- intersect(
      c(
        "Study_ID",
        "Sample_ID",
        "file_value",
        "file_name",
        "file_url",
        "file_class",
        "source"
      ),
      names(file_inventory)
    )

    file_evidence <-
      file_inventory %>%
      dplyr::select(
        dplyr::all_of(
          file_columns
        )
      )

    if (nrow(file_evidence) > 0L) {
      file_evidence$evidence_class <-
        "File evidence"

      evidence <- dplyr::bind_rows(
        evidence,
        file_evidence
      )
    }
  }

  evidence
}


# -----------------------------------------------------------------------------
# 6. Study-level NMD summary
# -----------------------------------------------------------------------------

ae_nmd_make_study_summary <- function(
    accession,
    project_title,
    sample_summary,
    idf_long,
    file_inventory,
    download_method,
    prescreen_summary = NULL) {

  study_text <-
    ae_core_collapse_unique(
      paste(
        idf_long$field_raw,
        idf_long$value,
        sep = ": "
      ),
      max_values = 500L
    )

  n_samples <-
    dplyr::n_distinct(
      sample_summary$Sample_ID
    )

  n_nmd_perturbation_samples <-
    sum(
      sample_summary$nmd_perturbation_flag,
      na.rm = TRUE
    )

  n_control_samples <-
    sum(
      sample_summary$nmd_control_sample_flag,
      na.rm = TRUE
    )

  n_eligible_samples <-
    sum(
      sample_summary$eligible_nmd_sample_flag,
      na.rm = TRUE
    )

  n_raw_sample_links <-
    sum(
      sample_summary$raw_sequence_flag,
      na.rm = TRUE
    )

  file_text <- NA_character_

  if (
    !is.null(file_inventory) &&
      nrow(file_inventory) > 0L
  ) {
    file_text <-
      ae_core_collapse_unique(
        unlist(
          file_inventory,
          recursive = TRUE,
          use.names = FALSE
        ),
        max_values = 1000L
      )
  }

  study_raw_flag <-
    ae_nmd_has(
      file_text,
      AE_NMD_RAW_SEQUENCE_TERMS
    ) ||
    n_raw_sample_links > 0L

  study_targets <-
    ae_core_collapse_unique(
      c(
        sample_summary$nmd_target_or_study_target,
        ae_nmd_hits(
          study_text,
          AE_NMD_TARGET_TERMS
        )
      )
    )

  perturbation_methods <-
    ae_core_collapse_unique(
      sample_summary$perturbation_type[
        sample_summary$nmd_perturbation_flag
      ]
    )

  expected_nmd_effects <-
    ae_core_collapse_unique(
      sample_summary$expected_nmd_effect[
        sample_summary$nmd_perturbation_flag
      ]
    )

  cell_models <-
    ae_core_collapse_unique(
      sample_summary$source_model_values[
        sample_summary$cell_model_flag
      ],
      max_values = 40L
    )

  control_types <-
    ae_core_collapse_unique(
      sample_summary$control_type[
        sample_summary$nmd_control_sample_flag
      ]
    )

  library_strategy <-
    ae_core_collapse_unique(
      sample_summary$library_strategy
    )

  library_source <-
    ae_core_collapse_unique(
      sample_summary$library_source
    )

  library_selection <-
    ae_core_collapse_unique(
      sample_summary$library_selection
    )

  library_layout <-
    ae_core_collapse_unique(
      sample_summary$library_layout
    )

  sequencing_instrument <-
    ae_core_collapse_unique(
      sample_summary$sequencing_instrument
    )

  strandedness <-
    ae_core_collapse_unique(
      sample_summary$strandedness
    )

  screening_tier <-
    dplyr::case_when(
      n_nmd_perturbation_samples > 0L &&
        n_control_samples > 0L &&
        study_raw_flag ~
        "Strong NMD RNA-seq candidate",

      n_nmd_perturbation_samples > 0L &&
        study_raw_flag ~
        "NMD perturbation with raw RNA-seq; control needs review",

      n_nmd_perturbation_samples > 0L ~
        "NMD perturbation identified; raw availability needs review",

      TRUE ~
        "NMD perturbation not confirmed after full parsing"
    )

  tibble::tibble(
    Study_ID = accession,
    project_title =
      project_title,
    screening_tier =
      screening_tier,
    download_method =
      download_method,

    nmd_targets =
      study_targets,
    perturbation_methods =
      perturbation_methods,
    expected_nmd_effects =
      expected_nmd_effects,
    cell_models =
      cell_models,
    control_types =
      control_types,

    n_samples =
      n_samples,
    n_nmd_perturbation_samples =
      n_nmd_perturbation_samples,
    n_control_samples =
      n_control_samples,
    n_eligible_nmd_samples =
      n_eligible_samples,

    raw_rnaseq_available =
      study_raw_flag,

    library_strategy =
      library_strategy,
    library_source =
      library_source,
    library_selection =
      library_selection,
    library_layout =
      library_layout,
    sequencing_instrument =
      sequencing_instrument,
    strandedness =
      strandedness,

    direct_nmd_terms =
      ae_nmd_hits(
        study_text,
        AE_NMD_DIRECT_TERMS
      ),

    idf_nmd_target_terms =
      ae_nmd_hits(
        study_text,
        AE_NMD_TARGET_TERMS
      ),

    idf_genetic_perturbation_terms =
      ae_nmd_hits(
        study_text,
        AE_NMD_GENETIC_PERTURBATION_TERMS
      ),

    idf_direct_inhibitor_terms =
      ae_nmd_hits(
        study_text,
        AE_NMD_DIRECT_INHIBITOR_TERMS
      ),

    idf_translation_inhibitor_terms =
      ae_nmd_hits(
        study_text,
        AE_NMD_TRANSLATION_INHIBITOR_TERMS
      ),

    idf_stress_terms =
      ae_nmd_hits(
        study_text,
        AE_NMD_STRESS_TERMS
      )
  )
}


# -----------------------------------------------------------------------------
# 7. Parse one study
# -----------------------------------------------------------------------------

ae_screen_one_nmd_study <- function(
    accession,
    base_dir =
      "ArrayExpress_NMD_Screen",
    overwrite = FALSE,
    use_api_fallback = TRUE,
    query_biostudies_inventory = TRUE,
    verbose = TRUE,
    biostudies_json = NULL,
    prescreen_summary = NULL) {

  accession <- toupper(
    trimws(accession)
  )

  downloaded <-
    ae_core_download_magetab(
      accession = accession,
      base_dir = base_dir,
      overwrite = overwrite,
      use_api_fallback =
        use_api_fallback,
      verbose = verbose,
      biostudies_json =
        biostudies_json
    )

  study_dir <- downloaded$path

  sdrf_files <- list.files(
    study_dir,
    pattern = "sdrf(?:\\.txt)?$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  idf_files <- list.files(
    study_dir,
    pattern = "idf(?:\\.txt)?$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(sdrf_files) == 0L) {
    stop(
      "No SDRF file was found for ",
      accession,
      "."
    )
  }

  if (length(idf_files) == 0L) {
    stop(
      "No IDF file was found for ",
      accession,
      "."
    )
  }

  sdrf_parsed <- purrr::map(
    sdrf_files,
    ae_core_read_sdrf,
    accession = accession
  )

  sdrf_wide <- purrr::map_dfr(
    sdrf_parsed,
    "wide"
  )

  sdrf_long <- purrr::map_dfr(
    sdrf_parsed,
    "long"
  )

  idf_long <- purrr::map_dfr(
    idf_files,
    ae_core_read_idf,
    accession = accession
  )

  sample_summary <-
    ae_nmd_make_sample_summary(
      accession = accession,
      sdrf_long = sdrf_long,
      idf_long = idf_long,
      prescreen_summary =
        prescreen_summary
    )

  sdrf_inventory <-
    ae_core_sdrf_file_inventory(
      sdrf_long
    )

  api_inventory <-
    tibble::tibble()

  if (isTRUE(
    query_biostudies_inventory
  )) {
    api_inventory <-
      ae_core_biostudies_file_inventory(
        accession = accession,
        verbose = verbose,
        biostudies_json =
          biostudies_json
      )
  }

  file_inventory <-
    dplyr::bind_rows(
      sdrf_inventory,
      api_inventory
    ) %>%
    dplyr::distinct()

  evidence <-
    ae_nmd_make_evidence(
      idf_long = idf_long,
      sdrf_long = sdrf_long,
      file_inventory =
        file_inventory
    )

  project_title <-
    ae_core_extract_project_title(
      idf_long = idf_long,
      prescreen_summary =
        prescreen_summary
    )

  study_summary <-
    ae_nmd_make_study_summary(
      accession = accession,
      project_title =
        project_title,
      sample_summary =
        sample_summary,
      idf_long = idf_long,
      file_inventory =
        file_inventory,
      download_method =
        downloaded$download_method,
      prescreen_summary =
        prescreen_summary
    )

  eligible_sample_summary <-
    sample_summary %>%
    dplyr::filter(
      .data$eligible_nmd_sample_flag
    )

  perturbation_sample_summary <-
    sample_summary %>%
    dplyr::filter(
      .data$nmd_perturbation_flag
    )

  control_sample_summary <-
    sample_summary %>%
    dplyr::filter(
      .data$nmd_control_sample_flag
    )

  result <- list(
    Study_ID = accession,
    study_summary =
      study_summary,
    sample_summary =
      sample_summary,
    eligible_sample_summary =
      eligible_sample_summary,
    perturbation_sample_summary =
      perturbation_sample_summary,
    control_sample_summary =
      control_sample_summary,
    evidence = evidence,
    file_inventory =
      file_inventory,
    sdrf_wide =
      sdrf_wide,
    sdrf_long =
      sdrf_long,
    idf_long =
      idf_long,
    prescreen_summary =
      prescreen_summary
  )

  saveRDS(
    result,
    file = file.path(
      study_dir,
      paste0(
        accession,
        "_NMD_screening.rds"
      )
    )
  )

  result
}


# -----------------------------------------------------------------------------
# 8. Batch NMD screen
# -----------------------------------------------------------------------------

ae_screen_nmd_accessions <- function(
    accessions,
    base_dir =
      "ArrayExpress_NMD_Screen",
    overwrite = FALSE,
    use_api_fallback = TRUE,
    query_biostudies_inventory = TRUE,
    verbose = TRUE,
    prescreen = TRUE,
    prescreen_require_raw_evidence = TRUE,
    prescreen_require_perturbation_evidence = TRUE,
    prescreen_require_human_evidence = FALSE,
    prescreen_exclude_clear_nonhuman = TRUE,
    prescreen_require_cell_model_evidence = FALSE,
    prescreen_include_stress_perturbations = FALSE,
    prescreen_on_error =
      c("continue", "skip")) {

  prescreen_on_error <-
    match.arg(
      prescreen_on_error
    )

  accessions <- unique(
    toupper(
      trimws(
        as.character(
          accessions
        )
      )
    )
  )

  accessions <- accessions[
    !is.na(accessions) &
      nzchar(accessions)
  ]

  results <- stats::setNames(
    vector(
      "list",
      length(accessions)
    ),
    accessions
  )

  errors <- tibble::tibble(
    Study_ID = character(),
    stage = character(),
    error = character()
  )

  prescreen_summary <-
    tibble::tibble()

  prescreen_evidence <-
    tibble::tibble()

  skipped_prescreen <-
    tibble::tibble()

  for (accession in accessions) {

    if (verbose) {
      message(
        "\n========== ",
        accession,
        " =========="
      )
    }

    quick <- NULL
    preloaded_json <- NULL
    this_prescreen_summary <-
      NULL

    if (isTRUE(prescreen)) {

      quick <- tryCatch(
        ae_nmd_quick_prescreen(
          accession = accession,
          require_raw_evidence =
            prescreen_require_raw_evidence,
          require_perturbation_evidence =
            prescreen_require_perturbation_evidence,
          require_human_evidence =
            prescreen_require_human_evidence,
          exclude_clear_nonhuman =
            prescreen_exclude_clear_nonhuman,
          require_cell_model_evidence =
            prescreen_require_cell_model_evidence,
          include_stress_perturbations =
            prescreen_include_stress_perturbations,
          verbose = verbose
        ),
        error = function(e) e
      )

      if (inherits(
        quick,
        "error"
      )) {

        fallback_continue <-
          identical(
            prescreen_on_error,
            "continue"
          )

        errors <-
          dplyr::bind_rows(
            errors,
            tibble::tibble(
              Study_ID = accession,
              stage =
                "BioStudies NMD pre-screen",
              error =
                conditionMessage(quick)
            )
          )

        this_prescreen_summary <-
          tibble::tibble(
            Study_ID = accession,
            project_title =
              NA_character_,
            prescreen_status =
              ifelse(
                fallback_continue,
                "Continue: NMD pre-screen failed; full parser used as fallback",
                "Skip: NMD pre-screen failed"
              ),
            continue_full_parse =
              fallback_continue,
            nmd_flag = NA,
            perturbation_flag = NA,
            human_flag = NA,
            nonhuman_flag = NA,
            cell_model_flag = NA,
            rnaseq_flag = NA,
            raw_sequence_flag = NA
          )

        prescreen_summary <-
          dplyr::bind_rows(
            prescreen_summary,
            this_prescreen_summary
          )

        if (!fallback_continue) {
          skipped_prescreen <-
            dplyr::bind_rows(
              skipped_prescreen,
              this_prescreen_summary
            )
          next
        }

      } else {

        this_prescreen_summary <-
          quick$summary

        preloaded_json <-
          quick$json

        prescreen_summary <-
          dplyr::bind_rows(
            prescreen_summary,
            quick$summary
          )

        prescreen_evidence <-
          dplyr::bind_rows(
            prescreen_evidence,
            quick$evidence
          )

        if (verbose) {
          message(
            quick$summary$prescreen_status[1]
          )
        }

        if (!isTRUE(
          quick$summary$continue_full_parse[1]
        )) {

          skipped_prescreen <-
            dplyr::bind_rows(
              skipped_prescreen,
              quick$summary
            )

          next
        }
      }
    }

    result <- tryCatch(
      ae_screen_one_nmd_study(
        accession = accession,
        base_dir = base_dir,
        overwrite = overwrite,
        use_api_fallback =
          use_api_fallback,
        query_biostudies_inventory =
          query_biostudies_inventory,
        verbose = verbose,
        biostudies_json =
          preloaded_json,
        prescreen_summary =
          this_prescreen_summary
      ),
      error = function(e) e
    )

    if (inherits(
      result,
      "error"
    )) {

      errors <-
        dplyr::bind_rows(
          errors,
          tibble::tibble(
            Study_ID = accession,
            stage =
              "Full NMD MAGE-TAB parsing",
            error =
              conditionMessage(result)
          )
        )

      if (verbose) {
        message(
          "Failed ",
          accession,
          ": ",
          conditionMessage(result)
        )
      }

      next
    }

    results[[accession]] <-
      result
  }

  valid_results <- results[
    vapply(
      results,
      function(x) {
        is.list(x) &&
          !is.null(
            x$study_summary
          )
      },
      logical(1)
    )
  ]

  combined <- list(
    prescreen_summary =
      prescreen_summary,

    prescreen_evidence =
      prescreen_evidence,

    skipped_prescreen =
      skipped_prescreen,

    study_summary =
      purrr::map_dfr(
        valid_results,
        "study_summary"
      ),

    sample_summary =
      purrr::map_dfr(
        valid_results,
        "sample_summary"
      ),

    eligible_sample_summary =
      purrr::map_dfr(
        valid_results,
        "eligible_sample_summary"
      ),

    perturbation_sample_summary =
      purrr::map_dfr(
        valid_results,
        "perturbation_sample_summary"
      ),

    control_sample_summary =
      purrr::map_dfr(
        valid_results,
        "control_sample_summary"
      ),

    evidence =
      purrr::map_dfr(
        valid_results,
        "evidence"
      ),

    file_inventory =
      purrr::map_dfr(
        valid_results,
        "file_inventory"
      ),

    sdrf_wide =
      purrr::map_dfr(
        valid_results,
        "sdrf_wide"
      ),

    sdrf_long =
      purrr::map_dfr(
        valid_results,
        "sdrf_long"
      ),

    idf_long =
      purrr::map_dfr(
        valid_results,
        "idf_long"
      ),

    errors =
      errors,

    per_study =
      valid_results
  )

  if (
    nrow(
      combined$study_summary
    ) > 0L
  ) {
    combined$study_summary <-
      combined$study_summary %>%
      dplyr::arrange(
        factor(
          .data$screening_tier,
          levels = c(
            "Strong NMD RNA-seq candidate",
            "NMD perturbation with raw RNA-seq; control needs review",
            "NMD perturbation identified; raw availability needs review",
            "NMD perturbation not confirmed after full parsing"
          )
        ),
        .data$Study_ID
      )
  }

  combined
}


# -----------------------------------------------------------------------------
# 9. Write NMD output tables
# -----------------------------------------------------------------------------

ae_write_nmd_screening_results <- function(
    screening_results,
    out_dir =
      "ArrayExpress_NMD_Screening_Tables",
    prefix =
      "ArrayExpress_NMD") {

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  tables <- c(
    "prescreen_summary",
    "prescreen_evidence",
    "skipped_prescreen",
    "study_summary",
    "sample_summary",
    "eligible_sample_summary",
    "perturbation_sample_summary",
    "control_sample_summary",
    "evidence",
    "file_inventory",
    "sdrf_wide",
    "sdrf_long",
    "idf_long",
    "errors"
  )

  output_files <-
    stats::setNames(
      file.path(
        out_dir,
        paste0(
          prefix,
          "_",
          tables,
          ".tsv"
        )
      ),
      tables
    )

  for (table_name in tables) {

    table <-
      screening_results[[table_name]]

    if (is.null(table)) {
      table <- data.frame()
    }

    data.table::fwrite(
      table,
      file =
        output_files[[table_name]],
      sep = "\t",
      quote = TRUE,
      na = ""
    )
  }

  invisible(
    output_files
  )
}


# -----------------------------------------------------------------------------
# 10. Example workflow
# -----------------------------------------------------------------------------

# source(
#   "arrayexpress_workflow.R"
# )
#
# source(
#   "arrayexpress_nmd_perturbation_screening.R"
# )
#
# nmd_screen <- ae_screen_nmd_accessions(
#   accessions = nmd_search$accessions,
#
#   base_dir =
#     "ArrayExpress_NMD_Metadata",
#
#   prescreen = TRUE,
#
#   # Recommended for the first pass:
#   prescreen_require_raw_evidence = TRUE,
#   prescreen_require_perturbation_evidence = TRUE,
#
#   # Do not require the phrase "cell line" during the cheap API stage because
#   # it may only appear in the SDRF. Full sample parsing does require a cell
#   # model before a sample is marked eligible.
#   prescreen_require_cell_model_evidence = FALSE,
#
#   prescreen_require_human_evidence = FALSE,
#   prescreen_exclude_clear_nonhuman = TRUE,
#
#   # Stress perturbations are noisy, so leave them off unless specifically
#   # exploring stress-mediated NMD regulation.
#   prescreen_include_stress_perturbations = FALSE,
#
#   prescreen_on_error = "continue",
#   overwrite = FALSE,
#   verbose = TRUE
# )
#
# Study-level overview:
#
# View(
#   nmd_screen$study_summary
# )
#
# NMD perturbation samples:
#
# View(
#   nmd_screen$perturbation_sample_summary
# )
#
# Matching controls:
#
# View(
#   nmd_screen$control_sample_summary
# )
#
# Eligible perturbation + control samples:
#
# View(
#   nmd_screen$eligible_sample_summary
# )
#
# Quick summary:
#
# nmd_screen$eligible_sample_summary %>%
#   dplyr::count(
#     .data$Study_ID,
#     .data$sample_role,
#     .data$nmd_target_or_study_target,
#     .data$perturbation_type
#   )
#
# Export:
#
# ae_write_nmd_screening_results(
#   nmd_screen,
#   out_dir = "ArrayExpress_NMD_Tables"
# )


# =============================================================================
# D. Clean discovery wrapper used by run_nmd_workflow.R
# =============================================================================

ae_search_nmd <- function(
    context = nmd_default_context(),
    include_broad = TRUE,
    include_indirect = TRUE,
    include_stress = FALSE,
    page_size = 100L,
    pause_seconds = 0.15,
    verbose = TRUE) {

  queries <- build_arrayexpress_nmd_queries(
    context = context,
    include_broad = include_broad,
    include_indirect = include_indirect,
    include_stress = include_stress
  )

  biostudies_search_many(
    queries = queries,
    page_size = page_size,
    pause_seconds = pause_seconds,
    verbose = verbose
  )
}

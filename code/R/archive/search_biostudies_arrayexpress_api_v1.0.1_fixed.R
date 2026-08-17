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


# -----------------------------------------------------------------------------
# 7. Suggested AML / normal hematopoietic searches
# -----------------------------------------------------------------------------

arrayexpress_aml_queries <- c(
  AML_RNAseq = paste(
    '(',
    '"acute myeloid leukemia"',
    'OR "acute myeloid leukaemia"',
    'OR AML',
    'OR "acute promyelocytic leukemia"',
    'OR "acute promyelocytic leukaemia"',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR "transcriptome sequencing"',
    'OR "expression profiling by high throughput sequencing"',
    ')'
  ),

  AML_methylation = paste(
    '(',
    '"acute myeloid leukemia"',
    'OR "acute myeloid leukaemia"',
    'OR AML',
    'OR "acute promyelocytic leukemia"',
    'OR "acute promyelocytic leukaemia"',
    ')',
    'AND',
    '(',
    'methylation',
    'OR methylome',
    'OR bisulfite',
    'OR bisulphite',
    'OR RRBS',
    'OR ERRBS',
    'OR WGBS',
    'OR IDAT',
    ')'
  ),

  Normal_blood_RNAseq = paste(
    '(',
    '"healthy donor"',
    'OR "healthy control"',
    'OR "normal donor"',
    'OR "normal blood"',
    'OR "normal bone marrow"',
    'OR PBMC',
    'OR "cord blood"',
    'OR "hematopoietic stem cells"',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR "transcriptome sequencing"',
    ')'
  ),

  Normal_blood_methylation = paste(
    '(',
    '"healthy donor"',
    'OR "healthy control"',
    'OR "normal donor"',
    'OR "normal blood"',
    'OR "normal bone marrow"',
    'OR PBMC',
    'OR "cord blood"',
    'OR "hematopoietic stem cells"',
    ')',
    'AND',
    '(',
    'methylation',
    'OR methylome',
    'OR bisulfite',
    'OR bisulphite',
    'OR RRBS',
    'OR ERRBS',
    'OR WGBS',
    'OR IDAT',
    ')'
  )
)


# -----------------------------------------------------------------------------
# 8. Example
# -----------------------------------------------------------------------------

# Search every page from the four suggested queries:
#
# ae_search <- biostudies_search_many(
#   queries = arrayexpress_aml_queries,
#   page_size = 100,
#   pause_seconds = 0.15,
#   verbose = TRUE
# )
#
# View the deduplicated study table:
#
# View(ae_search$unique_studies)
#
# Plain accession vector:
#
# arrayexpress_ids <- ae_search$accessions
#
# Save TSV, TXT, and RDS outputs:
#
# write_biostudies_search_results(
#   search_result = ae_search,
#   output_prefix =
#     "AML_Normal_Blood_ArrayExpress"
# )
#
# Feed the accessions into the previously created metadata parser:
#
# source(
#   "arrayexpress_aml_normal_blood_screening_prefiltered.R"
# )
#
# aml_screen <- ae_screen_aml_accessions(
#   accessions = arrayexpress_ids,
#   prescreen = TRUE,
#   prescreen_require_raw_evidence = TRUE,
#   prescreen_include_normal_controls = TRUE,
#   prescreen_on_error = "continue"
# )
#
# Search the entire ArrayExpress collection.
# This can return tens of thousands of studies and many API pages:
#
# all_arrayexpress <- biostudies_search_all(
#   query = NULL,
#   query_name = "All ArrayExpress",
#   page_size = 100,
#   sort_by = "release_date",
#   sort_order = "descending"
# )

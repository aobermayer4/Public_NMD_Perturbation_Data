# =============================================================================
# GEO discovery, metadata retrieval, and NMD screening
# =============================================================================
# Public functions used by the driver:
#   geo_search_nmd()
#   geo_fetch_many_nmd()
#   parse_nmd_geo_metadata()
# =============================================================================

# ----------------------------------------------------------------------------
# GEO DataSets discovery through NCBI Entrez E-utilities
# ----------------------------------------------------------------------------

geo_eutils_json <- function(
  endpoint,
  query,
  retries = 4L,
  timeout_seconds = 90,
  verbose = FALSE
) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for GEO discovery.")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for GEO discovery.")
  }

  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
    endpoint,
    ".fcgi"
  )

  response <- httr::RETRY(
    "GET",
    url,
    query = query,
    httr::accept_json(),
    httr::user_agent("NMD perturbation GEO metadata screening workflow"),
    httr::timeout(timeout_seconds),
    times = retries,
    pause_min = 1,
    pause_cap = 10,
    quiet = !verbose
  )

  httr::stop_for_status(response)

  jsonlite::fromJSON(
    httr::content(response, as = "text", encoding = "UTF-8"),
    simplifyVector = TRUE,
    flatten = FALSE
  )
}


geo_first_nonempty <- function(x, default = NA_character_) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) default else x[[1]]
}


geo_extract_gse_accession <- function(summary_record) {
  preferred <- c("accession", "Accession", "gse", "GSE")

  for (nm in preferred) {
    if (!is.null(summary_record[[nm]])) {
      x <- toupper(as.character(summary_record[[nm]]))
      hit <- regmatches(x, regexpr("GSE[0-9]+", x, perl = TRUE))
      hit <- hit[nzchar(hit)]
      if (length(hit) > 0L) return(hit[[1]])
    }
  }

  flat <- toupper(
    paste(
      as.character(unlist(summary_record, recursive = TRUE, use.names = FALSE)),
      collapse = " "
    )
  )

  hit <- regmatches(flat, regexpr("GSE[0-9]+", flat, perl = TRUE))
  if (length(hit) == 0L || !nzchar(hit)) NA_character_ else hit
}


geo_summary_records <- function(payload) {
  result <- payload$result

  if (is.null(result)) {
    return(list())
  }

  uids <- result$uids
  if (is.null(uids)) {
    keep <- setdiff(names(result), "uids")
  } else {
    keep <- intersect(as.character(uids), names(result))
  }

  result[keep]
}


geo_search_one <- function(
  query,
  query_name = NULL,
  retmax = 10000L,
  summary_batch_size = 200L,
  api_key = NULL,
  email = NULL,
  pause_seconds = 0.4,
  verbose = TRUE
) {
  if (is.null(query_name) || !nzchar(query_name)) {
    query_name <- query
  }

  # NCBI's Entrez database used for searching GEO datasets is: 'gds'
  esearch_query <- list(
    db = "gds",
    term = query,
    retmax = as.integer(retmax),
    retmode = "json"
  )

  if (!is.null(api_key) && nzchar(api_key)) {
    esearch_query$api_key <- api_key
  }
  if (!is.null(email) && nzchar(email)) {
    esearch_query$email <- email
  }
  esearch_query$tool <- "NMDPerturbationWorkflow"

  if (verbose) {
    message("GEO search: ", query_name)
  }

  search_payload <- geo_eutils_json(
    "esearch",
    esearch_query,
    verbose = FALSE
  )

  # returns internal Entrez IDs for GSE projects
  ids <- as.character(search_payload$esearchresult$idlist)
  ids <- ids[!is.na(ids) & nzchar(ids)]

  if (length(ids) == 0L) {
    return(list(
      query = query,
      query_name = query_name,
      total_hits_reported = as.numeric(search_payload$esearchresult$count),
      hits = data.frame()
    ))
  }

  batches <- split(
    ids,
    ceiling(seq_along(ids) / as.integer(summary_batch_size))
  )

  hit_tables <- lapply(seq_along(batches), function(i) {
    if (i > 1L && pause_seconds > 0) {
      Sys.sleep(pause_seconds)
    }

    q <- list(
      db = "gds",
      id = paste(batches[[i]], collapse = ","),
      retmode = "json",
      version = "2.0",
      tool = "NMDPerturbationWorkflow"
    )

    if (!is.null(api_key) && nzchar(api_key)) {
      q$api_key <- api_key
    }
    if (!is.null(email) && nzchar(email)) {
      q$email <- email
    }

    payload <- geo_eutils_json("esummary", q, verbose = FALSE)
    records <- geo_summary_records(payload) # loops and unlist records that were actually retrieved

    if (length(records) == 0L) {
      return(data.frame())
    }

    # Format extracted GEO data as table
    rows <- lapply(records, function(rec) {
      accession <- geo_extract_gse_accession(rec)

      data.frame(
        accession = accession,
        uid = geo_first_nonempty(rec$uid),
        title = geo_first_nonempty(c(rec$title, rec$Title)),
        summary = geo_first_nonempty(c(rec$summary, rec$Summary)),
        organism = geo_first_nonempty(c(rec$taxon, rec$organism, rec$Organism)),
        gds_type = geo_first_nonempty(c(rec$gdstype, rec$GDS_Type, rec$type)),
        n_samples = suppressWarnings(as.numeric(
          geo_first_nonempty(c(rec$n_samples, rec$nsamples, rec$Samples))
        )),
        pub_date = geo_first_nonempty(c(rec$pdat, rec$PDAT, rec$pubdate)),
        matched_query_name = query_name,
        matched_query = query,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, rows)
  })

  hits <- do.call(rbind, hit_tables)
  if (is.null(hits)) {
    hits <- data.frame()
  }

  if (nrow(hits) > 0L) {
    hits <- hits[
      !is.na(hits$accession) & grepl("^GSE[0-9]+$", hits$accession),
      ,
      drop = FALSE
    ]
  }

  list(
    query = query,
    query_name = query_name,
    total_hits_reported = as.numeric(search_payload$esearchresult$count),
    hits = hits
  )
}


geo_search_many <- function(
  queries,
  retmax_per_query = 10000L, # for each query retrieve up to 10,000 NCBI records
  summary_batch_size = 200L, # retrieve NCBI IDs descriptive metadata 200 records at a time
  api_key = NULL, # optional
  email = NULL, # optional
  pause_seconds = 0.4, # brief pause between requests
  verbose = TRUE
) {
  queries <- as.character(queries)
  if (is.null(names(queries))) {
    names(queries) <- paste0("query_", seq_along(queries))
  }

  # Search and retain results of each query
  per_query <- lapply(seq_along(queries), function(i) {
    if (i > 1L && pause_seconds > 0) {
      Sys.sleep(pause_seconds)
    }

    geo_search_one(
      query = queries[[i]],
      query_name = names(queries)[[i]],
      retmax = retmax_per_query,
      summary_batch_size = summary_batch_size,
      api_key = api_key,
      email = email,
      pause_seconds = pause_seconds,
      verbose = verbose
    )
  })
  names(per_query) <- names(queries)

  # Combine query results
  all_hits <- do.call(
    rbind,
    lapply(per_query, function(x) x$hits)
  )

  if (is.null(all_hits) || nrow(all_hits) == 0L) {
    all_hits <- data.frame()
    unique_studies <- data.frame()
    accessions <- character()
  } else {
    split_hits <- split(all_hits, all_hits$accession)

    unique_studies <- do.call(
      rbind,
      lapply(split_hits, function(df) {
        data.frame(
          accession = df$accession[[1]],
          title = geo_first_nonempty(df$title),
          summary = geo_first_nonempty(df$summary),
          organism = geo_first_nonempty(df$organism),
          gds_type = geo_first_nonempty(df$gds_type),
          n_samples = suppressWarnings(max(df$n_samples, na.rm = TRUE)),
          pub_date = geo_first_nonempty(df$pub_date),
          matched_query_names = paste(
            unique(df$matched_query_name),
            collapse = "; "
          ),
          n_query_matches = length(unique(df$matched_query_name)),
          stringsAsFactors = FALSE
        )
      })
    )

    unique_studies$n_samples[is.infinite(unique_studies$n_samples)] <- NA_real_
    unique_studies <- unique_studies[
      order(-unique_studies$n_query_matches, unique_studies$accession),
      ,
      drop = FALSE
    ]
    rownames(unique_studies) <- NULL
    accessions <- unique_studies$accession
  }

  list(
    queries = queries,
    per_query = per_query,
    all_query_hits = all_hits,
    unique_studies = unique_studies,
    accessions = accessions
  )
}

# Build the GEO queries and hand them off to be searched
geo_search_nmd <- function(
  context = nmd_default_context(),
  include_broad = TRUE,
  include_indirect = TRUE,
  include_stress = FALSE,
  api_key = NULL,
  email = NULL,
  verbose = TRUE
) {
  queries <- build_geo_nmd_queries(
    context = context,
    include_broad = include_broad,
    include_indirect = include_indirect,
    include_stress = include_stress
  )

  geo_search_many(
    queries = queries,
    api_key = api_key,
    email = email,
    verbose = verbose
  )
}


# =============================================================================
# Existing GEO metadata parser, retained from the prior workflow
# =============================================================================
# NMD GEO metadata parser
# Designed for a named list of sample-level GEO metadata tables such as:
#   gse_metas[["GSE12345"]] -> one row per GSM/sample
#
# Required packages:
#   install.packages(c("dplyr", "tidyr", "stringr", "purrr", "tibble"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})


# -------------------------------------------------------------------------
# General helpers
# -------------------------------------------------------------------------

normalize_geo_field <- function(x) {
  x |>
    as.character() |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_replace_all("^_+|_+$", "")
}


collapse_unique <- function(x, sep = "; ", max_values = Inf) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x <- x[
    !is.na(x) &
      nzchar(x) &
      !tolower(x) %in% c("na", "n/a", "null", "none", "not applicable")
  ]

  x <- unique(x)

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


as_character_column <- function(x) {
  if (is.list(x)) {
    return(vapply(
      x,
      function(z) {
        collapse_unique(unlist(z, recursive = TRUE, use.names = FALSE))
      },
      character(1)
    ))
  }

  as.character(x)
}


dictionary_regex <- function(dictionary, selected = names(dictionary)) {
  selected <- intersect(selected, names(dictionary))

  if (length(selected) == 0L) {
    return("(?!)")
  }

  paste(unlist(dictionary[selected], use.names = FALSE), collapse = "|")
}


match_dictionary_one <- function(text, dictionary) {
  if (length(text) == 0L || is.na(text) || !nzchar(text)) {
    return(NA_character_)
  }

  hits <- names(dictionary)[vapply(
    dictionary,
    function(patterns) {
      stringr::str_detect(
        text,
        stringr::regex(paste(patterns, collapse = "|"), ignore_case = TRUE)
      )
    },
    logical(1)
  )]

  collapse_unique(hits)
}


match_dictionary <- function(text, dictionary) {
  vapply(text, match_dictionary_one, character(1), dictionary = dictionary)
}


has_dictionary_hit <- function(text, dictionary, selected = names(dictionary)) {
  pattern <- dictionary_regex(dictionary, selected)

  stringr::str_detect(
    dplyr::coalesce(text, ""),
    stringr::regex(pattern, ignore_case = TRUE)
  )
}


has_nearby_dictionary_pair <- function(
  text,
  dictionary_a,
  dictionary_b,
  max_chars = 150L
) {
  pattern_a <- dictionary_regex(dictionary_a)
  pattern_b <- dictionary_regex(dictionary_b)

  pair_pattern <- paste0(
    "(?:(?:",
    pattern_a,
    ").{0,",
    max_chars,
    "}(?:",
    pattern_b,
    ")",
    "|(?:",
    pattern_b,
    ").{0,",
    max_chars,
    "}(?:",
    pattern_a,
    "))"
  )

  stringr::str_detect(
    dplyr::coalesce(text, ""),
    stringr::regex(
      pair_pattern,
      ignore_case = TRUE,
      dotall = TRUE
    )
  )
}


# -------------------------------------------------------------------------
# Default dictionaries
#
# These are deliberately editable. GEO submitters use creative spelling,
# aliases, abbreviations, and occasionally metadata generated by chaos.
# -------------------------------------------------------------------------

NMD_TARGET_PATTERNS <- list(
  UPF1 = c("\\bUPF[-_ ]?1\\b", "\\bRENT[-_ ]?1\\b", "\\bSMG[-_ ]?2\\b"),
  UPF2 = c("\\bUPF[-_ ]?2\\b", "\\bRENT[-_ ]?2\\b", "\\bSMG[-_ ]?3\\b"),
  UPF3A = c("\\bUPF[-_ ]?3A\\b"),
  UPF3B = c("\\bUPF[-_ ]?3B\\b", "\\bUPF3X\\b"),
  SMG1 = c("\\bSMG[-_ ]?1\\b"),
  SMG5 = c("\\bSMG[-_ ]?5\\b"),
  SMG6 = c("\\bSMG[-_ ]?6\\b", "\\bEST1A\\b"),
  SMG7 = c("\\bSMG[-_ ]?7\\b"),
  SMG8 = c("\\bSMG[-_ ]?8\\b"),
  SMG9 = c("\\bSMG[-_ ]?9\\b"),
  DHX34 = c("\\bDHX[-_ ]?34\\b"),
  NBAS = c("\\bNBAS\\b"),
  SEC13 = c("\\bSEC[-_ ]?13\\b"),
  EIF4A3 = c("\\bEIF4A[-_ ]?3\\b"),
  RBM8A = c("\\bRBM[-_ ]?8A\\b", "\\bY14\\b"),
  MAGOH = c("\\bMAGOH\\b"),
  CASC3 = c("\\bCASC[-_ ]?3\\b", "\\bMLN51\\b"),
  RNPS1 = c("\\bRNPS[-_ ]?1\\b")
)


NMD_MECHANISM_PATTERNS <- list(
  "RNAi/knockdown" = c(
    "\\bsi\\s*RNA\\b",
    "\\bsh\\s*RNA\\b",
    "small interfering RNA",
    "short hairpin RNA",
    "knock\\s*-?down",
    "\\bKD\\b",
    "deplet(?:e|ed|ion)",
    "silenc(?:e|ed|ing)"
  ),
  "CRISPR/knockout" = c(
    "\\bCRISPR(?:[-/ ]?Cas9)?\\b",
    "knock\\s*-?out",
    "\\bKO\\b",
    "gene disruption",
    "null mutant"
  ),
  "degron/acute depletion" = c(
    "\\bdegron\\b",
    "\\bdTAG\\b",
    "auxin[- ]inducible",
    "induced degradation",
    "acute depletion"
  ),
  "small-molecule NMD inhibition" = c(
    "\\bNMDI[-_ ]?1\\b",
    "\\bNMDI[-_ ]?14\\b",
    "\\bSMG1i\\b",
    "SMG[-_ ]?1 inhibitor",
    "NMD inhibitor",
    "inhibition of nonsense[- ]mediated"
  ),
  "translation inhibition" = c(
    "\\bcycloheximide\\b",
    "\\bCHX\\b",
    "\\bemetine\\b",
    "\\banisomycin\\b",
    "\\bharringtonine\\b",
    "\\bpuromycin\\b"
  ),
  "overexpression/induction" = c(
    "over[- ]?express(?:ed|ion)?",
    "ectopic expression",
    "inducible expression",
    "transfect(?:ed|ion)",
    "transduc(?:ed|tion)"
  ),
  "rescue/complementation" = c(
    "\\brescue\\b",
    "complement(?:ed|ation)",
    "re[- ]?expression",
    "add[- ]back"
  ),
  "mutation/dominant negative" = c(
    "dominant[- ]negative",
    "\\bmutant\\b",
    "loss[- ]of[- ]function",
    "\\bLOF\\b"
  ),
  "readthrough treatment" = c(
    "\\bataluren\\b",
    "\\bPTC[-_ ]?124\\b",
    "\\bG418\\b",
    "\\bgeneticin\\b",
    "\\bgentamicin\\b",
    "read[- ]?through"
  )
)


NMD_AGENT_PATTERNS <- list(
  NMDI1 = c("\\bNMDI[-_ ]?1\\b"),
  NMDI14 = c("\\bNMDI[-_ ]?14\\b"),
  SMG1_inhibitor = c("\\bSMG1i\\b", "SMG[-_ ]?1 inhibitor"),
  cycloheximide = c("\\bcycloheximide\\b", "\\bCHX\\b"),
  emetine = c("\\bemetine\\b"),
  anisomycin = c("\\banisomycin\\b"),
  harringtonine = c("\\bharringtonine\\b"),
  puromycin = c("\\bpuromycin\\b"),
  ataluren_PTC124 = c("\\bataluren\\b", "\\bPTC[-_ ]?124\\b"),
  G418_geneticin = c("\\bG418\\b", "\\bgeneticin\\b"),
  gentamicin = c("\\bgentamicin\\b")
)


CONTROL_PATTERNS <- list(
  "untreated control" = c("\\buntreated\\b", "\\bno treatment\\b"),
  "vehicle control" = c("\\bvehicle\\b", "\\bDMSO\\b"),
  "mock control" = c("\\bmock\\b"),
  "non-targeting RNA control" = c(
    "non[- ]targeting",
    "scrambl(?:e|ed)",
    "\\bsiCTRL\\b",
    "\\bsiCON\\b",
    "\\bshCTRL\\b",
    "control siRNA",
    "control shRNA"
  ),
  "wild type" = c("\\bwild[- ]?type\\b", "\\bWT\\b"),
  "empty vector" = c("empty vector", "vector control")
)


NMD_EXPLICIT_PATTERN <- paste(
  c(
    "\\bnonsense[- ]mediated (?:mRNA |RNA )?decay\\b",
    "\\bnonsense[- ]mediated decay\\b",
    "\\bNMD\\b",
    "premature termination codon",
    "\\bPTC[- ]containing\\b"
  ),
  collapse = "|"
)


RNA_SEQ_PATTERN <- paste(
  c(
    "\\bRNA[- ]?seq\\b",
    "expression profiling by high throughput sequencing",
    "\\bIllumina\\b",
    "\\bSRA\\b",
    "\\bSRR[0-9]+\\b"
  ),
  collapse = "|"
)


PROTOCOL_FIELD_PATTERN <- paste(
  c(
    "treatment",
    "protocol",
    "description",
    "characteristics",
    "title",
    "source_name",
    "genotype",
    "transfection",
    "transduction",
    "growth",
    "extract",
    "library",
    "agent",
    "compound",
    "drug",
    "dose",
    "time"
  ),
  collapse = "|"
)


DISEASE_FIELD_PATTERN <- paste(
  c(
    "disease",
    "diagnosis",
    "disease_state",
    "health_status",
    "pathology",
    "tumou?r_type",
    "cancer_type",
    "clinical_status",
    "phenotype",
    "condition"
  ),
  collapse = "|"
)


TISSUE_FIELD_PATTERN <- paste(
  c(
    "(^|_)tissue($|_)",
    "(^|_)organ($|_)",
    "source_name",
    "cell_type",
    "cell_line",
    "sample_type",
    "anatom",
    "biospecimen",
    "primary_site"
  ),
  collapse = "|"
)


# Optional broad disease/context labels. These do not replace manual review.
DISEASE_KEYWORD_PATTERNS <- list(
  cancer = c(
    "\\bcancer\\b",
    "\\btumou?r\\b",
    "\\bcarcinoma\\b",
    "\\bleukemia\\b",
    "\\blymphoma\\b",
    "\\bmelanoma\\b",
    "\\bsarcoma\\b",
    "\\bglioma\\b"
  ),
  cystic_fibrosis = c("\\bcystic fibrosis\\b", "\\bCFTR\\b"),
  beta_thalassemia = c("\\bthalass?emia\\b", "\\bHBB\\b"),
  muscular_dystrophy = c("\\bmuscular dystrophy\\b", "\\bDMD\\b"),
  neurodegenerative_disease = c(
    "\\bneurodegenerat",
    "\\bALS\\b",
    "\\bAlzheimer",
    "\\bParkinson"
  ),
  developmental_disorder = c(
    "\\bdevelopmental disorder\\b",
    "\\bintellectual disability\\b",
    "\\bneurodevelopment"
  ),
  infection = c(
    "\\binfection\\b",
    "\\binfected\\b",
    "\\bviral\\b",
    "\\bbacterial\\b"
  )
)


# -------------------------------------------------------------------------
# Convert the named list of heterogeneous tables into one long table
# -------------------------------------------------------------------------

geo_metadata_to_long <- function(gse_metas) {
  if (!is.list(gse_metas)) {
    stop("gse_metas must be a list of data frames.")
  }

  if (is.null(names(gse_metas)) || any(!nzchar(names(gse_metas)))) {
    stop("gse_metas must be named with GSE accession IDs.")
  }

  # For each study metadata and it's GSE ID
  purrr::imap_dfr(gse_metas, function(df, gse_id) {
    if (length(df) == 1L && all(is.na(df))) {
      return(tibble())
    }

    if (!is.data.frame(df)) {
      df <- as.data.frame(df, stringsAsFactors = FALSE)
    }

    original_rownames <- rownames(df)

    # coerce all data to character
    df <- tibble::as_tibble(
      setNames(lapply(df, as_character_column), names(df)),
      .name_repair = "unique"
    )

    # identify candidate variables that define sample ID
    sample_candidates <- c(
      "Sample_ID",
      "sample_id",
      "geo_accession",
      "accession",
      "gsm",
      "GSM"
    )
    sample_col <- sample_candidates[sample_candidates %in% names(df)][1]

    # if no or unlikely sample ID column found, mutate into proper sample ID
    if (length(sample_col) == 0L || is.na(sample_col)) {
      sample_id <- original_rownames
      if (
        is.null(sample_id) ||
          length(sample_id) != nrow(df) ||
          all(sample_id %in% as.character(seq_len(nrow(df))))
      ) {
        sample_id <- paste0(gse_id, "_sample_", seq_len(nrow(df)))
      }
      df$Sample_ID <- sample_id
    } else if (sample_col != "Sample_ID") {
      df$Sample_ID <- df[[sample_col]]
    }

    df |>
      # Append the GSE accession to all sample rows
      mutate(
        GSE_ID = toupper(gse_id),
        Sample_ID = dplyr::coalesce(
          as.character(.data$Sample_ID),
          paste0(toupper(gse_id), "_sample_", row_number())
        )
      ) |>
      # Transform to long table to standardize
      pivot_longer(
        cols = -c(GSE_ID, Sample_ID),
        names_to = "field_raw",
        values_to = "value"
      ) |>
      # Preserve raw and normalized field names
      mutate(
        field_raw = as.character(.data$field_raw),
        field = normalize_geo_field(.data$field_raw),
        value = stringr::str_squish(as.character(.data$value)),
        field_value = paste0(.data$field_raw, ": ", .data$value)
      ) |>
      filter(
        !is.na(.data$value),
        nzchar(.data$value),
        !tolower(.data$value) %in%
          c("na", "n/a", "null", "none", "not applicable")
      )
  })
}


summarize_geo_fields <- function(long_meta, example_n = 5L) {
  long_meta |>
    group_by(.data$field, .data$field_raw) |>
    summarise(
      n_studies = n_distinct(.data$GSE_ID),
      n_samples = n_distinct(paste(.data$GSE_ID, .data$Sample_ID, sep = "::")),
      n_nonempty_values = n(),
      n_unique_values = n_distinct(.data$value),
      example_values = collapse_unique(.data$value, max_values = example_n),
      .groups = "drop"
    ) |>
    arrange(desc(.data$n_studies), desc(.data$n_samples), .data$field)
}


# -------------------------------------------------------------------------
# Build one searchable record per GSM/sample
# -------------------------------------------------------------------------

# trying to extract and format the main fields of interest frome each sample
## Using search term patterns to find fields that relate to specific categories, such as study title, tissue, disease type and others
build_geo_sample_corpus <- function(long_meta) {
  long_meta |>
    # Group by study and sample to format
    group_by(.data$GSE_ID, .data$Sample_ID) |>
    # Summarise the fields that may be represented by different field names
    summarise(
      study_title = collapse_unique(
        .data$value[.data$field %in% c("gse_title", "series_title")]
      ),
      sample_title = collapse_unique(
        .data$value[.data$field %in% c("title", "sample_title")]
      ),
      organism = collapse_unique(
        .data$value[stringr::str_detect(.data$field, "organism")]
      ),
      disease_field_values = collapse_unique(
        paste0(.data$field_raw, "=", .data$value)[
          stringr::str_detect(
            .data$field,
            stringr::regex(DISEASE_FIELD_PATTERN, ignore_case = TRUE)
          )
        ],
        max_values = 12
      ),
      tissue_field_values = collapse_unique(
        paste0(.data$field_raw, "=", .data$value)[
          stringr::str_detect(
            .data$field,
            stringr::regex(TISSUE_FIELD_PATTERN, ignore_case = TRUE)
          ) &
            !stringr::str_detect(.data$field, "organism")
        ],
        max_values = 12
      ),
      protocol_metadata = collapse_unique(
        .data$field_value[
          stringr::str_detect(
            .data$field,
            stringr::regex(PROTOCOL_FIELD_PATTERN, ignore_case = TRUE)
          )
        ],
        sep = " | ",
        max_values = 40
      ),
      all_metadata = collapse_unique(
        .data$field_value,
        sep = " | ",
        max_values = 100
      ),
      .groups = "drop"
    )
}


# -------------------------------------------------------------------------
# Sample-level annotation
# -------------------------------------------------------------------------

annotate_nmd_samples <- function(
  sample_corpus,
  target_patterns = NMD_TARGET_PATTERNS,
  mechanism_patterns = NMD_MECHANISM_PATTERNS,
  agent_patterns = NMD_AGENT_PATTERNS,
  control_patterns = CONTROL_PATTERNS,
  disease_patterns = DISEASE_KEYWORD_PATTERNS
) {
  genetic_mechanisms <- c(
    "RNAi/knockdown",
    "CRISPR/knockout",
    "degron/acute depletion",
    "overexpression/induction",
    "rescue/complementation",
    "mutation/dominant negative"
  )

  sample_corpus |>
    mutate(
      # format one giant searchable string per sample
      search_text = paste(
        dplyr::coalesce(.data$sample_title, ""),
        dplyr::coalesce(.data$protocol_metadata, ""),
        dplyr::coalesce(.data$all_metadata, ""),
        sep = " | "
      ),

      # looks for explicit NMD terminology, in case the NMD factors are not found, as it may still relate, could be new mech
      explicit_nmd_term = stringr::str_detect(
        .data$search_text,
        stringr::regex(NMD_EXPLICIT_PATTERN, ignore_case = TRUE)
      ),

      # match dictionary terms to sample search text and return the dictionary terms that matched
      ## Extract the matched terms raw value
      nmd_targets = match_dictionary(.data$search_text, target_patterns),
      perturbation_mechanisms = match_dictionary(
        .data$search_text,
        mechanism_patterns
      ),
      agents = match_dictionary(.data$search_text, agent_patterns),
      control_terms = match_dictionary(.data$search_text, control_patterns),
      disease_keyword_hits = match_dictionary(
        .data$search_text,
        disease_patterns
      ),

      # Annotate as boolean if a hit or not
      has_nmd_target = !is.na(.data$nmd_targets),
      # Is an NMD factor + perturbation mechanism within 150 characters of each other in search text?
      has_target_mechanism_pair = has_nearby_dictionary_pair(
        .data$search_text,
        target_patterns,
        mechanism_patterns,
        max_chars = 150L
      ),
      has_genetic_perturbation = has_dictionary_hit(
        .data$search_text,
        mechanism_patterns,
        genetic_mechanisms
      ),
      has_direct_nmd_inhibitor = has_dictionary_hit(
        .data$search_text,
        mechanism_patterns,
        "small-molecule NMD inhibition"
      ),
      has_translation_inhibitor = has_dictionary_hit(
        .data$search_text,
        mechanism_patterns,
        "translation inhibition"
      ),
      has_readthrough_treatment = has_dictionary_hit(
        .data$search_text,
        mechanism_patterns,
        "readthrough treatment"
      ),
      has_control_term = !is.na(.data$control_terms),

      rna_seq_evidence = stringr::str_detect(
        .data$search_text,
        stringr::regex(RNA_SEQ_PATTERN, ignore_case = TRUE)
      ),

      nmd_relevance = case_when(
        .data$has_target_mechanism_pair &
          .data$has_genetic_perturbation ~
          "direct NMD-factor perturbation candidate",

        .data$has_nmd_target & .data$has_genetic_perturbation ~
          "possible NMD-factor perturbation; target and mechanism separated",

        .data$has_direct_nmd_inhibitor ~
          "pharmacologic NMD-inhibition candidate",

        .data$has_translation_inhibitor & .data$explicit_nmd_term ~
          "indirect NMD inhibition by translation blockade",

        .data$has_readthrough_treatment ~
          "readthrough treatment; not necessarily direct NMD inhibition",

        .data$explicit_nmd_term & !is.na(.data$perturbation_mechanisms) ~
          "NMD perturbation candidate; target unclear",

        .data$has_nmd_target | .data$explicit_nmd_term ~
          "NMD-related sample; perturbation unclear",

        TRUE ~ "no clear NMD perturbation evidence"
      ),

      sample_role = case_when(
        .data$has_control_term &
          (.data$has_nmd_target |
            .data$has_direct_nmd_inhibitor |
            .data$has_translation_inhibitor |
            .data$has_genetic_perturbation) ~ "ambiguous: control and perturbation language both present",

        .data$has_control_term ~ "candidate control",

        .data$nmd_relevance %in%
          c(
            "direct NMD-factor perturbation candidate",
            "possible NMD-factor perturbation; target and mechanism separated",
            "pharmacologic NMD-inhibition candidate",
            "indirect NMD inhibition by translation blockade",
            "NMD perturbation candidate; target unclear"
          ) ~ "candidate perturbation",

        TRUE ~ "unclear"
      ),

      priority_score = 5L *
        as.integer(
          .data$nmd_relevance == "direct NMD-factor perturbation candidate"
        ) +
        4L *
          as.integer(
            .data$nmd_relevance ==
              "possible NMD-factor perturbation; target and mechanism separated"
          ) +
        5L *
          as.integer(
            .data$nmd_relevance == "pharmacologic NMD-inhibition candidate"
          ) +
        3L *
          as.integer(
            .data$nmd_relevance ==
              "indirect NMD inhibition by translation blockade"
          ) +
        2L * as.integer(.data$explicit_nmd_term) +
        1L * as.integer(.data$rna_seq_evidence) +
        1L * as.integer(!is.na(.data$disease_field_values)) +
        1L * as.integer(!is.na(.data$tissue_field_values)),

      evidence_confidence = case_when(
        (.data$has_target_mechanism_pair &
          .data$has_genetic_perturbation) |
          .data$has_direct_nmd_inhibitor ~ "high",

        .data$has_nmd_target & .data$has_genetic_perturbation ~ "moderate",

        .data$explicit_nmd_term &
          (.data$has_translation_inhibitor |
            !is.na(.data$perturbation_mechanisms)) ~ "moderate",

        .data$has_nmd_target | .data$explicit_nmd_term ~ "low",

        TRUE ~ "none"
      ),

      manual_review_reason = case_when(
        .data$sample_role ==
          "ambiguous: control and perturbation language both present" ~
          "Control and perturbation terms appear in the same sample record.",

        .data$has_translation_inhibitor &
          stringr::str_detect(
            dplyr::coalesce(.data$agents, ""),
            stringr::regex("puromycin", ignore_case = TRUE)
          ) ~
          "Puromycin may represent selection rather than intentional NMD inhibition.",

        .data$has_readthrough_treatment ~
          "Readthrough therapy changes PTC recognition but is not equivalent to direct NMD-factor inhibition.",

        .data$evidence_confidence == "low" ~
          "NMD terminology or an NMD factor is present, but the perturbation is not explicit.",

        TRUE ~ NA_character_
      )
    ) |>
    select(-.data$search_text) |>
    arrange(desc(.data$priority_score), .data$GSE_ID, .data$Sample_ID)
}


# -------------------------------------------------------------------------
# Evidence table: retains exact GEO fields and values that caused hits
# -------------------------------------------------------------------------

build_nmd_evidence_table <- function(
  long_meta,
  target_patterns = NMD_TARGET_PATTERNS,
  mechanism_patterns = NMD_MECHANISM_PATTERNS,
  agent_patterns = NMD_AGENT_PATTERNS,
  control_patterns = CONTROL_PATTERNS
) {
  target_regex <- dictionary_regex(target_patterns)
  mechanism_regex <- dictionary_regex(mechanism_patterns)
  agent_regex <- dictionary_regex(agent_patterns)
  control_regex <- dictionary_regex(control_patterns)

  long_meta |>
    mutate(
      explicit_nmd_term = stringr::str_detect(
        .data$field_value,
        stringr::regex(NMD_EXPLICIT_PATTERN, ignore_case = TRUE)
      ),
      nmd_target = stringr::str_detect(
        .data$field_value,
        stringr::regex(target_regex, ignore_case = TRUE)
      ),
      perturbation_mechanism = stringr::str_detect(
        .data$field_value,
        stringr::regex(mechanism_regex, ignore_case = TRUE)
      ),
      named_agent = stringr::str_detect(
        .data$field_value,
        stringr::regex(agent_regex, ignore_case = TRUE)
      ),
      control_term = stringr::str_detect(
        .data$field_value,
        stringr::regex(control_regex, ignore_case = TRUE)
      ),
      disease_field = stringr::str_detect(
        .data$field,
        stringr::regex(DISEASE_FIELD_PATTERN, ignore_case = TRUE)
      ),
      tissue_field = stringr::str_detect(
        .data$field,
        stringr::regex(TISSUE_FIELD_PATTERN, ignore_case = TRUE)
      ) &
        !stringr::str_detect(.data$field, "organism")
    ) |>
    pivot_longer(
      cols = c(
        explicit_nmd_term,
        nmd_target,
        perturbation_mechanism,
        named_agent,
        control_term,
        disease_field,
        tissue_field
      ),
      names_to = "evidence_type",
      values_to = "matched"
    ) |>
    filter(.data$matched) |>
    select(
      .data$GSE_ID,
      .data$Sample_ID,
      .data$evidence_type,
      .data$field_raw,
      .data$value
    ) |>
    distinct() |>
    arrange(.data$GSE_ID, .data$Sample_ID, .data$evidence_type)
}


# -------------------------------------------------------------------------
# Study-level summary
# -------------------------------------------------------------------------

summarize_nmd_studies <- function(sample_annotations) {
  sample_annotations |>
    group_by(.data$GSE_ID) |>
    summarise(
      study_title = collapse_unique(.data$study_title),
      organism = collapse_unique(.data$organism),

      n_samples = n(),
      n_candidate_perturbation = sum(
        .data$sample_role == "candidate perturbation",
        na.rm = TRUE
      ),
      n_candidate_controls = sum(
        .data$sample_role == "candidate control",
        na.rm = TRUE
      ),
      n_ambiguous_samples = sum(
        stringr::str_detect(
          .data$sample_role,
          stringr::fixed("ambiguous")
        ),
        na.rm = TRUE
      ),

      nmd_targets = collapse_unique(.data$nmd_targets),
      perturbation_mechanisms = collapse_unique(.data$perturbation_mechanisms),
      agents = collapse_unique(.data$agents),

      nmd_relevance = collapse_unique(.data$nmd_relevance),
      disease_field_values = collapse_unique(
        .data$disease_field_values,
        max_values = 15
      ),
      disease_keyword_hits = collapse_unique(.data$disease_keyword_hits),
      tissue_field_values = collapse_unique(
        .data$tissue_field_values,
        max_values = 15
      ),

      has_rna_seq_evidence = any(.data$rna_seq_evidence, na.rm = TRUE),
      max_sample_priority = max(.data$priority_score, na.rm = TRUE),

      study_priority_score = max(.data$priority_score, na.rm = TRUE) +
        2L *
          as.integer(
            any(.data$sample_role == "candidate perturbation", na.rm = TRUE)
          ) +
        2L *
          as.integer(
            any(.data$sample_role == "candidate control", na.rm = TRUE)
          ) +
        1L *
          as.integer(
            any(.data$rna_seq_evidence, na.rm = TRUE)
          ),

      review_category = case_when(
        .data$n_candidate_perturbation > 0L &
          .data$n_candidate_controls > 0L ~
          "highest priority: candidate perturbation and control arms",

        .data$n_candidate_perturbation > 0L ~
          "candidate perturbation study; control arm not identified",

        any(
          .data$evidence_confidence %in% c("high", "moderate"),
          na.rm = TRUE
        ) ~
          "NMD-related study requiring sample-role review",

        TRUE ~ "low-priority or unclear"
      ),

      .groups = "drop"
    ) |>
    arrange(
      factor(
        .data$review_category,
        levels = c(
          "highest priority: candidate perturbation and control arms",
          "candidate perturbation study; control arm not identified",
          "NMD-related study requiring sample-role review",
          "low-priority or unclear"
        )
      ),
      desc(.data$study_priority_score),
      .data$GSE_ID
    )
}


# -------------------------------------------------------------------------
# Interactive metadata search
# -------------------------------------------------------------------------

search_geo_metadata <- function(
  long_meta,
  terms,
  fields = NULL,
  require_all_terms = FALSE,
  ignore_case = TRUE
) {
  terms <- terms[!is.na(terms) & nzchar(terms)]

  if (length(terms) == 0L) {
    stop("Provide at least one search term or regular expression.")
  }

  out <- long_meta

  if (!is.null(fields)) {
    field_pattern <- paste(fields, collapse = "|")
    out <- out |>
      filter(
        stringr::str_detect(
          .data$field,
          stringr::regex(field_pattern, ignore_case = ignore_case)
        )
      )
  }

  term_hits <- lapply(terms, function(term) {
    stringr::str_detect(
      out$field_value,
      stringr::regex(term, ignore_case = ignore_case)
    )
  })

  keep <- if (require_all_terms) {
    Reduce(`&`, term_hits)
  } else {
    Reduce(`|`, term_hits)
  }

  out |>
    filter(keep) |>
    arrange(.data$GSE_ID, .data$Sample_ID, .data$field)
}


# -------------------------------------------------------------------------
# Main wrapper and export
# -------------------------------------------------------------------------

parse_nmd_geo_metadata <- function(
  gse_metas,
  target_patterns = NMD_TARGET_PATTERNS,
  mechanism_patterns = NMD_MECHANISM_PATTERNS,
  agent_patterns = NMD_AGENT_PATTERNS,
  control_patterns = CONTROL_PATTERNS,
  disease_patterns = DISEASE_KEYWORD_PATTERNS
) {
  long_meta <- geo_metadata_to_long(gse_metas)
  field_summary <- summarize_geo_fields(long_meta)
  sample_corpus <- build_geo_sample_corpus(long_meta)

  # Now start to apply biological logic
  sample_annotations <- annotate_nmd_samples(
    sample_corpus = sample_corpus,
    target_patterns = target_patterns,
    mechanism_patterns = mechanism_patterns,
    agent_patterns = agent_patterns,
    control_patterns = control_patterns,
    disease_patterns = disease_patterns
  )

  evidence <- build_nmd_evidence_table(
    long_meta = long_meta,
    target_patterns = target_patterns,
    mechanism_patterns = mechanism_patterns,
    agent_patterns = agent_patterns,
    control_patterns = control_patterns
  )

  study_summary <- summarize_nmd_studies(sample_annotations)

  list(
    studies = study_summary,
    samples = sample_annotations,
    evidence = evidence,
    fields = field_summary,
    long_metadata = long_meta
  )
}


write_nmd_geo_results <- function(
  parsed,
  out_dir = "NMD_GEO_Metadata_Review"
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  outputs <- list(
    studies = file.path(out_dir, "NMD_GEO_study_summary.tsv"),
    samples = file.path(out_dir, "NMD_GEO_sample_annotations.tsv"),
    evidence = file.path(out_dir, "NMD_GEO_evidence.tsv"),
    fields = file.path(out_dir, "NMD_GEO_field_summary.tsv"),
    long_metadata = file.path(out_dir, "NMD_GEO_metadata_long.tsv")
  )

  for (nm in names(outputs)) {
    utils::write.table(
      parsed[[nm]],
      file = outputs[[nm]],
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE,
      na = ""
    )
  }

  invisible(outputs)
}


# -------------------------------------------------------------------------
# Example usage with your existing objects
# -------------------------------------------------------------------------

# source("nmd_geo_metadata_parser.R")
#
# parsed_nmd <- parse_nmd_geo_metadata(gse_metas)
#
# View(parsed_nmd$studies)
# View(parsed_nmd$samples)
# View(parsed_nmd$evidence)
# View(parsed_nmd$fields)
#
# write_nmd_geo_results(
#   parsed_nmd,
#   out_dir = "NMD_GEO_Metadata_Review_20260712"
# )
#
# High-priority candidate studies:
# parsed_nmd$studies |>
#   filter(review_category != "low-priority or unclear")
#
# Inspect all evidence for one GSE:
# parsed_nmd$evidence |>
#   filter(GSE_ID == "GSE12345")
#
# Search arbitrary metadata terms:
# search_geo_metadata(
#   parsed_nmd$long_metadata,
#   terms = c("UPF1", "knockdown")
# )
#
# Add an additional compound or alias:
# my_agents <- NMD_AGENT_PATTERNS
# my_agents$custom_compound <- c("\\bCompoundName\\b", "\\balias\\b")
#
# parsed_nmd <- parse_nmd_geo_metadata(
#   gse_metas,
#   agent_patterns = my_agents
# )

# -------------------------------------------------------------------------
# Optional improved GEO extraction
#
# Your existing saved gse_metas can be parsed without re-downloading.
# For future downloads, this version also retains the GSE summary and
# overall design, and it parses all characteristics fields using only the
# first colon as the key/value separator.
# -------------------------------------------------------------------------

geo_retry <- function(expr, n = 3L, wait = 10, verbose = TRUE) {
  for (i in seq_len(n)) {
    result <- tryCatch(expr, error = function(e) e)

    if (!inherits(result, "error")) {
      return(result)
    }

    if (verbose) {
      message(
        sprintf("Attempt %d/%d failed: %s", i, n, result$message)
      )
    }

    if (i < n) {
      Sys.sleep(wait)
    }
  }

  stop(result)
}


parse_characteristics_vector <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x <- x[!is.na(x) & nzchar(x)]

  if (length(x) == 0L) {
    return(list())
  }

  has_colon <- stringr::str_detect(x, ":")

  keys <- ifelse(
    has_colon,
    stringr::str_replace(x, "^\\s*([^:]+):.*$", "\\1"),
    "unparsed_characteristic"
  )

  values <- ifelse(
    has_colon,
    stringr::str_replace(x, "^\\s*[^:]+:\\s*", ""),
    x
  )

  keys <- normalize_geo_field(keys)

  split_values <- split(values, keys)

  lapply(split_values, collapse_unique)
}

# Extracts GEO study sample data -> outputs wide sample level table
extract_gse_meta_v2 <- function(
  gse_id,
  verbose = FALSE,
  retries = 3L,
  retry_wait = 10
) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop(
      "The GEOquery package is required. Install it with ",
      "BiocManager::install('GEOquery')."
    )
  }

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = 300)

  if (verbose) {
    message("Fetching ", gse_id)
  }

  gse <- geo_retry(
    GEOquery::getGEO(gse_id, GSEMatrix = FALSE),
    n = retries,
    wait = retry_wait,
    verbose = verbose
  )

  series_meta <- GEOquery::Meta(gse)

  series_fields <- list(
    GSE_Title = collapse_unique(series_meta$title),
    GSE_Summary = collapse_unique(series_meta$summary),
    GSE_Overall_Design = collapse_unique(series_meta$overall_design),
    GSE_Type = collapse_unique(series_meta$type),
    GSE_PubMed_ID = collapse_unique(series_meta$pubmed_id),
    GSE_Organism = collapse_unique(series_meta$organism)
  )

  samples <- GEOquery::GSMList(gse)

  # Loop over the list of samples for each single sample object and the sample name
  # _dfr will end up row-binding the results
  purrr::imap_dfr(samples, function(sample_obj, gsm_id) {
    sample_meta <- GEOquery::Meta(sample_obj)

    characteristic_idx <- grep(
      "^characteristics_ch",
      names(sample_meta),
      ignore.case = TRUE
    )

    parsed_characteristics <- list()

    if (length(characteristic_idx) > 0L) {
      parsed_characteristics <- unlist(
        lapply(
          sample_meta[characteristic_idx],
          parse_characteristics_vector
        ),
        recursive = FALSE,
        use.names = TRUE
      )

      if (length(parsed_characteristics) > 0L) {
        names(parsed_characteristics) <- make.unique(
          names(parsed_characteristics),
          sep = "_duplicate_"
        )
      }
    }

    remaining_meta <- sample_meta[-characteristic_idx]

    remaining_meta <- lapply(
      remaining_meta,
      function(x) {
        collapse_unique(unlist(x, recursive = TRUE, use.names = FALSE))
      }
    )

    names(remaining_meta) <- make.unique(
      names(remaining_meta),
      sep = "_duplicate_"
    )

    row_list <- c(
      list(Sample_ID = gsm_id),
      series_fields,
      remaining_meta,
      parsed_characteristics
    )

    tibble::as_tibble_row(row_list, .name_repair = "unique")
  })
}


# Example future re-extraction:
#
# gse_metas_v2 <- lapply(GSE_IDs, extract_gse_meta_v2, verbose = TRUE)
# gse_metas_v2 <- setNames(gse_metas_v2, GSE_IDs)
# saveRDS(
#   gse_metas_v2,
#   file = "NMD_Perturbation_GSE_Metas_v2_20260712.rds"
# )
#
# parsed_nmd_v2 <- parse_nmd_geo_metadata(gse_metas_v2)

# =============================================================================
# Simple batch retrieval wrapper for the cleaned workflow
# =============================================================================

geo_fetch_many_nmd <- function(
  accessions,
  verbose = TRUE,
  cache_file = NULL
) {
  accessions <- unique(toupper(trimws(as.character(accessions))))
  accessions <- accessions[grepl("^GSE[0-9]+$", accessions)]

  if (length(accessions) == 0L) {
    return(list(
      metadata = list(),
      errors = data.frame(
        GSE_ID = character(),
        error = character(),
        stringsAsFactors = FALSE
      )
    ))
  }

  metadata <- vector("list", length(accessions))
  names(metadata) <- accessions
  errors <- list()

  for (i in seq_along(accessions)) {
    gse_id <- accessions[[i]]

    if (verbose) {
      message("[", i, "/", length(accessions), "] GEO metadata: ", gse_id)
    }

    ans <- tryCatch(
      extract_gse_meta_v2(gse_id, verbose = FALSE),
      error = function(e) e
    )

    # If extracting the GSE data resulted in error, write to error file
    if (inherits(ans, "error")) {
      errors[[length(errors) + 1L]] <- data.frame(
        GSE_ID = gse_id,
        error = conditionMessage(ans),
        stringsAsFactors = FALSE
      )
      metadata[[gse_id]] <- NULL
    } else {
      metadata[[gse_id]] <- ans
    }
  }

  metadata <- metadata[!vapply(metadata, is.null, logical(1))]

  error_df <- if (length(errors) == 0L) {
    data.frame(
      GSE_ID = character(),
      error = character(),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, errors)
  }

  if (!is.null(cache_file) && length(metadata) > 0L) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(metadata, cache_file)
  }

  list(metadata = metadata, errors = error_df)
}


geo_run_nmd_screen <- function(
  accessions,
  verbose = TRUE,
  cache_file = NULL
) {
  fetched <- geo_fetch_many_nmd(
    accessions = accessions,
    verbose = verbose,
    cache_file = cache_file
  )

  parsed <- if (length(fetched$metadata) == 0L) {
    list(
      studies = data.frame(),
      samples = data.frame(),
      evidence = data.frame(),
      fields = data.frame(),
      long_metadata = data.frame()
    )
  } else {
    parse_nmd_geo_metadata(fetched$metadata)
  }

  parsed$errors <- fetched$errors
  parsed$raw_metadata <- fetched$metadata
  parsed
}

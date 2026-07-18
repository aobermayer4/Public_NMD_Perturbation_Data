# ArrayExpress MAGE-TAB metadata downloader and parser
#
# Purpose:
#   - Download ArrayExpress study metadata for a vector of accessions
#   - Read SDRF sample metadata
#   - Read IDF study and protocol metadata
#   - Create long-format tables that can be searched like GEO metadata
#
# The function first tries ArrayExpress::getAE(type = "mage"), which is
# supported by newer ArrayExpress releases. If an older package such as
# version 1.32.0 fails, it can fall back to the current BioStudies API and
# download only the IDF and SDRF files.
#
# Required packages:
#   ArrayExpress, data.table, dplyr, tidyr, purrr, stringr, tibble, jsonlite
#
# Install:
#   BiocManager::install("ArrayExpress")
#   install.packages(c(
#     "data.table", "dplyr", "tidyr", "purrr",
#     "stringr", "tibble", "jsonlite"
#   ))

suppressPackageStartupMessages({
  library(ArrayExpress)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(jsonlite)
})


# -------------------------------------------------------------------------
# General helpers
# -------------------------------------------------------------------------

ae_normalize_field <- function(x) {
  x |>
    as.character() |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_replace_all("^_+|_+$", "")
}


ae_collapse_unique <- function(x, sep = "; ", max_values = Inf) {
  x <- as.character(x)
  x <- stringr::str_squish(x)

  x <- x[
    !is.na(x) &
      nzchar(x) &
      !tolower(x) %in% c(
        "na", "n/a", "null", "none", "not applicable"
      )
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


ae_first_nonempty <- function(values, default = NA_character_) {
  values <- as.character(values)
  values <- stringr::str_squish(values)

  keep <- !is.na(values) & nzchar(values)

  if (!any(keep)) {
    return(default)
  }

  values[which(keep)[1]]
}


ae_retry <- function(expr_fun, n = 3L, wait = 5, verbose = TRUE) {
  last_error <- NULL

  for (i in seq_len(n)) {
    ans <- tryCatch(
      expr_fun(),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )

    if (!is.null(ans)) {
      return(ans)
    }

    if (verbose && !is.null(last_error)) {
      message(
        sprintf(
          "Attempt %d/%d failed: %s",
          i,
          n,
          conditionMessage(last_error)
        )
      )
    }

    if (i < n) {
      Sys.sleep(wait)
    }
  }

  if (!is.null(last_error)) {
    stop(last_error)
  }

  stop("Operation failed without returning an error object.")
}


ae_resolve_paths <- function(paths, study_dir) {
  paths <- unique(as.character(unlist(paths, use.names = FALSE)))
  paths <- paths[!is.na(paths) & nzchar(paths)]

  if (length(paths) == 0L) {
    return(character())
  }

  resolved <- vapply(
    paths,
    function(x) {
      candidates <- unique(c(
        x,
        file.path(study_dir, x),
        file.path(study_dir, basename(x))
      ))

      existing <- candidates[file.exists(candidates)]

      if (length(existing) == 0L) {
        return(NA_character_)
      }

      normalizePath(existing[1], winslash = "/", mustWork = TRUE)
    },
    character(1)
  )

  unique(resolved[!is.na(resolved)])
}


# -------------------------------------------------------------------------
# Metadata-only download
# -------------------------------------------------------------------------

ae_download_with_package <- function(
    accession,
    study_dir,
    overwrite = FALSE,
    verbose = TRUE) {

  getae_formals <- names(formals(ArrayExpress::getAE))

  args <- list(
    accession = accession,
    path = study_dir,
    type = "mage",
    extract = TRUE
  )

  if ("overwrite" %in% getae_formals) {
    args$overwrite <- overwrite
  }

  if (verbose) {
    message(
      "Trying ArrayExpress::getAE(type = \"mage\") for ",
      accession
    )
  }

  result <- do.call(ArrayExpress::getAE, args)

  if (is.null(result)) {
    stop(
      "ArrayExpress::getAE() returned NULL. ",
      "The installed package may not support type = \"mage\"."
    )
  }

  result
}


ae_json_from_url <- function(url, retries = 3L, wait = 5) {
  ae_retry(
    function() {
      jsonlite::fromJSON(
        url,
        simplifyVector = TRUE
      )
    },
    n = retries,
    wait = wait,
    verbose = TRUE
  )
}


ae_download_file <- function(
    url,
    destination,
    overwrite = FALSE,
    retries = 3L,
    wait = 5) {

  if (file.exists(destination) && !overwrite) {
    return(normalizePath(
      destination,
      winslash = "/",
      mustWork = TRUE
    ))
  }

  dir.create(
    dirname(destination),
    recursive = TRUE,
    showWarnings = FALSE
  )

  ae_retry(
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
    },
    n = retries,
    wait = wait,
    verbose = TRUE
  )
}


ae_download_from_biostudies <- function(
    accession,
    study_dir,
    overwrite = FALSE,
    verbose = TRUE) {

  if (verbose) {
    message(
      "Using BioStudies metadata-only fallback for ",
      accession
    )
  }

  study_url <- paste0(
    "https://www.ebi.ac.uk/biostudies/api/v1/studies/",
    accession
  )

  info_url <- paste0(study_url, "/info")

  study_json <- ae_json_from_url(study_url)
  info_json <- ae_json_from_url(info_url)

  # The study JSON contains the paths of all files. Restrict this to IDF
  # and SDRF files so processed matrices and raw sequence files are not
  # downloaded merely because humans enjoy filling disks.
  all_strings <- unique(as.character(
    unlist(study_json, recursive = TRUE, use.names = FALSE)
  ))

  magetab_paths <- all_strings[
    stringr::str_detect(
      basename(all_strings),
      stringr::regex(
        "(?:^|\\.)(?:idf|sdrf)(?:\\.txt)?$|(?:idf|sdrf)\\.txt$",
        ignore_case = TRUE
      )
    )
  ]

  magetab_paths <- unique(magetab_paths[
    !is.na(magetab_paths) & nzchar(magetab_paths)
  ])

  if (length(magetab_paths) == 0L) {
    stop(
      "No IDF or SDRF paths were found in the BioStudies record for ",
      accession,
      "."
    )
  }

  ftp_link <- NULL

  if (!is.null(info_json$ftpLink)) {
    ftp_link <- as.character(info_json$ftpLink)[1]
  }

  if (is.null(ftp_link) || is.na(ftp_link) || !nzchar(ftp_link)) {
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
        message("Downloading metadata file: ", basename(relative_path))
      }

      ae_download_file(
        url = file_url,
        destination = destination,
        overwrite = overwrite
      )
    },
    character(1)
  )

  sdrf <- downloaded[
    stringr::str_detect(
      basename(downloaded),
      stringr::regex("sdrf(?:\\.txt)?$", ignore_case = TRUE)
    )
  ]

  idf <- downloaded[
    stringr::str_detect(
      basename(downloaded),
      stringr::regex("idf(?:\\.txt)?$", ignore_case = TRUE)
    )
  ]

  list(
    path = study_dir,
    rawFiles = NULL,
    rawArchive = NULL,
    processedFiles = NULL,
    processedArchive = NULL,
    mageTabFiles = downloaded,
    sdrf = sdrf,
    idf = idf,
    adf = NULL,
    dataFiles = data.frame(
      type = "MAGE-TAB Files",
      file = basename(downloaded),
      url = NA_character_,
      stringsAsFactors = FALSE
    )
  )
}


ae_download_magetab <- function(
    accession,
    base_dir = "ArrayExpress_Metadata",
    overwrite = FALSE,
    use_api_fallback = TRUE,
    verbose = TRUE) {

  accession <- toupper(trimws(accession))

  if (!stringr::str_detect(
    accession,
    "^E-[A-Z0-9]+-[0-9]+$"
  )) {
    stop("Invalid-looking ArrayExpress accession: ", accession)
  }

  study_dir <- file.path(base_dir, accession)

  dir.create(
    study_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  package_result <- tryCatch(
    ae_download_with_package(
      accession = accession,
      study_dir = study_dir,
      overwrite = overwrite,
      verbose = verbose
    ),
    error = function(e) e
  )

  if (!inherits(package_result, "error")) {
    package_result$download_method <-
      "ArrayExpress::getAE(type = 'mage')"
    return(package_result)
  }

  if (!use_api_fallback) {
    stop(
      "ArrayExpress package download failed: ",
      conditionMessage(package_result),
      "\nThe installed package may predate the BioStudies migration."
    )
  }

  if (verbose) {
    message(
      "Package download failed: ",
      conditionMessage(package_result)
    )
  }

  api_result <- ae_download_from_biostudies(
    accession = accession,
    study_dir = study_dir,
    overwrite = overwrite,
    verbose = verbose
  )

  api_result$download_method <- "BioStudies API fallback"
  api_result
}


# -------------------------------------------------------------------------
# SDRF parsing
# -------------------------------------------------------------------------

ae_make_sample_id <- function(df, original_names) {
  preferred_fields <- c(
    "Sample Name",
    "Source Name",
    "Assay Name",
    "Extract Name",
    "Scan Name"
  )

  sample_id <- rep(NA_character_, nrow(df))

  for (field in preferred_fields) {
    matching_columns <- names(df)[
      tolower(original_names) == tolower(field)
    ]

    for (column in matching_columns) {
      candidate <- stringr::str_squish(as.character(df[[column]]))
      valid <- !is.na(candidate) & nzchar(candidate)
      replace <- (is.na(sample_id) | !nzchar(sample_id)) & valid
      sample_id[replace] <- candidate[replace]
    }
  }

  missing <- is.na(sample_id) | !nzchar(sample_id)
  sample_id[missing] <- paste0("SDRF_row_", which(missing))

  sample_id
}


ae_read_sdrf_file <- function(file, accession) {
  sdrf <- data.table::fread(
    file,
    sep = "\t",
    header = TRUE,
    fill = TRUE,
    quote = "",
    na.strings = c("", "NA", "N/A", "null"),
    data.table = FALSE,
    check.names = FALSE,
    encoding = "UTF-8"
  )

  original_names <- names(sdrf)
  unique_names <- make.unique(
    original_names,
    sep = "__duplicate_"
  )

  names(sdrf) <- unique_names

  sample_id <- ae_make_sample_id(
    sdrf,
    original_names = original_names
  )

  wide <- tibble::as_tibble(sdrf) |>
    mutate(
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

  long <- wide |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(unique_names),
      names_to = "field_unique",
      values_to = "value"
    ) |>
    mutate(
      field_raw = unname(name_map[.data$field_unique]),
      field = ae_normalize_field(.data$field_raw),
      value = stringr::str_squish(as.character(.data$value))
    ) |>
    filter(
      !is.na(.data$value),
      nzchar(.data$value)
    ) |>
    select(
      .data$Study_ID,
      .data$Sample_ID,
      .data$SDRF_Row,
      .data$SDRF_File,
      .data$field_raw,
      .data$field,
      .data$value
    )

  protocol_refs <- long |>
    filter(
      stringr::str_detect(
        .data$field_raw,
        stringr::regex("^Protocol REF", ignore_case = TRUE)
      )
    ) |>
    rename(Protocol_REF = .data$value)

  list(
    wide = wide,
    long = long,
    protocol_refs = protocol_refs
  )
}


# -------------------------------------------------------------------------
# IDF and protocol parsing
# -------------------------------------------------------------------------

ae_read_idf_file <- function(file, accession) {
  idf_raw <- data.table::fread(
    file,
    sep = "\t",
    header = FALSE,
    fill = TRUE,
    quote = "",
    blank.lines.skip = FALSE,
    na.strings = NULL,
    data.table = FALSE,
    encoding = "UTF-8"
  )

  if (ncol(idf_raw) < 2L) {
    stop("IDF file has fewer than two columns: ", file)
  }

  # Remove columns that are completely empty.
  keep_column <- vapply(
    idf_raw,
    function(x) {
      x <- stringr::str_squish(as.character(x))
      any(!is.na(x) & nzchar(x))
    },
    logical(1)
  )

  idf_raw <- idf_raw[, keep_column, drop = FALSE]
  names(idf_raw) <- paste0("V", seq_len(ncol(idf_raw)))

  idf_raw[] <- lapply(
    idf_raw,
    function(x) stringr::str_squish(as.character(x))
  )

  field_values <- idf_raw[[1]]
  field_values[
    is.na(field_values) | !nzchar(field_values)
  ] <- "Unlabeled IDF field"

  idf_long <- purrr::map_dfr(
    seq_len(nrow(idf_raw)),
    function(i) {
      values <- as.character(
        unlist(idf_raw[i, -1, drop = FALSE], use.names = FALSE)
      )

      tibble(
        Study_ID = accession,
        IDF_File = basename(file),
        IDF_Row = i,
        field_raw = field_values[i],
        field = ae_normalize_field(field_values[i]),
        value_index = seq_along(values),
        value = values
      )
    }
  ) |>
    mutate(
      value = stringr::str_squish(as.character(.data$value))
    ) |>
    filter(
      !is.na(.data$value),
      nzchar(.data$value)
    )

  protocol_rows <- which(
    stringr::str_detect(
      field_values,
      stringr::regex("^Protocol\\b", ignore_case = TRUE)
    )
  )

  protocols <- tibble()

  if (length(protocol_rows) > 0L) {
    protocol_block <- idf_raw[
      protocol_rows,
      ,
      drop = FALSE
    ]

    protocol_fields <- ae_normalize_field(
      protocol_block[[1]]
    )

    protocol_fields <- make.unique(
      protocol_fields,
      sep = "__duplicate_"
    )

    protocol_values <- protocol_block[
      ,
      -1,
      drop = FALSE
    ]

    protocols <- purrr::map_dfr(
      seq_len(ncol(protocol_values)),
      function(j) {
        values <- as.character(protocol_values[[j]])

        names(values) <- protocol_fields

        row <- as.list(values)
        row$Study_ID <- accession
        row$IDF_File <- basename(file)
        row$Protocol_Index <- j

        tibble::as_tibble_row(row)
      }
    ) |>
      relocate(
        .data$Study_ID,
        .data$IDF_File,
        .data$Protocol_Index
      )

    data_columns <- setdiff(
      names(protocols),
      c("Study_ID", "IDF_File", "Protocol_Index")
    )

    keep_protocol <- apply(
      protocols[, data_columns, drop = FALSE],
      1,
      function(x) {
        x <- stringr::str_squish(as.character(x))
        any(!is.na(x) & nzchar(x))
      }
    )

    protocols <- protocols[keep_protocol, , drop = FALSE]
  }

  idf_summary <- idf_long |>
    group_by(
      .data$Study_ID,
      .data$IDF_File,
      .data$field_raw,
      .data$field
    ) |>
    summarise(
      value = ae_collapse_unique(.data$value),
      .groups = "drop"
    )

  list(
    raw = idf_raw,
    long = idf_long,
    summary = idf_summary,
    protocols = protocols
  )
}


# -------------------------------------------------------------------------
# Study-level summary
# -------------------------------------------------------------------------

ae_extract_idf_value <- function(
    idf_long,
    patterns,
    max_values = 10L) {

  if (nrow(idf_long) == 0L) {
    return(NA_character_)
  }

  matched <- idf_long |>
    filter(
      stringr::str_detect(
        .data$field_raw,
        stringr::regex(
          paste(patterns, collapse = "|"),
          ignore_case = TRUE
        )
      )
    ) |>
    pull(.data$value)

  ae_collapse_unique(
    matched,
    max_values = max_values
  )
}


ae_make_study_summary <- function(
    accession,
    sdrf_long,
    idf_long,
    protocols,
    download_method) {

  organism_values <- sdrf_long |>
    filter(
      stringr::str_detect(
        .data$field_raw,
        stringr::regex(
          "Characteristics\\[organism\\]|organism",
          ignore_case = TRUE
        )
      )
    ) |>
    pull(.data$value)

  factor_fields <- sdrf_long |>
    filter(
      stringr::str_detect(
        .data$field_raw,
        stringr::regex(
          "^Factor Value\\[",
          ignore_case = TRUE
        )
      )
    ) |>
    distinct(.data$field_raw) |>
    pull(.data$field_raw)

  protocol_names <- character()

  if (nrow(protocols) > 0L) {
    protocol_name_columns <- names(protocols)[
      stringr::str_detect(
        names(protocols),
        "^protocol_name"
      )
    ]

    if (length(protocol_name_columns) > 0L) {
      protocol_names <- unlist(
        protocols[protocol_name_columns],
        use.names = FALSE
      )
    }
  }

  tibble(
    Study_ID = accession,
    package_version = as.character(
      utils::packageVersion("ArrayExpress")
    ),
    download_method = download_method,
    title = ae_extract_idf_value(
      idf_long,
      "^Investigation Title$"
    ),
    description = ae_extract_idf_value(
      idf_long,
      c(
        "^Experiment Description$",
        "^Experimental Design$",
        "^Investigation Description$"
      )
    ),
    experiment_type = ae_extract_idf_value(
      idf_long,
      c(
        "^Experimental Design$",
        "^Experiment Type$"
      )
    ),
    pubmed_id = ae_extract_idf_value(
      idf_long,
      "Publication PubMed ID"
    ),
    release_date = ae_extract_idf_value(
      idf_long,
      "Public Release Date|Release Date"
    ),
    organism = ae_collapse_unique(
      organism_values,
      max_values = 10L
    ),
    n_sdrf_rows = dplyr::n_distinct(
      paste(
        sdrf_long$SDRF_File,
        sdrf_long$SDRF_Row,
        sep = "::"
      )
    ),
    n_sample_ids = dplyr::n_distinct(
      sdrf_long$Sample_ID
    ),
    factor_fields = ae_collapse_unique(
      factor_fields,
      max_values = 20L
    ),
    protocol_names = ae_collapse_unique(
      protocol_names,
      max_values = 20L
    )
  )
}


# -------------------------------------------------------------------------
# One-study and batch wrappers
# -------------------------------------------------------------------------

ae_fetch_study_metadata <- function(
    accession,
    base_dir = "ArrayExpress_Metadata",
    overwrite = FALSE,
    use_api_fallback = TRUE,
    verbose = TRUE) {

  accession <- toupper(trimws(accession))

  downloaded <- ae_download_magetab(
    accession = accession,
    base_dir = base_dir,
    overwrite = overwrite,
    use_api_fallback = use_api_fallback,
    verbose = verbose
  )

  study_dir <- downloaded$path

  reported_sdrf <- ae_resolve_paths(
    downloaded$sdrf,
    study_dir
  )

  reported_idf <- ae_resolve_paths(
    downloaded$idf,
    study_dir
  )

  # Directory search protects against package-version differences in the
  # exact structure of the list returned by getAE().
  found_sdrf <- list.files(
    study_dir,
    pattern = "sdrf(?:\\.txt)?$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  found_idf <- list.files(
    study_dir,
    pattern = "idf(?:\\.txt)?$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  sdrf_files <- unique(c(reported_sdrf, found_sdrf))
  idf_files <- unique(c(reported_idf, found_idf))

  if (length(sdrf_files) == 0L) {
    stop("No SDRF file was found for ", accession, ".")
  }

  if (length(idf_files) == 0L) {
    stop("No IDF file was found for ", accession, ".")
  }

  sdrf_parsed <- purrr::map(
    sdrf_files,
    ae_read_sdrf_file,
    accession = accession
  )

  idf_parsed <- purrr::map(
    idf_files,
    ae_read_idf_file,
    accession = accession
  )

  sdrf <- purrr::map_dfr(sdrf_parsed, "wide")
  sdrf_long <- purrr::map_dfr(sdrf_parsed, "long")
  protocol_refs <- purrr::map_dfr(
    sdrf_parsed,
    "protocol_refs"
  )

  idf_long <- purrr::map_dfr(idf_parsed, "long")
  idf_summary <- purrr::map_dfr(idf_parsed, "summary")
  protocols <- purrr::map_dfr(idf_parsed, "protocols")

  study_summary <- ae_make_study_summary(
    accession = accession,
    sdrf_long = sdrf_long,
    idf_long = idf_long,
    protocols = protocols,
    download_method = downloaded$download_method
  )

  result <- list(
    Study_ID = accession,
    files = downloaded,
    study_summary = study_summary,
    sdrf = sdrf,
    sdrf_long = sdrf_long,
    idf_long = idf_long,
    idf_summary = idf_summary,
    protocols = protocols,
    protocol_refs = protocol_refs
  )

  saveRDS(
    result,
    file = file.path(
      study_dir,
      paste0(accession, "_metadata.rds")
    )
  )

  result
}


ae_fetch_many <- function(
    accessions,
    base_dir = "ArrayExpress_Metadata",
    overwrite = FALSE,
    use_api_fallback = TRUE,
    verbose = TRUE) {

  accessions <- unique(toupper(trimws(accessions)))
  accessions <- accessions[
    !is.na(accessions) & nzchar(accessions)
  ]

  results <- stats::setNames(
    vector("list", length(accessions)),
    accessions
  )

  errors <- tibble(
    Study_ID = character(),
    error = character()
  )

  for (accession in accessions) {
    if (verbose) {
      message(
        "\n========== ",
        accession,
        " =========="
      )
    }

    result <- tryCatch(
      ae_fetch_study_metadata(
        accession = accession,
        base_dir = base_dir,
        overwrite = overwrite,
        use_api_fallback = use_api_fallback,
        verbose = verbose
      ),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      errors <- bind_rows(
        errors,
        tibble(
          Study_ID = accession,
          error = conditionMessage(result)
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
    } else {
      results[[accession]] <- result
    }
  }

  attr(results, "errors") <- errors

  saveRDS(
    results,
    file = file.path(
      base_dir,
      "ArrayExpress_metadata_results.rds"
    )
  )

  results
}


ae_combine_results <- function(results) {
  valid <- results[
    !vapply(results, is.null, logical(1))
  ]

  list(
    studies = purrr::map_dfr(
      valid,
      "study_summary"
    ),
    samples = purrr::map_dfr(
      valid,
      "sdrf"
    ),
    sample_metadata_long = purrr::map_dfr(
      valid,
      "sdrf_long"
    ),
    idf_metadata_long = purrr::map_dfr(
      valid,
      "idf_long"
    ),
    idf_summary = purrr::map_dfr(
      valid,
      "idf_summary"
    ),
    protocols = purrr::map_dfr(
      valid,
      "protocols"
    ),
    protocol_refs = purrr::map_dfr(
      valid,
      "protocol_refs"
    ),
    errors = attr(results, "errors")
  )
}


ae_write_results <- function(
    combined,
    out_dir = "ArrayExpress_Metadata_Tables") {

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  output_files <- list(
    studies = file.path(
      out_dir,
      "ArrayExpress_study_summary.tsv"
    ),
    samples = file.path(
      out_dir,
      "ArrayExpress_SDRF_samples.tsv"
    ),
    sample_metadata_long = file.path(
      out_dir,
      "ArrayExpress_SDRF_long.tsv"
    ),
    idf_metadata_long = file.path(
      out_dir,
      "ArrayExpress_IDF_long.tsv"
    ),
    idf_summary = file.path(
      out_dir,
      "ArrayExpress_IDF_summary.tsv"
    ),
    protocols = file.path(
      out_dir,
      "ArrayExpress_protocols.tsv"
    ),
    protocol_refs = file.path(
      out_dir,
      "ArrayExpress_sample_protocol_refs.tsv"
    ),
    errors = file.path(
      out_dir,
      "ArrayExpress_download_errors.tsv"
    )
  )

  for (name in names(output_files)) {
    table <- combined[[name]]

    if (is.null(table)) {
      next
    }

    data.table::fwrite(
      table,
      file = output_files[[name]],
      sep = "\t",
      quote = FALSE,
      na = ""
    )
  }

  invisible(output_files)
}


# -------------------------------------------------------------------------
# Search helpers
# -------------------------------------------------------------------------

ae_search_metadata <- function(
    combined,
    terms,
    scope = c("samples", "idf", "protocols"),
    require_all = FALSE,
    ignore_case = TRUE) {

  scope <- match.arg(scope)

  terms <- terms[
    !is.na(terms) & nzchar(terms)
  ]

  if (length(terms) == 0L) {
    stop("Provide at least one search term or regular expression.")
  }

  if (scope == "samples") {
    table <- combined$sample_metadata_long |>
      mutate(
        search_text = paste(
          .data$field_raw,
          .data$value,
          sep = ": "
        )
      )
  } else if (scope == "idf") {
    table <- combined$idf_metadata_long |>
      mutate(
        search_text = paste(
          .data$field_raw,
          .data$value,
          sep = ": "
        )
      )
  } else {
    table <- combined$protocols |>
      tidyr::pivot_longer(
        cols = -c(
          .data$Study_ID,
          .data$IDF_File,
          .data$Protocol_Index
        ),
        names_to = "field_raw",
        values_to = "value"
      ) |>
      mutate(
        search_text = paste(
          .data$field_raw,
          .data$value,
          sep = ": "
        )
      )
  }

  hit_list <- lapply(
    terms,
    function(term) {
      stringr::str_detect(
        table$search_text,
        stringr::regex(
          term,
          ignore_case = ignore_case
        )
      )
    }
  )

  keep <- if (require_all) {
    Reduce(`&`, hit_list)
  } else {
    Reduce(`|`, hit_list)
  }

  table[keep, , drop = FALSE]
}


# -------------------------------------------------------------------------
# Example usage
# -------------------------------------------------------------------------

# study_ids <- c(
#   "E-MTAB-9330",
#   "E-MTAB-1234"
# )
#
# ae_results <- ae_fetch_many(
#   accessions = study_ids,
#   base_dir = "D:/R/NMD_Perturbation_Data/ArrayExpress_Metadata",
#   overwrite = FALSE,
#   use_api_fallback = TRUE,
#   verbose = TRUE
# )
#
# ae_tables <- ae_combine_results(ae_results)
#
# View(ae_tables$studies)
# View(ae_tables$samples)
# View(ae_tables$protocols)
# View(ae_tables$protocol_refs)
#
# ae_write_results(
#   ae_tables,
#   out_dir = "D:/R/NMD_Perturbation_Data/ArrayExpress_Metadata_Tables"
# )
#
# Search sample-level SDRF metadata:
# ae_search_metadata(
#   ae_tables,
#   terms = c("UPF1", "knockdown"),
#   scope = "samples"
# )
#
# Search study-level IDF descriptions:
# ae_search_metadata(
#   ae_tables,
#   terms = c("nonsense-mediated decay", "\\bNMD\\b"),
#   scope = "idf"
# )
#
# Search protocol names and descriptions:
# ae_search_metadata(
#   ae_tables,
#   terms = c("siRNA", "shRNA", "CRISPR", "inhibitor"),
#   scope = "protocols"
# )

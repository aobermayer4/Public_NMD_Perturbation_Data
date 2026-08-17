# ArrayExpress one-row-per-sample patch
#
# Source the main parser first, then this file:
#
#   source("arrayexpress_metadata_parser.R")
#   source("arrayexpress_one_row_samples_patch.R")
#
# The function can be run directly on the existing combined samples table:
#
#   ae_samples_one_row <- ae_collapse_sdrf_samples(
#     ae_tables$samples
#   )
#
# To integrate it permanently into the parser, use the edits documented
# at the bottom of this file.

# -------------------------------------------------------------------------
# Collapse SDRF file-level rows to one row per sequenced sample
# -------------------------------------------------------------------------
#
# ArrayExpress SDRF files commonly contain one row per FASTQ file rather
# than one row per biological/sequenced sample. Paired-end samples may
# therefore have separate R1 and R2 rows, and samples sequenced across
# multiple lanes or runs may have more than two rows.
#
# This function:
#   1. Removes exact duplicate SDRF records.
#   2. Detects R1/R2 using submitted filenames first, then scan names,
#      and finally FASTQ URIs.
#   3. Collapses all metadata to one row per Study_ID + Sample_ID.
#   4. Creates separate R1 and R2 columns for filenames and FASTQ URIs.
#   5. Preserves multiple lanes/runs as semicolon-separated values.
#
# By default, one output row represents one Study_ID + Sample_ID.
# To retain separate assays or libraries, supply additional group_cols,
# for example:
#   group_cols = c("Study_ID", "Sample_ID", "Assay Name")
# -------------------------------------------------------------------------

ae_collapse_sdrf_samples <- function(
    sdrf,
    group_cols = c("Study_ID", "Sample_ID"),
    read_column_map = c(
      "Scan Name" = "Scan_Name",
      "Comment[SUBMITTED_FILE_NAME]" = "Submitted_File",
      "Comment[FASTQ_URI]" = "FASTQ_URI"
    ),
    read_label_cols = c(
      "Comment[SUBMITTED_FILE_NAME]",
      "Scan Name",
      "Comment[FASTQ_URI]"
    ),
    layout_col = "Comment[LIBRARY_LAYOUT]",
    deduplicate = TRUE,
    collapse_sep = "; ",
    verbose = TRUE) {

  if (!is.data.frame(sdrf)) {
    stop("sdrf must be a data.frame or tibble.")
  }

  missing_group_cols <- setdiff(group_cols, names(sdrf))

  if (length(missing_group_cols) > 0L) {
    stop(
      "The following grouping columns are missing from sdrf: ",
      paste(missing_group_cols, collapse = ", ")
    )
  }

  # SDRF values are metadata. Treat all fields as character even when a
  # column contains numeric-looking replicate or time values.
  x_raw <- as.data.frame(sdrf, stringsAsFactors = FALSE)
  x_raw[] <- lapply(x_raw, as.character)

  is_nonempty <- function(x) {
    !is.na(x) & nzchar(trimws(as.character(x)))
  }

  collapse_values <- function(x) {
    x <- trimws(as.character(x))
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

    paste(x, collapse = collapse_sep)
  }

  count_unique_values <- function(x) {
    x <- trimws(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    length(unique(x))
  }

  # Count the rows before duplicate removal so the output records how much
  # the original SDRF was expanded.
  raw_counts <- x_raw %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      n_sdrf_rows_raw = dplyr::n(),
      .groups = "drop"
    )

  x <- x_raw

  if (isTRUE(deduplicate)) {
    # SDRF_Row and SDRF_File describe where a record was found, not its
    # scientific identity. Ignoring them also removes repeated reads of
    # the same SDRF file from differently represented local paths.
    dedupe_cols <- setdiff(
      names(x),
      c("SDRF_Row", "SDRF_File")
    )

    x <- x %>%
      dplyr::distinct(
        dplyr::across(dplyr::all_of(dedupe_cols)),
        .keep_all = TRUE
      )
  }

  detect_mate_in_vector <- function(values) {
    values <- as.character(values)

    # These patterns intentionally require a separator around the mate
    # number so ordinary sample numbers are not mistaken for read mates.
    r1 <- stringr::str_detect(
      values,
      stringr::regex(
        "(^|[/._-])R?1(?=([._-]|\\.f(?:ast)?q|$))",
        ignore_case = TRUE
      )
    )

    r2 <- stringr::str_detect(
      values,
      stringr::regex(
        "(^|[/._-])R?2(?=([._-]|\\.f(?:ast)?q|$))",
        ignore_case = TRUE
      )
    )

    r1[is.na(r1)] <- FALSE
    r2[is.na(r2)] <- FALSE

    out <- rep(NA_character_, length(values))
    out[r1 & !r2] <- "R1"
    out[r2 & !r1] <- "R2"
    out[r1 & r2] <- "UNRESOLVED"
    out
  }

  # Detect the mate from the most informative field available. Submitted
  # names are preferred over ENA URI numbering because those labels can
  # occasionally disagree after archival processing.
  x$.Read_Mate <- rep(NA_character_, nrow(x))

  available_label_cols <- intersect(read_label_cols, names(x))

  for (column in available_label_cols) {
    mate_here <- detect_mate_in_vector(x[[column]])
    fill <- is.na(x$.Read_Mate) & !is.na(mate_here)
    x$.Read_Mate[fill] <- mate_here[fill]
  }

  if (layout_col %in% names(x)) {
    layout <- trimws(as.character(x[[layout_col]]))

    is_single <- stringr::str_detect(
      layout,
      stringr::regex(
        "^SINGLE$|single[- ]?end",
        ignore_case = TRUE
      )
    )

    x$.Read_Mate[is.na(x$.Read_Mate) & is_single] <- "SE"
  }

  x$.Read_Mate[is.na(x$.Read_Mate)] <- "UNRESOLVED"

  # Choose one file identifier per row for counts. The URI is preferred
  # because it is usually unique across runs and lanes.
  x$.Primary_File <- rep(NA_character_, nrow(x))

  primary_candidates <- intersect(
    c(
      "Comment[FASTQ_URI]",
      "Comment[SUBMITTED_FILE_NAME]",
      "Scan Name"
    ),
    names(x)
  )

  for (column in primary_candidates) {
    candidate <- trimws(as.character(x[[column]]))
    fill <- !is_nonempty(x$.Primary_File) & is_nonempty(candidate)
    x$.Primary_File[fill] <- candidate[fill]
  }

  read_cols <- intersect(
    names(read_column_map),
    names(x)
  )

  # Metadata columns remain unsuffixed. File-level columns are widened to
  # R1/R2/SE/UNRESOLVED columns below.
  metadata_cols <- setdiff(
    names(x),
    c(
      group_cols,
      read_cols,
      ".Read_Mate",
      ".Primary_File",
      "SDRF_Row"
    )
  )

  metadata_one_row <- x %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(metadata_cols),
        collapse_values
      ),
      .groups = "drop"
    )

  record_counts <- x %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      n_sdrf_records_unique = dplyr::n(),

      n_read1_records = sum(.data$.Read_Mate == "R1", na.rm = TRUE),
      n_read2_records = sum(.data$.Read_Mate == "R2", na.rm = TRUE),
      n_single_end_records = sum(.data$.Read_Mate == "SE", na.rm = TRUE),
      n_unresolved_read_records = sum(
        .data$.Read_Mate == "UNRESOLVED",
        na.rm = TRUE
      ),

      n_read1_files = count_unique_values(
        .data$.Primary_File[.data$.Read_Mate == "R1"]
      ),
      n_read2_files = count_unique_values(
        .data$.Primary_File[.data$.Read_Mate == "R2"]
      ),
      n_single_end_files = count_unique_values(
        .data$.Primary_File[.data$.Read_Mate == "SE"]
      ),
      n_unresolved_files = count_unique_values(
        .data$.Primary_File[
          .data$.Read_Mate == "UNRESOLVED"
        ]
      ),

      .groups = "drop"
    )

  # Count run, assay, and experiment accessions when those fields exist.
  count_specs <- c(
    n_ena_runs = "Comment[ENA_RUN]",
    n_assays = "Assay Name",
    n_ena_experiments = "Comment[ENA_EXPERIMENT]"
  )

  for (output_name in names(count_specs)) {
    source_column <- unname(count_specs[[output_name]])

    if (!source_column %in% names(x)) {
      next
    }

    additional_count <- x %>%
      dplyr::group_by(
        dplyr::across(dplyr::all_of(group_cols))
      ) %>%
      dplyr::summarise(
        .count_value = count_unique_values(
          .data[[source_column]]
        ),
        .groups = "drop"
      )

    names(additional_count)[
      names(additional_count) == ".count_value"
    ] <- output_name

    record_counts <- dplyr::left_join(
      record_counts,
      additional_count,
      by = group_cols
    )
  }

  record_counts <- raw_counts %>%
    dplyr::left_join(
      record_counts,
      by = group_cols
    ) %>%
    dplyr::mutate(
      n_exact_duplicate_rows_removed =
        .data$n_sdrf_rows_raw -
        .data$n_sdrf_records_unique,

      read_structure = dplyr::case_when(
        .data$n_read1_files > 0L &
          .data$n_read2_files > 0L ~
          "paired-end",

        .data$n_single_end_files > 0L &
          .data$n_read1_files == 0L &
          .data$n_read2_files == 0L ~
          "single-end",

        .data$n_unresolved_files > 0L &
          .data$n_read1_files == 0L &
          .data$n_read2_files == 0L &
          .data$n_single_end_files == 0L ~
          "unresolved",

        TRUE ~ "mixed or incomplete"
      ),

      paired_files_complete =
        .data$n_read1_files > 0L &
        .data$n_read2_files > 0L
    )

  # Create separate columns for each read mate. Multiple files for the same
  # mate, such as multiple sequencing lanes, are retained in one cell and
  # separated by collapse_sep.
  files_one_row <- NULL

  if (length(read_cols) > 0L) {
    files_long <- x %>%
      dplyr::select(
        dplyr::all_of(group_cols),
        .data$.Read_Mate,
        dplyr::all_of(read_cols)
      ) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(read_cols),
        names_to = ".source_column",
        values_to = ".file_value",
        values_transform = list(
          .file_value = as.character
        )
      ) %>%
      dplyr::mutate(
        .file_value = trimws(
          as.character(.data$.file_value)
        ),
        .output_prefix = unname(
          read_column_map[.data$.source_column]
        ),
        .output_column = paste0(
          .data$.output_prefix,
          "_",
          .data$.Read_Mate
        )
      ) %>%
      dplyr::filter(
        !is.na(.data$.file_value),
        nzchar(.data$.file_value)
      ) %>%
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(
            c(group_cols, ".output_column")
          )
        )
      ) %>%
      dplyr::summarise(
        value = collapse_values(.data$.file_value),
        .groups = "drop"
      )

    if (nrow(files_long) > 0L) {
      files_one_row <- files_long %>%
        tidyr::pivot_wider(
          names_from = ".output_column",
          values_from = "value"
        )
    }
  }

  output <- metadata_one_row %>%
    dplyr::left_join(
      record_counts,
      by = group_cols
    )

  if (!is.null(files_one_row)) {
    output <- output %>%
      dplyr::left_join(
        files_one_row,
        by = group_cols
      )
  }

  # Put identifiers and QC/count columns first.
  count_columns <- grep(
    "^n_|^read_structure$|^paired_files_complete$",
    names(output),
    value = TRUE
  )

  ordered_columns <- unique(c(
    group_cols,
    count_columns,
    setdiff(
      names(output),
      c(group_cols, count_columns)
    )
  ))

  output <- output[, ordered_columns, drop = FALSE]

  expected_rows <- nrow(
    unique(x[, group_cols, drop = FALSE])
  )

  if (nrow(output) != expected_rows) {
    stop(
      "Unexpected collapse result: expected ",
      expected_rows,
      " sample rows but produced ",
      nrow(output),
      "."
    )
  }

  if (verbose) {
    message(
      "Collapsed ",
      nrow(x_raw),
      " SDRF row(s) to ",
      nrow(output),
      " unique sequenced sample row(s)."
    )

    unresolved_samples <- sum(
      output$read_structure %in%
        c("unresolved", "mixed or incomplete"),
      na.rm = TRUE
    )

    if (unresolved_samples > 0L) {
      message(
        unresolved_samples,
        " sample(s) have unresolved or incomplete read-pair labels; ",
        "review the *_UNRESOLVED columns."
      )
    }
  }

  tibble::as_tibble(output)
}


# -------------------------------------------------------------------------
# Permanent integration edits
# -------------------------------------------------------------------------
#
# 1. In ae_fetch_study_metadata(), immediately after:
#
#      sdrf <- purrr::map_dfr(sdrf_parsed, "wide")
#
#    add:
#
#      sdrf_samples_one_row <- ae_collapse_sdrf_samples(sdrf)
#
# 2. In the result <- list(...) returned by ae_fetch_study_metadata(), add:
#
#      sdrf_samples_one_row = sdrf_samples_one_row,
#
# 3. In ae_combine_results(), immediately after the existing samples item,
#    add:
#
#      samples_one_row = purrr::map_dfr(
#        valid,
#        "sdrf_samples_one_row"
#      ),
#
# 4. In ae_write_results(), add this entry to output_files:
#
#      samples_one_row = file.path(
#        out_dir,
#        "ArrayExpress_SDRF_samples_one_row.tsv"
#      ),
#
# Then ae_write_results() will export the additional table automatically.

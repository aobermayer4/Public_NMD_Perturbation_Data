# Patch for ArrayExpress SDRF mixed-column-type pivot_longer error
#
# Usage:
#   source("arrayexpress_metadata_parser.R")
#   source("arrayexpress_sdrf_character_patch.R")
#
# This redefines ae_read_sdrf_file() so every SDRF metadata column is
# character before tidyr::pivot_longer() is called.

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
    colClasses = "character",
    encoding = "UTF-8"
  )

  # Redundant by design: protects against unexpected fread behavior and
  # guarantees that pivot_longer() combines only character columns.
  sdrf[] <- lapply(sdrf, as.character)

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

  long <- wide |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(unique_names),
      names_to = "field_unique",
      values_to = "value",
      values_transform = list(value = as.character)
    ) |>
    dplyr::mutate(
      field_raw = unname(name_map[.data$field_unique]),
      field = ae_normalize_field(.data$field_raw),
      value = stringr::str_squish(as.character(.data$value))
    ) |>
    dplyr::filter(
      !is.na(.data$value),
      nzchar(.data$value)
    ) |>
    dplyr::select(
      .data$Study_ID,
      .data$Sample_ID,
      .data$SDRF_Row,
      .data$SDRF_File,
      .data$field_raw,
      .data$field,
      .data$value
    )

  protocol_refs <- long |>
    dplyr::filter(
      stringr::str_detect(
        .data$field_raw,
        stringr::regex("^Protocol REF", ignore_case = TRUE)
      )
    ) |>
    dplyr::rename(Protocol_REF = .data$value)

  list(
    wide = wide,
    long = long,
    protocol_refs = protocol_refs
  )
}

message(
  "Patched ae_read_sdrf_file(): SDRF columns will be coerced to character ",
  "before pivot_longer()."
)

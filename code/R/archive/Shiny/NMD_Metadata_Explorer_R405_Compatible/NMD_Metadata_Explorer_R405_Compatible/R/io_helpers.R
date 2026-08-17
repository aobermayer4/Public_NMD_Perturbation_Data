GEO_FILES <- c(
  studies = "NMD_GEO_study_summary.tsv",
  samples = "NMD_GEO_sample_annotations.tsv",
  evidence = "NMD_GEO_evidence.tsv",
  fields = "NMD_GEO_field_summary.tsv",
  metadata = "NMD_GEO_metadata_long.tsv"
)

AE_FILES <- c(
  studies = "ArrayExpress_study_summary.tsv",
  samples = "ArrayExpress_SDRF_samples.tsv",
  sdrf_long = "ArrayExpress_SDRF_long.tsv",
  idf_long = "ArrayExpress_IDF_long.tsv",
  idf_summary = "ArrayExpress_IDF_summary.tsv",
  protocols = "ArrayExpress_protocols.tsv",
  protocol_refs = "ArrayExpress_sample_protocol_refs.tsv",
  errors = "ArrayExpress_download_errors.tsv"
)

empty_dt <- function() data.table::data.table()

read_tsv_character <- function(path) {
  if (is.null(path) || length(path) == 0L || is.na(path) || !file.exists(path)) {
    return(empty_dt())
  }

  data.table::fread(
    path,
    sep = "\t",
    header = TRUE,
    fill = TRUE,
    quote = "",
    na.strings = c("", "NA", "N/A", "null", "NULL"),
    colClasses = "character",
    data.table = TRUE,
    encoding = "UTF-8",
    showProgress = FALSE
  )
}

paths_from_folder <- function(folder) {
  if (is.null(folder) || !nzchar(trimws(folder)) || !dir.exists(folder)) {
    stop("The selected folder does not exist: ", folder)
  }

  files <- list.files(folder, full.names = TRUE, recursive = FALSE)
  stats::setNames(files, basename(files))
}

paths_from_upload <- function(upload) {
  if (is.null(upload) || nrow(upload) == 0L) {
    stop("No files were uploaded.")
  }
  stats::setNames(upload$datapath, upload$name)
}

find_expected_path <- function(paths, expected_name) {
  if (length(paths) == 0L) return(NA_character_)
  hit <- which(tolower(names(paths)) == tolower(expected_name))
  if (length(hit) == 0L) return(NA_character_)
  unname(paths[hit[1]])
}

load_platform_tables <- function(paths, expected_files) {
  out <- setNames(vector("list", length(expected_files)), names(expected_files))
  manifest <- data.table::data.table(
    table_key = names(expected_files),
    expected_file = unname(expected_files),
    path = NA_character_,
    loaded = FALSE,
    rows = 0L,
    columns = 0L,
    error = NA_character_
  )

  for (i in seq_along(expected_files)) {
    key <- names(expected_files)[i]
    expected <- unname(expected_files[i])
    path <- find_expected_path(paths, expected)
    manifest$path[i] <- path

    if (is.na(path)) {
      out[[key]] <- empty_dt()
      manifest$error[i] <- "File not found"
      next
    }

    value <- tryCatch(read_tsv_character(path), error = function(e) e)
    if (inherits(value, "error")) {
      out[[key]] <- empty_dt()
      manifest$error[i] <- conditionMessage(value)
    } else {
      out[[key]] <- value
      manifest$loaded[i] <- TRUE
      manifest$rows[i] <- nrow(value)
      manifest$columns[i] <- ncol(value)
    }
  }

  list(tables = out, manifest = manifest)
}

load_metadata_bundle <- function(paths, source_label = "Selected files") {
  geo <- load_platform_tables(paths, GEO_FILES)
  ae <- load_platform_tables(paths, AE_FILES)

  manifest <- data.table::rbindlist(list(
    data.table::copy(geo$manifest)[, platform := "GEO"],
    data.table::copy(ae$manifest)[, platform := "ArrayExpress"]
  ), fill = TRUE)
  data.table::setcolorder(manifest, c("platform", setdiff(names(manifest), "platform")))

  list(
    geo = geo$tables,
    ae = ae$tables,
    manifest = manifest,
    source = source_label,
    loaded_at = Sys.time()
  )
}

validate_geo_bundle <- function(geo) {
  required <- c("studies", "samples", "metadata")
  missing <- required[vapply(required, function(x) is.null(geo[[x]]) || nrow(geo[[x]]) == 0L, logical(1))]
  if (length(missing) > 0L) {
    paste0("Missing or empty GEO tables: ", paste(missing, collapse = ", "))
  } else {
    NULL
  }
}

validate_ae_bundle <- function(ae) {
  required <- c("studies", "samples", "sdrf_long", "idf_long")
  missing <- required[vapply(required, function(x) is.null(ae[[x]]) || nrow(ae[[x]]) == 0L, logical(1))]
  if (length(missing) > 0L) {
    paste0("Missing or empty ArrayExpress tables: ", paste(missing, collapse = ", "))
  } else {
    NULL
  }
}

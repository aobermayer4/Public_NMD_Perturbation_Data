# Backend smoke test. Run from the NMD_Metadata_Explorer project directory.
required <- c("data.table", "zip")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) stop("Install required packages first: ", paste(missing, collapse = ", "))

source(file.path("R", "dictionaries.R"))
source(file.path("R", "search_helpers.R"))
source(file.path("R", "io_helpers.R"))
source(file.path("R", "module_geo.R"))
source(file.path("R", "module_arrayexpress.R"))

bundle <- load_metadata_bundle(paths_from_folder("data"), "Smoke-test data")
stopifnot(nrow(bundle$geo$studies) > 0L)
stopifnot(nrow(bundle$geo$samples) > 0L)
stopifnot(nrow(bundle$ae$studies) > 0L)
stopifnot(nrow(bundle$ae$samples) > 0L)

geo_filters <- list(
  query = "UPF1 knockdown in HEK293 cells",
  match_mode = "all",
  use_synonyms = TRUE,
  include_study_context = TRUE,
  targets = character(),
  mechanisms = character(),
  agents = character(),
  model_presets = character(),
  model_text = "",
  organisms = character(),
  sample_roles = character(),
  confidence = character(),
  review_categories = character(),
  study_ids = character()
)
geo_result <- run_geo_search(bundle$geo, build_geo_index(bundle$geo), geo_filters)
stopifnot(nrow(geo_result$tables$GEO_studies) > 0L)
stopifnot(nrow(geo_result$tables$GEO_samples) > 0L)

#ae example based on E-MTAB-9330 in the bundled data
ae_filters <- list(
  query = "SMG7 knockout",
  match_mode = "all",
  use_synonyms = TRUE,
  include_study_context = TRUE,
  targets = character(),
  mechanisms = character(),
  agents = character(),
  model_presets = character(),
  model_text = "",
  organisms = character(),
  protocol_types = character(),
  characteristic_field = "",
  characteristic_values = character(),
  study_ids = character()
)
ae_result <- run_ae_search(bundle$ae, build_ae_index(bundle$ae), ae_filters)
stopifnot(nrow(ae_result$tables$ArrayExpress_studies) > 0L)
stopifnot(nrow(ae_result$tables$ArrayExpress_samples) > 0L)

zip_file <- tempfile(fileext = ".zip")
write_result_zip(
  geo_result,
  platform = "GEO",
  filters = geo_filters,
  interpretation = geo_result$interpretation,
  output_file = zip_file
)
stopifnot(file.exists(zip_file), file.info(zip_file)$size > 0L)
unlink(zip_file)

message("All backend smoke tests passed.")

cat("NMD Metadata Explorer server preflight\n")
cat("======================================\n")
cat("R version: ", R.version.string, "\n", sep = "")

required <- c("shiny", "DT", "data.table", "zip")
versions <- vapply(
  required,
  function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) return("NOT INSTALLED")
    as.character(utils::packageVersion(pkg))
  },
  character(1)
)
print(versions)

missing <- names(versions)[versions == "NOT INSTALLED"]
if (length(missing) > 0L) {
  stop("Missing required package(s): ", paste(missing, collapse = ", "))
}

if (utils::packageVersion("shiny") < "1.5.0") {
  stop("Shiny 1.5.0 or newer is required for moduleServer().")
}

app_files <- c(
  "app.R",
  "R/dictionaries.R",
  "R/search_helpers.R",
  "R/io_helpers.R",
  "R/module_geo.R",
  "R/module_arrayexpress.R",
  "www/app.css"
)
missing_files <- app_files[!file.exists(app_files)]
if (length(missing_files) > 0L) {
  stop("Missing application file(s): ", paste(missing_files, collapse = ", "))
}

# Confirm the deployment directory only needs to be readable.
readable <- file.access(app_files, mode = 4) == 0
if (!all(readable)) {
  stop("Unreadable application file(s): ", paste(app_files[!readable], collapse = ", "))
}

# Runtime exports are created under tempdir(), not in the app directory.
test_dir <- tempfile("nmd_preflight_")
dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(test_dir)) {
  stop("R could not create a temporary directory under: ", tempdir())
}
test_file <- file.path(test_dir, "write_test.txt")
writeLines("ok", test_file)
if (!file.exists(test_file)) {
  stop("R could not write to its temporary directory: ", tempdir())
}
unlink(test_dir, recursive = TRUE, force = TRUE)

app_text <- paste(readLines("app.R", warn = FALSE), collapse = "\n")
if (grepl("bslib::|bs_theme\\(", app_text, perl = TRUE)) {
  stop("This app.R still contains a bslib theme call. Use the compatibility build.")
}

cat("\nPreflight passed.\n")
cat("The app does not require bslib or Sass compilation.\n")
cat("Temporary directory: ", tempdir(), "\n", sep = "")

# NMD Metadata Explorer

A local R Shiny application for interactively searching the parsed GEO and ArrayExpress metadata tables produced by:

- `write_nmd_geo_results()`
- `ae_write_results()`

The two repositories are intentionally kept separate because their metadata models differ. GEO results are organized around GSE/GSM annotations and evidence records. ArrayExpress results are organized around SDRF sample rows, IDF descriptions, protocol definitions, and sample-to-protocol references.

## Features

- Loads a folder containing the expected TSV files or accepts multiple uploaded files.
- Includes the supplied GEO and ArrayExpress tables in `data/`, so the app runs immediately as an example.
- Plain-language free-text search with optional synonym expansion.
- Preset filters for NMD factors, perturbation mechanisms, treatments, common cell models, organisms, study IDs, and platform-specific metadata.
- Results update only after **Run search** is clicked.
- Separate tabbed results for study, sample, evidence, protocol, and long-format metadata tables.
- Interactive DT tables with column filters, sorting, copy, CSV, and column visibility controls.
- ZIP download containing all filtered tables plus a record of the search interpretation and selected filters.

## Required packages

```r
install.packages(c(
  "shiny",
  "bslib",
  "DT",
  "data.table",
  "zip"
))
```

Or run:

```r
source("install_packages.R")
```

## Launch the app

Set the R working directory to the project folder and run:

```r
shiny::runApp()
```

Alternatively:

```r
source("run_app.R")
```

In RStudio or Positron, opening `app.R` and selecting **Run App** should also work.

## Loading updated output tables

The easiest local workflow is:

1. Run your GEO and ArrayExpress parsing scripts.
2. Put their output TSV files in one folder.
3. Open the app's **Load data** tab.
4. Enter the folder path and click **Load folder**.

The recognized filenames are listed in the app's Help tab. Missing optional tables are shown in the loaded-table inventory.

## How smart search works

For a query such as:

```text
UPF1 knockdown in HEK293 cells
```

The app recognizes concept groups and expands them. For example:

- `UPF1` also searches `UPF-1`, `RENT1`, and `SMG2`.
- `knockdown` also searches RNAi, siRNA, shRNA, depletion, and related terms.
- `HEK293` also searches common spellings such as `HEK-293`, `293T`, and `Flp-In T-REx 293`.

With **Match all concepts**, every concept group must be present somewhere in the searchable sample, protocol, or study text. The exact interpretation appears in the **Search interpretation** results tab and in downloaded ZIP files.

## Important interpretation note

This app is a metadata screening and review tool. A keyword or synonym match is not proof that a treatment directly inhibited NMD. Translation inhibitors, readthrough compounds, selection antibiotics, and genetic manipulations can appear in metadata for different experimental reasons. The detailed evidence, protocol, SDRF, and IDF tabs are retained so promising matches can be manually verified.

## Project structure

```text
NMD_Metadata_Explorer/
├── app.R
├── install_packages.R
├── run_app.R
├── README.md
├── R/
│   ├── dictionaries.R
│   ├── io_helpers.R
│   ├── search_helpers.R
│   ├── module_geo.R
│   └── module_arrayexpress.R
├── data/
│   └── parsed metadata TSV files
└── www/
    └── app.css
```

required_packages <- c("shiny", "DT", "data.table", "zip")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "),
    "\nRun install_packages.R, then restart R and launch the app again."
  )
}

options(
  shiny.maxRequestSize = 250 * 1024^2,
  DT.options = list(pageLength = 15)
)

source(file.path("R", "dictionaries.R"), local = TRUE)
source(file.path("R", "search_helpers.R"), local = TRUE)
source(file.path("R", "io_helpers.R"), local = TRUE)
source(file.path("R", "module_geo.R"), local = TRUE)
source(file.path("R", "module_arrayexpress.R"), local = TRUE)

default_data_dir <- normalizePath(file.path(getwd(), "data"), winslash = "/", mustWork = FALSE)

# This compatibility build intentionally uses Shiny's bundled Bootstrap 3
# and static CSS. It does not invoke bslib or Sass at runtime.

ui <- shiny::tagList(
  shiny::tags$head(
    shiny::includeCSS(file.path("www", "app.css")),
    shiny::tags$title("NMD Metadata Explorer")
  ),
  shiny::navbarPage(
    title = "NMD Metadata Explorer",
    id = "main_navigation",
    collapsible = TRUE,

    shiny::tabPanel(
      "Load data",
      shiny::fluidPage(
        shiny::div(
          class = "page-intro",
          shiny::h2("Load parsed GEO and ArrayExpress metadata"),
          shiny::p(
            "The app recognizes the table names created by write_nmd_geo_results() and ae_write_results(). ",
            "You can load a local folder or upload the TSV files directly."
          )
        ),
        shiny::fluidRow(
          shiny::column(
            width = 6,
            shiny::wellPanel(
              shiny::h4("Load from a local folder"),
              shiny::p(
                "Best when running Shiny on your own computer. The folder may contain both GEO and ArrayExpress outputs."
              ),
              shiny::textInput(
                "data_folder",
                "Folder path",
                value = default_data_dir,
                width = "100%"
              ),
              shiny::actionButton("load_folder", "Load folder", class = "btn-primary"),
              shiny::actionButton("load_bundled", "Reload bundled example data")
            )
          ),
          shiny::column(
            width = 6,
            shiny::wellPanel(
              shiny::h4("Upload TSV files"),
              shiny::p(
                "Select any or all recognized output tables. Existing loaded data are replaced when you click Load uploaded files."
              ),
              shiny::fileInput(
                "uploaded_files",
                "Choose TSV files",
                multiple = TRUE,
                accept = c(".tsv", ".txt", "text/tab-separated-values", "text/plain")
              ),
              shiny::actionButton("load_upload", "Load uploaded files", class = "btn-primary")
            )
          )
        ),
        shiny::uiOutput("load_status"),
        shiny::h3("Loaded table inventory"),
        DT::DTOutput("manifest_table")
      )
    ),

    shiny::tabPanel("GEO Explorer", geo_explorer_ui("geo")),
    shiny::tabPanel("ArrayExpress Explorer", ae_explorer_ui("ae")),

    shiny::tabPanel(
      "Help",
      shiny::fluidPage(
        shiny::div(
          class = "page-intro",
          shiny::h2("How the search works"),
          shiny::p(
            "GEO and ArrayExpress remain separate because their metadata structures are not equivalent. ",
            "The app searches each platform through a platform-specific sample index and then returns related tables for the matching studies and samples."
          )
        ),
        shiny::h3("Suggested searches"),
        shiny::tags$ul(
          shiny::tags$li(shiny::strong("UPF1 knockdown in HEK293 cells"), " using Match all concepts."),
          shiny::tags$li(shiny::strong("SMG7 knockout with SMG6 knockdown"), " for compound genetic perturbations."),
          shiny::tags$li(shiny::strong("emetine"), " or select Emetine under Drug or treatment."),
          shiny::tags$li(shiny::strong("control samples in HeLa"), " using a model preset and the GEO sample classification filter."),
          shiny::tags$li(shiny::strong("nasopharyngeal carcinoma"), " in the custom biological model field or free-text query.")
        ),
        shiny::h3("Free-text modes"),
        shiny::tags$dl(
          shiny::tags$dt("Match all concepts"),
          shiny::tags$dd("Every interpreted concept must occur in the sample, protocol, or study text."),
          shiny::tags$dt("Match any concept"),
          shiny::tags$dd("At least one interpreted concept must occur."),
          shiny::tags$dt("Match the exact phrase"),
          shiny::tags$dd("The full phrase is searched literally without synonym expansion.")
        ),
        shiny::h3("Synonym expansion"),
        shiny::p(
          "The built-in dictionaries recognize common NMD-factor aliases, RNAi and CRISPR terminology, ",
          "several NMD or translation inhibitors, control terminology, and common cell-line spellings. ",
          "The Search interpretation tab shows exactly how the query was interpreted."
        ),
        shiny::h3("Downloads"),
        shiny::p(
          "Each explorer creates a ZIP file containing all currently filtered tables and a text file recording ",
          "the search interpretation, selected filters, generation time, and exported row counts."
        ),
        shiny::h3("Expected filenames"),
        shiny::fluidRow(
          shiny::column(
            6,
            shiny::h4("GEO"),
            shiny::tags$pre(paste(unname(GEO_FILES), collapse = "\n"))
          ),
          shiny::column(
            6,
            shiny::h4("ArrayExpress"),
            shiny::tags$pre(paste(unname(AE_FILES), collapse = "\n"))
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  initial_bundle <- tryCatch(
    load_metadata_bundle(paths_from_folder(default_data_dir), source_label = "Bundled example data"),
    error = function(e) NULL
  )
  data_state <- shiny::reactiveVal(initial_bundle)
  last_error <- shiny::reactiveVal(NULL)

  load_from_paths <- function(paths, label) {
    bundle <- tryCatch(
      load_metadata_bundle(paths, source_label = label),
      error = function(e) e
    )

    if (inherits(bundle, "error")) {
      last_error(conditionMessage(bundle))
      shiny::showNotification(conditionMessage(bundle), type = "error", duration = NULL)
    } else {
      data_state(bundle)
      last_error(NULL)
      shiny::showNotification(paste0("Loaded metadata from ", label), type = "message")
    }
  }

  shiny::observeEvent(input$load_folder, {
    paths <- tryCatch(paths_from_folder(input$data_folder), error = function(e) e)
    if (inherits(paths, "error")) {
      last_error(conditionMessage(paths))
      shiny::showNotification(conditionMessage(paths), type = "error")
    } else {
      load_from_paths(paths, paste0("folder: ", normalizePath(input$data_folder, winslash = "/", mustWork = FALSE)))
    }
  })

  shiny::observeEvent(input$load_bundled, {
    paths <- tryCatch(paths_from_folder(default_data_dir), error = function(e) e)
    if (inherits(paths, "error")) {
      last_error(conditionMessage(paths))
      shiny::showNotification(conditionMessage(paths), type = "error")
    } else {
      shiny::updateTextInput(session, "data_folder", value = default_data_dir)
      load_from_paths(paths, "Bundled example data")
    }
  })

  shiny::observeEvent(input$load_upload, {
    paths <- tryCatch(paths_from_upload(input$uploaded_files), error = function(e) e)
    if (inherits(paths, "error")) {
      last_error(conditionMessage(paths))
      shiny::showNotification(conditionMessage(paths), type = "error")
    } else {
      load_from_paths(paths, "uploaded files")
    }
  })

  output$load_status <- shiny::renderUI({
    bundle <- data_state()
    if (is.null(bundle)) {
      return(shiny::div(
        class = "status-panel status-error",
        shiny::strong("No data loaded."),
        shiny::p(last_error() %||% "Choose a folder or upload the expected TSV files.")
      ))
    }

    loaded <- bundle$manifest[loaded == TRUE]
    missing <- bundle$manifest[loaded == FALSE]
    shiny::div(
      class = "status-panel status-success",
      shiny::strong(paste0("Data source: ", bundle$source)),
      shiny::p(
        paste0(
          nrow(loaded), " tables loaded; ", nrow(missing), " expected tables missing or unreadable. ",
          "Loaded at ", format(bundle$loaded_at, "%Y-%m-%d %H:%M:%S"), "."
        )
      )
    )
  })

  output$manifest_table <- DT::renderDT({
    bundle <- data_state()
    shiny::req(bundle)
    result_table(bundle$manifest, page_length = 20L)
  }, server = TRUE)

  geo_explorer_server("geo", shiny::reactive({
    bundle <- data_state()
    if (is.null(bundle)) return(NULL)
    bundle$geo
  }))

  ae_explorer_server("ae", shiny::reactive({
    bundle <- data_state()
    if (is.null(bundle)) return(NULL)
    bundle$ae
  }))
}

shiny::shinyApp(ui = ui, server = server)

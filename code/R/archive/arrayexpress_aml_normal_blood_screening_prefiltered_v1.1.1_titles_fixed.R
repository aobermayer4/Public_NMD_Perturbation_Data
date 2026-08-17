# =============================================================================
# ArrayExpress / BioStudies AML and normal blood-control screening parser
#
# Purpose
# -------
# Screen a vector of ArrayExpress accessions for:
#   * Human AML relevance and/or normal/healthy hematopoietic controls
#   * Patient/donor-derived blood, PBMC, bone marrow, blast, or HSPC material
#   * RNA-seq and/or DNA methylation assays
#   * Raw-file availability (FASTQ/SRA/ENA runs, IDAT, BAM/CRAM)
#   * Approximate sample count
#   * Explicit patient/donor/subject count when such identifiers are present
#
# This is a screening tool, not a final biological annotation system.
# It retains the metadata evidence supporting every keyword-based flag.
#
# Main functions
# --------------
#   ae_aml_quick_prescreen()
#   ae_screen_aml_accessions()
#   ae_write_aml_screening_results()
#
# Two-stage behavior
# ------------------
# By default, ae_screen_aml_accessions() first queries the small BioStudies
# JSON records. It only downloads and parses MAGE-TAB metadata when the
# pre-screen supports AML and/or normal hematopoietic controls plus raw
# RNA-seq and/or raw methylation data.
#
# Required packages
# -----------------
#   ArrayExpress, data.table, dplyr, tidyr, purrr,
#   stringr, tibble, jsonlite
#
# R compatibility
# ---------------
# Written without the native |> pipe, so it can run under R 4.0.5.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------

required_packages <- c(
  "ArrayExpress",
  "data.table",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "tibble",
  "jsonlite"
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
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

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


# -----------------------------------------------------------------------------
# 1. Dictionaries used for loose screening
# -----------------------------------------------------------------------------

AE_AML_DISEASE_TERMS <- c(
  "Acute myeloid leukemia" =
    "acute[ -]+myeloid[ -]+leuk(?:e|ae)mia",
  "Acute myelogenous leukemia" =
    "acute[ -]+myelogenous[ -]+leuk(?:e|ae)mia",
  "Acute myeloblastic leukemia" =
    "acute[ -]+myeloblastic[ -]+leuk(?:e|ae)mia",
  "Acute myelocytic leukemia" =
    "acute[ -]+myelocytic[ -]+leuk(?:e|ae)mia",
  "AML" =
    "(^|[^A-Za-z0-9])AML([^A-Za-z0-9]|$)"
)

# Explicit normal/healthy-control phrases. "Normal" alone is intentionally
# excluded because phrases such as "normal-karyotype AML" are not controls.
AE_AML_NORMAL_CONTROL_TERMS <- c(
  "Healthy donor" =
    "healthy[ -]+donors?",
  "Healthy control" =
    "healthy[ -]+controls?",
  "Normal donor" =
    "normal[ -]+donors?",
  "Normal control" =
    "normal[ -]+controls?",
  "Control donor" =
    "control[ -]+donors?",
  "Unaffected donor or control" =
    "unaffected[ -]+(donors?|controls?|individuals?|subjects?)",
  "Non-leukemic control" =
    "non[ -]?leuk(?:e|ae)mic[ -]+(controls?|donors?|samples?)",
  "Healthy volunteer" =
    "healthy[ -]+volunteers?",
  "Normal blood" =
    "normal[ -]+(peripheral[ -]+)?blood",
  "Normal bone marrow" =
    "normal[ -]+bone[ -]+marrow",
  "Normal PBMC" =
    "normal[ -]+PBMCs?",
  "Normal CD34-positive cells" =
    "normal[ -]+CD34(?:\\+|[ -]+positive)?[ -]+cells?"
)

# Broader hematopoietic-source terms used when deciding whether a normal-only
# study is relevant. A normal skin, liver, or brain study should not pass merely
# because it contains the phrase "healthy donor."
AE_AML_HEMATOPOIETIC_SOURCE_TERMS <- c(
  "Peripheral blood mononuclear cells" =
    "peripheral[ -]+blood[ -]+mononuclear[ -]+cells?",
  "PBMC" =
    "(^|[^A-Za-z0-9])PBMCs?([^A-Za-z0-9]|$)",
  "Peripheral blood" =
    "peripheral[ -]+blood",
  "Whole blood" =
    "whole[ -]+blood",
  "Blood" =
    "(^|[^A-Za-z])blood([^A-Za-z]|$)",
  "Bone marrow aspirate" =
    "bone[ -]+marrow[ -]+aspirate",
  "Bone marrow" =
    "bone[ -]+marrow",
  "Cord blood" =
    "cord[ -]+blood",
  "Mononuclear cells" =
    "mononuclear[ -]+cells?",
  "White blood cells" =
    "white[ -]+blood[ -]+cells?",
  "Leukocytes" =
    "leu?kocytes?",
  "CD34-positive cells" =
    "CD34(?:\\+|[ -]+positive)?[ -]+cells?",
  "Hematopoietic stem or progenitor cells" =
    "h[ae]matopoietic[ -]+stem[ -]+(?:and|or)?[ -]*progenitor[ -]+cells?",
  "HSPC" =
    "(^|[^A-Za-z0-9])HSPCs?([^A-Za-z0-9]|$)",
  "Leukemic blasts" =
    "leuk(?:e|ae)mic[ -]+blasts?",
  "Blast cells" =
    "blast[ -]+cells?",
  "Myeloblasts" =
    "myeloblasts?"
)

AE_AML_HUMAN_TERMS <- c(
  "Homo sapiens" = "homo[ _-]*sapiens",
  "Human" = "(^|[^A-Za-z])human([^A-Za-z]|$)"
)

AE_AML_NONHUMAN_TERMS <- c(
  "Mus musculus" = "mus[ _-]*musculus",
  "Mouse or murine" = "(^|[^A-Za-z])(mouse|mice|murine)([^A-Za-z]|$)",
  "Rattus norvegicus" = "rattus[ _-]*norvegicus",
  "Rat" = "(^|[^A-Za-z])rats?([^A-Za-z]|$)"
)

AE_AML_PATIENT_SOURCE_TERMS <- c(
  "Peripheral blood mononuclear cells" =
    "peripheral[ -]+blood[ -]+mononuclear[ -]+cells?",
  "PBMC" =
    "(^|[^A-Za-z0-9])PBMCs?([^A-Za-z0-9]|$)",
  "Peripheral blood" =
    "peripheral[ -]+blood",
  "Whole blood" =
    "whole[ -]+blood",
  "Blood" =
    "(^|[^A-Za-z])blood([^A-Za-z]|$)",
  "Bone marrow aspirate" =
    "bone[ -]+marrow[ -]+aspirate",
  "Bone marrow" =
    "bone[ -]+marrow",
  "Leukemic blasts" =
    "leuk(?:e|ae)mic[ -]+blasts?",
  "Blast cells" =
    "blast[ -]+cells?",
  "Myeloblasts" =
    "myeloblasts?",
  "Mononuclear cells" =
    "mononuclear[ -]+cells?",
  "Cord blood" =
    "cord[ -]+blood",
  "White blood cells" =
    "white[ -]+blood[ -]+cells?",
  "Leukocytes" =
    "leu?kocytes?",
  "CD34-positive cells" =
    "CD34(?:\\+|[ -]+positive)?[ -]+cells?",
  "Hematopoietic stem or progenitor cells" =
    "h[ae]matopoietic[ -]+stem[ -]+(?:and|or)?[ -]*progenitor[ -]+cells?",
  "HSPC" =
    "(^|[^A-Za-z0-9])HSPCs?([^A-Za-z0-9]|$)"
)

AE_AML_CELL_LINE_TERMS <- c(
  "Cell line" = "cell[ -]+line",
  "Xenograft" = "xenograft",
  "THP-1" = "THP[ -]?1",
  "HL-60" = "HL[ -]?60",
  "MOLM-13" = "MOLM[ -]?13",
  "MV4-11" = "MV4[ -]?11",
  "OCI-AML" = "OCI[ -]?AML",
  "KG-1" = "KG[ -]?1",
  "K562" = "(^|[^A-Za-z0-9])K562([^A-Za-z0-9]|$)",
  "NB4" = "(^|[^A-Za-z0-9])NB4([^A-Za-z0-9]|$)",
  "Kasumi" = "Kasumi"
)

AE_AML_RNASEQ_TERMS <- c(
  "RNA-seq" = "RNA[ -]?seq",
  "Transcriptome sequencing" =
    "transcriptome[ -]+sequenc",
  "Expression profiling by high-throughput sequencing" =
    "expression[ -]+profiling[ -]+by[ -]+high[ -]+throughput[ -]+sequenc",
  "mRNA sequencing" =
    "mRNA[ -]+sequenc",
  "Total RNA sequencing" =
    "total[ -]+RNA[ -]+sequenc"
)

AE_AML_METHYLATION_TERMS <- c(
  "DNA methylation" =
    "DNA[ -]+methylation|methylation[ -]+profil",
  "Methylome" =
    "methylome",
  "Bisulfite sequencing" =
    "bisulfite[ -]+sequenc|bisulphite[ -]+sequenc",
  "Whole-genome bisulfite sequencing" =
    "(^|[^A-Za-z0-9])WGBS([^A-Za-z0-9]|$)|whole[ -]+genome[ -]+bisulfite",
  "Reduced representation bisulfite sequencing" =
    "(^|[^A-Za-z0-9])RRBS([^A-Za-z0-9]|$)|reduced[ -]+representation[ -]+bisulfite",
  "Illumina methylation array" =
    "HumanMethylation|MethylationEPIC|Infinium|450K|850K|EPIC[ -]+array",
  "Methyl-seq" =
    "methyl[ -]?seq"
)

AE_AML_RAW_SEQUENCE_TERMS <- c(
  "FASTQ" = "\\.(fastq|fq)(\\.gz)?([?#].*)?$",
  "SRA file" = "\\.sra([?#].*)?$",
  "ENA/SRA run accession" =
    "(^|[^A-Za-z0-9])(ERR|SRR|DRR)[0-9]+([^A-Za-z0-9]|$)",
  "FASTQ URI field" =
    "FASTQ_URI|FASTQ[ -]+URI"
)

AE_AML_IDAT_TERMS <- c(
  "IDAT" = "\\.idat(\\.gz)?([?#].*)?$|(^|[^A-Za-z0-9])IDAT([^A-Za-z0-9]|$)"
)

# Additional raw methylation-array formats occur in older imported studies.
# These are only used for the quick gate and are retained as review evidence.
AE_AML_RAW_METHYLATION_ARRAY_TERMS <- c(
  "IDAT" =
    "\\.idat(\\.gz)?([?#].*)?$|(^|[^A-Za-z0-9])IDAT([^A-Za-z0-9]|$)",
  "Affymetrix CEL" =
    "\\.cel(\\.gz)?([?#].*)?$",
  "NimbleGen pair file" =
    "\\.pair(\\.gz)?([?#].*)?$",
  "GenePix result" =
    "\\.gpr(\\.gz)?([?#].*)?$",
  "Illumina red/green channel file" =
    "(_Red|_Grn)\\.(idat|txt)(\\.gz)?([?#].*)?$"
)

AE_AML_GENERIC_RAW_TERMS <- c(
  "Raw data wording" =
    "(^|[^A-Za-z])raw[ -]+(data|file|files|read|reads)([^A-Za-z]|$)",
  "ENA" =
    "(^|[^A-Za-z0-9])ENA([^A-Za-z0-9]|$)|European[ -]+Nucleotide[ -]+Archive",
  "SRA" =
    "(^|[^A-Za-z0-9])SRA([^A-Za-z0-9]|$)|Sequence[ -]+Read[ -]+Archive",
  "BioProject accession" =
    "(^|[^A-Za-z0-9])PRJ(E|N|D)[A-Z]?[0-9]+([^A-Za-z0-9]|$)",
  "Sequencing project accession" =
    "(^|[^A-Za-z0-9])(ERP|SRP|DRP)[0-9]+([^A-Za-z0-9]|$)"
)

AE_AML_METHYLATION_ARRAY_TERMS <- c(
  "Illumina methylation array" =
    "HumanMethylation|MethylationEPIC|Infinium|450K|850K|EPIC[ -]+array",
  "Methylation profiling by array" =
    "methylation[ -]+profiling[ -]+by[ -]+array"
)

AE_AML_METHYLATION_SEQUENCING_TERMS <- c(
  "Bisulfite sequencing" =
    "bisulfite[ -]+sequenc|bisulphite[ -]+sequenc",
  "WGBS" =
    "(^|[^A-Za-z0-9])WGBS([^A-Za-z0-9]|$)|whole[ -]+genome[ -]+bisulfite",
  "RRBS or ERRBS" =
    "(^|[^A-Za-z0-9])E?RRBS([^A-Za-z0-9]|$)|reduced[ -]+representation[ -]+bisulfite",
  "MeDIP-seq" =
    "MeDIP[ -]?seq",
  "Methyl-seq" =
    "methyl[ -]?seq",
  "Methylation high-throughput sequencing" =
    "methylation[ -]+profiling[ -]+by[ -]+high[ -]+throughput[ -]+sequenc"
)

AE_AML_ALIGNED_TERMS <- c(
  "BAM" = "\\.bam([?#].*)?$",
  "CRAM" = "\\.cram([?#].*)?$"
)

AE_AML_PROCESSED_FILE_TERMS <- c(
  "Expression or count matrix" =
    "count[s]?[ _-]*matrix|expression[ _-]*matrix|normalized[ _-]*matrix",
  "Methylation beta values" =
    "beta[ _-]*value|methylation[ _-]*matrix",
  "Tabular file" =
    "\\.(csv|tsv|txt)(\\.gz)?([?#].*)?$",
  "Spreadsheet" =
    "\\.(xls|xlsx)([?#].*)?$",
  "Genomic interval or signal" =
    "\\.(bed|bedgraph|bigwig|bw)(\\.gz)?([?#].*)?$"
)


# -----------------------------------------------------------------------------
# 2. General text and matching helpers
# -----------------------------------------------------------------------------

ae_aml_safe_utf8 <- function(x) {
  x <- as.character(x)
  original_na <- is.na(x)
  out <- rep(NA_character_, length(x))

  keep <- !original_na

  if (any(keep)) {
    out[keep] <- suppressWarnings(
      iconv(
        x[keep],
        from = "",
        to = "UTF-8",
        sub = NA
      )
    )

    bad <- keep & is.na(out)

    if (any(bad)) {
      out[bad] <- suppressWarnings(
        iconv(
          x[bad],
          from = "CP1252",
          to = "UTF-8",
          sub = NA
        )
      )
    }

    bad <- keep & is.na(out)

    if (any(bad)) {
      out[bad] <- suppressWarnings(
        iconv(
          x[bad],
          from = "latin1",
          to = "UTF-8",
          sub = NA
        )
      )
    }

    bad <- keep & is.na(out)

    if (any(bad)) {
      out[bad] <- suppressWarnings(
        iconv(
          x[bad],
          from = "",
          to = "UTF-8",
          sub = "?"
        )
      )
    }
  }

  out[original_na] <- NA_character_
  out
}


ae_aml_clean_text <- function(x) {
  x <- ae_aml_safe_utf8(x)
  x <- stringr::str_squish(x)

  x[
    !is.na(x) &
      tolower(x) %in% c(
        "", "na", "n/a", "null", "none",
        "not applicable", "unknown"
      )
  ] <- NA_character_

  x
}


ae_aml_normalize_field <- function(x) {
  x <- ae_aml_clean_text(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}


ae_aml_collapse_unique <- function(
    x,
    sep = "; ",
    max_values = Inf) {

  x <- ae_aml_clean_text(x)
  x <- unique(x[!is.na(x) & nzchar(x)])

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


ae_aml_first_nonempty <- function(x) {
  x <- ae_aml_clean_text(x)
  x <- x[!is.na(x) & nzchar(x)]

  if (length(x) == 0L) {
    return(NA_character_)
  }

  x[1]
}


ae_aml_dictionary_hits <- function(text, dictionary) {
  text <- ae_aml_clean_text(text)
  text <- paste(text[!is.na(text)], collapse = " ; ")

  if (!nzchar(text)) {
    return(character())
  }

  hit <- vapply(
    dictionary,
    function(pattern) {
      stringr::str_detect(
        text,
        stringr::regex(pattern, ignore_case = TRUE)
      )
    },
    logical(1)
  )

  names(dictionary)[hit]
}


ae_aml_dictionary_hit_text <- function(text, dictionary) {
  hits <- ae_aml_dictionary_hits(text, dictionary)

  if (length(hits) == 0L) {
    return(NA_character_)
  }

  paste(hits, collapse = "; ")
}


ae_aml_has_dictionary_hit <- function(text, dictionary) {
  length(ae_aml_dictionary_hits(text, dictionary)) > 0L
}


ae_aml_retry <- function(
    expression_function,
    n = 3L,
    wait_seconds = 4,
    verbose = TRUE) {

  last_error <- NULL

  for (attempt in seq_len(n)) {
    answer <- tryCatch(
      expression_function(),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )

    if (!is.null(answer)) {
      return(answer)
    }

    if (verbose && !is.null(last_error)) {
      message(
        "Attempt ",
        attempt,
        "/",
        n,
        " failed: ",
        conditionMessage(last_error)
      )
    }

    if (attempt < n) {
      Sys.sleep(wait_seconds)
    }
  }

  if (!is.null(last_error)) {
    stop(last_error)
  }

  stop("Operation failed without an error object.")
}


# -----------------------------------------------------------------------------
# 3. MAGE-TAB metadata download
# -----------------------------------------------------------------------------

ae_aml_download_file <- function(
    url,
    destination,
    overwrite = FALSE) {

  if (file.exists(destination) && !overwrite) {
    return(
      normalizePath(
        destination,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  dir.create(
    dirname(destination),
    recursive = TRUE,
    showWarnings = FALSE
  )

  ae_aml_retry(
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
    }
  )
}


ae_aml_download_with_arrayexpress <- function(
    accession,
    study_dir,
    overwrite = FALSE,
    verbose = TRUE) {

  getae_formals <- names(formals(ArrayExpress::getAE))

  args <- list(
    accession,
    path = study_dir,
    type = "mage",
    extract = TRUE
  )

  # The first argument has been called "input" in some older releases and
  # "accession" in newer releases. Passing it positionally avoids that drama.
  if ("overwrite" %in% getae_formals) {
    args$overwrite <- overwrite
  }

  if (verbose) {
    message(
      "Trying ArrayExpress::getAE(type = 'mage') for ",
      accession
    )
  }

  result <- do.call(
    ArrayExpress::getAE,
    args
  )

  if (is.null(result)) {
    stop("ArrayExpress::getAE() returned NULL.")
  }

  result
}


ae_aml_biostudies_json <- function(accession) {
  study_url <- paste0(
    "https://www.ebi.ac.uk/biostudies/api/v1/studies/",
    accession
  )

  list(
    study = ae_aml_retry(
      function() {
        jsonlite::fromJSON(
          study_url,
          simplifyVector = TRUE
        )
      }
    ),
    info = ae_aml_retry(
      function() {
        jsonlite::fromJSON(
          paste0(study_url, "/info"),
          simplifyVector = TRUE
        )
      }
    )
  )
}


# Convert the PageTab study JSON and /info JSON into a small field/value table.
# This is much cheaper than downloading and expanding the MAGE-TAB archive.
ae_aml_flatten_biostudies_json <- function(
    biostudies_json,
    accession) {

  flatten_one <- function(x, source_name) {
    flat <- unlist(
      x,
      recursive = TRUE,
      use.names = TRUE
    )

    if (length(flat) == 0L) {
      return(
        tibble::tibble(
          Study_ID = character(),
          json_source = character(),
          field = character(),
          value = character()
        )
      )
    }

    field_names <- names(flat)

    if (is.null(field_names)) {
      field_names <- paste0(
        source_name,
        "_value_",
        seq_along(flat)
      )
    }

    tibble::tibble(
      Study_ID = accession,
      json_source = source_name,
      field = ae_aml_safe_utf8(field_names),
      value = ae_aml_clean_text(
        as.character(flat)
      )
    ) %>%
      dplyr::filter(
        !is.na(.data$value),
        nzchar(.data$value)
      )
  }

  dplyr::bind_rows(
    flatten_one(
      biostudies_json$study,
      "study"
    ),
    flatten_one(
      biostudies_json$info,
      "info"
    )
  ) %>%
    dplyr::distinct()
}


ae_aml_extract_prescreen_title <- function(flat_json) {
  if (nrow(flat_json) == 0L) {
    return(NA_character_)
  }

  title_rows <- stringr::str_detect(
    flat_json$field,
    stringr::regex(
      "(^|\\.)title($|\\.)|section.*title|attributes.*title",
      ignore_case = TRUE
    )
  )

  title <- ae_aml_first_nonempty(
    flat_json$value[title_rows]
  )

  if (is.na(title)) {
    title <- ae_aml_first_nonempty(
      flat_json$value
    )
  }

  title
}


# Fast gate used before any MAGE-TAB download.
#
# The gate intentionally does not require PBMC/blood/bone-marrow terms because
# these are often only present in the SDRF. It can, however, reject studies that
# are clearly non-human when no human evidence is present.
ae_aml_quick_prescreen <- function(
    accession,
    biostudies_json = NULL,
    require_raw_evidence = TRUE,
    require_human_evidence = FALSE,
    exclude_clear_nonhuman = TRUE,
    include_normal_controls = TRUE,
    verbose = TRUE) {

  accession <- toupper(trimws(accession))

  if (is.null(biostudies_json)) {
    if (verbose) {
      message(
        "Quick BioStudies pre-screen: ",
        accession
      )
    }

    biostudies_json <- ae_aml_biostudies_json(
      accession
    )
  }

  flat_json <- ae_aml_flatten_biostudies_json(
    biostudies_json = biostudies_json,
    accession = accession
  )

  # Exclude obvious citation/reference fields from disease screening. This
  # reduces false positives where AML appears only in a cited paper title.
  screening_rows <- !stringr::str_detect(
    flat_json$field,
    stringr::regex(
      "reference|publication|pubmed|citation|bibliograph",
      ignore_case = TRUE
    )
  )

  screening_text <- ae_aml_collapse_unique(
    paste(
      flat_json$field[screening_rows],
      flat_json$value[screening_rows],
      sep = ": "
    ),
    max_values = 2000L
  )

  all_text <- ae_aml_collapse_unique(
    paste(
      flat_json$field,
      flat_json$value,
      sep = ": "
    ),
    max_values = 2500L
  )

  aml_terms <- ae_aml_dictionary_hit_text(
    screening_text,
    AE_AML_DISEASE_TERMS
  )

  normal_control_terms <-
    ae_aml_dictionary_hit_text(
      screening_text,
      AE_AML_NORMAL_CONTROL_TERMS
    )

  hematopoietic_source_terms <-
    ae_aml_dictionary_hit_text(
      screening_text,
      AE_AML_HEMATOPOIETIC_SOURCE_TERMS
    )

  human_terms <- ae_aml_dictionary_hit_text(
    screening_text,
    AE_AML_HUMAN_TERMS
  )

  nonhuman_terms <- ae_aml_dictionary_hit_text(
    screening_text,
    AE_AML_NONHUMAN_TERMS
  )

  rnaseq_terms <- ae_aml_dictionary_hit_text(
    screening_text,
    AE_AML_RNASEQ_TERMS
  )

  methylation_terms <- ae_aml_dictionary_hit_text(
    screening_text,
    AE_AML_METHYLATION_TERMS
  )

  methylation_array_terms <-
    ae_aml_dictionary_hit_text(
      screening_text,
      AE_AML_METHYLATION_ARRAY_TERMS
    )

  methylation_sequencing_terms <-
    ae_aml_dictionary_hit_text(
      screening_text,
      AE_AML_METHYLATION_SEQUENCING_TERMS
    )

  # Raw files and archive accessions can live anywhere in the JSON, including
  # file sections that are excluded from disease-context screening.
  raw_sequence_terms <-
    ae_aml_dictionary_hit_text(
      all_text,
      AE_AML_RAW_SEQUENCE_TERMS
    )

  raw_methylation_array_terms <-
    ae_aml_dictionary_hit_text(
      all_text,
      AE_AML_RAW_METHYLATION_ARRAY_TERMS
    )

  generic_raw_terms <-
    ae_aml_dictionary_hit_text(
      all_text,
      AE_AML_GENERIC_RAW_TERMS
    )

  aml_flag <- !is.na(aml_terms)
  normal_control_term_flag <-
    !is.na(normal_control_terms)
  hematopoietic_source_flag <-
    !is.na(hematopoietic_source_terms)

  # A normal-only study must also mention a hematopoietic source. This keeps
  # healthy non-blood tissues out of the full parser.
  normal_hematopoietic_control_flag <-
    isTRUE(include_normal_controls) &&
    normal_control_term_flag &&
    hematopoietic_source_flag

  disease_or_control_flag <-
    aml_flag ||
    normal_hematopoietic_control_flag

  human_flag <- !is.na(human_terms)
  nonhuman_flag <- !is.na(nonhuman_terms)
  rnaseq_flag <- !is.na(rnaseq_terms)
  methylation_flag <- !is.na(methylation_terms)

  raw_sequence_flag <- !is.na(
    raw_sequence_terms
  )

  raw_methylation_array_flag <- !is.na(
    raw_methylation_array_terms
  )

  generic_raw_flag <- !is.na(
    generic_raw_terms
  )

  methylation_array_flag <- !is.na(
    methylation_array_terms
  )

  methylation_sequencing_flag <- !is.na(
    methylation_sequencing_terms
  )

  human_ok <- TRUE

  if (isTRUE(require_human_evidence)) {
    human_ok <- human_flag
  } else if (
    isTRUE(exclude_clear_nonhuman) &&
      nonhuman_flag &&
      !human_flag
  ) {
    human_ok <- FALSE
  }

  rnaseq_raw_candidate <-
    disease_or_control_flag &&
    human_ok &&
    rnaseq_flag &&
    (
      raw_sequence_flag ||
        generic_raw_flag
    )

  methylation_raw_candidate <-
    disease_or_control_flag &&
    human_ok &&
    methylation_flag &&
    (
      raw_methylation_array_flag ||
        generic_raw_flag ||
        (
          methylation_sequencing_flag &&
            raw_sequence_flag
        ) ||
        (
          !methylation_array_flag &&
            raw_sequence_flag
        )
    )

  relevant_assay_flag <-
    rnaseq_flag ||
    methylation_flag

  raw_relevant_flag <-
    rnaseq_raw_candidate ||
    methylation_raw_candidate

  study_scope <- dplyr::case_when(
    aml_flag &&
      normal_hematopoietic_control_flag ~
      "AML and normal hematopoietic controls",

    aml_flag ~
      "AML",

    normal_hematopoietic_control_flag ~
      "Normal hematopoietic controls",

    normal_control_term_flag &&
      !hematopoietic_source_flag ~
      "Normal controls mentioned; hematopoietic source unclear",

    TRUE ~
      "Neither AML nor normal hematopoietic controls identified"
  )

  continue_full_parse <- dplyr::case_when(
    !disease_or_control_flag ~ FALSE,
    !human_ok ~ FALSE,
    !relevant_assay_flag ~ FALSE,
    isTRUE(require_raw_evidence) ~
      raw_relevant_flag,
    TRUE ~ TRUE
  )

  prescreen_status <- dplyr::case_when(
    !disease_or_control_flag ~
      "Skip: neither AML nor normal hematopoietic controls identified",

    !human_ok ~
      "Skip: clearly non-human or human evidence required",

    !relevant_assay_flag ~
      "Skip: RNA-seq or methylation assay not identified",

    raw_relevant_flag ~
      paste0(
        "Continue: ",
        study_scope,
        " with raw RNA-seq and/or methylation evidence"
      ),

    !isTRUE(require_raw_evidence) ~
      paste0(
        "Continue: ",
        study_scope,
        "; relevant assay found but raw availability unclear"
      ),

    TRUE ~
      "Skip: raw RNA-seq or methylation evidence not identified"
  )

  summary <- tibble::tibble(
    Study_ID = accession,
    title_or_first_value =
      ae_aml_extract_prescreen_title(
        flat_json
      ),
    study_scope = study_scope,
    prescreen_status = prescreen_status,
    continue_full_parse =
      continue_full_parse,

    aml_flag = aml_flag,
    normal_control_term_flag =
      normal_control_term_flag,
    hematopoietic_source_flag =
      hematopoietic_source_flag,
    normal_hematopoietic_control_flag =
      normal_hematopoietic_control_flag,
    disease_or_control_flag =
      disease_or_control_flag,

    human_flag = human_flag,
    nonhuman_flag = nonhuman_flag,

    rnaseq_flag = rnaseq_flag,
    methylation_flag =
      methylation_flag,
    methylation_array_flag =
      methylation_array_flag,
    methylation_sequencing_flag =
      methylation_sequencing_flag,

    raw_sequence_flag =
      raw_sequence_flag,
    raw_methylation_array_flag =
      raw_methylation_array_flag,
    generic_raw_flag =
      generic_raw_flag,

    rnaseq_raw_candidate =
      rnaseq_raw_candidate,
    methylation_raw_candidate =
      methylation_raw_candidate,

    aml_terms = aml_terms,
    normal_control_terms =
      normal_control_terms,
    hematopoietic_source_terms =
      hematopoietic_source_terms,
    human_terms = human_terms,
    nonhuman_terms = nonhuman_terms,
    rnaseq_terms = rnaseq_terms,
    methylation_terms =
      methylation_terms,
    raw_sequence_terms =
      raw_sequence_terms,
    raw_methylation_array_terms =
      raw_methylation_array_terms,
    generic_raw_terms =
      generic_raw_terms
  )

  evidence <- flat_json %>%
    dplyr::mutate(
      row_text = paste(
        .data$field,
        .data$value,
        sep = ": "
      ),
      evidence_class = dplyr::case_when(
        vapply(
          .data$row_text,
          ae_aml_has_dictionary_hit,
          logical(1),
          dictionary = AE_AML_DISEASE_TERMS
        ) ~ "AML disease",

        vapply(
          .data$row_text,
          ae_aml_has_dictionary_hit,
          logical(1),
          dictionary =
            AE_AML_NORMAL_CONTROL_TERMS
        ) ~ "Normal/healthy control",

        vapply(
          .data$row_text,
          ae_aml_has_dictionary_hit,
          logical(1),
          dictionary =
            AE_AML_HEMATOPOIETIC_SOURCE_TERMS
        ) ~ "Hematopoietic source",

        vapply(
          .data$row_text,
          ae_aml_has_dictionary_hit,
          logical(1),
          dictionary = AE_AML_RNASEQ_TERMS
        ) ~ "RNA-seq assay",

        vapply(
          .data$row_text,
          ae_aml_has_dictionary_hit,
          logical(1),
          dictionary = AE_AML_METHYLATION_TERMS
        ) ~ "Methylation assay",

        vapply(
          .data$row_text,
          ae_aml_has_dictionary_hit,
          logical(1),
          dictionary = AE_AML_RAW_SEQUENCE_TERMS
        ) ~ "Raw sequencing evidence",

        vapply(
          .data$row_text,
          ae_aml_has_dictionary_hit,
          logical(1),
          dictionary =
            AE_AML_RAW_METHYLATION_ARRAY_TERMS
        ) ~ "Raw methylation-array evidence",

        vapply(
          .data$row_text,
          ae_aml_has_dictionary_hit,
          logical(1),
          dictionary = AE_AML_GENERIC_RAW_TERMS
        ) ~ "Generic raw-data evidence",

        vapply(
          .data$row_text,
          ae_aml_has_dictionary_hit,
          logical(1),
          dictionary = AE_AML_HUMAN_TERMS
        ) ~ "Human organism",

        vapply(
          .data$row_text,
          ae_aml_has_dictionary_hit,
          logical(1),
          dictionary = AE_AML_NONHUMAN_TERMS
        ) ~ "Non-human organism",

        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(
      !is.na(.data$evidence_class)
    ) %>%
    dplyr::select(
      .data$Study_ID,
      .data$json_source,
      .data$evidence_class,
      .data$field,
      .data$value
    ) %>%
    dplyr::distinct()

  list(
    summary = summary,
    evidence = evidence,
    json = biostudies_json
  )
}


ae_aml_download_from_biostudies <- function(
    accession,
    study_dir,
    overwrite = FALSE,
    verbose = TRUE,
    biostudies_json = NULL) {

  if (verbose) {
    message(
      "Using BioStudies metadata fallback for ",
      accession
    )
  }

  json <- biostudies_json

  if (is.null(json)) {
    json <- ae_aml_biostudies_json(accession)
  }

  all_strings <- unique(
    ae_aml_safe_utf8(
      as.character(
        unlist(
          c(json$study, json$info),
          recursive = TRUE,
          use.names = FALSE
        )
      )
    )
  )

  magetab_paths <- all_strings[
    !is.na(all_strings) &
      stringr::str_detect(
        basename(all_strings),
        stringr::regex(
          "(idf|sdrf)(\\.txt)?$",
          ignore_case = TRUE
        )
      )
  ]

  magetab_paths <- unique(magetab_paths)

  if (length(magetab_paths) == 0L) {
    stop(
      "No IDF or SDRF paths were found for ",
      accession,
      "."
    )
  }

  ftp_link <- NA_character_

  if (!is.null(json$info$ftpLink)) {
    ftp_link <- as.character(json$info$ftpLink)[1]
  }

  if (is.na(ftp_link) || !nzchar(ftp_link)) {
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
      relative_path <- sub(
        "^/?Files/",
        "",
        relative_path,
        ignore.case = TRUE
      )

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
        message(
          "Downloading metadata file: ",
          basename(relative_path)
        )
      }

      ae_aml_download_file(
        url = file_url,
        destination = destination,
        overwrite = overwrite
      )
    },
    character(1)
  )

  list(
    path = study_dir,
    mageTabFiles = downloaded,
    download_method = "BioStudies API fallback"
  )
}


ae_aml_download_magetab <- function(
    accession,
    base_dir = "ArrayExpress_AML_Screen",
    overwrite = FALSE,
    use_api_fallback = TRUE,
    verbose = TRUE,
    biostudies_json = NULL) {

  accession <- toupper(trimws(accession))

  if (!stringr::str_detect(
    accession,
    "^E-[A-Z0-9]+-[0-9]+$"
  )) {
    stop(
      "Invalid-looking ArrayExpress accession: ",
      accession
    )
  }

  study_dir <- file.path(base_dir, accession)

  dir.create(
    study_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  package_result <- tryCatch(
    ae_aml_download_with_arrayexpress(
      accession = accession,
      study_dir = study_dir,
      overwrite = overwrite,
      verbose = verbose
    ),
    error = function(e) e
  )

  if (!inherits(package_result, "error")) {
    package_result$path <- study_dir
    package_result$download_method <-
      "ArrayExpress::getAE(type = 'mage')"
    return(package_result)
  }

  if (!use_api_fallback) {
    stop(
      "ArrayExpress package download failed: ",
      conditionMessage(package_result)
    )
  }

  if (verbose) {
    message(
      "Package download failed: ",
      conditionMessage(package_result)
    )
  }

  ae_aml_download_from_biostudies(
    accession = accession,
    study_dir = study_dir,
    overwrite = overwrite,
    verbose = verbose,
    biostudies_json = biostudies_json
  )
}


# -----------------------------------------------------------------------------
# 4. IDF and SDRF parsing
# -----------------------------------------------------------------------------

ae_aml_make_sample_id <- function(
    data,
    original_names) {

  preferred_fields <- c(
    "Sample Name",
    "Source Name",
    "Assay Name",
    "Extract Name",
    "Scan Name",
    "Comment[ENA_SAMPLE]",
    "Comment[BIOSAMPLE]"
  )

  sample_id <- rep(
    NA_character_,
    nrow(data)
  )

  for (field in preferred_fields) {
    matching_columns <- names(data)[
      tolower(original_names) == tolower(field)
    ]

    for (column in matching_columns) {
      candidate <- ae_aml_clean_text(
        data[[column]]
      )

      replace <- (
        is.na(sample_id) |
          !nzchar(sample_id)
      ) &
        !is.na(candidate) &
        nzchar(candidate)

      sample_id[replace] <- candidate[replace]
    }
  }

  missing <- is.na(sample_id) | !nzchar(sample_id)

  sample_id[missing] <- paste0(
    "SDRF_row_",
    which(missing)
  )

  sample_id
}


ae_aml_read_sdrf <- function(
    file,
    accession) {

  sdrf <- data.table::fread(
    file,
    sep = "\t",
    header = TRUE,
    fill = TRUE,
    quote = "",
    na.strings = c(
      "", "NA", "N/A", "null", "NULL"
    ),
    data.table = FALSE,
    check.names = FALSE,
    colClasses = "character",
    encoding = "UTF-8",
    showProgress = FALSE
  )

  sdrf[] <- lapply(
    sdrf,
    ae_aml_clean_text
  )

  original_names <- names(sdrf)

  unique_names <- make.unique(
    original_names,
    sep = "__duplicate_"
  )

  names(sdrf) <- unique_names

  sample_id <- ae_aml_make_sample_id(
    data = sdrf,
    original_names = original_names
  )

  wide <- tibble::as_tibble(sdrf) %>%
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

  long <- wide %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(unique_names),
      names_to = "field_unique",
      values_to = "value",
      values_transform = list(
        value = as.character
      )
    ) %>%
    dplyr::mutate(
      field_raw = unname(
        name_map[.data$field_unique]
      ),
      field = ae_aml_normalize_field(
        .data$field_raw
      ),
      value = ae_aml_clean_text(
        .data$value
      )
    ) %>%
    dplyr::filter(
      !is.na(.data$value),
      nzchar(.data$value)
    ) %>%
    dplyr::select(
      .data$Study_ID,
      .data$Sample_ID,
      .data$SDRF_Row,
      .data$SDRF_File,
      .data$field_raw,
      .data$field,
      .data$value
    )

  list(
    wide = wide,
    long = long
  )
}


ae_aml_read_idf <- function(
    file,
    accession) {

  idf <- data.table::fread(
    file,
    sep = "\t",
    header = FALSE,
    fill = TRUE,
    quote = "",
    blank.lines.skip = FALSE,
    na.strings = NULL,
    data.table = FALSE,
    colClasses = "character",
    encoding = "UTF-8",
    showProgress = FALSE
  )

  idf[] <- lapply(
    idf,
    ae_aml_clean_text
  )

  keep_column <- vapply(
    idf,
    function(x) {
      any(!is.na(x) & nzchar(x))
    },
    logical(1)
  )

  idf <- idf[, keep_column, drop = FALSE]

  if (ncol(idf) < 2L) {
    stop(
      "IDF file has fewer than two populated columns: ",
      file
    )
  }

  names(idf) <- paste0(
    "V",
    seq_len(ncol(idf))
  )

  field_values <- idf[[1]]
  field_values[
    is.na(field_values) |
      !nzchar(field_values)
  ] <- "Unlabeled IDF field"

  idf_long <- purrr::map_dfr(
    seq_len(nrow(idf)),
    function(i) {
      values <- as.character(
        unlist(
          idf[i, -1, drop = FALSE],
          use.names = FALSE
        )
      )

      tibble::tibble(
        Study_ID = accession,
        IDF_File = basename(file),
        IDF_Row = i,
        field_raw = field_values[i],
        field = ae_aml_normalize_field(
          field_values[i]
        ),
        value_index = seq_along(values),
        value = ae_aml_clean_text(values)
      )
    }
  ) %>%
    dplyr::filter(
      !is.na(.data$value),
      nzchar(.data$value)
    )

  idf_long
}


# -----------------------------------------------------------------------------
# 5. File inventory
# -----------------------------------------------------------------------------

ae_aml_classify_file_value <- function(value) {
  value <- ae_aml_clean_text(value)

  if (is.na(value) || !nzchar(value)) {
    return("Other metadata")
  }

  if (ae_aml_has_dictionary_hit(
    value,
    AE_AML_IDAT_TERMS
  )) {
    return("Raw methylation array file")
  }

  if (ae_aml_has_dictionary_hit(
    value,
    AE_AML_RAW_SEQUENCE_TERMS
  )) {
    return("Raw sequencing file or run")
  }

  if (ae_aml_has_dictionary_hit(
    value,
    AE_AML_ALIGNED_TERMS
  )) {
    return("Aligned sequencing file")
  }

  if (ae_aml_has_dictionary_hit(
    value,
    AE_AML_PROCESSED_FILE_TERMS
  )) {
    return("Processed or derived file")
  }

  "Other file or accession"
}


ae_aml_sdrf_file_inventory <- function(sdrf_long) {
  if (nrow(sdrf_long) == 0L) {
    return(
      tibble::tibble(
        Study_ID = character(),
        Sample_ID = character(),
        inventory_source = character(),
        field_raw = character(),
        value = character(),
        file_class = character()
      )
    )
  }

  possible_file_row <- stringr::str_detect(
    sdrf_long$field,
    stringr::regex(
      "file|uri|url|scan_name|ena_run|sra_run|fastq|raw_data|array_data|submitted_file|bam|cram|idat",
      ignore_case = TRUE
    )
  ) |
    stringr::str_detect(
      sdrf_long$value,
      stringr::regex(
        "\\.(fastq|fq|sra|bam|cram|idat|csv|tsv|txt|xls|xlsx|bed|bigwig|bw)(\\.gz)?([?#].*)?$|(^|[^A-Za-z0-9])(ERR|SRR|DRR)[0-9]+([^A-Za-z0-9]|$)",
        ignore_case = TRUE
      )
    )

  sdrf_long[possible_file_row, , drop = FALSE] %>%
    dplyr::transmute(
      Study_ID = .data$Study_ID,
      Sample_ID = .data$Sample_ID,
      inventory_source = "SDRF",
      field_raw = .data$field_raw,
      value = .data$value,
      file_class = vapply(
        .data$value,
        ae_aml_classify_file_value,
        character(1)
      )
    ) %>%
    dplyr::distinct()
}


ae_aml_biostudies_file_inventory <- function(
    accession,
    verbose = TRUE,
    biostudies_json = NULL) {

  answer <- biostudies_json

  if (is.null(answer)) {
    answer <- tryCatch(
      ae_aml_biostudies_json(accession),
      error = function(e) e
    )
  }

  if (inherits(answer, "error")) {
    if (verbose) {
      message(
        "Could not retrieve BioStudies inventory for ",
        accession,
        ": ",
        conditionMessage(answer)
      )
    }

    return(
      tibble::tibble(
        Study_ID = character(),
        Sample_ID = character(),
        inventory_source = character(),
        field_raw = character(),
        value = character(),
        file_class = character()
      )
    )
  }

  strings <- unique(
    ae_aml_clean_text(
      as.character(
        unlist(
          c(answer$study, answer$info),
          recursive = TRUE,
          use.names = FALSE
        )
      )
    )
  )

  strings <- strings[
    !is.na(strings) &
      (
        stringr::str_detect(
          strings,
          stringr::regex(
            "\\.(fastq|fq|sra|bam|cram|idat|csv|tsv|txt|xls|xlsx|bed|bedgraph|bigwig|bw|zip|tar|gz)([?#].*)?$",
            ignore_case = TRUE
          )
        ) |
          stringr::str_detect(
            strings,
            stringr::regex(
              "(^|[^A-Za-z0-9])(ERR|SRR|DRR)[0-9]+([^A-Za-z0-9]|$)",
              ignore_case = TRUE
            )
          )
      )
  ]

  tibble::tibble(
    Study_ID = accession,
    Sample_ID = NA_character_,
    inventory_source = "BioStudies API",
    field_raw = "BioStudies file or accession inventory",
    value = strings,
    file_class = vapply(
      strings,
      ae_aml_classify_file_value,
      character(1)
    )
  ) %>%
    dplyr::distinct()
}


# -----------------------------------------------------------------------------
# 6. Sample-level screening
# -----------------------------------------------------------------------------

ae_aml_collapse_where <- function(
    values,
    keep,
    max_values = 20L) {

  ae_aml_collapse_unique(
    values[keep],
    max_values = max_values
  )
}


ae_aml_make_sample_summary <- function(
    accession,
    sdrf_long,
    idf_long) {

  if (nrow(sdrf_long) == 0L) {
    return(tibble::tibble())
  }

  study_context <- ae_aml_collapse_unique(
    paste(
      idf_long$field_raw,
      idf_long$value,
      sep = ": "
    ),
    max_values = 200L
  )

  study_aml_terms <- ae_aml_dictionary_hit_text(
    study_context,
    AE_AML_DISEASE_TERMS
  )

  tagged <- sdrf_long %>%
    dplyr::mutate(
      field = ae_aml_normalize_field(
        .data$field_raw
      ),

      is_organism_field = stringr::str_detect(
        .data$field,
        "organism|species"
      ),

      is_disease_field = stringr::str_detect(
        .data$field,
        "disease|diagnos|phenotype|health_status|clinical_status|leukemia|leukaemia|case_control"
      ),

      is_source_field = stringr::str_detect(
        .data$field,
        "organism_part|tissue|source|specimen|material|cell_type|cell_line|biosource|body_site|sample_type"
      ),

      is_patient_id_field = (
        stringr::str_detect(
          .data$field,
          "(^|_)(patient|donor|subject|individual|participant)(_|$)"
        ) &
          !stringr::str_detect(
            .data$field,
            "age|sex|gender|disease|diagnos|status|treatment|group"
          )
      ),

      is_assay_field = stringr::str_detect(
        .data$field,
        "assay|technology|library|sequenc|methyl|rna|experimental_factor|factor_value"
      ),

      is_file_field = stringr::str_detect(
        .data$field,
        "file|uri|url|scan_name|ena_run|sra_run|fastq|raw_data|array_data|submitted_file|bam|cram|idat"
      )
    )

  sample_summary <- tagged %>%
    dplyr::group_by(
      .data$Study_ID,
      .data$Sample_ID
    ) %>%
    dplyr::summarise(
      all_sample_metadata = ae_aml_collapse_unique(
        paste(
          .data$field_raw,
          .data$value,
          sep = ": "
        ),
        max_values = 150L
      ),

      organism_values = ae_aml_collapse_where(
        .data$value,
        .data$is_organism_field
      ),

      disease_values = ae_aml_collapse_where(
        .data$value,
        .data$is_disease_field
      ),

      source_values = ae_aml_collapse_where(
        .data$value,
        .data$is_source_field
      ),

      explicit_patient_id = ae_aml_collapse_where(
        .data$value,
        .data$is_patient_id_field,
        max_values = 10L
      ),

      assay_values = ae_aml_collapse_where(
        .data$value,
        .data$is_assay_field
      ),

      file_values = ae_aml_collapse_where(
        .data$value,
        .data$is_file_field,
        max_values = 50L
      ),

      .groups = "drop"
    ) %>%
    dplyr::mutate(
      study_context = study_context,

      sample_search_text = coalesce(
        .data$all_sample_metadata,
        ""
      ),

      sample_disease_search_text = paste(
        coalesce(.data$disease_values, ""),
        .data$sample_search_text
      ),

      organism_search_text = paste(
        coalesce(.data$organism_values, ""),
        .data$sample_search_text,
        coalesce(.data$study_context, "")
      ),

      source_search_text = dplyr::if_else(
        !is.na(.data$source_values) &
          nzchar(.data$source_values),
        .data$source_values,
        .data$sample_search_text
      ),

      assay_sample_search_text = paste(
        coalesce(.data$assay_values, ""),
        coalesce(.data$file_values, ""),
        .data$sample_search_text
      ),

      sample_aml_terms = vapply(
        .data$sample_disease_search_text,
        ae_aml_dictionary_hit_text,
        character(1),
        dictionary = AE_AML_DISEASE_TERMS
      ),

      normal_control_terms = vapply(
        .data$sample_disease_search_text,
        ae_aml_dictionary_hit_text,
        character(1),
        dictionary = AE_AML_NORMAL_CONTROL_TERMS
      ),

      human_terms = vapply(
        .data$organism_search_text,
        ae_aml_dictionary_hit_text,
        character(1),
        dictionary = AE_AML_HUMAN_TERMS
      ),

      patient_source_terms = vapply(
        .data$source_search_text,
        ae_aml_dictionary_hit_text,
        character(1),
        dictionary = AE_AML_PATIENT_SOURCE_TERMS
      ),

      hematopoietic_source_terms = vapply(
        .data$source_search_text,
        ae_aml_dictionary_hit_text,
        character(1),
        dictionary =
          AE_AML_HEMATOPOIETIC_SOURCE_TERMS
      ),

      cell_line_terms = vapply(
        .data$source_search_text,
        ae_aml_dictionary_hit_text,
        character(1),
        dictionary = AE_AML_CELL_LINE_TERMS
      ),

      rnaseq_sample_terms = vapply(
        .data$assay_sample_search_text,
        ae_aml_dictionary_hit_text,
        character(1),
        dictionary = AE_AML_RNASEQ_TERMS
      ),

      rnaseq_study_terms = vapply(
        coalesce(.data$study_context, ""),
        ae_aml_dictionary_hit_text,
        character(1),
        dictionary = AE_AML_RNASEQ_TERMS
      ),

      rnaseq_terms = dplyr::coalesce(
        .data$rnaseq_sample_terms,
        .data$rnaseq_study_terms
      ),

      methylation_sample_terms = vapply(
        .data$assay_sample_search_text,
        ae_aml_dictionary_hit_text,
        character(1),
        dictionary = AE_AML_METHYLATION_TERMS
      ),

      methylation_study_terms = vapply(
        coalesce(.data$study_context, ""),
        ae_aml_dictionary_hit_text,
        character(1),
        dictionary = AE_AML_METHYLATION_TERMS
      ),

      methylation_terms = dplyr::coalesce(
        .data$methylation_sample_terms,
        .data$methylation_study_terms
      ),

      assay_assignment_source = dplyr::case_when(
        !is.na(.data$rnaseq_sample_terms) |
          !is.na(.data$methylation_sample_terms) ~
          "Sample-level assay metadata",

        !is.na(.data$rnaseq_study_terms) |
          !is.na(.data$methylation_study_terms) ~
          "Study-level assay metadata fallback",

        TRUE ~
          "Assay not identified"
      ),

      raw_sequence_terms = vapply(
        coalesce(.data$file_values, ""),
        ae_aml_dictionary_hit_text,
        character(1),
        dictionary = AE_AML_RAW_SEQUENCE_TERMS
      ),

      idat_terms = vapply(
        coalesce(.data$file_values, ""),
        ae_aml_dictionary_hit_text,
        character(1),
        dictionary = AE_AML_IDAT_TERMS
      ),

      aligned_file_terms = vapply(
        coalesce(.data$file_values, ""),
        ae_aml_dictionary_hit_text,
        character(1),
        dictionary = AE_AML_ALIGNED_TERMS
      ),

      sample_aml_flag =
        !is.na(.data$sample_aml_terms),

      normal_control_term_flag =
        !is.na(.data$normal_control_terms),

      human_flag = !is.na(.data$human_terms),

      patient_source_term_flag =
        !is.na(.data$patient_source_terms),

      hematopoietic_source_term_flag =
        !is.na(.data$hematopoietic_source_terms),

      cell_line_term_flag =
        !is.na(.data$cell_line_terms),

      rnaseq_flag = !is.na(.data$rnaseq_terms),

      methylation_flag =
        !is.na(.data$methylation_terms),

      likely_patient_source_flag =
        (
          .data$patient_source_term_flag |
            .data$hematopoietic_source_term_flag
        ) &
        !.data$cell_line_term_flag,

      normal_control_flag =
        .data$normal_control_term_flag &
        .data$likely_patient_source_flag,

      # Only use the study-level AML context when the sample has no explicit
      # disease/control annotation. This prevents healthy controls in an AML
      # study from inheriting the AML label.
      aml_context_fallback_flag =
        !.data$sample_aml_flag &
        !.data$normal_control_term_flag &
        (
          is.na(.data$disease_values) |
            !nzchar(.data$disease_values)
        ) &
        !is.na(study_aml_terms),

      aml_flag =
        .data$sample_aml_flag |
        .data$aml_context_fallback_flag,

      disease_group = dplyr::case_when(
        .data$sample_aml_flag &
          .data$normal_control_term_flag ~
          "Mixed or ambiguous AML/control annotation",

        .data$normal_control_flag ~
          "Normal/healthy hematopoietic control",

        .data$sample_aml_flag ~
          "AML",

        .data$aml_context_fallback_flag ~
          "AML study context; sample status unconfirmed",

        .data$normal_control_term_flag &
          !.data$likely_patient_source_flag ~
          "Normal/healthy control; source unclear",

        TRUE ~
          "Other or unclear"
      ),

      patient_source_status = dplyr::case_when(
        .data$likely_patient_source_flag ~
          "Likely patient/donor-derived hematopoietic sample",

        (
          .data$patient_source_term_flag |
            .data$hematopoietic_source_term_flag
        ) &
          .data$cell_line_term_flag ~
          "Mixed or ambiguous source metadata",

        .data$cell_line_term_flag ~
          "Likely cell-line/model sample",

        TRUE ~
          "Patient/donor source not confirmed"
      ),

      eligible_patient_sample_flag =
        .data$human_flag &
        .data$likely_patient_source_flag &
        (
          .data$aml_flag |
            .data$normal_control_flag
        ),

      high_confidence_eligible_sample_flag =
        .data$human_flag &
        .data$likely_patient_source_flag &
        (
          .data$sample_aml_flag |
            .data$normal_control_flag
        ),

      sample_raw_sequence_flag =
        !is.na(.data$raw_sequence_terms),

      sample_idat_flag =
        !is.na(.data$idat_terms),

      sample_aligned_file_flag =
        !is.na(.data$aligned_file_terms),

      rnaseq_biological_candidate =
        .data$eligible_patient_sample_flag &
        .data$rnaseq_flag,

      methylation_biological_candidate =
        .data$eligible_patient_sample_flag &
        .data$methylation_flag,

      aml_rnaseq_candidate =
        .data$aml_flag &
        .data$human_flag &
        .data$likely_patient_source_flag &
        .data$rnaseq_flag,

      normal_control_rnaseq_candidate =
        .data$normal_control_flag &
        .data$human_flag &
        .data$rnaseq_flag,

      aml_methylation_candidate =
        .data$aml_flag &
        .data$human_flag &
        .data$likely_patient_source_flag &
        .data$methylation_flag,

      normal_control_methylation_candidate =
        .data$normal_control_flag &
        .data$human_flag &
        .data$methylation_flag
    ) %>%
    dplyr::select(
      .data$Study_ID,
      .data$Sample_ID,
      .data$explicit_patient_id,
      .data$organism_values,
      .data$disease_values,
      .data$disease_group,
      .data$source_values,
      .data$assay_values,
      .data$assay_assignment_source,
      .data$file_values,

      .data$sample_aml_flag,
      .data$aml_context_fallback_flag,
      .data$aml_flag,
      .data$normal_control_term_flag,
      .data$normal_control_flag,

      .data$human_flag,
      .data$likely_patient_source_flag,
      .data$patient_source_status,
      .data$cell_line_term_flag,

      .data$eligible_patient_sample_flag,
      .data$high_confidence_eligible_sample_flag,

      .data$rnaseq_flag,
      .data$methylation_flag,
      .data$sample_raw_sequence_flag,
      .data$sample_idat_flag,
      .data$sample_aligned_file_flag,

      .data$rnaseq_biological_candidate,
      .data$methylation_biological_candidate,
      .data$aml_rnaseq_candidate,
      .data$normal_control_rnaseq_candidate,
      .data$aml_methylation_candidate,
      .data$normal_control_methylation_candidate,

      .data$sample_aml_terms,
      .data$normal_control_terms,
      .data$human_terms,
      .data$patient_source_terms,
      .data$hematopoietic_source_terms,
      .data$cell_line_terms,
      .data$rnaseq_terms,
      .data$methylation_terms,
      .data$raw_sequence_terms,
      .data$idat_terms,
      .data$aligned_file_terms
    )

  sample_summary
}


# -----------------------------------------------------------------------------
# 7. Evidence table
# -----------------------------------------------------------------------------

ae_aml_tag_evidence_rows <- function(
    metadata,
    source_name,
    dictionary,
    criterion,
    sample_column = NULL) {

  if (nrow(metadata) == 0L) {
    return(tibble::tibble())
  }

  row_text <- paste(
    metadata$field_raw,
    metadata$value,
    sep = ": "
  )

  matched_terms <- vapply(
    row_text,
    ae_aml_dictionary_hit_text,
    character(1),
    dictionary = dictionary
  )

  keep <- !is.na(matched_terms)

  if (!any(keep)) {
    return(tibble::tibble())
  }

  sample_id <- rep(
    NA_character_,
    nrow(metadata)
  )

  if (!is.null(sample_column) &&
      sample_column %in% names(metadata)) {
    sample_id <- as.character(
      metadata[[sample_column]]
    )
  }

  tibble::tibble(
    Study_ID = metadata$Study_ID[keep],
    Sample_ID = sample_id[keep],
    evidence_source = source_name,
    criterion = criterion,
    field_raw = metadata$field_raw[keep],
    value = metadata$value[keep],
    matched_terms = matched_terms[keep]
  )
}


ae_aml_make_evidence <- function(
    idf_long,
    sdrf_long,
    file_inventory) {

  evidence <- dplyr::bind_rows(
    ae_aml_tag_evidence_rows(
      idf_long,
      "IDF",
      AE_AML_DISEASE_TERMS,
      "AML disease"
    ),
    ae_aml_tag_evidence_rows(
      sdrf_long,
      "SDRF",
      AE_AML_DISEASE_TERMS,
      "AML disease",
      sample_column = "Sample_ID"
    ),
    ae_aml_tag_evidence_rows(
      idf_long,
      "IDF",
      AE_AML_NORMAL_CONTROL_TERMS,
      "Normal/healthy control"
    ),
    ae_aml_tag_evidence_rows(
      sdrf_long,
      "SDRF",
      AE_AML_NORMAL_CONTROL_TERMS,
      "Normal/healthy control",
      sample_column = "Sample_ID"
    ),
    ae_aml_tag_evidence_rows(
      idf_long,
      "IDF",
      AE_AML_HEMATOPOIETIC_SOURCE_TERMS,
      "Hematopoietic source"
    ),
    ae_aml_tag_evidence_rows(
      sdrf_long,
      "SDRF",
      AE_AML_HEMATOPOIETIC_SOURCE_TERMS,
      "Hematopoietic source",
      sample_column = "Sample_ID"
    ),
    ae_aml_tag_evidence_rows(
      idf_long,
      "IDF",
      AE_AML_HUMAN_TERMS,
      "Human organism"
    ),
    ae_aml_tag_evidence_rows(
      sdrf_long,
      "SDRF",
      AE_AML_HUMAN_TERMS,
      "Human organism",
      sample_column = "Sample_ID"
    ),
    ae_aml_tag_evidence_rows(
      idf_long,
      "IDF",
      AE_AML_PATIENT_SOURCE_TERMS,
      "Patient-derived source"
    ),
    ae_aml_tag_evidence_rows(
      sdrf_long,
      "SDRF",
      AE_AML_PATIENT_SOURCE_TERMS,
      "Patient-derived source",
      sample_column = "Sample_ID"
    ),
    ae_aml_tag_evidence_rows(
      idf_long,
      "IDF",
      AE_AML_RNASEQ_TERMS,
      "RNA-seq assay"
    ),
    ae_aml_tag_evidence_rows(
      sdrf_long,
      "SDRF",
      AE_AML_RNASEQ_TERMS,
      "RNA-seq assay",
      sample_column = "Sample_ID"
    ),
    ae_aml_tag_evidence_rows(
      idf_long,
      "IDF",
      AE_AML_METHYLATION_TERMS,
      "Methylation assay"
    ),
    ae_aml_tag_evidence_rows(
      sdrf_long,
      "SDRF",
      AE_AML_METHYLATION_TERMS,
      "Methylation assay",
      sample_column = "Sample_ID"
    )
  )

  file_evidence <- file_inventory %>%
    dplyr::filter(
      .data$file_class %in% c(
        "Raw sequencing file or run",
        "Raw methylation array file",
        "Aligned sequencing file"
      )
    ) %>%
    dplyr::transmute(
      Study_ID = .data$Study_ID,
      Sample_ID = .data$Sample_ID,
      evidence_source = .data$inventory_source,
      criterion = .data$file_class,
      field_raw = .data$field_raw,
      value = .data$value,
      matched_terms = .data$file_class
    )

  dplyr::bind_rows(
    evidence,
    file_evidence
  ) %>%
    dplyr::distinct()
}


# -----------------------------------------------------------------------------
# 8. Study-title helpers
# -----------------------------------------------------------------------------

# Extract a project title from the quick pre-screen when available, otherwise
# fall back to title-like fields in the IDF metadata.
ae_aml_extract_project_title <- function(
    idf_long = NULL,
    prescreen_summary = NULL) {

  if (
    !is.null(prescreen_summary) &&
      nrow(prescreen_summary) > 0L &&
      "title_or_first_value" %in%
        names(prescreen_summary)
  ) {
    title <- ae_aml_first_nonempty(
      prescreen_summary$title_or_first_value
    )

    if (!is.na(title) && nzchar(title)) {
      return(title)
    }
  }

  if (
    is.null(idf_long) ||
      nrow(idf_long) == 0L
  ) {
    return(NA_character_)
  }

  field_normalized <- ae_aml_normalize_field(
    idf_long$field_raw
  )

  preferred_fields <- c(
    "investigation_title",
    "study_title",
    "experiment_title",
    "project_title",
    "title"
  )

  for (preferred_field in preferred_fields) {
    values <- idf_long$value[
      field_normalized == preferred_field
    ]

    title <- ae_aml_first_nonempty(
      values
    )

    if (!is.na(title) && nzchar(title)) {
      return(title)
    }
  }

  title_rows <- stringr::str_detect(
    field_normalized,
    "(^|_)(investigation|study|experiment|project)?_?title($|_)"
  )

  ae_aml_first_nonempty(
    idf_long$value[title_rows]
  )
}


# Retrieve study titles for the project IDs in a data frame.
#
# Parameters
# ----------
# data:
#   A data.frame or tibble containing ArrayExpress/BioStudies accessions.
#
# id_col:
#   Character string naming the accession column, for example "Study_ID".
#
# search_result:
#   Optional object returned by biostudies_search_many(). Titles found in
#   search_result$unique_studies are used before making any API requests.
#
# query_api:
#   Query the BioStudies detailed-study API for accessions not found in
#   search_result.
#
# Returns
# -------
# A character vector aligned to the rows of data. The vector is named with the
# accession values.
ae_get_study_titles <- function(
    data,
    id_col,
    search_result = NULL,
    query_api = TRUE,
    pause_seconds = 0.15,
    verbose = TRUE) {

  if (!is.data.frame(data)) {
    stop("data must be a data.frame or tibble.")
  }

  if (
    length(id_col) != 1L ||
      !is.character(id_col) ||
      !id_col %in% names(data)
  ) {
    stop(
      "id_col must be the name of a column in data."
    )
  }

  ids <- toupper(
    trimws(
      as.character(
        data[[id_col]]
      )
    )
  )

  unique_ids <- unique(
    ids[
      !is.na(ids) &
        nzchar(ids)
    ]
  )

  title_lookup <- stats::setNames(
    rep(
      NA_character_,
      length(unique_ids)
    ),
    unique_ids
  )

  # First use the BioStudies search results when supplied. This avoids
  # unnecessary API calls because the search table already contains titles.
  if (!is.null(search_result)) {

    search_table <- NULL

    if (
      is.list(search_result) &&
        "unique_studies" %in%
          names(search_result)
    ) {
      search_table <-
        search_result$unique_studies
    } else if (
      is.data.frame(search_result)
    ) {
      search_table <- search_result
    }

    if (!is.null(search_table)) {

      accession_candidates <- c(
        "accession",
        "Study_ID",
        "study_id",
        "project_id"
      )

      title_candidates <- c(
        "title",
        "project_title",
        "title_or_first_value"
      )

      accession_col <- accession_candidates[
        accession_candidates %in%
          names(search_table)
      ]

      title_col <- title_candidates[
        title_candidates %in%
          names(search_table)
      ]

      if (
        length(accession_col) > 0L &&
          length(title_col) > 0L
      ) {
        search_ids <- toupper(
          trimws(
            as.character(
              search_table[[accession_col[1]]]
            )
          )
        )

        search_titles <- as.character(
          search_table[[title_col[1]]]
        )

        keep <- !is.na(search_ids) &
          nzchar(search_ids) &
          !is.na(search_titles) &
          nzchar(trimws(search_titles))

        search_map <- stats::setNames(
          search_titles[keep],
          search_ids[keep]
        )

        matched <- intersect(
          names(title_lookup),
          names(search_map)
        )

        title_lookup[matched] <-
          search_map[matched]
      }
    }
  }

  missing_ids <- names(title_lookup)[
    is.na(title_lookup) |
      !nzchar(title_lookup)
  ]

  if (
    isTRUE(query_api) &&
      length(missing_ids) > 0L
  ) {

    for (i in seq_along(missing_ids)) {
      accession <- missing_ids[i]

      if (verbose) {
        message(
          "Retrieving title for ",
          accession,
          " (",
          i,
          "/",
          length(missing_ids),
          ")"
        )
      }

      title <- tryCatch(
        {
          json <- ae_aml_biostudies_json(
            accession
          )

          flat_json <-
            ae_aml_flatten_biostudies_json(
              biostudies_json = json,
              accession = accession
            )

          ae_aml_extract_prescreen_title(
            flat_json
          )
        },
        error = function(e) {
          if (verbose) {
            message(
              "Title lookup failed for ",
              accession,
              ": ",
              conditionMessage(e)
            )
          }

          NA_character_
        }
      )

      title_lookup[accession] <- title

      if (
        pause_seconds > 0 &&
          i < length(missing_ids)
      ) {
        Sys.sleep(pause_seconds)
      }
    }
  }

  output <- unname(
    title_lookup[ids]
  )

  names(output) <- ids
  output
}


# Add a title column directly to a data frame.
ae_add_study_titles <- function(
    data,
    id_col,
    title_col = "project_title",
    search_result = NULL,
    query_api = TRUE,
    pause_seconds = 0.15,
    verbose = TRUE) {

  titles <- ae_get_study_titles(
    data = data,
    id_col = id_col,
    search_result = search_result,
    query_api = query_api,
    pause_seconds = pause_seconds,
    verbose = verbose
  )

  output <- data

  output[[title_col]] <- unname(
    titles
  )

  # Move the title immediately after the project accession.
  id_position <- match(
    id_col,
    names(output)
  )

  title_position <- match(
    title_col,
    names(output)
  )

  new_order <- names(output)

  new_order <- new_order[
    new_order != title_col
  ]

  new_order <- append(
    new_order,
    title_col,
    after = id_position
  )

  output[
    new_order
  ]
}


# Add titles already present in prescreen_summary to an existing aml_screen
# result. No API queries are needed.
ae_add_titles_to_screening_result <- function(
    screening_result,
    overwrite_existing = FALSE) {

  if (!is.list(screening_result)) {
    stop(
      "screening_result must be the list returned by ",
      "ae_screen_aml_accessions()."
    )
  }

  if (
    is.null(
      screening_result$study_summary
    ) ||
      is.null(
        screening_result$prescreen_summary
      )
  ) {
    stop(
      "screening_result must contain study_summary ",
      "and prescreen_summary."
    )
  }

  title_lookup <-
    screening_result$prescreen_summary %>%
    dplyr::select(
      .data$Study_ID,
      project_title =
        .data$title_or_first_value
    ) %>%
    dplyr::filter(
      !is.na(.data$Study_ID),
      nzchar(.data$Study_ID)
    ) %>%
    dplyr::distinct(
      .data$Study_ID,
      .keep_all = TRUE
    )

  study_summary <-
    screening_result$study_summary

  if (
    "project_title" %in%
      names(study_summary)
  ) {
    if (isTRUE(overwrite_existing)) {
      study_summary <-
        study_summary %>%
        dplyr::select(
          -.data$project_title
        )
    } else {
      missing_title <- is.na(
        study_summary$project_title
      ) |
        !nzchar(
          study_summary$project_title
        )

      replacement <- title_lookup$project_title[
        match(
          study_summary$Study_ID,
          title_lookup$Study_ID
        )
      ]

      study_summary$project_title[
        missing_title
      ] <- replacement[
        missing_title
      ]

      screening_result$study_summary <-
        study_summary

      return(
        screening_result
      )
    }
  }

  study_summary <-
    study_summary %>%
    dplyr::left_join(
      title_lookup,
      by = "Study_ID"
    ) %>%
    dplyr::relocate(
      .data$project_title,
      .after = .data$Study_ID
    )

  screening_result$study_summary <-
    study_summary

  screening_result
}


# -----------------------------------------------------------------------------
# 9. Study-level summary
# -----------------------------------------------------------------------------

ae_aml_make_study_summary <- function(
    accession,
    project_title = NA_character_,
    sample_summary,
    sdrf_long,
    idf_long,
    file_inventory,
    download_method) {

  study_text <- ae_aml_collapse_unique(
    c(
      paste(
        idf_long$field_raw,
        idf_long$value,
        sep = ": "
      ),
      paste(
        sdrf_long$field_raw,
        sdrf_long$value,
        sep = ": "
      )
    ),
    max_values = 500L
  )

  explicit_patient_values <- sdrf_long %>%
    dplyr::mutate(
      field = ae_aml_normalize_field(
        .data$field_raw
      ),
      is_patient_id_field = (
        stringr::str_detect(
          .data$field,
          "(^|_)(patient|donor|subject|individual|participant)(_|$)"
        ) &
          !stringr::str_detect(
            .data$field,
            "age|sex|gender|disease|diagnos|status|treatment|group"
          )
      )
    ) %>%
    dplyr::filter(
      .data$is_patient_id_field,
      !is.na(.data$value),
      nzchar(.data$value)
    ) %>%
    dplyr::distinct(.data$value)

  n_explicit_patients <- nrow(
    explicit_patient_values
  )

  if (n_explicit_patients == 0L) {
    n_explicit_patients <- NA_integer_
    patient_count_method <-
      "No explicit patient/donor/subject identifier field found"
  } else {
    patient_count_method <-
      "Distinct values in explicit patient/donor/subject fields"
  }

  file_counts <- file_inventory %>%
    dplyr::count(
      .data$file_class,
      name = "n"
    )

  count_class <- function(class_name) {
    value <- file_counts$n[
      file_counts$file_class == class_name
    ]

    if (length(value) == 0L) {
      return(0L)
    }

    as.integer(sum(value))
  }

  has_aml <- any(
    sample_summary$aml_flag,
    na.rm = TRUE
  ) ||
    ae_aml_has_dictionary_hit(
      study_text,
      AE_AML_DISEASE_TERMS
    )

  has_normal_control <- any(
    sample_summary$normal_control_flag,
    na.rm = TRUE
  )

  has_human <- any(
    sample_summary$human_flag,
    na.rm = TRUE
  ) ||
    ae_aml_has_dictionary_hit(
      study_text,
      AE_AML_HUMAN_TERMS
    )

  has_patient_source <- any(
    sample_summary$likely_patient_source_flag,
    na.rm = TRUE
  ) ||
    ae_aml_has_dictionary_hit(
      study_text,
      AE_AML_PATIENT_SOURCE_TERMS
    ) ||
    ae_aml_has_dictionary_hit(
      study_text,
      AE_AML_HEMATOPOIETIC_SOURCE_TERMS
    )

  has_cell_line_terms <- any(
    sample_summary$cell_line_term_flag,
    na.rm = TRUE
  )

  has_rnaseq <- any(
    sample_summary$rnaseq_flag,
    na.rm = TRUE
  ) ||
    ae_aml_has_dictionary_hit(
      study_text,
      AE_AML_RNASEQ_TERMS
    )

  has_methylation <- any(
    sample_summary$methylation_flag,
    na.rm = TRUE
  ) ||
    ae_aml_has_dictionary_hit(
      study_text,
      AE_AML_METHYLATION_TERMS
    )

  has_eligible_patient_samples <- any(
    sample_summary$eligible_patient_sample_flag,
    na.rm = TRUE
  )

  n_raw_sequence_records <- count_class(
    "Raw sequencing file or run"
  )

  n_idat_records <- count_class(
    "Raw methylation array file"
  )

  n_aligned_records <- count_class(
    "Aligned sequencing file"
  )

  raw_sequence_available <-
    n_raw_sequence_records > 0L

  idat_available <-
    n_idat_records > 0L

  aligned_sequence_available <-
    n_aligned_records > 0L

  rnaseq_biological_candidate <-
    has_eligible_patient_samples &&
    has_rnaseq

  methylation_biological_candidate <-
    has_eligible_patient_samples &&
    has_methylation

  rnaseq_raw_candidate <-
    rnaseq_biological_candidate &&
    raw_sequence_available

  methylation_raw_candidate <-
    methylation_biological_candidate &&
    (
      idat_available ||
        raw_sequence_available
    )

  study_scope <- dplyr::case_when(
    has_aml && has_normal_control ~
      "AML with normal/healthy hematopoietic controls",

    has_aml ~
      "AML samples identified",

    has_normal_control ~
      "Normal/healthy hematopoietic controls identified",

    TRUE ~
      "Eligible AML/control samples not confirmed"
  )

  rnaseq_raw_status <- dplyr::case_when(
    !has_rnaseq ~
      "RNA-seq not identified",

    raw_sequence_available ~
      "Raw sequencing files or run accessions identified",

    aligned_sequence_available ~
      "Only aligned BAM/CRAM evidence identified",

    TRUE ~
      "RNA-seq identified; raw files not confirmed from metadata"
  )

  methylation_raw_status <- dplyr::case_when(
    !has_methylation ~
      "Methylation assay not identified",

    idat_available ~
      "Raw IDAT files identified",

    raw_sequence_available ~
      "Raw sequencing files identified; review for WGBS/RRBS",

    aligned_sequence_available ~
      "Only aligned BAM/CRAM evidence identified",

    TRUE ~
      "Methylation assay identified; raw files not confirmed from metadata"
  )

  screening_tier <- dplyr::case_when(
    rnaseq_raw_candidate ||
      methylation_raw_candidate ~
      "Strong candidate",

    rnaseq_biological_candidate ||
      methylation_biological_candidate ~
      "Biologically relevant; raw availability needs review",

    (
      has_aml ||
        has_normal_control
    ) &&
      has_human &&
      (
        has_rnaseq ||
          has_methylation
      ) ~
      "Potential candidate; sample source or disease group unclear",

    has_aml ||
      has_normal_control ||
      has_human ||
      has_rnaseq ||
      has_methylation ~
      "Partial match; manual review needed",

    TRUE ~
      "Criteria not supported by retrieved metadata"
  )

  tibble::tibble(
    Study_ID = accession,
    project_title = project_title,
    study_scope = study_scope,
    screening_tier = screening_tier,
    download_method = download_method,

    aml_flag = has_aml,
    normal_control_flag =
      has_normal_control,
    human_flag = has_human,
    patient_source_flag =
      has_patient_source,
    eligible_patient_samples_flag =
      has_eligible_patient_samples,
    cell_line_terms_present =
      has_cell_line_terms,

    rnaseq_flag = has_rnaseq,
    methylation_flag = has_methylation,

    raw_sequence_available =
      raw_sequence_available,
    idat_available = idat_available,
    aligned_sequence_available =
      aligned_sequence_available,

    rnaseq_biological_candidate =
      rnaseq_biological_candidate,
    rnaseq_raw_candidate =
      rnaseq_raw_candidate,

    methylation_biological_candidate =
      methylation_biological_candidate,
    methylation_raw_candidate =
      methylation_raw_candidate,

    rnaseq_raw_status = rnaseq_raw_status,
    methylation_raw_status =
      methylation_raw_status,

    n_samples = dplyr::n_distinct(
      sample_summary$Sample_ID
    ),

    n_explicit_patients =
      n_explicit_patients,

    patient_count_method =
      patient_count_method,

    n_likely_patient_source_samples =
      sum(
        sample_summary$likely_patient_source_flag,
        na.rm = TRUE
      ),

    n_eligible_patient_samples =
      sum(
        sample_summary$eligible_patient_sample_flag,
        na.rm = TRUE
      ),

    n_high_confidence_eligible_samples =
      sum(
        sample_summary$high_confidence_eligible_sample_flag,
        na.rm = TRUE
      ),

    n_aml_samples =
      sum(
        sample_summary$aml_flag &
          sample_summary$human_flag &
          sample_summary$likely_patient_source_flag,
        na.rm = TRUE
      ),

    n_normal_control_samples =
      sum(
        sample_summary$normal_control_flag &
          sample_summary$human_flag,
        na.rm = TRUE
      ),

    n_rnaseq_candidate_samples =
      sum(
        sample_summary$rnaseq_biological_candidate,
        na.rm = TRUE
      ),

    n_methylation_candidate_samples =
      sum(
        sample_summary$methylation_biological_candidate,
        na.rm = TRUE
      ),

    n_aml_rnaseq_samples =
      sum(
        sample_summary$aml_rnaseq_candidate,
        na.rm = TRUE
      ),

    n_normal_control_rnaseq_samples =
      sum(
        sample_summary$normal_control_rnaseq_candidate,
        na.rm = TRUE
      ),

    n_aml_methylation_samples =
      sum(
        sample_summary$aml_methylation_candidate,
        na.rm = TRUE
      ),

    n_normal_control_methylation_samples =
      sum(
        sample_summary$normal_control_methylation_candidate,
        na.rm = TRUE
      ),

    n_raw_sequence_records =
      n_raw_sequence_records,
    n_idat_records =
      n_idat_records,
    n_aligned_records =
      n_aligned_records,

    aml_terms = ae_aml_dictionary_hit_text(
      study_text,
      AE_AML_DISEASE_TERMS
    ),

    normal_control_terms =
      ae_aml_dictionary_hit_text(
        study_text,
        AE_AML_NORMAL_CONTROL_TERMS
      ),

    human_terms = ae_aml_dictionary_hit_text(
      study_text,
      AE_AML_HUMAN_TERMS
    ),

    patient_source_terms =
      ae_aml_dictionary_hit_text(
        study_text,
        AE_AML_PATIENT_SOURCE_TERMS
      ),

    hematopoietic_source_terms =
      ae_aml_dictionary_hit_text(
        study_text,
        AE_AML_HEMATOPOIETIC_SOURCE_TERMS
      ),

    cell_line_terms =
      ae_aml_dictionary_hit_text(
        study_text,
        AE_AML_CELL_LINE_TERMS
      ),

    rnaseq_terms =
      ae_aml_dictionary_hit_text(
        study_text,
        AE_AML_RNASEQ_TERMS
      ),

    methylation_terms =
      ae_aml_dictionary_hit_text(
        study_text,
        AE_AML_METHYLATION_TERMS
      )
  )
}


# -----------------------------------------------------------------------------
# 10. One-study and batch wrappers
# -----------------------------------------------------------------------------

ae_screen_one_aml_study <- function(
    accession,
    base_dir = "ArrayExpress_AML_Screen",
    overwrite = FALSE,
    use_api_fallback = TRUE,
    query_biostudies_inventory = TRUE,
    verbose = TRUE,
    biostudies_json = NULL,
    prescreen_summary = NULL) {

  accession <- toupper(trimws(accession))

  downloaded <- ae_aml_download_magetab(
    accession = accession,
    base_dir = base_dir,
    overwrite = overwrite,
    use_api_fallback = use_api_fallback,
    verbose = verbose,
    biostudies_json = biostudies_json
  )

  study_dir <- downloaded$path

  sdrf_files <- list.files(
    study_dir,
    pattern = "sdrf(?:\\.txt)?$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  idf_files <- list.files(
    study_dir,
    pattern = "idf(?:\\.txt)?$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(sdrf_files) == 0L) {
    stop(
      "No SDRF file was found for ",
      accession,
      "."
    )
  }

  if (length(idf_files) == 0L) {
    stop(
      "No IDF file was found for ",
      accession,
      "."
    )
  }

  sdrf_parsed <- purrr::map(
    sdrf_files,
    ae_aml_read_sdrf,
    accession = accession
  )

  sdrf_wide <- purrr::map_dfr(
    sdrf_parsed,
    "wide"
  )

  sdrf_long <- purrr::map_dfr(
    sdrf_parsed,
    "long"
  )

  idf_long <- purrr::map_dfr(
    idf_files,
    ae_aml_read_idf,
    accession = accession
  )

  sample_summary <- ae_aml_make_sample_summary(
    accession = accession,
    sdrf_long = sdrf_long,
    idf_long = idf_long
  )

  sdrf_inventory <- ae_aml_sdrf_file_inventory(
    sdrf_long
  )

  api_inventory <- tibble::tibble()

  if (isTRUE(query_biostudies_inventory)) {
    api_inventory <-
      ae_aml_biostudies_file_inventory(
        accession = accession,
        verbose = verbose,
        biostudies_json = biostudies_json
      )
  }

  file_inventory <- dplyr::bind_rows(
    sdrf_inventory,
    api_inventory
  ) %>%
    dplyr::distinct()

  evidence <- ae_aml_make_evidence(
    idf_long = idf_long,
    sdrf_long = sdrf_long,
    file_inventory = file_inventory
  )

  project_title <- ae_aml_extract_project_title(
    idf_long = idf_long,
    prescreen_summary = prescreen_summary
  )

  study_summary <- ae_aml_make_study_summary(
    accession = accession,
    project_title = project_title,
    sample_summary = sample_summary,
    sdrf_long = sdrf_long,
    idf_long = idf_long,
    file_inventory = file_inventory,
    download_method = downloaded$download_method
  )

  eligible_sample_summary <- sample_summary %>%
    dplyr::filter(
      .data$eligible_patient_sample_flag
    )

  result <- list(
    Study_ID = accession,
    study_summary = study_summary,
    sample_summary = sample_summary,
    eligible_sample_summary =
      eligible_sample_summary,
    evidence = evidence,
    file_inventory = file_inventory,
    sdrf_wide = sdrf_wide,
    sdrf_long = sdrf_long,
    idf_long = idf_long,
    prescreen_summary = prescreen_summary
  )

  saveRDS(
    result,
    file = file.path(
      study_dir,
      paste0(
        accession,
        "_AML_screening.rds"
      )
    )
  )

  result
}


ae_screen_aml_accessions <- function(
    accessions,
    base_dir = "ArrayExpress_AML_Screen",
    overwrite = FALSE,
    use_api_fallback = TRUE,
    query_biostudies_inventory = TRUE,
    verbose = TRUE,
    prescreen = TRUE,
    prescreen_require_raw_evidence = TRUE,
    prescreen_require_human_evidence = FALSE,
    prescreen_exclude_clear_nonhuman = TRUE,
    prescreen_include_normal_controls = TRUE,
    prescreen_on_error = c("continue", "skip")) {

  prescreen_on_error <- match.arg(
    prescreen_on_error
  )

  accessions <- unique(
    toupper(
      trimws(
        as.character(accessions)
      )
    )
  )

  accessions <- accessions[
    !is.na(accessions) &
      nzchar(accessions)
  ]

  results <- stats::setNames(
    vector(
      "list",
      length(accessions)
    ),
    accessions
  )

  errors <- tibble::tibble(
    Study_ID = character(),
    stage = character(),
    error = character()
  )

  prescreen_summary <- tibble::tibble()
  prescreen_evidence <- tibble::tibble()
  skipped_prescreen <- tibble::tibble()

  for (accession in accessions) {
    if (verbose) {
      message(
        "\n========== ",
        accession,
        " =========="
      )
    }

    quick <- NULL
    preloaded_json <- NULL
    this_prescreen_summary <- NULL

    if (isTRUE(prescreen)) {
      quick <- tryCatch(
        ae_aml_quick_prescreen(
          accession = accession,
          require_raw_evidence =
            prescreen_require_raw_evidence,
          require_human_evidence =
            prescreen_require_human_evidence,
          exclude_clear_nonhuman =
            prescreen_exclude_clear_nonhuman,
          include_normal_controls =
            prescreen_include_normal_controls,
          verbose = verbose
        ),
        error = function(e) e
      )

      if (inherits(quick, "error")) {
        errors <- dplyr::bind_rows(
          errors,
          tibble::tibble(
            Study_ID = accession,
            stage = "BioStudies pre-screen",
            error = conditionMessage(quick)
          )
        )

        fallback_continue <-
          identical(
            prescreen_on_error,
            "continue"
          )

        this_prescreen_summary <- tibble::tibble(
          Study_ID = accession,
          title_or_first_value = NA_character_,
          prescreen_status = ifelse(
            fallback_continue,
            "Continue: pre-screen failed; full parsing used as fallback",
            "Skip: pre-screen failed"
          ),
          continue_full_parse =
            fallback_continue,
          aml_flag = NA,
          human_flag = NA,
          nonhuman_flag = NA,
          rnaseq_flag = NA,
          methylation_flag = NA,
          methylation_array_flag = NA,
          methylation_sequencing_flag = NA,
          raw_sequence_flag = NA,
          raw_methylation_array_flag = NA,
          generic_raw_flag = NA,
          rnaseq_raw_candidate = NA,
          methylation_raw_candidate = NA,
          aml_terms = NA_character_,
          human_terms = NA_character_,
          nonhuman_terms = NA_character_,
          rnaseq_terms = NA_character_,
          methylation_terms = NA_character_,
          raw_sequence_terms = NA_character_,
          raw_methylation_array_terms =
            NA_character_,
          generic_raw_terms = NA_character_
        )

        prescreen_summary <- dplyr::bind_rows(
          prescreen_summary,
          this_prescreen_summary
        )

        if (!fallback_continue) {
          skipped_prescreen <- dplyr::bind_rows(
            skipped_prescreen,
            this_prescreen_summary
          )
          next
        }
      } else {
        this_prescreen_summary <- quick$summary
        preloaded_json <- quick$json

        prescreen_summary <- dplyr::bind_rows(
          prescreen_summary,
          quick$summary
        )

        prescreen_evidence <- dplyr::bind_rows(
          prescreen_evidence,
          quick$evidence
        )

        if (verbose) {
          message(
            quick$summary$prescreen_status[1]
          )
        }

        if (!isTRUE(
          quick$summary$continue_full_parse[1]
        )) {
          skipped_prescreen <- dplyr::bind_rows(
            skipped_prescreen,
            quick$summary
          )
          next
        }
      }
    }

    result <- tryCatch(
      ae_screen_one_aml_study(
        accession = accession,
        base_dir = base_dir,
        overwrite = overwrite,
        use_api_fallback =
          use_api_fallback,
        query_biostudies_inventory =
          query_biostudies_inventory,
        verbose = verbose,
        biostudies_json =
          preloaded_json,
        prescreen_summary =
          this_prescreen_summary
      ),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      errors <- dplyr::bind_rows(
        errors,
        tibble::tibble(
          Study_ID = accession,
          stage = "Full MAGE-TAB parsing",
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

      next
    }

    results[[accession]] <- result
  }

  valid_results <- results[
    vapply(
      results,
      function(x) {
        is.list(x) &&
          !is.null(x$study_summary)
      },
      logical(1)
    )
  ]

  combined <- list(
    prescreen_summary =
      prescreen_summary,

    prescreen_evidence =
      prescreen_evidence,

    skipped_prescreen =
      skipped_prescreen,

    study_summary = purrr::map_dfr(
      valid_results,
      "study_summary"
    ),

    sample_summary = purrr::map_dfr(
      valid_results,
      "sample_summary"
    ),

    eligible_sample_summary = purrr::map_dfr(
      valid_results,
      "eligible_sample_summary"
    ),

    evidence = purrr::map_dfr(
      valid_results,
      "evidence"
    ),

    file_inventory = purrr::map_dfr(
      valid_results,
      "file_inventory"
    ),

    sdrf_wide = purrr::map_dfr(
      valid_results,
      "sdrf_wide"
    ),

    sdrf_long = purrr::map_dfr(
      valid_results,
      "sdrf_long"
    ),

    idf_long = purrr::map_dfr(
      valid_results,
      "idf_long"
    ),

    errors = errors,
    per_study = valid_results
  )

  if (nrow(combined$study_summary) > 0L) {
    combined$study_summary <-
      combined$study_summary %>%
      dplyr::arrange(
        factor(
          .data$screening_tier,
          levels = c(
            "Strong candidate",
            "Biologically relevant; raw availability needs review",
            "Potential candidate; patient source unclear",
            "Partial match; manual review needed",
            "Criteria not supported by retrieved metadata"
          )
        ),
        .data$Study_ID
      )
  }

  combined
}


# -----------------------------------------------------------------------------
# 11. Write output tables
# -----------------------------------------------------------------------------

ae_write_aml_screening_results <- function(
    screening_results,
    out_dir = "ArrayExpress_AML_Screening_Tables",
    prefix = "AML_ArrayExpress") {

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  tables <- c(
    "prescreen_summary",
    "prescreen_evidence",
    "skipped_prescreen",
    "study_summary",
    "sample_summary",
    "eligible_sample_summary",
    "evidence",
    "file_inventory",
    "sdrf_wide",
    "sdrf_long",
    "idf_long",
    "errors"
  )

  output_files <- stats::setNames(
    file.path(
      out_dir,
      paste0(
        prefix,
        "_",
        tables,
        ".tsv"
      )
    ),
    tables
  )

  for (table_name in tables) {
    table <- screening_results[[table_name]]

    if (is.null(table)) {
      table <- data.frame()
    }

    data.table::fwrite(
      table,
      file = output_files[[table_name]],
      sep = "\t",
      quote = TRUE,
      na = ""
    )
  }

  invisible(output_files)
}


# -----------------------------------------------------------------------------
# 12. Example
# -----------------------------------------------------------------------------

# aml_screen <- ae_screen_aml_accessions(
#   accessions = arrayexpress_ids,
#   base_dir = "D:/R/AML_Normal_Blood_ArrayExpress_Metadata",
#
#   prescreen = TRUE,
#   prescreen_require_raw_evidence = TRUE,
#
#   # Continue studies containing AML and/or normal healthy blood, PBMC,
#   # bone-marrow, cord-blood, or HSPC controls.
#   prescreen_include_normal_controls = TRUE,
#
#   # Human evidence may only occur in the SDRF. Clearly non-human-only
#   # studies are still excluded.
#   prescreen_require_human_evidence = FALSE,
#   prescreen_exclude_clear_nonhuman = TRUE,
#
#   # Do not silently lose a study if the small API pre-screen fails.
#   prescreen_on_error = "continue",
#
#   overwrite = FALSE,
#   use_api_fallback = TRUE,
#   query_biostudies_inventory = TRUE,
#   verbose = TRUE
# )
#
# Quick study decisions:
#
# View(aml_screen$prescreen_summary)
#
# Fully parsed studies:
#
# View(aml_screen$study_summary)
#
# All samples from the retained studies:
#
# View(aml_screen$sample_summary)
#
# Only human AML and normal/healthy hematopoietic samples:
#
# View(aml_screen$eligible_sample_summary)
#
# Useful sample grouping:
#
# aml_screen$eligible_sample_summary %>%
#   dplyr::count(
#     .data$Study_ID,
#     .data$disease_group,
#     .data$rnaseq_flag,
#     .data$methylation_flag
#   )
#
# RNA-seq AML and normal-control samples:
#
# aml_screen$eligible_sample_summary %>%
#   dplyr::filter(
#     .data$rnaseq_biological_candidate
#   )
#
# Methylation AML and normal-control samples:
#
# aml_screen$eligible_sample_summary %>%
#   dplyr::filter(
#     .data$methylation_biological_candidate
#   )
#
# Export all screening and parsed tables:
#
# ae_write_aml_screening_results(
#   screening_results = aml_screen,
#   out_dir = "D:/R/AML_Normal_Blood_ArrayExpress_Tables"
# )
#
# To restore AML-only pre-screening:
#
# aml_only_screen <- ae_screen_aml_accessions(
#   accessions = arrayexpress_ids,
#   prescreen = TRUE,
#   prescreen_include_normal_controls = FALSE
# )


# -----------------------------------------------------------------------------
# Title examples
# -----------------------------------------------------------------------------

# Add titles to an existing screening result without re-running the parser:
#
# aml_screen <- ae_add_titles_to_screening_result(
#   aml_screen
# )
#
# View(
#   aml_screen$study_summary
# )
#
# Retrieve a title vector for an arbitrary data frame:
#
# study_titles <- ae_get_study_titles(
#   data = my_study_table,
#   id_col = "Study_ID",
#
#   # Optional: use titles already returned by the search before querying API.
#   search_result = ae_search
# )
#
# Add the titles as a new column:
#
# my_study_table <- ae_add_study_titles(
#   data = my_study_table,
#   id_col = "Study_ID",
#   title_col = "project_title",
#   search_result = ae_search
# )

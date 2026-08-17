# =============================================================================
# ArrayExpress / BioStudies NMD perturbation screening extension
#
# PURPOSE
# -------
# Screen ArrayExpress/BioStudies studies for human cell-line RNA-seq experiments
# that perturb nonsense-mediated mRNA decay (NMD), then parse study- and
# sample-level experimental and technical metadata.
#
# This script intentionally REUSES the generic BioStudies/MAGE-TAB functions
# from:
#
#   arrayexpress_aml_normal_blood_screening_prefiltered_v1.1.1_titles_fixed.R
#
# The underlying download / IDF / SDRF functions are generic despite their
# historical ae_aml_* names. The biological screening below is NMD-specific.
#
# Recommended source order
# ------------------------
# source("arrayexpress_aml_normal_blood_screening_prefiltered_v1.1.1_titles_fixed.R")
# source("arrayexpress_nmd_perturbation_screening.R")
#
# Main functions
# --------------
# ae_nmd_quick_prescreen()
# ae_screen_one_nmd_study()
# ae_screen_nmd_accessions()
# ae_write_nmd_screening_results()
#
# Example
# -------
# nmd_screen <- ae_screen_nmd_accessions(
#   accessions = nmd_search$accessions,
#   prescreen = TRUE,
#   prescreen_require_raw_evidence = TRUE,
#   prescreen_require_perturbation_evidence = TRUE,
#   prescreen_require_cell_model_evidence = FALSE,
#   prescreen_on_error = "continue"
# )
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Confirm that the generic ArrayExpress helpers have already been sourced
# -----------------------------------------------------------------------------

ae_nmd_required_base_functions <- c(
  "ae_aml_biostudies_json",
  "ae_aml_flatten_biostudies_json",
  "ae_aml_extract_prescreen_title",
  "ae_aml_dictionary_hit_text",
  "ae_aml_has_dictionary_hit",
  "ae_aml_collapse_unique",
  "ae_aml_first_nonempty",
  "ae_aml_clean_text",
  "ae_aml_normalize_field",
  "ae_aml_download_magetab",
  "ae_aml_read_sdrf",
  "ae_aml_read_idf",
  "ae_aml_sdrf_file_inventory",
  "ae_aml_biostudies_file_inventory",
  "ae_aml_extract_project_title"
)

ae_nmd_missing_base_functions <- ae_nmd_required_base_functions[
  !vapply(
    ae_nmd_required_base_functions,
    exists,
    logical(1),
    mode = "function",
    inherits = TRUE
  )
]

if (length(ae_nmd_missing_base_functions) > 0L) {
  stop(
    paste0(
      "The generic ArrayExpress parser must be sourced before this NMD ",
      "extension.\n\nSource:\n",
      "  arrayexpress_aml_normal_blood_screening_prefiltered_v1.1.1_titles_fixed.R",
      "\n\nMissing functions:\n  ",
      paste(
        ae_nmd_missing_base_functions,
        collapse = "\n  "
      )
    )
  )
}


# -----------------------------------------------------------------------------
# 1. NMD dictionaries
# -----------------------------------------------------------------------------

AE_NMD_DIRECT_TERMS <- c(
  "Nonsense-mediated decay" =
    "nonsense[ -]+mediated[ -]+(?:mRNA[ -]+)?decay",
  "NMD" =
    "(^|[^A-Za-z0-9])NMD([^A-Za-z0-9]|$)"
)


AE_NMD_TARGET_TERMS <- c(
  "UPF1" =
    "(^|[^A-Za-z0-9])(UPF1|hUPF1|RENT1)([^A-Za-z0-9]|$)",
  "UPF2" =
    "(^|[^A-Za-z0-9])(UPF2|hUPF2|RENT2)([^A-Za-z0-9]|$)",
  "UPF3A" =
    "(^|[^A-Za-z0-9])UPF3A([^A-Za-z0-9]|$)",
  "UPF3B" =
    "(^|[^A-Za-z0-9])(UPF3B|UPF3X)([^A-Za-z0-9]|$)",
  "SMG1" =
    "(^|[^A-Za-z0-9])SMG[ -]?1([^A-Za-z0-9]|$)",
  "SMG5" =
    "(^|[^A-Za-z0-9])SMG[ -]?5([^A-Za-z0-9]|$)",
  "SMG6" =
    "(^|[^A-Za-z0-9])SMG[ -]?6([^A-Za-z0-9]|$)",
  "SMG7" =
    "(^|[^A-Za-z0-9])SMG[ -]?7([^A-Za-z0-9]|$)",
  "SMG8" =
    "(^|[^A-Za-z0-9])SMG[ -]?8([^A-Za-z0-9]|$)",
  "SMG9" =
    "(^|[^A-Za-z0-9])SMG[ -]?9([^A-Za-z0-9]|$)",
  "DHX34" =
    "(^|[^A-Za-z0-9])DHX34([^A-Za-z0-9]|$)",
  "NBAS" =
    "(^|[^A-Za-z0-9])NBAS([^A-Za-z0-9]|$)"
)


AE_NMD_GENETIC_PERTURBATION_TERMS <- c(
  "siRNA" =
    "(^|[^A-Za-z0-9])siRNA([^A-Za-z0-9]|$)|small[ -]+interfering[ -]+RNA",
  "shRNA" =
    "(^|[^A-Za-z0-9])shRNA([^A-Za-z0-9]|$)|short[ -]+hairpin[ -]+RNA",
  "RNAi" =
    "(^|[^A-Za-z0-9])RNAi([^A-Za-z0-9]|$)",
  "Knockdown" =
    "knock[ -]?down|deplet(?:e|ed|ion)|silenc(?:e|ed|ing)",
  "Knockout" =
    "knock[ -]?out|(^|[^A-Za-z0-9])KO([^A-Za-z0-9]|$)",
  "CRISPR" =
    "CRISPR|Cas9",
  "Degron" =
    "degron|auxin[ -]+inducible|(^|[^A-Za-z0-9])AID(?:[0-9]+)?([^A-Za-z0-9]|$)",
  "Overexpression" =
    "over[ -]?express(?:ion|ed|ing)?|(^|[^A-Za-z0-9])OE([^A-Za-z0-9]|$)",
  "Inducible expression" =
    "inducible|doxycycline|(^|[^A-Za-z0-9])dox([^A-Za-z0-9]|$)"
)


AE_NMD_DIRECT_INHIBITOR_TERMS <- c(
  "SMG1 inhibitor" =
    "SMG[ -]?1[ -]+inhibitor|SMG1i",
  "SMG1 inhibitor 11j" =
    "(^|[^A-Za-z0-9])11j([^A-Za-z0-9]|$)",
  "NMDI-1" =
    "(^|[^A-Za-z0-9])NMDI[ -]?1([^A-Za-z0-9]|$)",
  "VG1" =
    "(^|[^A-Za-z0-9])VG1([^A-Za-z0-9]|$)"
)


AE_NMD_TRANSLATION_INHIBITOR_TERMS <- c(
  "Emetine" =
    "(^|[^A-Za-z0-9])emetine([^A-Za-z0-9]|$)",
  "Cycloheximide" =
    "cycloheximide|(^|[^A-Za-z0-9])CHX([^A-Za-z0-9]|$)",
  "Puromycin" =
    "(^|[^A-Za-z0-9])puromycin([^A-Za-z0-9]|$)",
  "Caffeine" =
    "(^|[^A-Za-z0-9])caffeine([^A-Za-z0-9]|$)"
)


AE_NMD_STRESS_TERMS <- c(
  "Thapsigargin" =
    "(^|[^A-Za-z0-9])thapsigargin([^A-Za-z0-9]|$)",
  "Tunicamycin" =
    "(^|[^A-Za-z0-9])tunicamycin([^A-Za-z0-9]|$)",
  "Arsenite" =
    "sodium[ -]+arsenite|(^|[^A-Za-z0-9])arsenite([^A-Za-z0-9]|$)",
  "ER stress" =
    "ER[ -]+stress|endoplasmic[ -]+reticulum[ -]+stress",
  "Integrated stress response" =
    "integrated[ -]+stress[ -]+response|(^|[^A-Za-z0-9])ISR([^A-Za-z0-9]|$)"
)


AE_NMD_CONTROL_TERMS <- c(
  "Untreated control" =
    "untreated|no[ -]+treatment|mock[ -]+treated",
  "Control" =
    "(^|[^A-Za-z0-9])controls?([^A-Za-z0-9]|$)",
  "Negative control" =
    "negative[ -]+control",
  "Scramble control" =
    "scrambl(?:e|ed)|non[ -]?targeting|nontargeting",
  "siRNA control" =
    "si(?:RNA)?[ -]?(?:control|ctrl)|siCTRL|siNC",
  "shRNA control" =
    "sh(?:RNA)?[ -]?(?:control|ctrl)",
  "DMSO vehicle" =
    "(^|[^A-Za-z0-9])DMSO([^A-Za-z0-9]|$)|vehicle",
  "Wild type" =
    "wild[ -]?type|(^|[^A-Za-z0-9])WT([^A-Za-z0-9]|$)",
  "Parental" =
    "(^|[^A-Za-z0-9])parental([^A-Za-z0-9]|$)",
  "Empty vector" =
    "empty[ -]+vector|vector[ -]+control",
  "Non-induced" =
    "non[ -]?induced|uninduced|without[ -]+dox(?:ycycline)?"
)


AE_NMD_HUMAN_TERMS <- c(
  "Homo sapiens" =
    "homo[ _-]*sapiens",
  "Human" =
    "(^|[^A-Za-z])human([^A-Za-z]|$)"
)


AE_NMD_NONHUMAN_TERMS <- c(
  "Mus musculus" =
    "mus[ _-]*musculus",
  "Mouse or murine" =
    "(^|[^A-Za-z])(mouse|mice|murine)([^A-Za-z]|$)",
  "Rattus norvegicus" =
    "rattus[ _-]*norvegicus",
  "Rat" =
    "(^|[^A-Za-z])rats?([^A-Za-z]|$)"
)


AE_NMD_CELL_MODEL_TERMS <- c(
  "Cell line" =
    "cell[ -]+line|cellline",
  "Cultured cells" =
    "cultured[ -]+cells?",
  "HEK293" =
    "HEK[ -]?293(?:T|TO)?|293T|293TO|Flp[ -]?In[ -]?T[ -]?REx[ -]?293",
  "HeLa" =
    "(^|[^A-Za-z0-9])HeLa([^A-Za-z0-9]|$)",
  "K562" =
    "(^|[^A-Za-z0-9])K562([^A-Za-z0-9]|$)",
  "Huh7" =
    "(^|[^A-Za-z0-9])Huh[ -]?7([^A-Za-z0-9]|$)",
  "HT1080" =
    "(^|[^A-Za-z0-9])HT[ -]?1080([^A-Za-z0-9]|$)",
  "A498" =
    "(^|[^A-Za-z0-9])A[ -]?498([^A-Za-z0-9]|$)",
  "A704" =
    "(^|[^A-Za-z0-9])A[ -]?704([^A-Za-z0-9]|$)",
  "HK1" =
    "(^|[^A-Za-z0-9])HK[ -]?1([^A-Za-z0-9]|$)",
  "SUNE1" =
    "(^|[^A-Za-z0-9])SUNE[ -]?1([^A-Za-z0-9]|$)"
)


AE_NMD_RNASEQ_TERMS <- c(
  "RNA-seq" =
    "RNA[ -]?seq|RNASeq",
  "RNA sequencing" =
    "RNA[ -]+sequencing",
  "Transcriptome sequencing" =
    "transcriptome[ -]+sequenc",
  "Expression profiling by high-throughput sequencing" =
    "expression[ -]+profiling[ -]+by[ -]+high[ -]+throughput[ -]+sequenc",
  "mRNA sequencing" =
    "mRNA[ -]+sequenc",
  "Total RNA sequencing" =
    "total[ -]+RNA[ -]+sequenc"
)


AE_NMD_RAW_SEQUENCE_TERMS <- c(
  "FASTQ" =
    "\\.(fastq|fq)(\\.gz)?([?#].*)?$",
  "SRA file" =
    "\\.sra([?#].*)?$",
  "ENA/SRA run accession" =
    "(^|[^A-Za-z0-9])(ERR|SRR|DRR)[0-9]+([^A-Za-z0-9]|$)",
  "SRA experiment accession" =
    "(^|[^A-Za-z0-9])(ERX|SRX|DRX)[0-9]+([^A-Za-z0-9]|$)",
  "BAM" =
    "\\.bam([?#].*)?$",
  "CRAM" =
    "\\.cram([?#].*)?$"
)


AE_NMD_GENERIC_RAW_TERMS <- c(
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


# -----------------------------------------------------------------------------
# 2. Classification helpers
# -----------------------------------------------------------------------------

ae_nmd_hits <- function(text, dictionary) {
  ae_aml_dictionary_hit_text(
    text,
    dictionary
  )
}


ae_nmd_has <- function(text, dictionary) {
  ae_aml_has_dictionary_hit(
    text,
    dictionary
  )
}


ae_nmd_classify_perturbation <- function(text) {

  text <- ae_aml_clean_text(text)

  if (is.na(text) || !nzchar(text)) {
    return("Unclear")
  }

  if (ae_nmd_has(text, AE_NMD_CONTROL_TERMS)) {
    return("Control")
  }

  if (
    grepl(
      "degron|auxin[ -]+inducible|(^|[^A-Za-z0-9])AID(?:[0-9]+)?([^A-Za-z0-9]|$)",
      text,
      ignore.case = TRUE,
      perl = TRUE
    )
  ) {
    return("Inducible degradation")
  }

  if (
    grepl(
      "knock[ -]?out|CRISPR|Cas9|(^|[^A-Za-z0-9])KO([^A-Za-z0-9]|$)",
      text,
      ignore.case = TRUE,
      perl = TRUE
    )
  ) {
    return("Knockout / CRISPR")
  }

  if (
    grepl(
      "(^|[^A-Za-z0-9])siRNA([^A-Za-z0-9]|$)|small[ -]+interfering[ -]+RNA",
      text,
      ignore.case = TRUE,
      perl = TRUE
    )
  ) {
    return("siRNA knockdown")
  }

  if (
    grepl(
      "(^|[^A-Za-z0-9])shRNA([^A-Za-z0-9]|$)|short[ -]+hairpin[ -]+RNA",
      text,
      ignore.case = TRUE,
      perl = TRUE
    )
  ) {
    return("shRNA knockdown")
  }

  if (
    grepl(
      "knock[ -]?down|deplet(?:e|ed|ion)|silenc(?:e|ed|ing)|(^|[^A-Za-z0-9])RNAi([^A-Za-z0-9]|$)",
      text,
      ignore.case = TRUE,
      perl = TRUE
    )
  ) {
    return("Other knockdown / depletion")
  }

  if (ae_nmd_has(text, AE_NMD_DIRECT_INHIBITOR_TERMS)) {
    return("Direct NMD / SMG1 inhibitor")
  }

  if (ae_nmd_has(text, AE_NMD_TRANSLATION_INHIBITOR_TERMS)) {
    return("Translation inhibitor")
  }

  if (ae_nmd_has(text, AE_NMD_STRESS_TERMS)) {
    return("Stress-mediated NMD modulation")
  }

  if (
    grepl(
      "over[ -]?express(?:ion|ed|ing)?|(^|[^A-Za-z0-9])OE([^A-Za-z0-9]|$)",
      text,
      ignore.case = TRUE,
      perl = TRUE
    )
  ) {
    return("NMD-factor overexpression")
  }

  if (
    ae_nmd_has(text, AE_NMD_TARGET_TERMS) ||
      ae_nmd_has(text, AE_NMD_DIRECT_TERMS)
  ) {
    return("NMD perturbation; method unclear")
  }

  "Other / unclear"
}


ae_nmd_expected_effect <- function(
    perturbation_type,
    text = NA_character_) {

  if (is.na(perturbation_type)) {
    return("Unclear")
  }

  if (perturbation_type == "Control") {
    return("Control / baseline NMD")
  }

  if (
    perturbation_type %in% c(
      "Inducible degradation",
      "Knockout / CRISPR",
      "siRNA knockdown",
      "shRNA knockdown",
      "Other knockdown / depletion",
      "Direct NMD / SMG1 inhibitor"
    )
  ) {
    return("Likely NMD inhibition")
  }

  if (
    perturbation_type ==
      "Translation inhibitor"
  ) {
    return("Indirect NMD inhibition via translation inhibition")
  }

  if (
    perturbation_type ==
      "Stress-mediated NMD modulation"
  ) {
    return("Potential stress-mediated NMD inhibition")
  }

  if (
    perturbation_type ==
      "NMD-factor overexpression"
  ) {
    return("NMD-factor overexpression / altered NMD activity")
  }

  "NMD effect unclear"
}


ae_nmd_control_type <- function(text) {

  hits <- ae_nmd_hits(
    text,
    AE_NMD_CONTROL_TERMS
  )

  if (is.na(hits)) {
    return(NA_character_)
  }

  hits
}


# -----------------------------------------------------------------------------
# 3. Fast BioStudies pre-screen
# -----------------------------------------------------------------------------

ae_nmd_quick_prescreen <- function(
    accession,
    biostudies_json = NULL,
    require_raw_evidence = TRUE,
    require_perturbation_evidence = TRUE,
    require_human_evidence = FALSE,
    exclude_clear_nonhuman = TRUE,
    require_cell_model_evidence = FALSE,
    include_stress_perturbations = FALSE,
    verbose = TRUE) {

  accession <- toupper(
    trimws(accession)
  )

  if (is.null(biostudies_json)) {
    if (verbose) {
      message(
        "Quick NMD BioStudies pre-screen: ",
        accession
      )
    }

    biostudies_json <-
      ae_aml_biostudies_json(
        accession
      )
  }

  flat_json <-
    ae_aml_flatten_biostudies_json(
      biostudies_json =
        biostudies_json,
      accession = accession
    )

  # Avoid using publication/citation titles as biological evidence.
  screening_rows <-
    !stringr::str_detect(
      flat_json$field,
      stringr::regex(
        "reference|publication|pubmed|citation|bibliograph",
        ignore_case = TRUE
      )
    )

  screening_text <-
    ae_aml_collapse_unique(
      paste(
        flat_json$field[screening_rows],
        flat_json$value[screening_rows],
        sep = ": "
      ),
      max_values = 2500L
    )

  all_text <-
    ae_aml_collapse_unique(
      paste(
        flat_json$field,
        flat_json$value,
        sep = ": "
      ),
      max_values = 3500L
    )

  direct_nmd_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_DIRECT_TERMS
    )

  target_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_TARGET_TERMS
    )

  genetic_perturbation_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_GENETIC_PERTURBATION_TERMS
    )

  direct_inhibitor_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_DIRECT_INHIBITOR_TERMS
    )

  translation_inhibitor_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_TRANSLATION_INHIBITOR_TERMS
    )

  stress_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_STRESS_TERMS
    )

  cell_model_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_CELL_MODEL_TERMS
    )

  human_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_HUMAN_TERMS
    )

  nonhuman_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_NONHUMAN_TERMS
    )

  rnaseq_terms <-
    ae_nmd_hits(
      screening_text,
      AE_NMD_RNASEQ_TERMS
    )

  raw_sequence_terms <-
    ae_nmd_hits(
      all_text,
      AE_NMD_RAW_SEQUENCE_TERMS
    )

  generic_raw_terms <-
    ae_nmd_hits(
      all_text,
      AE_NMD_GENERIC_RAW_TERMS
    )

  nmd_flag <-
    !is.na(direct_nmd_terms) ||
    !is.na(target_terms) ||
    !is.na(direct_inhibitor_terms)

  perturbation_flag <-
    !is.na(genetic_perturbation_terms) ||
    !is.na(direct_inhibitor_terms) ||
    !is.na(translation_inhibitor_terms) ||
    (
      isTRUE(include_stress_perturbations) &&
        !is.na(stress_terms)
    )

  human_flag <- !is.na(human_terms)
  nonhuman_flag <- !is.na(nonhuman_terms)
  cell_model_flag <- !is.na(cell_model_terms)
  rnaseq_flag <- !is.na(rnaseq_terms)

  raw_sequence_flag <-
    !is.na(raw_sequence_terms) ||
    !is.na(generic_raw_terms)

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

  cell_model_ok <-
    !isTRUE(require_cell_model_evidence) ||
    cell_model_flag

  perturbation_ok <-
    !isTRUE(require_perturbation_evidence) ||
    perturbation_flag

  raw_ok <-
    !isTRUE(require_raw_evidence) ||
    raw_sequence_flag

  continue_full_parse <-
    nmd_flag &&
    perturbation_ok &&
    human_ok &&
    cell_model_ok &&
    rnaseq_flag &&
    raw_ok

  prescreen_status <- dplyr::case_when(
    !nmd_flag ~
      "Skip: NMD pathway or NMD factor not identified",

    !perturbation_ok ~
      "Skip: NMD perturbation evidence not identified",

    !human_ok ~
      "Skip: clearly non-human or human evidence required",

    !cell_model_ok ~
      "Skip: cell-line/model evidence not identified",

    !rnaseq_flag ~
      "Skip: RNA-seq not identified",

    !raw_ok ~
      "Skip: raw RNA-seq evidence not identified",

    TRUE ~
      "Continue: NMD perturbation RNA-seq candidate"
  )

  summary <- tibble::tibble(
    Study_ID = accession,
    project_title =
      ae_aml_extract_prescreen_title(
        flat_json
      ),
    prescreen_status =
      prescreen_status,
    continue_full_parse =
      continue_full_parse,

    nmd_flag = nmd_flag,
    perturbation_flag =
      perturbation_flag,
    human_flag = human_flag,
    nonhuman_flag =
      nonhuman_flag,
    cell_model_flag =
      cell_model_flag,
    rnaseq_flag =
      rnaseq_flag,
    raw_sequence_flag =
      raw_sequence_flag,

    direct_nmd_terms =
      direct_nmd_terms,
    nmd_target_terms =
      target_terms,
    genetic_perturbation_terms =
      genetic_perturbation_terms,
    direct_inhibitor_terms =
      direct_inhibitor_terms,
    translation_inhibitor_terms =
      translation_inhibitor_terms,
    stress_terms =
      stress_terms,
    cell_model_terms =
      cell_model_terms,
    human_terms =
      human_terms,
    nonhuman_terms =
      nonhuman_terms,
    rnaseq_terms =
      rnaseq_terms,
    raw_sequence_terms =
      raw_sequence_terms,
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
      evidence_class =
        dplyr::case_when(
          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_DIRECT_TERMS
          ) ~ "NMD terminology",

          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_TARGET_TERMS
          ) ~ "NMD factor",

          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_GENETIC_PERTURBATION_TERMS
          ) ~ "Genetic perturbation",

          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_DIRECT_INHIBITOR_TERMS
          ) ~ "Direct NMD inhibitor",

          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_TRANSLATION_INHIBITOR_TERMS
          ) ~ "Translation inhibitor",

          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_STRESS_TERMS
          ) ~ "Stress perturbation",

          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_CELL_MODEL_TERMS
          ) ~ "Cell model",

          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_RNASEQ_TERMS
          ) ~ "RNA-seq",

          vapply(
            .data$row_text,
            ae_nmd_has,
            logical(1),
            dictionary =
              AE_NMD_RAW_SEQUENCE_TERMS
          ) ~ "Raw sequencing evidence",

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


# -----------------------------------------------------------------------------
# 4. Sample-level NMD metadata parser
# -----------------------------------------------------------------------------

ae_nmd_collapse_where <- function(
    values,
    keep,
    max_values = 30L) {

  ae_aml_collapse_unique(
    values[keep],
    max_values = max_values
  )
}


ae_nmd_make_sample_summary <- function(
    accession,
    sdrf_long,
    idf_long,
    prescreen_summary = NULL) {

  if (nrow(sdrf_long) == 0L) {
    return(tibble::tibble())
  }

  study_context <-
    ae_aml_collapse_unique(
      paste(
        idf_long$field_raw,
        idf_long$value,
        sep = ": "
      ),
      max_values = 300L
    )

  study_nmd_terms <-
    ae_nmd_hits(
      study_context,
      AE_NMD_DIRECT_TERMS
    )

  study_target_terms <-
    ae_nmd_hits(
      study_context,
      AE_NMD_TARGET_TERMS
    )

  study_rnaseq_terms <-
    ae_nmd_hits(
      study_context,
      AE_NMD_RNASEQ_TERMS
    )

  tagged <- sdrf_long %>%
    dplyr::mutate(
      field =
        ae_aml_normalize_field(
          .data$field_raw
        ),

      is_organism_field =
        stringr::str_detect(
          .data$field,
          "organism|species"
        ),

      is_source_field =
        stringr::str_detect(
          .data$field,
          "source_name|cell_line|cell_type|model|biosource|material|organism_part|tissue"
        ),

      is_characteristic_field =
        stringr::str_detect(
          .data$field,
          "characteristic|factor_value|genotype|phenotype|disease|cell"
        ),

      is_perturbation_field =
        stringr::str_detect(
          .data$field,
          "factor_value|treatment|compound|drug|perturb|transfect|transduc|RNAi|siRNA|shRNA|CRISPR|knock|deplet|genotype|construct|vector|induc|dox|auxin|protocol"
        ),

      is_control_field =
        stringr::str_detect(
          .data$field,
          "factor_value|treatment|control|group|condition|genotype|construct|vector"
        ),

      is_assay_field =
        stringr::str_detect(
          .data$field,
          "assay|technology|sequenc|library|rna"
        ),

      is_library_strategy_field =
        stringr::str_detect(
          .data$field,
          "library_strategy|sequencing_strategy"
        ),

      is_library_source_field =
        stringr::str_detect(
          .data$field,
          "library_source"
        ),

      is_library_selection_field =
        stringr::str_detect(
          .data$field,
          "library_selection"
        ),

      is_library_layout_field =
        stringr::str_detect(
          .data$field,
          "library_layout|paired|single_end|singleend"
        ),

      is_library_prep_field =
        stringr::str_detect(
          .data$field,
          "library_prep|library_preparation|library_construction|protocol_ref"
        ),

      is_instrument_field =
        stringr::str_detect(
          .data$field,
          "instrument|sequencer|sequencing_platform|platform"
        ),

      is_stranded_field =
        stringr::str_detect(
          .data$field,
          "stranded|strand_specific|strand_specificity"
        ),

      is_read_length_field =
        stringr::str_detect(
          .data$field,
          "read_length|readlength|cycle"
        ),

      is_file_field =
        stringr::str_detect(
          .data$field,
          "file|uri|url|fastq|fq|bam|cram|ena|sra|run_accession|raw"
        )
    )

  sample_summary <- tagged %>%
    dplyr::group_by(
      .data$Study_ID,
      .data$Sample_ID
    ) %>%
    dplyr::summarise(
      all_sample_metadata =
        ae_aml_collapse_unique(
          paste(
            .data$field_raw,
            .data$value,
            sep = ": "
          ),
          max_values = 250L
        ),

      organism_values =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_organism_field
        ),

      source_model_values =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_source_field
        ),

      characteristic_values =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_characteristic_field,
          max_values = 60L
        ),

      perturbation_values =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_perturbation_field,
          max_values = 80L
        ),

      control_condition_values =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_control_field,
          max_values = 60L
        ),

      assay_values =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_assay_field
        ),

      library_strategy =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_library_strategy_field
        ),

      library_source =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_library_source_field
        ),

      library_selection =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_library_selection_field
        ),

      library_layout =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_library_layout_field
        ),

      library_preparation =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_library_prep_field,
          max_values = 40L
        ),

      sequencing_instrument =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_instrument_field
        ),

      strandedness =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_stranded_field
        ),

      read_length =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_read_length_field
        ),

      file_values =
        ae_nmd_collapse_where(
          .data$value,
          .data$is_file_field,
          max_values = 80L
        ),

      .groups = "drop"
    ) %>%
    dplyr::mutate(
      study_context =
        study_context,

      sample_search_text =
        paste(
          dplyr::coalesce(
            .data$all_sample_metadata,
            ""
          ),
          dplyr::coalesce(
            .data$source_model_values,
            ""
          ),
          dplyr::coalesce(
            .data$characteristic_values,
            ""
          ),
          dplyr::coalesce(
            .data$perturbation_values,
            ""
          )
        ),

      perturbation_search_text =
        paste(
          dplyr::coalesce(
            .data$perturbation_values,
            ""
          ),
          dplyr::coalesce(
            .data$control_condition_values,
            ""
          ),
          dplyr::coalesce(
            .data$characteristic_values,
            ""
          )
        ),

      assay_search_text =
        paste(
          dplyr::coalesce(
            .data$assay_values,
            ""
          ),
          dplyr::coalesce(
            .data$library_strategy,
            ""
          ),
          dplyr::coalesce(
            .data$file_values,
            ""
          ),
          dplyr::coalesce(
            .data$study_context,
            ""
          )
        ),

      human_terms =
        vapply(
          .data$sample_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_HUMAN_TERMS
        ),

      nonhuman_terms =
        vapply(
          .data$sample_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_NONHUMAN_TERMS
        ),

      cell_model_terms =
        vapply(
          .data$sample_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_CELL_MODEL_TERMS
        ),

      nmd_direct_terms =
        vapply(
          .data$sample_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_DIRECT_TERMS
        ),

      nmd_target_terms =
        vapply(
          .data$sample_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_TARGET_TERMS
        ),

      genetic_perturbation_terms =
        vapply(
          .data$perturbation_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_GENETIC_PERTURBATION_TERMS
        ),

      direct_inhibitor_terms =
        vapply(
          .data$perturbation_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_DIRECT_INHIBITOR_TERMS
        ),

      translation_inhibitor_terms =
        vapply(
          .data$perturbation_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_TRANSLATION_INHIBITOR_TERMS
        ),

      stress_terms =
        vapply(
          .data$perturbation_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_STRESS_TERMS
        ),

      control_terms =
        vapply(
          .data$perturbation_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_CONTROL_TERMS
        ),

      rnaseq_sample_terms =
        vapply(
          .data$assay_search_text,
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_RNASEQ_TERMS
        ),

      raw_sequence_terms =
        vapply(
          dplyr::coalesce(
            .data$file_values,
            ""
          ),
          ae_nmd_hits,
          character(1),
          dictionary =
            AE_NMD_RAW_SEQUENCE_TERMS
        ),

      human_flag =
        !is.na(.data$human_terms),

      nonhuman_flag =
        !is.na(.data$nonhuman_terms),

      cell_model_flag =
        !is.na(.data$cell_model_terms) |
        (
          !is.na(.data$source_model_values) &
            nzchar(.data$source_model_values)
        ),

      rnaseq_flag =
        !is.na(.data$rnaseq_sample_terms) |
        !is.na(study_rnaseq_terms),

      raw_sequence_flag =
        !is.na(.data$raw_sequence_terms),

      perturbation_type =
        vapply(
          .data$perturbation_search_text,
          ae_nmd_classify_perturbation,
          character(1)
        ),

      control_type =
        vapply(
          .data$perturbation_search_text,
          ae_nmd_control_type,
          character(1)
        ),

      control_flag =
        .data$perturbation_type ==
          "Control" |
        !is.na(.data$control_terms),

      direct_nmd_evidence_flag =
        !is.na(.data$nmd_direct_terms) |
        !is.na(.data$nmd_target_terms) |
        !is.na(.data$direct_inhibitor_terms),

      nmd_perturbation_flag =
        !.data$control_flag &
        (
          !is.na(.data$nmd_target_terms) |
          !is.na(.data$direct_inhibitor_terms) |
          (
            !is.na(.data$nmd_direct_terms) &
              .data$perturbation_type !=
                "Other / unclear"
          )
        ),

      expected_nmd_effect =
        mapply(
          ae_nmd_expected_effect,
          perturbation_type =
            .data$perturbation_type,
          text =
            .data$perturbation_search_text,
          USE.NAMES = FALSE
        )
    )

  # Controls in an NMD study are useful even though they do not contain the
  # NMD-factor name themselves.
  study_has_perturbation <- any(
    sample_summary$nmd_perturbation_flag,
    na.rm = TRUE
  )

  study_human_from_prescreen <- FALSE

  if (
    !is.null(prescreen_summary) &&
      nrow(prescreen_summary) > 0L &&
      "human_flag" %in%
        names(prescreen_summary)
  ) {
    study_human_from_prescreen <-
      isTRUE(
        prescreen_summary$human_flag[1]
      )
  }

  sample_summary <- sample_summary %>%
    dplyr::mutate(
      human_flag =
        .data$human_flag |
        study_human_from_prescreen,

      study_has_nmd_perturbation =
        study_has_perturbation,

      nmd_control_sample_flag =
        .data$control_flag &
        study_has_perturbation,

      eligible_nmd_sample_flag =
        .data$human_flag &
        .data$cell_model_flag &
        .data$rnaseq_flag &
        (
          .data$nmd_perturbation_flag |
          .data$nmd_control_sample_flag
        ),

      sample_role =
        dplyr::case_when(
          .data$nmd_perturbation_flag ~
            "NMD perturbation",

          .data$nmd_control_sample_flag ~
            "Control",

          TRUE ~
            "Other / unclear"
        ),

      nmd_target_or_study_target =
        dplyr::coalesce(
          .data$nmd_target_terms,
          study_target_terms
        )
    ) %>%
    dplyr::select(
      .data$Study_ID,
      .data$Sample_ID,
      .data$sample_role,

      .data$nmd_target_or_study_target,
      .data$nmd_target_terms,
      .data$perturbation_type,
      .data$expected_nmd_effect,
      .data$control_type,

      .data$perturbation_values,
      .data$control_condition_values,

      .data$organism_values,
      .data$source_model_values,
      .data$characteristic_values,

      .data$human_flag,
      .data$nonhuman_flag,
      .data$cell_model_flag,

      .data$nmd_perturbation_flag,
      .data$nmd_control_sample_flag,
      .data$eligible_nmd_sample_flag,

      .data$rnaseq_flag,
      .data$raw_sequence_flag,

      .data$assay_values,
      .data$library_strategy,
      .data$library_source,
      .data$library_selection,
      .data$library_layout,
      .data$library_preparation,
      .data$sequencing_instrument,
      .data$strandedness,
      .data$read_length,
      .data$file_values,

      .data$nmd_direct_terms,
      .data$genetic_perturbation_terms,
      .data$direct_inhibitor_terms,
      .data$translation_inhibitor_terms,
      .data$stress_terms,
      .data$control_terms
    )

  sample_summary
}


# -----------------------------------------------------------------------------
# 5. NMD evidence table
# -----------------------------------------------------------------------------

ae_nmd_make_evidence <- function(
    idf_long,
    sdrf_long,
    file_inventory) {

  tag_rows <- function(
      table,
      source,
      dictionary,
      evidence_class,
      sample_column = NULL) {

    if (nrow(table) == 0L) {
      return(tibble::tibble())
    }

    row_text <- paste(
      table$field_raw,
      table$value,
      sep = ": "
    )

    matched <- vapply(
      row_text,
      ae_nmd_hits,
      character(1),
      dictionary = dictionary
    )

    keep <- !is.na(matched)

    if (!any(keep)) {
      return(tibble::tibble())
    }

    sample_id <- rep(
      NA_character_,
      sum(keep)
    )

    if (
      !is.null(sample_column) &&
        sample_column %in% names(table)
    ) {
      sample_id <-
        as.character(
          table[[sample_column]][keep]
        )
    }

    tibble::tibble(
      Study_ID =
        table$Study_ID[keep],
      Sample_ID =
        sample_id,
      source = source,
      evidence_class =
        evidence_class,
      field =
        table$field_raw[keep],
      value =
        table$value[keep],
      matched_terms =
        matched[keep]
    )
  }

  evidence <- dplyr::bind_rows(
    tag_rows(
      idf_long,
      "IDF",
      AE_NMD_DIRECT_TERMS,
      "NMD terminology"
    ),
    tag_rows(
      sdrf_long,
      "SDRF",
      AE_NMD_DIRECT_TERMS,
      "NMD terminology",
      "Sample_ID"
    ),
    tag_rows(
      idf_long,
      "IDF",
      AE_NMD_TARGET_TERMS,
      "NMD factor"
    ),
    tag_rows(
      sdrf_long,
      "SDRF",
      AE_NMD_TARGET_TERMS,
      "NMD factor",
      "Sample_ID"
    ),
    tag_rows(
      idf_long,
      "IDF",
      AE_NMD_GENETIC_PERTURBATION_TERMS,
      "Genetic perturbation"
    ),
    tag_rows(
      sdrf_long,
      "SDRF",
      AE_NMD_GENETIC_PERTURBATION_TERMS,
      "Genetic perturbation",
      "Sample_ID"
    ),
    tag_rows(
      idf_long,
      "IDF",
      AE_NMD_DIRECT_INHIBITOR_TERMS,
      "Direct NMD inhibitor"
    ),
    tag_rows(
      sdrf_long,
      "SDRF",
      AE_NMD_DIRECT_INHIBITOR_TERMS,
      "Direct NMD inhibitor",
      "Sample_ID"
    ),
    tag_rows(
      sdrf_long,
      "SDRF",
      AE_NMD_CONTROL_TERMS,
      "Control condition",
      "Sample_ID"
    ),
    tag_rows(
      sdrf_long,
      "SDRF",
      AE_NMD_CELL_MODEL_TERMS,
      "Cell model",
      "Sample_ID"
    ),
    tag_rows(
      idf_long,
      "IDF",
      AE_NMD_RNASEQ_TERMS,
      "RNA-seq"
    ),
    tag_rows(
      sdrf_long,
      "SDRF",
      AE_NMD_RNASEQ_TERMS,
      "RNA-seq",
      "Sample_ID"
    )
  ) %>%
    dplyr::distinct()

  if (
    !is.null(file_inventory) &&
      nrow(file_inventory) > 0L
  ) {
    file_columns <- intersect(
      c(
        "Study_ID",
        "Sample_ID",
        "file_value",
        "file_name",
        "file_url",
        "file_class",
        "source"
      ),
      names(file_inventory)
    )

    file_evidence <-
      file_inventory %>%
      dplyr::select(
        dplyr::all_of(
          file_columns
        )
      )

    if (nrow(file_evidence) > 0L) {
      file_evidence$evidence_class <-
        "File evidence"

      evidence <- dplyr::bind_rows(
        evidence,
        file_evidence
      )
    }
  }

  evidence
}


# -----------------------------------------------------------------------------
# 6. Study-level NMD summary
# -----------------------------------------------------------------------------

ae_nmd_make_study_summary <- function(
    accession,
    project_title,
    sample_summary,
    idf_long,
    file_inventory,
    download_method,
    prescreen_summary = NULL) {

  study_text <-
    ae_aml_collapse_unique(
      paste(
        idf_long$field_raw,
        idf_long$value,
        sep = ": "
      ),
      max_values = 500L
    )

  n_samples <-
    dplyr::n_distinct(
      sample_summary$Sample_ID
    )

  n_nmd_perturbation_samples <-
    sum(
      sample_summary$nmd_perturbation_flag,
      na.rm = TRUE
    )

  n_control_samples <-
    sum(
      sample_summary$nmd_control_sample_flag,
      na.rm = TRUE
    )

  n_eligible_samples <-
    sum(
      sample_summary$eligible_nmd_sample_flag,
      na.rm = TRUE
    )

  n_raw_sample_links <-
    sum(
      sample_summary$raw_sequence_flag,
      na.rm = TRUE
    )

  file_text <- NA_character_

  if (
    !is.null(file_inventory) &&
      nrow(file_inventory) > 0L
  ) {
    file_text <-
      ae_aml_collapse_unique(
        unlist(
          file_inventory,
          recursive = TRUE,
          use.names = FALSE
        ),
        max_values = 1000L
      )
  }

  study_raw_flag <-
    ae_nmd_has(
      file_text,
      AE_NMD_RAW_SEQUENCE_TERMS
    ) ||
    n_raw_sample_links > 0L

  study_targets <-
    ae_aml_collapse_unique(
      c(
        sample_summary$nmd_target_or_study_target,
        ae_nmd_hits(
          study_text,
          AE_NMD_TARGET_TERMS
        )
      )
    )

  perturbation_methods <-
    ae_aml_collapse_unique(
      sample_summary$perturbation_type[
        sample_summary$nmd_perturbation_flag
      ]
    )

  expected_nmd_effects <-
    ae_aml_collapse_unique(
      sample_summary$expected_nmd_effect[
        sample_summary$nmd_perturbation_flag
      ]
    )

  cell_models <-
    ae_aml_collapse_unique(
      sample_summary$source_model_values[
        sample_summary$cell_model_flag
      ],
      max_values = 40L
    )

  control_types <-
    ae_aml_collapse_unique(
      sample_summary$control_type[
        sample_summary$nmd_control_sample_flag
      ]
    )

  library_strategy <-
    ae_aml_collapse_unique(
      sample_summary$library_strategy
    )

  library_source <-
    ae_aml_collapse_unique(
      sample_summary$library_source
    )

  library_selection <-
    ae_aml_collapse_unique(
      sample_summary$library_selection
    )

  library_layout <-
    ae_aml_collapse_unique(
      sample_summary$library_layout
    )

  sequencing_instrument <-
    ae_aml_collapse_unique(
      sample_summary$sequencing_instrument
    )

  strandedness <-
    ae_aml_collapse_unique(
      sample_summary$strandedness
    )

  screening_tier <-
    dplyr::case_when(
      n_nmd_perturbation_samples > 0L &&
        n_control_samples > 0L &&
        study_raw_flag ~
        "Strong NMD RNA-seq candidate",

      n_nmd_perturbation_samples > 0L &&
        study_raw_flag ~
        "NMD perturbation with raw RNA-seq; control needs review",

      n_nmd_perturbation_samples > 0L ~
        "NMD perturbation identified; raw availability needs review",

      TRUE ~
        "NMD perturbation not confirmed after full parsing"
    )

  tibble::tibble(
    Study_ID = accession,
    project_title =
      project_title,
    screening_tier =
      screening_tier,
    download_method =
      download_method,

    nmd_targets =
      study_targets,
    perturbation_methods =
      perturbation_methods,
    expected_nmd_effects =
      expected_nmd_effects,
    cell_models =
      cell_models,
    control_types =
      control_types,

    n_samples =
      n_samples,
    n_nmd_perturbation_samples =
      n_nmd_perturbation_samples,
    n_control_samples =
      n_control_samples,
    n_eligible_nmd_samples =
      n_eligible_samples,

    raw_rnaseq_available =
      study_raw_flag,

    library_strategy =
      library_strategy,
    library_source =
      library_source,
    library_selection =
      library_selection,
    library_layout =
      library_layout,
    sequencing_instrument =
      sequencing_instrument,
    strandedness =
      strandedness,

    direct_nmd_terms =
      ae_nmd_hits(
        study_text,
        AE_NMD_DIRECT_TERMS
      ),

    idf_nmd_target_terms =
      ae_nmd_hits(
        study_text,
        AE_NMD_TARGET_TERMS
      ),

    idf_genetic_perturbation_terms =
      ae_nmd_hits(
        study_text,
        AE_NMD_GENETIC_PERTURBATION_TERMS
      ),

    idf_direct_inhibitor_terms =
      ae_nmd_hits(
        study_text,
        AE_NMD_DIRECT_INHIBITOR_TERMS
      ),

    idf_translation_inhibitor_terms =
      ae_nmd_hits(
        study_text,
        AE_NMD_TRANSLATION_INHIBITOR_TERMS
      ),

    idf_stress_terms =
      ae_nmd_hits(
        study_text,
        AE_NMD_STRESS_TERMS
      )
  )
}


# -----------------------------------------------------------------------------
# 7. Parse one study
# -----------------------------------------------------------------------------

ae_screen_one_nmd_study <- function(
    accession,
    base_dir =
      "ArrayExpress_NMD_Screen",
    overwrite = FALSE,
    use_api_fallback = TRUE,
    query_biostudies_inventory = TRUE,
    verbose = TRUE,
    biostudies_json = NULL,
    prescreen_summary = NULL) {

  accession <- toupper(
    trimws(accession)
  )

  downloaded <-
    ae_aml_download_magetab(
      accession = accession,
      base_dir = base_dir,
      overwrite = overwrite,
      use_api_fallback =
        use_api_fallback,
      verbose = verbose,
      biostudies_json =
        biostudies_json
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

  sample_summary <-
    ae_nmd_make_sample_summary(
      accession = accession,
      sdrf_long = sdrf_long,
      idf_long = idf_long,
      prescreen_summary =
        prescreen_summary
    )

  sdrf_inventory <-
    ae_aml_sdrf_file_inventory(
      sdrf_long
    )

  api_inventory <-
    tibble::tibble()

  if (isTRUE(
    query_biostudies_inventory
  )) {
    api_inventory <-
      ae_aml_biostudies_file_inventory(
        accession = accession,
        verbose = verbose,
        biostudies_json =
          biostudies_json
      )
  }

  file_inventory <-
    dplyr::bind_rows(
      sdrf_inventory,
      api_inventory
    ) %>%
    dplyr::distinct()

  evidence <-
    ae_nmd_make_evidence(
      idf_long = idf_long,
      sdrf_long = sdrf_long,
      file_inventory =
        file_inventory
    )

  project_title <-
    ae_aml_extract_project_title(
      idf_long = idf_long,
      prescreen_summary =
        prescreen_summary
    )

  study_summary <-
    ae_nmd_make_study_summary(
      accession = accession,
      project_title =
        project_title,
      sample_summary =
        sample_summary,
      idf_long = idf_long,
      file_inventory =
        file_inventory,
      download_method =
        downloaded$download_method,
      prescreen_summary =
        prescreen_summary
    )

  eligible_sample_summary <-
    sample_summary %>%
    dplyr::filter(
      .data$eligible_nmd_sample_flag
    )

  perturbation_sample_summary <-
    sample_summary %>%
    dplyr::filter(
      .data$nmd_perturbation_flag
    )

  control_sample_summary <-
    sample_summary %>%
    dplyr::filter(
      .data$nmd_control_sample_flag
    )

  result <- list(
    Study_ID = accession,
    study_summary =
      study_summary,
    sample_summary =
      sample_summary,
    eligible_sample_summary =
      eligible_sample_summary,
    perturbation_sample_summary =
      perturbation_sample_summary,
    control_sample_summary =
      control_sample_summary,
    evidence = evidence,
    file_inventory =
      file_inventory,
    sdrf_wide =
      sdrf_wide,
    sdrf_long =
      sdrf_long,
    idf_long =
      idf_long,
    prescreen_summary =
      prescreen_summary
  )

  saveRDS(
    result,
    file = file.path(
      study_dir,
      paste0(
        accession,
        "_NMD_screening.rds"
      )
    )
  )

  result
}


# -----------------------------------------------------------------------------
# 8. Batch NMD screen
# -----------------------------------------------------------------------------

ae_screen_nmd_accessions <- function(
    accessions,
    base_dir =
      "ArrayExpress_NMD_Screen",
    overwrite = FALSE,
    use_api_fallback = TRUE,
    query_biostudies_inventory = TRUE,
    verbose = TRUE,
    prescreen = TRUE,
    prescreen_require_raw_evidence = TRUE,
    prescreen_require_perturbation_evidence = TRUE,
    prescreen_require_human_evidence = FALSE,
    prescreen_exclude_clear_nonhuman = TRUE,
    prescreen_require_cell_model_evidence = FALSE,
    prescreen_include_stress_perturbations = FALSE,
    prescreen_on_error =
      c("continue", "skip")) {

  prescreen_on_error <-
    match.arg(
      prescreen_on_error
    )

  accessions <- unique(
    toupper(
      trimws(
        as.character(
          accessions
        )
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

  prescreen_summary <-
    tibble::tibble()

  prescreen_evidence <-
    tibble::tibble()

  skipped_prescreen <-
    tibble::tibble()

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
    this_prescreen_summary <-
      NULL

    if (isTRUE(prescreen)) {

      quick <- tryCatch(
        ae_nmd_quick_prescreen(
          accession = accession,
          require_raw_evidence =
            prescreen_require_raw_evidence,
          require_perturbation_evidence =
            prescreen_require_perturbation_evidence,
          require_human_evidence =
            prescreen_require_human_evidence,
          exclude_clear_nonhuman =
            prescreen_exclude_clear_nonhuman,
          require_cell_model_evidence =
            prescreen_require_cell_model_evidence,
          include_stress_perturbations =
            prescreen_include_stress_perturbations,
          verbose = verbose
        ),
        error = function(e) e
      )

      if (inherits(
        quick,
        "error"
      )) {

        fallback_continue <-
          identical(
            prescreen_on_error,
            "continue"
          )

        errors <-
          dplyr::bind_rows(
            errors,
            tibble::tibble(
              Study_ID = accession,
              stage =
                "BioStudies NMD pre-screen",
              error =
                conditionMessage(quick)
            )
          )

        this_prescreen_summary <-
          tibble::tibble(
            Study_ID = accession,
            project_title =
              NA_character_,
            prescreen_status =
              ifelse(
                fallback_continue,
                "Continue: NMD pre-screen failed; full parser used as fallback",
                "Skip: NMD pre-screen failed"
              ),
            continue_full_parse =
              fallback_continue,
            nmd_flag = NA,
            perturbation_flag = NA,
            human_flag = NA,
            nonhuman_flag = NA,
            cell_model_flag = NA,
            rnaseq_flag = NA,
            raw_sequence_flag = NA
          )

        prescreen_summary <-
          dplyr::bind_rows(
            prescreen_summary,
            this_prescreen_summary
          )

        if (!fallback_continue) {
          skipped_prescreen <-
            dplyr::bind_rows(
              skipped_prescreen,
              this_prescreen_summary
            )
          next
        }

      } else {

        this_prescreen_summary <-
          quick$summary

        preloaded_json <-
          quick$json

        prescreen_summary <-
          dplyr::bind_rows(
            prescreen_summary,
            quick$summary
          )

        prescreen_evidence <-
          dplyr::bind_rows(
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

          skipped_prescreen <-
            dplyr::bind_rows(
              skipped_prescreen,
              quick$summary
            )

          next
        }
      }
    }

    result <- tryCatch(
      ae_screen_one_nmd_study(
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

    if (inherits(
      result,
      "error"
    )) {

      errors <-
        dplyr::bind_rows(
          errors,
          tibble::tibble(
            Study_ID = accession,
            stage =
              "Full NMD MAGE-TAB parsing",
            error =
              conditionMessage(result)
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

    results[[accession]] <-
      result
  }

  valid_results <- results[
    vapply(
      results,
      function(x) {
        is.list(x) &&
          !is.null(
            x$study_summary
          )
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

    study_summary =
      purrr::map_dfr(
        valid_results,
        "study_summary"
      ),

    sample_summary =
      purrr::map_dfr(
        valid_results,
        "sample_summary"
      ),

    eligible_sample_summary =
      purrr::map_dfr(
        valid_results,
        "eligible_sample_summary"
      ),

    perturbation_sample_summary =
      purrr::map_dfr(
        valid_results,
        "perturbation_sample_summary"
      ),

    control_sample_summary =
      purrr::map_dfr(
        valid_results,
        "control_sample_summary"
      ),

    evidence =
      purrr::map_dfr(
        valid_results,
        "evidence"
      ),

    file_inventory =
      purrr::map_dfr(
        valid_results,
        "file_inventory"
      ),

    sdrf_wide =
      purrr::map_dfr(
        valid_results,
        "sdrf_wide"
      ),

    sdrf_long =
      purrr::map_dfr(
        valid_results,
        "sdrf_long"
      ),

    idf_long =
      purrr::map_dfr(
        valid_results,
        "idf_long"
      ),

    errors =
      errors,

    per_study =
      valid_results
  )

  if (
    nrow(
      combined$study_summary
    ) > 0L
  ) {
    combined$study_summary <-
      combined$study_summary %>%
      dplyr::arrange(
        factor(
          .data$screening_tier,
          levels = c(
            "Strong NMD RNA-seq candidate",
            "NMD perturbation with raw RNA-seq; control needs review",
            "NMD perturbation identified; raw availability needs review",
            "NMD perturbation not confirmed after full parsing"
          )
        ),
        .data$Study_ID
      )
  }

  combined
}


# -----------------------------------------------------------------------------
# 9. Write NMD output tables
# -----------------------------------------------------------------------------

ae_write_nmd_screening_results <- function(
    screening_results,
    out_dir =
      "ArrayExpress_NMD_Screening_Tables",
    prefix =
      "ArrayExpress_NMD") {

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
    "perturbation_sample_summary",
    "control_sample_summary",
    "evidence",
    "file_inventory",
    "sdrf_wide",
    "sdrf_long",
    "idf_long",
    "errors"
  )

  output_files <-
    stats::setNames(
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

    table <-
      screening_results[[table_name]]

    if (is.null(table)) {
      table <- data.frame()
    }

    data.table::fwrite(
      table,
      file =
        output_files[[table_name]],
      sep = "\t",
      quote = TRUE,
      na = ""
    )
  }

  invisible(
    output_files
  )
}


# -----------------------------------------------------------------------------
# 10. Example workflow
# -----------------------------------------------------------------------------

# source(
#   "arrayexpress_aml_normal_blood_screening_prefiltered_v1.1.1_titles_fixed.R"
# )
#
# source(
#   "arrayexpress_nmd_perturbation_screening.R"
# )
#
# nmd_screen <- ae_screen_nmd_accessions(
#   accessions = nmd_search$accessions,
#
#   base_dir =
#     "ArrayExpress_NMD_Metadata",
#
#   prescreen = TRUE,
#
#   # Recommended for the first pass:
#   prescreen_require_raw_evidence = TRUE,
#   prescreen_require_perturbation_evidence = TRUE,
#
#   # Do not require the phrase "cell line" during the cheap API stage because
#   # it may only appear in the SDRF. Full sample parsing does require a cell
#   # model before a sample is marked eligible.
#   prescreen_require_cell_model_evidence = FALSE,
#
#   prescreen_require_human_evidence = FALSE,
#   prescreen_exclude_clear_nonhuman = TRUE,
#
#   # Stress perturbations are noisy, so leave them off unless specifically
#   # exploring stress-mediated NMD regulation.
#   prescreen_include_stress_perturbations = FALSE,
#
#   prescreen_on_error = "continue",
#   overwrite = FALSE,
#   verbose = TRUE
# )
#
# Study-level overview:
#
# View(
#   nmd_screen$study_summary
# )
#
# NMD perturbation samples:
#
# View(
#   nmd_screen$perturbation_sample_summary
# )
#
# Matching controls:
#
# View(
#   nmd_screen$control_sample_summary
# )
#
# Eligible perturbation + control samples:
#
# View(
#   nmd_screen$eligible_sample_summary
# )
#
# Quick summary:
#
# nmd_screen$eligible_sample_summary %>%
#   dplyr::count(
#     .data$Study_ID,
#     .data$sample_role,
#     .data$nmd_target_or_study_target,
#     .data$perturbation_type
#   )
#
# Export:
#
# ae_write_nmd_screening_results(
#   nmd_screen,
#   out_dir = "ArrayExpress_NMD_Tables"
# )

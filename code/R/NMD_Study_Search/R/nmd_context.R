# =============================================================================
# NMD search context and repository-query builders
# =============================================================================
# This is the main file to edit when the biological search definition changes.
# It is intentionally separate from repository-specific parsing code.
# =============================================================================

nmd_default_context <- function() {
  list(
    explicit_nmd = c(
      "nonsense-mediated decay",
      "nonsense mediated decay",
      "nonsense-mediated mRNA decay",
      "nonsense mediated mRNA decay",
      "premature termination codon",
      "NMD"
    ),

    core_factors = c(
      "UPF1",
      "UPF2",
      "UPF3A",
      "UPF3B",
      "SMG1",
      "SMG5",
      "SMG6",
      "SMG7"
    ),

    extended_factors = c(
      "SMG8",
      "SMG9",
      "DHX34",
      "NBAS"
    ),

    factor_aliases = c(
      "hUPF1",
      "hUPF2",
      "UPF3X",
      "RENT1",
      "RENT2",
      "SMG-1",
      "SMG-5",
      "SMG-6",
      "SMG-7",
      "SMG-8",
      "SMG-9"
    ),

    genetic_perturbation = c(
      "knockdown",
      "knock-down",
      "knock down",
      "depletion",
      "depleted",
      "silencing",
      "silenced",
      "siRNA",
      "shRNA",
      "RNAi",
      "knockout",
      "knock-out",
      "knock out",
      "KO",
      "CRISPR",
      "CRISPR-Cas9",
      "degron",
      "AID",
      "auxin-inducible degron",
      "overexpression",
      "over-expression",
      "overexpressed",
      "OE",
      "inducible",
      "doxycycline",
      "dox"
    ),

    direct_inhibitors = c(
      "SMG1 inhibitor",
      "SMG-1 inhibitor",
      "SMG1i",
      "11j",
      "NMDI-1",
      "NMDI1",
      "VG1"
    ),

    indirect_inhibitors = c(
      "emetine",
      "cycloheximide",
      "CHX",
      "caffeine"
    ),

    stress_modulators = c(
      "tunicamycin",
      "thapsigargin",
      "sodium arsenite",
      "arsenite",
      "ER stress",
      "integrated stress response",
      "eIF2alpha phosphorylation",
      "eIF2α phosphorylation"
    ),

    controls = c(
      "control",
      "negative control",
      "scramble",
      "scrambled",
      "non-targeting",
      "nontargeting",
      "siCTRL",
      "siNC",
      "DMSO",
      "vehicle",
      "untreated",
      "mock",
      "wild type",
      "WT",
      "parental",
      "empty vector",
      "non-induced",
      "uninduced"
    ),

    rnaseq = c(
      "RNA-seq",
      "RNA seq",
      "RNA sequencing",
      "RNASeq",
      "transcriptome sequencing",
      "transcriptome profiling",
      "expression profiling by high throughput sequencing",
      "next gen sequencing",
      "next generation sequencing"
    ),

    raw_rnaseq = c(
      "FASTQ",
      ".fastq",
      ".fastq.gz",
      ".fq",
      ".fq.gz",
      "BAM",
      ".bam",
      "CRAM",
      ".cram",
      "SRA",
      "ENA",
      "SRR",
      "ERR",
      "DRR",
      "SRX",
      "ERX",
      "DRX",
      "PRJNA",
      "PRJEB"
    ),

    cell_models = c(
      "cell line",
      "cell culture",
      "cultured cells",
      "cancer cells",
      "tumor cells",
      "tumour cells",
      "HEK293",
      "HEK-293",
      "HEK293T",
      "HEK293TO",
      "HEK-293T",
      "Flp-In T-REx 293",
      "HeLa",
      "K562",
      "Huh7",
      "HT1080",
      "A498",
      "A-498",
      "A704",
      "A-704",
      "HK1",
      "HK-1",
      "SUNE1"
    )
  )
}


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

nmd_context_table <- function(context = nmd_default_context()) {
  context_df <- stack(context)
  colnames(context_df) <- c("term", "category")
  context_df <- context_df[, c("category", "term")]
  return(context_df)
  # pieces <- lapply(names(context), function(category) {
  #   data.frame(
  #     category = category,
  #     term = as.character(context[[category]]),
  #     stringsAsFactors = FALSE
  #   )
  # })
  # do.call(rbind, pieces)
}


nmd_quote_terms <- function(x) {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]

  # Quoting literal terms is safer for spaces, hyphens, and phrases.
  paste0('"', gsub('"', '\\\\"', x, fixed = TRUE), '"')
}


nmd_or_group <- function(x) {
  paste0("(", paste(nmd_quote_terms(x), collapse = " OR "), ")")
}


# ----------------------------------------------------------------------------
# ArrayExpress / BioStudies discovery queries
# ----------------------------------------------------------------------------
# These are deliberately fewer and easier to maintain than the older set of
# many near-duplicate factor-specific query strings.

build_arrayexpress_nmd_queries <- function(
  context = nmd_default_context(),
  include_broad = TRUE,
  include_indirect = TRUE,
  include_stress = FALSE
) {
  all_core <- unique(c(context$core_factors, context$factor_aliases))
  all_extended <- unique(context$extended_factors)

  # Generate 4 main queries
  strict <- c(
    NMD_RNAseq = paste(
      nmd_or_group(context$explicit_nmd),
      "AND",
      nmd_or_group(context$rnaseq)
    ),

    Core_factor_perturbation_RNAseq = paste(
      nmd_or_group(all_core),
      "AND",
      nmd_or_group(context$genetic_perturbation),
      "AND",
      nmd_or_group(context$rnaseq)
    ),

    Extended_factor_perturbation_RNAseq = paste(
      nmd_or_group(all_extended),
      "AND",
      nmd_or_group(context$genetic_perturbation),
      "AND",
      nmd_or_group(context$rnaseq)
    ),

    Direct_inhibitor_RNAseq = paste(
      nmd_or_group(context$direct_inhibitors),
      "AND",
      nmd_or_group(context$rnaseq)
    )
  )

  if (isTRUE(include_indirect)) {
    strict <- c(
      strict,
      Indirect_inhibitor_NMD_RNAseq = paste(
        nmd_or_group(context$indirect_inhibitors),
        "AND",
        nmd_or_group(c(context$explicit_nmd, all_core)),
        "AND",
        nmd_or_group(context$rnaseq)
      )
    )
  }

  if (isTRUE(include_stress)) {
    strict <- c(
      strict,
      Stress_NMD_RNAseq = paste(
        nmd_or_group(context$stress_modulators),
        "AND",
        nmd_or_group(c(context$explicit_nmd, all_core)),
        "AND",
        nmd_or_group(context$rnaseq)
      )
    )
  }

  if (!isTRUE(include_broad)) {
    return(strict)
  }

  broad <- c(
    Core_factor_RNAseq_broad = paste(
      nmd_or_group(all_core),
      "AND",
      nmd_or_group(context$rnaseq)
    ),

    Extended_factor_RNAseq_broad = paste(
      nmd_or_group(all_extended),
      "AND",
      nmd_or_group(context$rnaseq)
    )
  )

  c(strict, broad)
}


# ----------------------------------------------------------------------------
# GEO / Entrez discovery queries
# ----------------------------------------------------------------------------
# NCBI GEO DataSets supports Entrez field qualifiers. We explicitly limit
# searches to Series records (GSE) and high-throughput sequencing studies.

nmd_geo_literal_group <- function(x, field = "All Fields") {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]

  terms <- paste0('"', gsub('"', '\\\\"', x, fixed = TRUE), '"[', field, ']')
  paste0("(", paste(terms, collapse = " OR "), ")")
}


build_geo_nmd_queries <- function(
  context = nmd_default_context(),
  include_broad = TRUE,
  include_indirect = TRUE,
  include_stress = FALSE
) {
  all_core <- unique(c(context$core_factors, context$factor_aliases))
  all_extended <- unique(context$extended_factors)

  # The first term keeps results at the study/Series level. The assay terms use
  # official GEO DataSets indexed fields where possible, with an all-fields
  # RNA-seq fallback because submitter metadata is not exactly famous for
  # consistency.
  series_filter <- "gse[Entry Type]"
  rnaseq_filter <- paste0(
    "(",
    '"expression profiling by high throughput sequencing"[DataSet Type]',
    " OR ",
    '"high throughput sequencing"[Platform Technology Type]',
    " OR ",
    nmd_geo_literal_group(context$rnaseq, "All Fields"),
    ")"
  )

  strict <- c(
    NMD_RNAseq = paste(
      series_filter,
      "AND",
      rnaseq_filter,
      "AND",
      nmd_geo_literal_group(context$explicit_nmd)
    ),

    Core_factor_perturbation_RNAseq = paste(
      series_filter,
      "AND",
      rnaseq_filter,
      "AND",
      nmd_geo_literal_group(all_core),
      "AND",
      nmd_geo_literal_group(context$genetic_perturbation)
    ),

    Extended_factor_perturbation_RNAseq = paste(
      series_filter,
      "AND",
      rnaseq_filter,
      "AND",
      nmd_geo_literal_group(all_extended),
      "AND",
      nmd_geo_literal_group(context$genetic_perturbation)
    ),

    Direct_inhibitor_RNAseq = paste(
      series_filter,
      "AND",
      rnaseq_filter,
      "AND",
      nmd_geo_literal_group(context$direct_inhibitors)
    )
  )

  if (isTRUE(include_indirect)) {
    strict <- c(
      strict,
      Indirect_inhibitor_NMD_RNAseq = paste(
        series_filter,
        "AND",
        rnaseq_filter,
        "AND",
        nmd_geo_literal_group(context$indirect_inhibitors),
        "AND",
        nmd_geo_literal_group(c(context$explicit_nmd, all_core))
      )
    )
  }

  if (isTRUE(include_stress)) {
    strict <- c(
      strict,
      Stress_NMD_RNAseq = paste(
        series_filter,
        "AND",
        rnaseq_filter,
        "AND",
        nmd_geo_literal_group(context$stress_modulators),
        "AND",
        nmd_geo_literal_group(c(context$explicit_nmd, all_core))
      )
    )
  }

  if (!isTRUE(include_broad)) {
    return(strict)
  }

  broad <- c(
    Core_factor_RNAseq_broad = paste(
      series_filter,
      "AND",
      rnaseq_filter,
      "AND",
      nmd_geo_literal_group(all_core)
    ),

    Extended_factor_RNAseq_broad = paste(
      series_filter,
      "AND",
      rnaseq_filter,
      "AND",
      nmd_geo_literal_group(all_extended)
    )
  )

  c(strict, broad)
}


nmd_query_manifest <- function(
  context = nmd_default_context(),
  include_broad = TRUE,
  include_indirect = TRUE,
  include_stress = FALSE
) {
  ae <- build_arrayexpress_nmd_queries(
    context = context,
    include_broad = include_broad,
    include_indirect = include_indirect,
    include_stress = include_stress
  )

  geo <- build_geo_nmd_queries(
    context = context,
    include_broad = include_broad,
    include_indirect = include_indirect,
    include_stress = include_stress
  )

  rbind(
    data.frame(
      repository = "ArrayExpress",
      query_name = names(ae),
      query = unname(ae),
      stringsAsFactors = FALSE
    ),
    data.frame(
      repository = "GEO",
      query_name = names(geo),
      query = unname(geo),
      stringsAsFactors = FALSE
    )
  )
}

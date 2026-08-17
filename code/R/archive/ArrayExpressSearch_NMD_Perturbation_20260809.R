# =============================================================================
# NMD perturbation search vocabulary
# =============================================================================

# -----------------------------------------------------------------------------
# Core NMD factors
# -----------------------------------------------------------------------------

nmd_core_factors <- c(
  "UPF1",
  "UPF2",
  "UPF3A",
  "UPF3B",
  "SMG1",
  "SMG5",
  "SMG6",
  "SMG7"
)


# -----------------------------------------------------------------------------
# Additional / regulatory NMD factors
#
# These are useful for discovery, but I would keep them separate from the
# highest-confidence core-factor search.
# -----------------------------------------------------------------------------

nmd_extended_factors <- c(
  "SMG8",
  "SMG9",
  "DHX34",
  "NBAS"
)


# -----------------------------------------------------------------------------
# Older / alternative names that may occur in study descriptions
# -----------------------------------------------------------------------------

nmd_factor_aliases <- c(
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
)


# -----------------------------------------------------------------------------
# Genetic perturbation terminology
# -----------------------------------------------------------------------------

nmd_genetic_perturbation_terms <- c(
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
  "auxin inducible degron",

  "overexpression",
  "over-expression",
  "overexpressed",
  "overexpress",
  "OE",

  "inducible",
  "doxycycline",
  "dox"
)

# -----------------------------------------------------------------------------
# Direct / relatively specific pharmacologic NMD perturbations
# -----------------------------------------------------------------------------

nmd_pharmacologic_terms <- c(
  "SMG1 inhibitor",
  "SMG-1 inhibitor",
  "SMG1i",
  "11j",

  "NMDI-1",
  "NMDI1",
  "VG1"
)


# -----------------------------------------------------------------------------
# Broader / indirect NMD inhibition terms
#
# These should be treated as discovery terms, NOT sufficient evidence on their
# own that the experiment was specifically designed to perturb NMD.
# -----------------------------------------------------------------------------

nmd_indirect_inhibition_terms <- c(
  "caffeine",
  "emetine"
)


nmd_rnaseq_terms <- c(
  "RNA-seq",
  "RNA seq",
  "RNA sequencing",
  "RNASeq",
  "transcriptome sequencing",
  "transcriptome profiling",
  "expression profiling by high throughput sequencing"
)


nmd_raw_rnaseq_terms <- c(
  "FASTQ",
  "fastq",
  ".fastq",
  ".fastq.gz",
  "BAM",
  "bam",
  ".bam",
  "SRA",
  "ENA",
  "ERR",
  "SRR",
  "DRR"
)


arrayexpress_nmd_queries <- c(
  # ---------------------------------------------------------------------------
  # Direct pathway terminology
  # ---------------------------------------------------------------------------

  NMD_RNAseq = paste(
    '(',
    '"nonsense-mediated decay"',
    'OR "nonsense mediated decay"',
    'OR "nonsense-mediated mRNA decay"',
    'OR "nonsense mediated mRNA decay"',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR "RNASeq"',
    'OR "transcriptome sequencing"',
    ')'
  ),

  # ---------------------------------------------------------------------------
  # UPF1
  # ---------------------------------------------------------------------------

  UPF1_perturbation_RNAseq = paste(
    '(',
    'UPF1',
    'OR hUPF1',
    'OR RENT1',
    ')',
    'AND',
    '(',
    'knockdown',
    'OR "knock-down"',
    'OR depletion',
    'OR siRNA',
    'OR shRNA',
    'OR RNAi',
    'OR knockout',
    'OR "knock-out"',
    'OR CRISPR',
    'OR degron',
    'OR overexpression',
    'OR "over-expression"',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR RNASeq',
    'OR "transcriptome sequencing"',
    ')'
  ),

  # ---------------------------------------------------------------------------
  # UPF2 / UPF3
  # ---------------------------------------------------------------------------

  UPF2_UPF3_perturbation_RNAseq = paste(
    '(',
    'UPF2',
    'OR RENT2',
    'OR UPF3A',
    'OR UPF3B',
    'OR UPF3X',
    ')',
    'AND',
    '(',
    'knockdown',
    'OR depletion',
    'OR siRNA',
    'OR shRNA',
    'OR RNAi',
    'OR knockout',
    'OR CRISPR',
    'OR degron',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR RNASeq',
    'OR "transcriptome sequencing"',
    ')'
  ),

  # ---------------------------------------------------------------------------
  # SMG1
  # ---------------------------------------------------------------------------

  SMG1_perturbation_RNAseq = paste(
    '(',
    'SMG1',
    'OR "SMG-1"',
    ')',
    'AND',
    '(',
    'knockdown',
    'OR depletion',
    'OR siRNA',
    'OR shRNA',
    'OR knockout',
    'OR CRISPR',
    'OR inhibitor',
    'OR inhibition',
    'OR SMG1i',
    'OR 11j',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR RNASeq',
    'OR "transcriptome sequencing"',
    ')'
  ),

  # ---------------------------------------------------------------------------
  # Downstream SMG5 / SMG6 / SMG7
  # ---------------------------------------------------------------------------

  SMG5_6_7_perturbation_RNAseq = paste(
    '(',
    'SMG5',
    'OR "SMG-5"',
    'OR SMG6',
    'OR "SMG-6"',
    'OR SMG7',
    'OR "SMG-7"',
    ')',
    'AND',
    '(',
    'knockdown',
    'OR depletion',
    'OR siRNA',
    'OR shRNA',
    'OR RNAi',
    'OR knockout',
    'OR CRISPR',
    'OR degron',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR RNASeq',
    'OR "transcriptome sequencing"',
    ')'
  ),

  # ---------------------------------------------------------------------------
  # SMG8 / SMG9
  # ---------------------------------------------------------------------------

  SMG8_9_perturbation_RNAseq = paste(
    '(',
    'SMG8',
    'OR "SMG-8"',
    'OR SMG9',
    'OR "SMG-9"',
    ')',
    'AND',
    '(',
    'knockdown',
    'OR depletion',
    'OR knockout',
    'OR CRISPR',
    'OR inhibitor',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR RNASeq',
    ')'
  ),

  # ---------------------------------------------------------------------------
  # Auxiliary NMD factors
  # ---------------------------------------------------------------------------

  DHX34_NBAS_RNAseq = paste(
    '(',
    'DHX34',
    'OR NBAS',
    ')',
    'AND',
    '(',
    'knockdown',
    'OR depletion',
    'OR siRNA',
    'OR shRNA',
    'OR RNAi',
    'OR knockout',
    'OR CRISPR',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR RNASeq',
    'OR "transcriptome sequencing"',
    ')'
  ),

  # ---------------------------------------------------------------------------
  # SMG1 pharmacologic inhibition
  # ---------------------------------------------------------------------------

  SMG1_inhibitor_RNAseq = paste(
    '(',
    '"SMG1 inhibitor"',
    'OR "SMG-1 inhibitor"',
    'OR SMG1i',
    'OR 11j',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR RNASeq',
    'OR "transcriptome sequencing"',
    ')'
  ),

  # ---------------------------------------------------------------------------
  # Other NMD inhibitors
  # ---------------------------------------------------------------------------

  NMD_inhibitor_RNAseq = paste(
    '(',
    '"NMD inhibitor"',
    'OR "NMD inhibition"',
    'OR "nonsense-mediated decay inhibitor"',
    'OR NMDI1',
    'OR "NMDI-1"',
    'OR VG1',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR RNASeq',
    'OR "transcriptome sequencing"',
    ')'
  ),

  # ---------------------------------------------------------------------------
  # Translation-based / indirect perturbations
  #
  # Keep separate because these compounds have many effects outside NMD.
  # ---------------------------------------------------------------------------

  Translation_inhibitor_NMD_RNAseq = paste(
    '(',
    'emetine',
    'OR caffeine',
    ')',
    'AND',
    '(',
    '"nonsense-mediated"',
    'OR "nonsense mediated"',
    'OR UPF1',
    'OR SMG1',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR RNASeq',
    ')'
  )
)

arrayexpress_nmd_broad_queries <- c(
  UPF1_RNAseq = paste(
    '(',
    'UPF1',
    'OR hUPF1',
    'OR RENT1',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR RNASeq',
    'OR transcriptome',
    ')'
  ),

  UPF2_UPF3_RNAseq = paste(
    '(',
    'UPF2',
    'OR UPF3A',
    'OR UPF3B',
    'OR UPF3X',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR RNASeq',
    'OR transcriptome',
    ')'
  ),

  SMG1_RNAseq = paste(
    '(',
    'SMG1',
    'OR "SMG-1"',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR RNASeq',
    'OR transcriptome',
    ')'
  ),

  SMG5_6_7_RNAseq = paste(
    '(',
    'SMG5',
    'OR SMG6',
    'OR SMG7',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR RNASeq',
    'OR transcriptome',
    ')'
  ),

  Auxiliary_NMD_RNAseq = paste(
    '(',
    'DHX34',
    'OR NBAS',
    'OR SMG8',
    'OR SMG9',
    ')',
    'AND',
    '(',
    '"RNA-seq"',
    'OR "RNA sequencing"',
    'OR RNASeq',
    'OR transcriptome',
    ')'
  )
)

nmd_cell_model_terms <- c(
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

  "HeLa",
  "K562",
  "Huh7",
  "HT1080",

  "A498",
  "A-498",
  "A704",
  "A-704"
)

nmd_perturbation_classes <- c(
  "siRNA / RNAi knockdown",
  "shRNA knockdown",
  "CRISPR knockout",
  "inducible degradation",
  "SMG1 pharmacologic inhibition",
  "other NMD inhibitor",
  "translation inhibition",
  "NMD-factor overexpression",
  "stress-mediated NMD modulation",
  "other / unclear"
)

nmd_target_factors <- c(
  "UPF1",
  "UPF2",
  "UPF3A",
  "UPF3B",
  "SMG1",
  "SMG5",
  "SMG6",
  "SMG7",
  "SMG8",
  "SMG9",
  "DHX34",
  "NBAS"
)

nmd_stress_perturbation_terms <- c(
  "tunicamycin",
  "thapsigargin",
  "sodium arsenite",
  "arsenite",
  "ER stress",
  "integrated stress response",
  "eIF2alpha phosphorylation",
  "eIF2α phosphorylation"
)

arrayexpress_nmd_first_pass_queries <- c(
  arrayexpress_nmd_queries[
    c(
      "NMD_RNAseq",
      "UPF1_perturbation_RNAseq",
      "UPF2_UPF3_perturbation_RNAseq",
      "SMG1_perturbation_RNAseq",
      "SMG5_6_7_perturbation_RNAseq",
      "SMG8_9_perturbation_RNAseq",
      "DHX34_NBAS_RNAseq",
      "SMG1_inhibitor_RNAseq",
      "NMD_inhibitor_RNAseq",
      "Translation_inhibitor_NMD_RNAseq"
    )
  ],
  arrayexpress_nmd_broad_queries
)
source(
  "D:/R/NMD_Perturbation_Data/search_biostudies_arrayexpress_api_v1.0.1_fixed.R"
)
source(
  "D:/R/NMD_Perturbation_Data/arrayexpress_aml_normal_blood_screening_prefiltered_v1.1.1_titles_fixed.R"
)
source("D:/R/NMD_Perturbation_Data/arrayexpress_nmd_perturbation_screening.R")


nmd_search <- biostudies_search_many(
  queries = arrayexpress_nmd_queries,
  page_size = 100,
  pause_seconds = 0.15,
  verbose = TRUE
)

arrayexpress_nmd_all_queries <- c(
  arrayexpress_nmd_queries,
  arrayexpress_nmd_broad_queries
)

nmd_search <- biostudies_search_many(
  queries = arrayexpress_nmd_all_queries,
  page_size = 100,
  pause_seconds = 0.15,
  verbose = TRUE
)


names(arrayexpress_nmd_queries)

names(arrayexpress_nmd_broad_queries)

names(arrayexpress_nmd_all_queries)

names(arrayexpress_nmd_first_pass_queries)


nmd_search <- biostudies_search_many(
  queries = arrayexpress_nmd_first_pass_queries,
  page_size = 100,
  pause_seconds = 0.15,
  verbose = TRUE
)

nmd_ids <- nmd_search$accessions


write_biostudies_search_results(
  search_result = nmd_search,
  output_prefix = "NMD_Perturb_Search_ArrayExpress"
)


nmd_ids <- nmd_search$accessions

nmd_screen <- ae_screen_nmd_accessions(
  accessions = nmd_ids,
  base_dir = "ArrayExpress_NMD_Metadata",
  prescreen = TRUE,
  prescreen_require_raw_evidence = TRUE,
  prescreen_require_perturbation_evidence = TRUE,
  prescreen_require_cell_model_evidence = FALSE,
  prescreen_require_human_evidence = FALSE,
  prescreen_exclude_clear_nonhuman = TRUE,
  prescreen_include_stress_perturbations = FALSE,
  prescreen_on_error = "continue",
  overwrite = FALSE,
  verbose = TRUE
)


ae_write_nmd_screening_results(
  nmd_screen,
  out_dir = "ArrayExpress_NMD_Tables"
)

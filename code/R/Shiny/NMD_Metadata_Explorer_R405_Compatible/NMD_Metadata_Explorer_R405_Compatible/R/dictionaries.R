`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

NMD_TARGETS <- list(
  UPF1 = c("UPF1", "UPF-1", "RENT1", "RENT-1", "SMG2", "SMG-2"),
  UPF2 = c("UPF2", "UPF-2", "RENT2", "RENT-2", "SMG3", "SMG-3"),
  UPF3A = c("UPF3A", "UPF-3A"),
  UPF3B = c("UPF3B", "UPF-3B", "UPF3X"),
  SMG1 = c("SMG1", "SMG-1"),
  SMG5 = c("SMG5", "SMG-5"),
  SMG6 = c("SMG6", "SMG-6", "EST1A"),
  SMG7 = c("SMG7", "SMG-7"),
  SMG8 = c("SMG8", "SMG-8"),
  SMG9 = c("SMG9", "SMG-9"),
  DHX34 = c("DHX34", "DHX-34"),
  NBAS = c("NBAS"),
  SEC13 = c("SEC13", "SEC-13"),
  EIF4A3 = c("EIF4A3", "eIF4AIII", "eIF4A3"),
  RBM8A = c("RBM8A", "Y14"),
  MAGOH = c("MAGOH"),
  CASC3 = c("CASC3", "MLN51", "BTZ"),
  RNPS1 = c("RNPS1")
)

NMD_MECHANISMS <- list(
  "RNAi / knockdown" = c(
    "siRNA", "shRNA", "RNAi", "knockdown", "knock-down", "depletion",
    "depleted", "silencing", "silenced", "KD"
  ),
  "CRISPR / knockout" = c(
    "CRISPR", "CRISPR-Cas9", "Cas9", "knockout", "knock-out", "KO",
    "gene disruption", "null mutant"
  ),
  "Small-molecule NMD inhibition" = c(
    "NMD inhibitor", "NMD inhibition", "NMDI-1", "NMDI1", "NMDI-14",
    "NMDI14", "SMG1 inhibitor", "SMG-1 inhibitor", "SMG1i"
  ),
  "Translation inhibition" = c(
    "translation inhibitor", "translation inhibition", "cycloheximide", "CHX",
    "emetine", "anisomycin", "harringtonine", "puromycin"
  ),
  "Overexpression / induction" = c(
    "overexpression", "over-expression", "overexpressed", "ectopic expression",
    "inducible expression", "transfection", "transfected", "transduction",
    "transduced"
  ),
  "Rescue / complementation" = c(
    "rescue", "complementation", "complemented", "re-expression",
    "reexpression", "add-back", "addback"
  ),
  "Degron / acute depletion" = c(
    "degron", "dTAG", "auxin-inducible", "auxin inducible",
    "induced degradation", "acute depletion"
  ),
  "Mutation / dominant negative" = c(
    "dominant-negative", "dominant negative", "mutant", "loss-of-function",
    "loss of function", "LOF"
  ),
  "Readthrough treatment" = c(
    "readthrough", "read-through", "ataluren", "PTC124", "PTC-124", "G418",
    "geneticin", "gentamicin"
  )
)

NMD_AGENTS <- list(
  "NMDI-1" = c("NMDI-1", "NMDI1"),
  "NMDI-14" = c("NMDI-14", "NMDI14"),
  "SMG1 inhibitor" = c("SMG1 inhibitor", "SMG-1 inhibitor", "SMG1i"),
  Cycloheximide = c("cycloheximide", "CHX"),
  Emetine = c("emetine"),
  Anisomycin = c("anisomycin"),
  Harringtonine = c("harringtonine"),
  Puromycin = c("puromycin"),
  Ataluren = c("ataluren", "PTC124", "PTC-124"),
  G418 = c("G418", "geneticin"),
  Gentamicin = c("gentamicin")
)

CONTROL_CONCEPTS <- list(
  "Untreated control" = c("untreated", "no treatment"),
  "Vehicle control" = c("vehicle", "DMSO"),
  "Mock control" = c("mock"),
  "Non-targeting RNA control" = c(
    "non-targeting", "nontargeting", "scramble", "scrambled", "siCTRL",
    "siCON", "shCTRL", "control siRNA", "control shRNA", "luciferase siRNA"
  ),
  "Wild type" = c("wild type", "wild-type", "WT"),
  "Empty vector" = c("empty vector", "vector control")
)

CELL_MODEL_CONCEPTS <- list(
  HEK293 = c(
    "HEK293", "HEK-293", "HEK 293", "293T", "HEK293T", "HEK-293T",
    "HEK-293TO", "HEK293TO", "Flp-In T-REx 293", "Flp-In-T-REx-293"
  ),
  HeLa = c("HeLa", "HeLa-S3", "HeLa S3"),
  HCT116 = c("HCT116", "HCT-116"),
  U2OS = c("U2OS", "U-2 OS"),
  MCF7 = c("MCF7", "MCF-7"),
  Huh7 = c("Huh7", "Huh-7"),
  HepG2 = c("HepG2", "Hep-G2"),
  K562 = c("K562", "K-562"),
  A549 = c("A549", "A-549"),
  Jurkat = c("Jurkat"),
  "Primary cells" = c("primary cell", "primary cells"),
  "Patient tissue" = c("patient tissue", "patient sample", "tumor tissue", "tumour tissue")
)

SMART_QUERY_CONCEPTS <- c(
  NMD_TARGETS,
  NMD_MECHANISMS,
  NMD_AGENTS,
  CONTROL_CONCEPTS,
  CELL_MODEL_CONCEPTS,
  list(
    "Nonsense-mediated decay" = c(
      "nonsense-mediated decay", "nonsense mediated decay", "NMD"
    ),
    "Premature termination codon" = c(
      "premature termination codon", "premature stop codon", "PTC"
    ),
    "RNA sequencing" = c("RNA-seq", "RNA seq", "RNA sequencing", "transcriptome sequencing")
  )
)

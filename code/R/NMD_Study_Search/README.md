# Clean NMD Perturbation Metadata Workflow

This folder condenses the previous GEO and ArrayExpress scripts into one workflow for:

1. defining the NMD search context;
2. discovering candidate studies in GEO and/or ArrayExpress;
3. retrieving study/sample metadata;
4. screening samples for NMD perturbation, controls, RNA-seq, and raw-data evidence;
5. ranking studies for review; and
6. writing the useful tables into a single Excel workbook.

## Files you normally care about

```text
NMD_workflow_clean/
├── run_nmd_workflow.R       # main script: source this
├── install_packages.R       # one-time package helper
├── README.md
└── R/
    ├── nmd_context.R        # biological search terms + query builders
    ├── geo_workflow.R       # GEO discovery + your GEO NMD parser
    ├── arrayexpress_workflow.R  # AE search + MAGE-TAB helpers + NMD screen
    └── export_workflow.R    # unified ranking + Excel workbook
```

The point is that you should normally interact with **`run_nmd_workflow.R`** and **`R/nmd_context.R`**. The repository-specific files are implementation details unless you are debugging metadata parsing.

## Basic use

From the workflow folder:

```r
source("run_nmd_workflow.R")

result <- run_nmd_workflow()
```

This searches both GEO and ArrayExpress, screens the resulting studies, and creates:

```text
data/NMD_metadata_screen/NMD_Perturbation_Metadata_Screen_YYYYMMDD.xlsx
```

It also saves the complete R object as an RDS file so you do not have to reconstruct the combined result just to inspect it later.

## Inspect only new IDs

This is probably the most useful mode while curating studies manually:

```r
result <- run_nmd_workflow(
  discover = FALSE,
  geo_ids = c(
    "GSE305669"
  ),
  arrayexpress_ids = c(
    "E-MTAB-9330",
    "E-MTAB-10716"
  )
)
```

You can provide GEO IDs, ArrayExpress IDs, or both.

## Search only one repository

```r
result <- run_nmd_workflow(
  repositories = "GEO"
)
```

or:

```r
result <- run_nmd_workflow(
  repositories = "ArrayExpress"
)
```

## Change the search definition

The main biological vocabulary is in:

```text
R/nmd_context.R
```

The default context contains:

- explicit NMD terminology;
- core NMD factors: UPF1, UPF2, UPF3A/B, SMG1, SMG5/6/7;
- extended factors: SMG8, SMG9, DHX34, NBAS;
- aliases such as RENT1/2 and UPF3X;
- knockdown/knockout/degron/overexpression terms;
- direct NMD/SMG1 inhibitors;
- indirect translation inhibitors;
- optional stress-mediated terms;
- control terms;
- RNA-seq and raw-file terms; and
- common cell-model terms.

You can create a custom context without editing the defaults:

```r
my_context <- nmd_default_context()
my_context$direct_inhibitors <- c(
  my_context$direct_inhibitors,
  "My new compound"
)

result <- run_nmd_workflow(
  context = my_context
)
```

## Strict vs broad discovery

By default the workflow runs both:

- **strict searches**, which combine NMD factor/agent evidence with perturbation language and RNA-seq; and
- **broad searches**, which combine NMD factors with RNA-seq even if the repository search record does not explicitly mention the perturbation method.

This intentionally favors recall during discovery. The downloaded sample metadata is then used for the more specific screening step.

To turn off the broad searches:

```r
result <- run_nmd_workflow(
  include_broad_queries = FALSE
)
```

Stress-mediated NMD modulation is off by default because it is much noisier:

```r
result <- run_nmd_workflow(
  include_stress_modulators = TRUE
)
```

## Excel output

The workbook begins with four useful sheets:

- `README`
- `Search_Context`
- `Search_Queries`
- `Unified_Ranking`

`Unified_Ranking` puts GEO and ArrayExpress studies into one review table while preserving the repository-specific classification.

The standardized 0-100 metadata score is deliberately simple and visible:

| Evidence | Points |
|---|---:|
| NMD perturbation sample(s) | 30 |
| candidate control arm | 20 |
| raw RNA-seq evidence | 20 |
| RNA-seq evidence | 10 |
| NMD target/inhibitor identified | 10 |
| explicit NMD terminology | 5 |
| cell-model evidence | 5 |

The score is for triage, not biological truth. Use the evidence and sample sheets for manual verification.

The remaining sheets retain repository-specific detail, including search hits, study summaries, sample summaries, evidence, file inventory, and optionally long-format metadata.

If the workbook gets too large, omit the long tables:

```r
result <- run_nmd_workflow(
  include_long_metadata_in_excel = FALSE
)
```

## What the workflow does with raw RNA-seq evidence

ArrayExpress/BioStudies has explicit file inventories and run/file metadata in the current screening code, so raw-sequence evidence can be screened directly.

For GEO, the workflow looks for FASTQ/BAM/SRA/ENA-style file or accession evidence retained in the GEO sample metadata. This is useful for triage but should still be confirmed before downloading data.

## Suggested cleanup of the old folder

I would **archive rather than delete** the older scripts until this workflow has been exercised on your known studies.

### Replaced by the new main workflow

- `NMD_Perturbation_Meta_Parsing_20260712.R`
- `NMD_Perturbation_Meta_Parsing_20260712_718.R`
- `ArrayExpressSearch_NMD_Perturbation_20260809.R`
- `arrayexpress_work.R`
- `arrayexpress_work_v2.R`
- `look_over_array_express_res.R`
- `Agg_GSE_AE_Sample_metadata.R`

These are mostly driver/exploration scripts whose useful pieces are now represented in `run_nmd_workflow.R`.

### ArrayExpress implementation history that can move to an archive folder

- `arrayexpress_metadata_parser.R`
- `arrayexpress_sdrf_character_patch.R`
- `arrayexpress_one_row_samples_patch.R`
- `arrayexpress_metadata_parser_with_one_row_samples.R`
- `search_biostudies_arrayexpress_api_v1.0.1_fixed.R`
- `arrayexpress_aml_normal_blood_screening_prefiltered_v1.1.1_titles_fixed.R`
- `arrayexpress_nmd_perturbation_screening.R`

The cleaned `R/arrayexpress_workflow.R` contains the search code, only the generic MAGE-TAB helpers needed by NMD screening, and the NMD-specific screen. The AML biology is no longer a dependency.

### GEO implementation history that can move to an archive folder

- `annotate_gse_ids.R`
- `annotate_gse_ids_718.R`
- `nmd_geo_metadata_parser.R`

The improved `extract_gse_meta_v2()` and the NMD parser are retained in `R/geo_workflow.R`.

### Keep separately

- `Setup_renv.R`: environment/setup history. Clean it later, but it is not part of the biological workflow.
- `Organizing_Biostools_Paths_apply_keywords.R`: unrelated to this NMD repository workflow.

## Recommended project layout

```text
Public_NMD_Perturbation_Data/
├── code/
│   ├── nmd_workflow/          # this folder
│   └── archive/               # old scripts after validation
├── data/
│   ├── curated/               # your manually curated study/sample tables
│   ├── NMD_metadata_screen/   # workflow output + metadata cache
│   └── raw/                   # later FASTQ/BAM downloads, if desired
├── renv/
├── renv.lock
└── README.md
```

## Validation before retiring the old scripts

I would test the clean workflow on a handful of known studies that exercise different designs:

- a straightforward UPF1 siRNA study;
- a direct inhibitor study;
- a degron/time-course study;
- a study with paired-end/multiple-file metadata;
- a study where the control does not itself mention the NMD factor; and
- one known poor match.

Compare the new `GEO_Samples`, `AE_Samples`, `AE_Files`, evidence sheets, and study ranking against your manually curated tables. Once those agree well enough, move the older scripts into `code/archive/` instead of keeping five generations in the active code directory, where they will inevitably develop opinions about which one is canonical.

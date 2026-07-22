# Summer 2026 Research Project Evolution

> The project evolved from systematic dataset identification into metadata harmonization, study-quality assessment, and development of a searchable resource.
# Working title

**Systematic Identification and Harmonization of Public RNA-seq Datasets for Nonsense-Mediated mRNA Decay Perturbation Analysis**
## Central research question

**What public RNA-seq datasets involving perturbation of nonsense-mediated mRNA decay are available, and which studies are sufficiently well annotated, controlled, replicated, and comparable for downstream analysis of NMD-sensitive transcripts and alternative splicing?**
# Mini-project aims
## Aim 1: Systematically identify candidate NMD perturbation RNA-seq studies

Identify publicly available studies through literature review, GEO, ArrayExpress/BioStudies, and accession cross-referencing.

Capture studies involving:
- Genetic depletion or knockout of NMD factors
- Overexpression or rescue of NMD factors
- Pharmacologic NMD inhibition
- Translation inhibition used to suppress NMD
- Readthrough treatments
- Related perturbations of EJC or NMD-associated proteins
## Aim 2: Harmonize heterogeneous study and sample metadata

Develop reproducible R workflows to:
- Retrieve study- and sample-level metadata
- Convert heterogeneous metadata into long-format tables
- Extract experimental factors, biological models, treatments, controls, and protocols
- Identify candidate perturbation and control samples
- Retain the original metadata evidence supporting each classification
- Standardize terminology and aliases across repositories
## Aim 3: Assess study readiness and cross-study comparability

Evaluate each study according to:
- Clarity of NMD perturbation
- Presence of a matched control
- Number of biological replicates
- Availability of raw RNA-seq data
- Biological model
- Organism
- Perturbed NMD factor
- Perturbation mechanism
- Sequencing design
- Metadata completeness
- Similarity to other candidate studies

The output would be a shortlist of studies suitable for:
1. Individual reanalysis
2. Comparisons among biologically similar studies
3. Cross-study integration
4. Exclusion or manual review
## <u>Exploratory</u> Aim 4: Develop a searchable study resource

Build a Shiny application allowing users to search and filter GEO and ArrayExpress metadata by:
- NMD factor
- Inhibitor or treatment
- Perturbation mechanism
- Cell line or tissue
- Disease context
- Organism
- Study accession
- Control type
- Study-readiness category
# What has been completed

- Conducted literature and repository searches for NMD perturbation studies
- Compiled an initial set of candidate GEO and ArrayExpress accessions
- Retrieved GEO metadata using `GEOquery`
- Retrieved ArrayExpress/BioStudies MAGE-TAB metadata
- Identified substantial metadata heterogeneity across studies
- Expanded aggregated GEO characteristics into individual metadata fields
- Converted metadata into standardized long-format structures
- Developed keyword and synonym-based annotation approaches
- Classified candidate NMD perturbations, controls, agents, and biological models
- Preserved field-level evidence for classifications
- Generated study-, sample-, protocol-, and evidence-level output tables
- Developed a Shiny application for searching and exporting filtered subsets
- Established a framework for manually validating studies before expression-level analysis
# Analysis to add

## Study-readiness and comparability assessment
### Study-readiness variables

One manually reviewed row per study with fields such as:

| Variable                   | Possible values                                 |
| -------------------------- | ----------------------------------------------- |
| Repository                 | GEO / ArrayExpress                              |
| NMD factor                 | UPF1, UPF2, SMG6, etc.                          |
| Perturbation mechanism     | siRNA, shRNA, CRISPR, inhibitor, overexpression |
| Direct NMD perturbation    | Yes / No / Unclear                              |
| Organism                   | Human, mouse, other                             |
| Biological model           | Cell line, primary cells, tissue                |
| Cell line or tissue        | Standardized name                               |
| Disease context            | Cancer, genetic disease, normal model           |
| Matched control            | Yes / No / Unclear                              |
| Control type               | Scramble, vehicle, wild type, empty vector      |
| Biological replicates      | Number per group                                |
| Raw RNA-seq available      | Yes / No                                        |
| Paired-end status          | Yes / No / Unknown                              |
| Metadata confidence        | High / Moderate / Low                           |
| Manual validation complete | Yes / No                                        |
| Comparable study group     | Group label                                     |
| Final recommendation       | Include / Review / Exclude                      |

### Readiness score

|Criterion|Points|
|---|---|
|Direct NMD-factor perturbation clearly defined|2|
|Matched control clearly identified|2|
|At least three replicates per primary group|2|
|Two replicates per group|1|
|Raw RNA-seq files available|2|
|Sample-to-treatment assignments clear|1|
|Protocol sufficiently documented|1|
|Biological model overlaps another study|1|
|Sequencing design documented|1|

Possible categories:
- **Tier 1: 9–12 points**  
    Strong candidate for analysis or cross-study comparison.
- **Tier 2: 6–8 points**  
    Potentially useful but requires additional review or has limited comparability.
- **Tier 3: 0–5 points**  
    Incomplete metadata, weak controls, insufficient replication, or unclear relevance.
# Figures
## Figure 1: Dataset discovery and curation workflow

Flowchart:

```
Literature and repository search
              ↓
Candidate accessions collected
              ↓
GEO and ArrayExpress metadata retrieval
              ↓
Metadata normalization and long-format conversion
              ↓
Automated annotation and synonym matching
              ↓
Manual sample and protocol validation
              ↓
Study-readiness and comparability assessment
              ↓
Searchable Shiny resource and prioritized dataset shortlist
```
## Figure 2: Landscape of candidate studies

A bar chart showing the number of studies by:
- Perturbed NMD factor
- Perturbation mechanism
- Organism
- Repository

Descriptive result:

> Most candidate studies relied on UPF1 depletion, while fewer studies targeted downstream NMD factors or used direct pharmacologic inhibition.

The exact conclusion would depend on final validated table.
## Figure 3: Study-readiness matrix

Create a heatmap-like table:

|Study|Direct perturbation|Control|≥3 replicates|Raw data|Clear protocol|Tier|
|---|---|---|---|---|---|---|
|GSE...|✓|✓|✓|✓|✓|1|
|E-MTAB...|✓|✓|✗|✓|✓|2|
|GSE...|?|✗|?|✓|✗|3|

## Figure 4: Comparable study groups

A table or bubble plot grouping studies into categories such as:
- Human UPF1 knockdown in cancer cell lines
- Human UPF1 knockdown in nonmalignant cell lines
- Mouse NMD-factor perturbation
- Translation-inhibition experiments
- Direct NMD inhibitor experiments
- NMD-factor rescue or overexpression
- EJC-associated perturbations

The main conclusion could be:

> Although many candidate datasets were identified, only a subset had sufficiently similar perturbations, controls, biological models, and replicate structures for direct cross-study comparison.
# Remaining work for the next week
## Step 1: Validate the ArrayExpress studies

For each study:
- Check the IDF description
- Check SDRF factor values
- Confirm the perturbation
- Confirm the control
- Count samples and replicates
- Verify raw sequencing availability
- Record unclear cases
## Step 2: Produce a final master inventory

Keep GEO and ArrayExpress source tables separate internally, but create one compact manually reviewed summary table with common study-level fields.
## Step 3: Assign readiness tiers

Apply the rubric and manually override obvious errors.

Include:

```
automated_score
manual_final_tier
manual_review_notes
```

This makes it clear that automated parsing supports rather than replaces expert review.
## Step 4: Identify the strongest comparable subset

Examples might include:
- UPF1 knockdown with non-targeting controls
- Human cell lines with at least three replicates
- Paired-end bulk RNA-seq
- Raw FASTQ or aligned reads available

You do not have to complete a full cross-study analysis. The deliverable can be:

> Three studies were identified as sufficiently similar for an initial comparison of NMD-sensitive transcript and splice-junction responses.
## Step 5: Demonstrate the Shiny application

Show one or two searches during the presentation:

```
UPF1 knockdown human
```

and:

```
SMG1 inhibitor
```

Then show:
- Filtered study table
- Sample evidence
- Protocol metadata
- ZIP download
# Recommended final deliverable package

Your submission or project folder could contain:

```
NMD_Perturbation_Summer_Project/
├── README.md
├── report/
│   ├── NMD_dataset_curation_report.pdf
│   └── figures/
├── presentation/
│   └── NMD_summer_project_slides.pptx
├── metadata/
│   ├── GEO_study_summary_reviewed.tsv
│   ├── ArrayExpress_study_summary_reviewed.tsv
│   ├── combined_study_readiness.tsv
│   └── comparable_study_groups.tsv
├── scripts/
│   ├── GEO_metadata_parser.R
│   ├── ArrayExpress_metadata_parser.R
│   ├── study_readiness_analysis.R
│   └── generate_summary_figures.R
└── shiny_app/
    └── NMD_Metadata_Explorer/
```

# Presentation structure

A seven-slide presentation would be enough.
## Slide 1: Motivation

- NMD regulates transcript stability and interacts with alternative splicing.
- Public perturbation studies could support identification of NMD-sensitive transcript features.
- These studies are difficult to compare because metadata and experimental designs are inconsistent.
## Slide 2: Research objective

> Develop a systematic and reproducible framework to identify, annotate, and prioritize public NMD perturbation RNA-seq datasets for downstream comparative analysis.
## Slide 3: Dataset discovery

- Literature search
- GEO accessions
- ArrayExpress accessions
- Initial inclusion criteria
- Challenges encountered
## Slide 4: Metadata harmonization workflow

Show the pipeline from repository metadata to standardized study and sample tables.
## Slide 5: Searchable resource

Include screenshots of the Shiny application and example searches.
## Slide 6: Study-readiness results

Show:
- Number of candidate studies
- Perturbation mechanisms
- Readiness tiers
- Most comparable study cluster
## Slide 7: Conclusions and next steps

- A curated public NMD perturbation resource was created.
- Automated parsing reduced manual review burden but did not eliminate it.
- A subset of studies appears suitable for downstream transcript and splice-junction analysis.
- Future work will use these studies to quantify NMD-sensitive transcript architecture and splicing events.
# Reusable mini-project description
## Systematic Identification and Harmonization of Public NMD Perturbation RNA-seq Datasets

Nonsense-mediated mRNA decay (NMD) is an RNA surveillance and regulatory pathway that influences transcript stability, alternative splicing outcomes, and gene-expression programs. Public RNA-sequencing studies involving genetic or pharmacologic perturbation of NMD provide an opportunity to identify transcript and splicing features associated with NMD susceptibility. However, reuse of these studies is limited by inconsistent metadata, variable experimental designs, unclear sample-to-treatment assignments, and differences in terminology across repositories.

The objective of this project was to develop a reproducible framework for identifying, annotating, and prioritizing public NMD perturbation RNA-seq datasets for downstream comparative analysis. Candidate studies were identified through literature review and searches of the Gene Expression Omnibus and ArrayExpress/BioStudies repositories. Study- and sample-level metadata were retrieved using R-based workflows and converted into standardized wide- and long-format tables. Because metadata fields differed substantially among studies, annotation methods were developed to search across study descriptions, sample characteristics, experimental factors, treatment protocols, and protocol definitions.

Candidate samples were classified according to the perturbed NMD factor, perturbation mechanism, treatment agent, biological model, organism, disease context, and control type. Synonym-based searches were used to recognize alternative names for NMD factors, cell lines, inhibitors, and experimental methods. Field-level evidence was retained for each annotation so that automated classifications could be manually reviewed. Separate workflows were developed for GEO and ArrayExpress metadata because of differences in repository structure.

A Shiny application was developed to make the resulting resource searchable and accessible. The application allows users to filter studies by NMD factor, perturbation mechanism, treatment, organism, cell line, tissue, disease context, control type, and accession. Filtered study, sample, protocol, and evidence tables can be reviewed within the application and exported as a ZIP archive.

The final phase of the project evaluates study readiness and cross-study comparability. Studies are assessed according to perturbation clarity, matched-control availability, biological replication, raw RNA-seq availability, metadata completeness, sequencing design, and overlap with other experimental models. This assessment will identify a prioritized subset of studies suitable for downstream analysis of NMD-sensitive transcripts and splice-junction usage.

The project produced a curated study inventory, reproducible metadata-processing scripts, standardized study and sample tables, a searchable Shiny resource, and a framework for selecting datasets for future comparative analyses. This work establishes the data foundation needed to investigate how transcript architecture and alternative splicing contribute to NMD susceptibility across experimental systems.
# Concise project outcome statement

You can summarize the accomplishment as:

> This project transformed a heterogeneous collection of public NMD perturbation studies into a reproducible, searchable, and manually verifiable resource for selecting RNA-seq datasets suitable for comparative transcript and splicing analysis.

The strongest remaining task is the **study-readiness matrix and prioritized shortlist**. That converts everything you have built into a clear research conclusion: not just what data exist, but which data can credibly answer the next biological question.


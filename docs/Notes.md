# 2026.07.18
- Setup and synched local positron environment at work with git

## array_fetch_many()

```r
study_ids <- c(
  "E-MTAB-9330","E-MTAB-10716","E-MTAB-8461","E-MTAB-13839","E-MTAB-13788",
  "E-MTAB-13789","E-MTAB-13787","E-MTAB-13829","E-MTAB-13829","E-MTAB-13836",
  "E-MTAB-16399","E-MTAB-13949","E-MTAB-14755","E-MTAB-14755","E-MTAB-14725",
  "E-MTAB-14725","E-MTAB-13837"
)
ae_results <- ae_fetch_many(
  accessions = study_ids,
  base_dir = "data/ArrayExpress_Metadata_v2",
  overwrite = FALSE,
  use_api_fallback = TRUE,
  verbose = TRUE
)
```

- Clean input accession variables
- setup empty output formats
- For each accession ID in the list of accession IDs:
	- ae_fetch_study_metadata()
		- Download mege-tab meta data
		- https://www.ebi.ac.uk/biostudies/misc/MAGE-TABv1.1_2011_07_28.pdf
		- ae_download_magetab()
			- ae_download_with_package()
				- Download project data path information
			- If data cannot be downloaded with ArrayExpress package:
				- ae_download_from_biostudies()
					- ae_json_from_url()
						- ae_retry()
						- Get SDRF study URL info
						- Get IDF study URL info
					- Parse all strings of the json object and detect SDRF and IDF urls
					- Get FTP Link
					- Download identified files via url
						- connect FTP and relative file path
						- Set file destination
						- ae_download_file()
							- Normalize paths
							- download via ae_retry() and libcurl
					- Format list for project data with newly downloaded file names
	- Now have downloaded project metatdata SDRF and IDF files
	- Read in SDRF files
		- ae_read_sdrf_file()
			- ae_make_sample_id()
				- Check to see if usable sample ID columns
					- Default is "Samples #" 
					- **May nee to look back at this , might need to expand search terms**
					- If so, remap, if not keep
			- Add download information and main study info to SDRF table
			- Pivot table to a long format
				- Clean field names
					- ae_normalize_field()
				- Clean values
					- stringr::str_squish()
				- Filter NAs
				- Select needed columns for new df
			- Subset out new df for project protocol fields
			- Return list of 3 dfs: original, long, protocol
	- Read in IDF files
		- ae_read_idf_file()

Currently around line 950 in arrayexpress_metadata_parser.R
- trying to see if we can stick to always using a fully resolved path
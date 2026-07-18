# GSE Study Lookup: Extract Disease and Tissue Information from GEO
# Requires: GEOquery, dplyr
# Install if needed:
#   BiocManager::install("GEOquery")
#   install.packages("dplyr")

library(GEOquery)
library(dplyr)


# Retry wrapper: attempts `expr` up to `n` times with a pause between tries
retry <- function(expr, n = 3, wait = 10, verbose = TRUE) {
  for (i in seq_len(n)) {
    result <- tryCatch(expr, error = function(e) e)
    if (!inherits(result, "error")) return(result)
    if (verbose) message(sprintf("  Attempt %d/%d failed: %s", i, n, result$message))
    if (i < n) {
      if (verbose) message(sprintf("  Retrying in %ds...", wait))
      Sys.sleep(wait)
    }
  }
  stop(result)  # re-throw the last error after all attempts exhausted
}

# -------------------------------------------------------------------------
# Core lookup function for a single GSE accession
# -------------------------------------------------------------------------

extract_gse_meta <- function(gse_id, verbose = FALSE, retries = 3, retry_wait = 10) {
  char_df_out <- tryCatch({
    if (verbose) message("Fetching: ", gse_id)
    
    # Set a longer curl timeout (seconds) before the call
    old_timeout <- getOption("timeout")
    on.exit(options(timeout = old_timeout))  # restore on exit no matter what
    options(timeout = 300)                   # 5 minutes
    
    gse <- retry(
      getGEO(gse_id, GSEMatrix = FALSE),
      n = retries,
      wait = retry_wait,
      verbose = verbose
    )
    
    #gse <- getGEO(gse_id, GSEMatrix = FALSE)
    gse_title <- Meta(gse)$title
    
    # --- Extract disease from sample characteristics ---
    samples <- GSMList(gse)
    
    # Collect all characteristic fields across samples
    all_chars <- lapply(samples, function(s) {
      chars_col <- grep("characteristics_ch1",names(Meta(s)))[1]
      if (length(chars_col) == 0 | is.na(chars_col)) {
        chars_col <- grep("description",names(Meta(s)))[1]
        if (length(chars_col) == 0 | is.na(chars_col)) return(unlist(Meta(s)))
      }
      chars <- Meta(s)[[chars_col]]
      # Parse "key: value" pairs
      pairs <- strsplit(chars, ": ", fixed = TRUE)
      keys <- sapply(pairs, function(p) tolower(trimws(p[1])))
      values <- sapply(pairs, function(p) if (length(p) >= 2) trimws(p[2]) else NA)
      values <- setNames(values, keys)
      meta_chars <- unlist(Meta(s)[-chars_col])
      meta_out <- c(meta_chars,GSE_Title = gse_title,values)
      return(meta_out)
    })
    
    # Flatten into a data frame of key/value pairs
    char_df <- bind_rows(lapply(all_chars, as.data.frame.list,
                                stringsAsFactors = FALSE), .id = "Sample_ID")
    
    return(char_df)
    
  }, error = function(e) {
    warning("Failed to fetch ", gse_id, " metadata: ", conditionMessage(e))
    return(NA_character_)
  })
  return(char_df_out)
}

lookup_gse <- function(gse_id, verbose = FALSE) {
  gse_df_out <- tryCatch({
    if (verbose) message("Fetching: ", gse_id)
    
    ## --- Extract disease from sample characteristics ---
    char_df <- if (is.character(gse_id)) extract_gse_meta(gse_id) else gse_id
    
    # --- Disease: look for common field names ---
    disease_keys <- c("disease", "diagnosis", "condition", "phenotype")
    disease_col <- grep(paste0(disease_keys, collapse = "|"), colnames(char_df), value = TRUE, ignore.case = TRUE)
    if (length(disease_col) > 0) {
      n_samples <- nrow(char_df)
      disease_col <- Filter(function(x) {
        vals <- na.omit(char_df[[x]])
        n_unique <- length(unique(vals))
        if (n_unique < 1 | n_unique > max(2, floor(n_samples * 0.25))) return(FALSE)
        !all(grepl("^\\s*-?[0-9]*\\.?[0-9]+\\s*$|^NA$", vals))
      }, disease_col)
      disease <- if (length(disease_col) == 0) NA_character_ else
        paste(sapply(disease_col, function(x) {
          paste(sort(unique(na.omit(char_df[[x]]))), collapse = "; ")
        }), collapse = ";; ")
      if (length(disease_col) > 0) {
        disease_cols <- paste(disease_col, collapse = "; ")
      } else {
        disease_cols <- NA_character_
      }
    } else {
      disease <- NA_character_
      disease_cols <- NA_character_
    }
    
    
    # --- Tissue: look for common field names ---
    tissue_keys <- c("tissue", "organ", "source", "source[[:punct:]]name", "cell[[:punct:]]type",
                     "cell[[:punct:]]line", "organ", "sample[[:punct:]]type",
                     "source name", "cell type", "cell line", "sample type")
    tissue_col <- grep(paste0(tissue_keys, collapse = "|"), colnames(char_df), value = TRUE, ignore.case = TRUE)
    tissue_col <- grep("organism", tissue_col, invert = TRUE, value = TRUE)
    if (length(tissue_col) > 0) {
      n_samples <- nrow(char_df)
      tissue_col <- Filter(function(x) {
        vals <- na.omit(char_df[[x]])
        n_unique <- length(unique(vals))
        if (n_unique < 1 | n_unique > max(2, floor(n_samples * 0.25))) return(FALSE)
        !all(grepl("^\\s*-?[0-9]*\\.?[0-9]+\\s*$|^NA$", vals))
      }, tissue_col)
      tissue <- if (length(tissue_col) == 0) NA_character_ else
        paste(sapply(tissue_col, function(x) {
          paste(sort(unique(na.omit(char_df[[x]]))), collapse = "; ")
        }), collapse = ";; ")
      if (length(tissue_col) > 0) {
        tissue_cols <- paste(tissue_col, collapse = "; ")
      } else {
        tissue_cols <- NA_character_
      }
    } else {
      tissue <- NA_character_
      tissue_cols <- NA_character_
    }
    
    # Fallback: use series-level title / summary for context
    study_title <- paste(unique(char_df$GSE_Title), collapse = "; ") %||% NA_character_
    organism_col <- grep("organism",colnames(char_df), value = T)[1]
    organism <- paste(unique(char_df[[organism_col]]), collapse = "; ")
    
    gse_df <- data.frame(
      title = study_title,
      organism = organism,
      disease_col = disease_cols,
      disease = disease,
      tissue_col = tissue_cols,
      tissue = tissue,
      stringsAsFactors = FALSE
    )
    return(gse_df)
    
  }, error = function(e) {
    warning("Failed to fetch ", gse_id, ": ", conditionMessage(e))
    gse_df <- data.frame(
      title = NA_character_,
      organism = NA_character_,
      disease_col = NA_character_,
      disease = NA_character_,
      tissue_col = NA_character_,
      tissue = NA_character_,
      stringsAsFactors = FALSE
    )
    return(gse_df)
  })
  return(gse_df_out)
}

# Null-coalescing helper (like %||% in newer R / rlang)
`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0) x else y


# -------------------------------------------------------------------------
# Batch function: takes a vector/column of mixed strings, detects GSE IDs,
# and returns a summary table
# -------------------------------------------------------------------------
#
# Parameters:
#   ids      – character vector that may contain GSE accession numbers
#              (e.g. "GSE12345") mixed with other strings / NAs
#   verbose  – print progress messages?
#
# Returns:
#   data.frame with columns: gse_id, title, organism, disease, tissue
#
annotate_gse_ids <- function(ids, verbose = TRUE) {
  
  # Extract only valid-looking GSE accessions (case-insensitive)
  gse_pattern <- "GSE[0-9]+"
  raw_matches <- regmatches(as.character(ids),
                            gregexpr(gse_pattern, as.character(ids),
                                     ignore.case = TRUE))
  found       <- unique(toupper(unlist(raw_matches)))
  found       <- found[nzchar(found)]
  
  if (length(found) == 0) {
    message("No GSE accession numbers detected in the provided vector.")
    return(invisible(NULL))
  }
  
  message(sprintf("Found %d unique GSE accession(s): %s",
                  length(found), paste(found, collapse = ", ")))
  
  gse_metas <- lapply(found, extract_gse_meta, verbose = F)
  gse_metas <- setNames(gse_metas,found)
  results <- lapply(gse_metas, lookup_gse, verbose = verbose)
  res_df <- bind_rows(results, .id = "GSE_ID")
  return(res_df)
}


# -------------------------------------------------------------------------
# USAGE EXAMPLES
# -------------------------------------------------------------------------

# --- Example 1: annotate a column from your data frame ---
# Assuming your table is called `my_df` and the column is `study_id`:
#
# result <- annotate_gse_ids(my_df$study_id)
# print(result)
#
# To join the annotations back to your original table:
# my_df <- my_df %>%
#   mutate(gse_id_clean = toupper(regmatches(study_id,
#                                  regexpr("GSE[0-9]+", study_id,
#                                          ignore.case = TRUE)))) %>%
#   left_join(result, by = c("gse_id_clean" = "gse_id"))


# --- Example 2: quick test with known accessions ---
# test_ids <- c("GSE144735", "some_other_value", "GSE98411", NA, "not_a_gse")
# result   <- annotate_gse_ids(test_ids, verbose = TRUE)
# print(result)






# -------------------------------------------------------------------------
# Annotate a vector of strings using a named keyword list
# Returns a data frame with one or more columns (name_1, name_2, etc.)
# if multiple vector names match for a single entry
# -------------------------------------------------------------------------
annotate_from_keywords <- function(strings, keyword_list, col_name, ignore.case = TRUE) {
  
  # For each string, find all named vectors whose keywords match
  matches_per_entry <- lapply(strings, function(s) {
    if (is.na(s)) return(character(0))
    matched_names <- names(keyword_list)[sapply(keyword_list, function(keywords) {
      any(sapply(keywords, function(kw) grepl(kw, s, ignore.case = ignore.case, perl = TRUE)))
    })]
    matched_names
  })
  
  # Find the max number of matches across all entries to know how many columns to make
  max_matches <- max(lengths(matches_per_entry), na.rm = TRUE)
  max_matches <- max(max_matches, 1L)  # always at least 1 column
  
  # Build one column per match slot
  col_names <- if (max_matches == 1) {
    col_name
  } else {
    paste0(col_name, "_", seq_len(max_matches))
  }
  
  out <- as.data.frame(
    setNames(
      lapply(seq_len(max_matches), function(i) {
        sapply(matches_per_entry, function(m) if (length(m) >= i) m[[i]] else NA_character_)
      }),
      col_names
    ),
    stringsAsFactors = FALSE
  )
  out[,paste0(col_name,"_Hits")] <- sapply(matches_per_entry, function(m) length(m))
  
  out
}


# -------------------------------------------------------------------------
# Apply multiple annotation lists to a data frame
# -------------------------------------------------------------------------
# Parameters:
#   df            - your data frame of file paths/names
#   string_col    - name of the column in df containing the strings to match against
#   annotations   - a named list of keyword lists, e.g.:
#                   list(Source = source_list, Region = region_list, ...)
#   ignore.case   - passed through to grepl
#
# Returns:
#   df with additional annotation columns appended
# -------------------------------------------------------------------------
annotate_table <- function(df, string_col, annotations, ignore.case = TRUE) {
  
  strings <- df[[string_col]]
  
  new_cols <- lapply(names(annotations), function(ann_name) {
    annotate_from_keywords(
      strings      = strings,
      keyword_list = annotations[[ann_name]],
      col_name     = ann_name,
      ignore.case  = ignore.case
    )
  })
  
  cbind(df, do.call(cbind, new_cols))
}

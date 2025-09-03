library(dplyr)

# ============================================================================
# FUNCTION 1: HARMONIZE KIDNEY DATA WITH CUSTOM MERGE GROUPS
# ============================================================================

#' Harmonize Kidney Biomarker Data for Machine Learning (Manual Merge Groups)
#' 
#' This function performs comprehensive data harmonization on kidney biomarker data,
#' with the option to provide custom merge groups for duplicate/similar columns.
#' This preserves your original manual approach while automating other steps.
#' 
#' @param df A dataframe containing kidney biomarker data with STUDYID as first column
#' @param unwanted_cols Vector of column names to remove (default: body weight related columns)
#' @param zero_threshold Threshold for removing columns with high percentage of zeros (default: 0.9)
#' @param remove_unremarkable Logical, whether to remove "highest_score" and "UNREMARKABLE" columns
#' @param custom_merge_groups List of character vectors, each containing column names to merge
#'                           If NULL, uses default kidney biomarker merge groups
#' @param show_available_cols Logical, if TRUE prints available column names before merging
#' 
#' @return A harmonized dataframe ready for machine learning
#' 
harmonize_kidney_data <- function(df, 
                                  unwanted_cols = c("BWSTRESN", "BWSTRESN_Init", "finalbodyweight", "BWZSCORE"),
                                  zero_threshold = 0.9,
                                  remove_unremarkable = TRUE,
                                  custom_merge_groups = NULL,
                                  show_available_cols = FALSE) {
  
  # ============================================================================
  # STEP 1: INITIAL DATA CLEANING AND PREPARATION
  # ============================================================================
  
  cat("Step 1: Initial data cleaning...\n")
  
  # Identify numeric columns (excluding STUDYID which should be column 1)
  num_cols <- sapply(df, is.numeric)
  num_cols[1] <- FALSE  # Exclude STUDYID from numeric processing
  
  # Replace NA and NaN with 0 in numeric columns only
  df[num_cols] <- lapply(df[num_cols], function(x) {
    x[is.na(x) | is.nan(x)] <- 0
    return(x)
  })
  
  # Remove unwanted columns
  df <- df[, !(names(df) %in% unwanted_cols)]
  
  cat("  - Removed", length(intersect(names(df), unwanted_cols)), "unwanted columns\n")
  cat("  - Replaced NA/NaN values with 0 in numeric columns\n")
  
  # ============================================================================
  # STEP 2: BIOMARKER SCORE HARMONIZATION AND TRANSFORMATION
  # ============================================================================
  
  cat("Step 2: Applying score harmonization logic...\n")
  
  harmonization_step1 <- df %>%
    mutate(across(
      2:ncol(.),
      ~ case_when(
        . == 5 ~ 5,         # Preserve severe findings
        . > 3 ~ 3,          # Cap moderate-high findings
        . == 3 ~ 2,         # Downgrade exact 3 to moderate
        . > 0 ~ 1,          # Any positive finding becomes mild
        TRUE ~ 0            # Zero/negative becomes no finding
      )
    ))
  
  cat("  - Applied standardized 0-5 scoring transformation\n")
  
  # ============================================================================
  # STEP 3: HANDLE EXACT DUPLICATE COLUMN NAMES
  # ============================================================================
  
  cat("Step 3: Merging exact duplicate columns...\n")
  
  dup_colnames_logical <- duplicated(names(harmonization_step1)) | 
    duplicated(names(harmonization_step1), fromLast = TRUE)
  dup_colnames <- unique(names(harmonization_step1)[dup_colnames_logical])
  
  if (length(dup_colnames) > 0) {
    for (name in dup_colnames) {
      cols_with_name <- which(names(harmonization_step1) == name)
      harmonization_step1[[name]] <- apply(
        harmonization_step1[, cols_with_name, drop = FALSE], 
        1, max, na.rm = TRUE
      )
      harmonization_step1 <- harmonization_step1[, -cols_with_name[-1], drop = FALSE]
    }
    cat("  - Merged", length(dup_colnames), "sets of duplicate column names\n")
  } else {
    cat("  - No exact duplicate column names found\n")
  }
  
  # ============================================================================
  # STEP 4: SHOW AVAILABLE COLUMNS (FOR MANUAL MERGE GROUP CREATION)
  # ============================================================================
  
  if (show_available_cols) {
    cat("\n=== AVAILABLE COLUMN NAMES FOR MERGE GROUP CREATION ===\n")
    available_cols <- names(harmonization_step1)[-1]  # Exclude STUDYID
    cat("Total columns available:", length(available_cols), "\n")
    for (i in seq_along(available_cols)) {
      cat(sprintf("%3d: %s\n", i, available_cols[i]))
    }
    cat("========================================================\n\n")
  }
  
  # ============================================================================
  # STEP 5: MERGE SEMANTICALLY SIMILAR BIOMARKER GROUPS (CUSTOM OR DEFAULT)
  # ============================================================================
  
  cat("Step 4: Merging semantically similar biomarker columns...\n")
  
  # Use custom merge groups if provided, otherwise use default
  if (is.null(custom_merge_groups)) {
    # Default merge groups based on common kidney biomarker patterns
    merge_groups <- list(
      c("INFILTRATION, MONONUCLEAR CELL", "INFILTRATE"),
      c("INFLAMMATORY CELL FOCI", "INFLAMMATION"),
      c("DEGENERATION", "DEGENERATION/NECROSIS", "NECROSIS"),
      c("REGENERATION", "DEGENERATION/REGENERATION"),
      c("BASOPHILIA", "BASOPHILIA, TUBULE", "BASOPHILIC TUBULE"),
      c("CASTS, HYALINE/GRANULAR", "CAST"),
      c("ACCUMULATION, HYALINE DROPLETS", "ACCUMULATION"),
      c("PIGMENT, INCREASED", "PIGMENT"),
      c("PERL'S PRUSSIAN BLUE STAIN, EXAMINED", "PERL'S PRUSSIAN BLUE STAIN, POSITIVE"),
      c("SCHMORL'S STAIN, EXAMINED", "SCHMORL'S STAIN, POSITIVE"),
      c("CHRONIC PROGRESSIVE NEPHROPATHY", "NEPHROPATHY")
    )
    cat("  - Using default kidney biomarker merge groups\n")
  } else {
    merge_groups <- custom_merge_groups
    cat("  - Using custom merge groups (", length(merge_groups), "groups provided)\n")
  }
  
  merged_groups_count <- 0
  
  for (group in merge_groups) {
    present <- group[group %in% names(harmonization_step1)]
    if (length(present) <= 1) next
    
    shortest_name <- present[which.min(nchar(present))]
    harmonization_step1[[shortest_name]] <- apply(
      harmonization_step1[, present, drop = FALSE], 
      1, max, na.rm = TRUE
    )
    harmonization_step1 <- harmonization_step1[, 
                                               !(names(harmonization_step1) %in% setdiff(present, shortest_name)), 
                                               drop = FALSE
    ]
    merged_groups_count <- merged_groups_count + 1
  }
  
  cat("  - Successfully merged", merged_groups_count, "semantic biomarker groups\n")
  
  # ============================================================================
  # STEPS 6-7: REMAINING CLEANUP (same as before)
  # ============================================================================
  
  cat("Step 5: Removing non-informative columns...\n")
  
  if (remove_unremarkable) {
    cols_to_remove <- c("highest_score", "UNREMARKABLE")
    harmonization_step1 <- harmonization_step1[, 
                                               !names(harmonization_step1) %in% cols_to_remove, 
                                               drop = FALSE
    ]
    cat("  - Removed 'highest_score' and 'UNREMARKABLE' columns\n")
  }
  
  initial_cols <- ncol(harmonization_step1)
  harmonization_step1 <- harmonization_step1[, 
                                             colSums(harmonization_step1 != 0, na.rm = TRUE) > 0, 
                                             drop = FALSE
  ]
  removed_zero_cols <- initial_cols - ncol(harmonization_step1)
  
  if (removed_zero_cols > 0) {
    cat("  - Removed", removed_zero_cols, "columns with all zero values\n")
  }
  
  cat("Step 6: Applying sparsity-based feature selection...\n")
  
  initial_feature_cols <- ncol(harmonization_step1)
  zero_percentages <- colMeans(harmonization_step1 == 0, na.rm = TRUE)
  cols_to_keep <- (zero_percentages <= zero_threshold) | (names(harmonization_step1) == "STUDYID")
  harmonization_step1 <- harmonization_step1[, cols_to_keep, drop = FALSE]
  
  removed_sparse_cols <- initial_feature_cols - ncol(harmonization_step1)
  cat("  - Removed", removed_sparse_cols, "columns with >=", 
      round(zero_threshold * 100, 1), "% zero values\n")
  
  cat("Step 7: Cleaning column names...\n")
  
  original_names <- names(harmonization_step1)
  names(harmonization_step1) <- gsub("[^A-Za-z0-9_]", "", names(harmonization_step1))
  
  modified_names <- sum(original_names != names(harmonization_step1))
  if (modified_names > 0) {
    cat("  - Cleaned", modified_names, "column names (removed special characters)\n")
  }
  
  # Final summary
  cat("\n=== HARMONIZATION SUMMARY ===\n")
  cat("Final dataset dimensions:", nrow(harmonization_step1), "rows x", 
      ncol(harmonization_step1), "columns\n")
  cat("Biomarker features (excluding STUDYID):", ncol(harmonization_step1) - 1, "\n")
  
  non_zero_percentage <- mean(harmonization_step1[, -1] != 0, na.rm = TRUE) * 100
  cat("Data density (non-zero values):", round(non_zero_percentage, 1), "%\n")
  cat("Harmonization completed successfully!\n")
  cat("==============================\n\n")
  
  return(harmonization_step1)
}

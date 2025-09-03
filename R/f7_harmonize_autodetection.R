# ============================================================================
# FUNCTION 2: AUTO-DUPLICATE DETECTION HARMONIZATION
# ============================================================================

#' Harmonize Kidney Data with Automatic Duplicate Detection
#' 
#' This function automatically detects similar column names using string similarity
#' and suggests/applies merges without requiring manual merge group creation.
#' 
#' @param df A dataframe containing kidney biomarker data
#' @param unwanted_cols Vector of column names to remove
#' @param zero_threshold Threshold for removing sparse columns
#' @param remove_unremarkable Logical, remove specific non-informative columns
#' @param similarity_threshold Minimum string similarity (0-1) to consider columns similar (default: 0.6)
#' @param min_common_words Minimum number of common words required for merging (default: 1)
#' @param show_detected_groups Logical, if TRUE shows detected similar groups before merging
#' @param auto_merge Logical, if TRUE automatically merges detected groups, if FALSE only shows suggestions
#' 
#' @return A harmonized dataframe with automatically detected duplicates merged
#' 
harmonize_kidney_auto_detect <- function(df,
                                         unwanted_cols = c("BWSTRESN", "BWSTRESN_Init", "finalbodyweight", "BWZSCORE"),
                                         zero_threshold = 0.9,
                                         remove_unremarkable = TRUE,
                                         similarity_threshold = 0.6,
                                         min_common_words = 1,
                                         show_detected_groups = TRUE,
                                         auto_merge = TRUE) {
  
  # Helper function to calculate string similarity (Jaccard similarity)
  jaccard_similarity <- function(str1, str2) {
    # Convert to lowercase and split into words
    words1 <- unlist(strsplit(tolower(str1), "[^a-z0-9]+"))
    words2 <- unlist(strsplit(tolower(str2), "[^a-z0-9]+"))
    
    # Remove empty strings
    words1 <- words1[words1 != ""]
    words2 <- words2[words2 != ""]
    
    # Calculate Jaccard similarity
    intersection <- length(intersect(words1, words2))
    union <- length(union(words1, words2))
    
    if (union == 0) return(0)
    return(intersection / union)
  }
  
  # Helper function to count common words
  count_common_words <- function(str1, str2) {
    words1 <- unlist(strsplit(tolower(str1), "[^a-z0-9]+"))
    words2 <- unlist(strsplit(tolower(str2), "[^a-z0-9]+"))
    words1 <- words1[words1 != ""]
    words2 <- words2[words2 != ""]
    return(length(intersect(words1, words2)))
  }
  
  # ============================================================================
  # STEPS 1-3: Same initial processing as Function 1
  # ============================================================================
  
  cat("Step 1: Initial data cleaning...\n")
  
  num_cols <- sapply(df, is.numeric)
  num_cols[1] <- FALSE
  
  df[num_cols] <- lapply(df[num_cols], function(x) {
    x[is.na(x) | is.nan(x)] <- 0
    return(x)
  })
  
  df <- df[, !(names(df) %in% unwanted_cols)]
  cat("  - Removed", length(intersect(names(df), unwanted_cols)), "unwanted columns\n")
  
  cat("Step 2: Applying score harmonization logic...\n")
  
  harmonization_step1 <- df %>%
    mutate(across(
      2:ncol(.),
      ~ case_when(
        . == 5 ~ 5, . > 3 ~ 3, . == 3 ~ 2, . > 0 ~ 1, TRUE ~ 0
      )
    ))
  
  cat("Step 3: Merging exact duplicate columns...\n")
  
  dup_colnames_logical <- duplicated(names(harmonization_step1)) | 
    duplicated(names(harmonization_step1), fromLast = TRUE)
  dup_colnames <- unique(names(harmonization_step1)[dup_colnames_logical])
  
  if (length(dup_colnames) > 0) {
    for (name in dup_colnames) {
      cols_with_name <- which(names(harmonization_step1) == name)
      harmonization_step1[[name]] <- apply(
        harmonization_step1[, cols_with_name, drop = FALSE], 1, max, na.rm = TRUE
      )
      harmonization_step1 <- harmonization_step1[, -cols_with_name[-1], drop = FALSE]
    }
    cat("  - Merged", length(dup_colnames), "sets of duplicate column names\n")
  }
  
  # ============================================================================
  # STEP 4: AUTOMATIC DUPLICATE DETECTION AND MERGING
  # ============================================================================
  
  cat("Step 4: Automatic duplicate detection...\n")
  
  # Get column names excluding STUDYID
  col_names <- names(harmonization_step1)[-1]
  
  # Find similar column groups
  detected_groups <- list()
  processed_cols <- character(0)
  
  for (i in seq_along(col_names)) {
    if (col_names[i] %in% processed_cols) next
    
    current_group <- col_names[i]
    
    # Find similar columns
    for (j in (i+1):length(col_names)) {
      if (j > length(col_names)) break
      if (col_names[j] %in% processed_cols) next
      
      # Calculate similarity
      similarity <- jaccard_similarity(col_names[i], col_names[j])
      common_words <- count_common_words(col_names[i], col_names[j])
      
      # Check if columns are similar enough to merge
      if (similarity >= similarity_threshold && common_words >= min_common_words) {
        current_group <- c(current_group, col_names[j])
      }
    }
    
    # Only add groups with more than one column
    if (length(current_group) > 1) {
      detected_groups[[length(detected_groups) + 1]] <- current_group
      processed_cols <- c(processed_cols, current_group)
    }
  }
  
  cat("  - Detected", length(detected_groups), "potential merge groups\n")
  
  # Show detected groups if requested
  if (show_detected_groups && length(detected_groups) > 0) {
    cat("\n=== DETECTED SIMILAR COLUMN GROUPS ===\n")
    for (i in seq_along(detected_groups)) {
      cat("Group", i, ":\n")
      for (col in detected_groups[[i]]) {
        similarity_scores <- sapply(detected_groups[[i]], function(x) {
          if (x == col) return("")
          return(paste0("(sim:", round(jaccard_similarity(col, x), 2), ")"))
        })
        cat("  -", col, paste(similarity_scores[similarity_scores != ""], collapse = " "), "\n")
      }
      cat("  -> Will merge as:", detected_groups[[i]][which.min(nchar(detected_groups[[i]]))], "\n\n")
    }
    cat("=====================================\n")
  }
  
  # Apply merging if auto_merge is TRUE
  if (auto_merge && length(detected_groups) > 0) {
    merged_count <- 0
    
    for (group in detected_groups) {
      # Check if all columns in group still exist
      present <- group[group %in% names(harmonization_step1)]
      if (length(present) <= 1) next
      
      # Use shortest name as merged column name
      shortest_name <- present[which.min(nchar(present))]
      
      # Merge by taking row-wise maximum
      harmonization_step1[[shortest_name]] <- apply(
        harmonization_step1[, present, drop = FALSE], 1, max, na.rm = TRUE
      )
      
      # Remove other columns in the group
      harmonization_step1 <- harmonization_step1[, 
                                                 !(names(harmonization_step1) %in% setdiff(present, shortest_name)), 
                                                 drop = FALSE
      ]
      
      merged_count <- merged_count + 1
    }
    
    cat("  - Successfully merged", merged_count, "detected groups\n")
    
  } else if (!auto_merge && length(detected_groups) > 0) {
    cat("  - Auto-merge disabled. Use detected groups above to create custom_merge_groups\n")
  }
  
  # ============================================================================
  # STEPS 5-7: Same cleanup as Function 1
  # ============================================================================
  
  # [Rest of the cleanup steps - same as Function 1]
  cat("Step 5: Removing non-informative columns...\n")
  
  if (remove_unremarkable) {
    cols_to_remove <- c("highest_score", "UNREMARKABLE")
    harmonization_step1 <- harmonization_step1[, 
                                               !names(harmonization_step1) %in% cols_to_remove, drop = FALSE
    ]
  }
  
  initial_cols <- ncol(harmonization_step1)
  harmonization_step1 <- harmonization_step1[, 
                                             colSums(harmonization_step1 != 0, na.rm = TRUE) > 0, drop = FALSE
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
    cat("  - Cleaned", modified_names, "column names\n")
  }
  
  # Final summary
  cat("\n=== AUTO-DETECTION HARMONIZATION SUMMARY ===\n")
  cat("Final dataset dimensions:", nrow(harmonization_step1), "rows x", 
      ncol(harmonization_step1), "columns\n")
  cat("Biomarker features:", ncol(harmonization_step1) - 1, "\n")
  cat("Harmonization completed successfully!\n")
  cat("============================================\n\n")
  
  return(harmonization_step1)
}


# ============================================================================
# USAGE EXAMPLES
# ============================================================================

# FUNCTION 1 EXAMPLES:

# Basic usage (your original approach)
# harmonized_data <- harmonize_kidney_data(df)

# Show available columns to help create custom merge groups
# harmonized_data <- harmonize_kidney_data(df, show_available_cols = TRUE)

# With custom merge groups (your manual Step B approach)
# custom_groups <- list(
#   c("INFILTRATION, MONONUCLEAR CELL", "INFILTRATE"),
#   c("DEGENERATION", "NECROSIS"),
#   c("BASOPHILIA", "BASOPHILIA, TUBULE")
# )
# harmonized_data <- harmonize_kidney_data(df, custom_merge_groups = custom_groups)


# FUNCTION 2 EXAMPLES:

# Auto-detect and show suggestions only (don't merge yet)
# harmonized_data <- harmonize_kidney_auto_detect(df, auto_merge = FALSE, show_detected_groups = TRUE)

# Auto-detect and merge automatically
# harmonized_data <- harmonize_kidney_auto_detect(df)

# More sensitive detection (lower threshold)
# harmonized_data <- harmonize_kidney_auto_detect(df, similarity_threshold = 0.4, min_common_words = 1)
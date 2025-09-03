#' @title process_kidney_scores
#'
#' @description
#' This function processes kidney organ toxicity scores: orchestrates BW, Kidney-to-BW, LB, and MI score calculations for the kidney domain.
#' External function dependencies: Assumes external functions loaded in the package:
#'                        get_compile_data(), get_bw_score(), kidneytobw_score(), kidney_lb_score(), kidney_mi_score()
#' for a set of studies or XPT files. It can output individual scores, z-scores by USUBJID, or averaged scores
#' for multiple studies, and handles errors during the processing steps.
#'
#' @param studyid_or_studyids A character vector or a single study ID to process.
#' If multiple studies (for machine learning) are provided, the function processes each study sequentially. (Mandatory)
#'
#' @param path_db A character string specifying the path to the database or directory containing the data files. (Mandatory)
#'
#' @param fake_study A boolean flag indicating if the study data is simulated (TRUE) or real (FALSE). Default is FALSE. (Optional)
#'
#' @param use_xpt_file A boolean flag indicating whether to use an XPT file for the study data. Default is FALSE. (Optional)
#'
#' @param output_individual_scores A boolean flag indicating whether individual scores should be returned (TRUE) or averaged scores (FALSE). Default is FALSE. (Optional)
#'
#' @param output_zscore_by_USUBJID A boolean flag indicating whether to output z-scores by USUBJID (TRUE) or averaged scores (FALSE). Default is FALSE. (Optional)
#'
#' @return A data frame containing the calculated scores for each study. The type of result depends on the flags passed:
#' - If output_individual_scores is TRUE, a data frame with individual scores for each study is returned.
#' - If output_zscore_by_USUBJID is TRUE, a data frame with z-scores by USUBJID for each study is returned.
#' - If neither flag is set, the function returns a data frame with averaged scores for each study.
#'
#' @examples
#' \dontrun{
#' # Get averaged scores for a single study
#' result_default <- process_kidney_scores("Study123", "path/to/database")
#'
#' # Get individual scores for multiple studies
#' result_individual <- process_kidney_scores(
#'   studyid_or_studyids = c("Study123", "Study456"),
#'   path_db = "path/to/database",
#'   output_individual_scores = TRUE
#' )
#'
#' # Get z-scores by USUBJID
#' result_zscore <- process_kidney_scores(
#'   studyid_or_studyids = "Study123",
#'   path_db = "path/to/database", 
#'   output_zscore_by_USUBJID = TRUE
#' )
#' }
#'
#' @export
Revised_process_kidney_scores <- function(studyid_or_studyids, path_db, fake_study = FALSE,
                                  use_xpt_file = FALSE, output_individual_scores = FALSE,
                                  output_zscore_by_USUBJID = FALSE) {
  
  # Load required libraries
  library(dplyr)
  library(tibble)
  
  # Helper function: safe error handling without global assignment
  safe_execute <- function(expr, studyid, block, use_xpt_file, path_db) {
    tryCatch(
      list(result = suppressWarnings(expr), error = NULL),
      error = function(e) {
        error_record <- tibble(
          STUDYID = if (use_xpt_file) path_db else studyid,
          Block = block,
          ErrorMessage = e$message
        )
        list(result = NULL, error = error_record)
      }
    )
  }
  
  # Nested function: calculate Kidney:BW ratio
  calculate_kidneyToBW_zscore <- function(studyid, path_db, fake_study, use_xpt_file, master_compiledata,
                                          output_individual_scores, output_zscore_by_USUBJID,
                                          bwzscore_BW) {
    studyid_local <- if (use_xpt_file) NULL else studyid
    
    if (output_zscore_by_USUBJID) {
      exec_result <- safe_execute(
        kidneytobw_score(
          studyid = studyid_local,
          path_db = path_db,
          fake_study = fake_study,
          use_xpt_file = use_xpt_file,
          master_compiledata = master_compiledata,
          bwzscore_BW = bwzscore_BW,
          return_individual_scores = FALSE,
          return_zscore_by_USUBJID = TRUE
        ),
        studyid, 'KidneyToBW', use_xpt_file, path_db
      )
      
      if (!is.null(exec_result$result) && !('USUBJID' %in% names(exec_result$result))) {
        exec_result$result <- exec_result$result %>% mutate(USUBJID = NA_character_)
      }
      return(list(result = exec_result$result, error = exec_result$error))
      
    } else if (!output_individual_scores) {
      exec_result <- safe_execute({
        avg_df <- kidneytobw_score(
          studyid = studyid_local,
          path_db = path_db,
          fake_study = fake_study,
          use_xpt_file = use_xpt_file,
          master_compiledata = master_compiledata,
          bwzscore_BW = bwzscore_BW,
          return_individual_scores = FALSE,
          return_zscore_by_USUBJID = FALSE
        )
        dplyr::rename(avg_df, kidneyToBW_avg = avg_kidneyToBW_zscore)$kidneyToBW_avg[1]
      }, studyid, 'KidneyToBW', use_xpt_file, path_db)
      
      return(list(result = exec_result$result, error = exec_result$error))
    } else {
      return(list(result = NULL, error = NULL))
    }
  }
  
  # Input validation
  if (missing(studyid_or_studyids) || length(studyid_or_studyids) == 0) {
    stop("studyid_or_studyids is required and cannot be empty")
  }
  if (missing(path_db)) {
    stop("path_db is required")
  }
  if (output_individual_scores && output_zscore_by_USUBJID) {
    stop('Cannot request both individual and USUBJID-level scores.')
  }
  
  # Initialize containers
  if (output_individual_scores) {
    master_bw <- tibble()
    master_lb <- tibble()
    master_mi <- tibble()
  }
  if (output_zscore_by_USUBJID) {
    master_usubjid_bw <- tibble()
    master_usubjid_lb <- tibble()
    master_usubjid_mi <- tibble()
  }
  
  FOUR_Kidney_Score_avg <- tibble(
    STUDYID = character(),
    BWZSCORE_avg = numeric(),
    kidneyToBW_avg = numeric(),
    LB_score_avg = numeric(),
    MI_score_avg = numeric()
  )
  
  all_errors <- tibble(STUDYID = character(), Block = character(), ErrorMessage = character())
  
  # Process each study
  for (study in studyid_or_studyids) {
    if (use_xpt_file) path_db <- study
    message('Processing study: ', study)
    
    # Get compiled data
    compile_result <- safe_execute(
      get_compile_data(
        studyid = if (use_xpt_file) NULL else study,
        path_db = path_db,
        fake_study = fake_study,
        use_xpt_file = use_xpt_file
      ),
      study, 'compiledata', use_xpt_file, path_db
    )
    
    if (!is.null(compile_result$error)) {
      all_errors <- bind_rows(all_errors, compile_result$error)
    }
    
    if (is.null(compile_result$result)) {
      next
    }
    
    mc <- compile_result$result
    
    # Initialize summary row for averaged scores mode
    if (!output_individual_scores && !output_zscore_by_USUBJID) {
      FOUR_Kidney_Score_avg <- bind_rows(
        FOUR_Kidney_Score_avg,
        tibble(
          STUDYID = as.character(unique(mc$STUDYID)),
          BWZSCORE_avg = NA_real_,
          kidneyToBW_avg = NA_real_,
          LB_score_avg = NA_real_,
          MI_score_avg = NA_real_
        )
      )
    }
    
    # Process BW z-scores
    bw_result <- safe_execute({
      sid <- if (use_xpt_file) NULL else study
      get_bw_score(
        studyid = sid,
        path_db = path_db,
        fake_study = fake_study,
        use_xpt_file = use_xpt_file,
        master_compiledata = mc,
        return_individual_scores = output_individual_scores,
        return_zscore_by_USUBJID = output_zscore_by_USUBJID
      )
    }, study, 'BWZscore', use_xpt_file, path_db)
    
    if (!is.null(bw_result$error)) {
      all_errors <- bind_rows(all_errors, bw_result$error)
    }
    
    bw_res <- bw_result$result
    
    if (output_individual_scores && !is.null(bw_res)) {
      master_bw <- bind_rows(master_bw, bw_res)
    }
    if (output_zscore_by_USUBJID && !is.null(bw_res)) {
      if (!('USUBJID' %in% names(bw_res))) bw_res <- bw_res %>% mutate(USUBJID = NA_character_)
      master_usubjid_bw <- bind_rows(master_usubjid_bw, bw_res %>% mutate(USUBJID = as.character(USUBJID)))
    }
    if (!output_individual_scores && !output_zscore_by_USUBJID && !is.null(bw_res)) {
      FOUR_Kidney_Score_avg$BWZSCORE_avg[FOUR_Kidney_Score_avg$STUDYID == unique(mc$STUDYID)] <- bw_res$BWZSCORE_avg[1]
    }
    
    # Process LB scores
    lb_result <- safe_execute({
      sid <- if (use_xpt_file) NULL else study
      kidney_lb_score(
        studyid = sid,
        path_db = path_db,
        fake_study = fake_study,
        use_xpt_file = use_xpt_file,
        master_compiledata = mc,
        return_individual_scores = output_individual_scores,
        return_zscore_by_USUBJID = output_zscore_by_USUBJID
      )
    }, study, 'LB', use_xpt_file, path_db)
    
    if (!is.null(lb_result$error)) {
      all_errors <- bind_rows(all_errors, lb_result$error)
    }
    
    lb_res <- lb_result$result
    
    if (output_individual_scores && !is.null(lb_res)) {
      master_lb <- bind_rows(master_lb, lb_res)
    }
    if (output_zscore_by_USUBJID && !is.null(lb_res)) {
      if (!('USUBJID' %in% names(lb_res))) lb_res <- lb_res %>% mutate(USUBJID = NA_character_)
      master_usubjid_lb <- bind_rows(master_usubjid_lb, lb_res %>% mutate(USUBJID = as.character(USUBJID)))
    }
    if (!output_individual_scores && !output_zscore_by_USUBJID && !is.null(lb_res)) {
      FOUR_Kidney_Score_avg$LB_score_avg[FOUR_Kidney_Score_avg$STUDYID == unique(mc$STUDYID)] <- lb_res$LB_score_avg[1]
    }
    
    # Process MI scores
    # For individual mode, we need to get the detailed USUBJID data first, then aggregate it
    # This ensures we capture all MI biomarker columns before aggregation
    mi_result <- safe_execute({
      sid <- if (use_xpt_file) NULL else study
      if (output_individual_scores) {
        # Get detailed USUBJID data first to capture all MI columns
        detailed_mi <- kidney_mi_score(
          studyid = sid,
          path_db = path_db,
          fake_study = fake_study,
          use_xpt_file = use_xpt_file,
          master_compiledata = mc,
          return_individual_scores = FALSE,
          return_zscore_by_USUBJID = TRUE
        )
        
        # Then aggregate to study level while preserving all columns
        if (!is.null(detailed_mi) && nrow(detailed_mi) > 0) {
          detailed_mi %>%
            group_by(STUDYID) %>%
            summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')
        } else {
          detailed_mi
        }
      } else {
        # Use the original call for other modes
        kidney_mi_score(
          studyid = sid,
          path_db = path_db,
          fake_study = fake_study,
          use_xpt_file = use_xpt_file,
          master_compiledata = mc,
          return_individual_scores = output_individual_scores,
          return_zscore_by_USUBJID = output_zscore_by_USUBJID
        )
      }
    }, study, 'MI', use_xpt_file, path_db)
    
    if (!is.null(mi_result$error)) {
      all_errors <- bind_rows(all_errors, mi_result$error)
    }
    
    mi_res <- mi_result$result
    
    if (output_individual_scores && !is.null(mi_res)) {
      master_mi <- bind_rows(master_mi, mi_res)
    }
    if (output_zscore_by_USUBJID && !is.null(mi_res)) {
      if (!('USUBJID' %in% names(mi_res))) mi_res <- mi_res %>% mutate(USUBJID = NA_character_)
      master_usubjid_mi <- bind_rows(master_usubjid_mi, mi_res %>% mutate(USUBJID = as.character(USUBJID)))
    }
    if (!output_individual_scores && !output_zscore_by_USUBJID && !is.null(mi_res)) {
      FOUR_Kidney_Score_avg$MI_score_avg[FOUR_Kidney_Score_avg$STUDYID == unique(mc$STUDYID)] <- mi_res$MI_score_avg[1]
    }
    
    # Process Kidney:BW ratio for summary mode
    if (!output_individual_scores && !output_zscore_by_USUBJID) {
      ktbw_result <- calculate_kidneyToBW_zscore(
        study, path_db, fake_study, use_xpt_file, mc,
        output_individual_scores = FALSE,
        output_zscore_by_USUBJID = FALSE,
        bwzscore_BW = NULL
      )
      
      if (!is.null(ktbw_result$error)) {
        all_errors <- bind_rows(all_errors, ktbw_result$error)
      }
      
      if (!is.null(ktbw_result$result)) {
        FOUR_Kidney_Score_avg$kidneyToBW_avg[FOUR_Kidney_Score_avg$STUDYID == unique(mc$STUDYID)] <- ktbw_result$result
      }
    }
  }
  
  # Print error summary if any errors occurred
  if (nrow(all_errors) > 0) {
    warning(paste("Errors occurred in", nrow(all_errors), "processing steps. Check the returned error log."))
    attr(FOUR_Kidney_Score_avg, "errors") <- all_errors
  }
  
  # Return results based on requested mode
  if (output_individual_scores) {
    # Aggregate to preserve all biomarker columns instead of just taking first row
    # This ensures all kidney biomarkers are included, not lost due to sparse matrix structure
    master_bw1 <- if (nrow(master_bw) > 0 && "STUDYID" %in% names(master_bw)) {
      master_bw %>% 
        group_by(STUDYID) %>% 
        summarise(across(where(is.numeric), ~ first(na.omit(.x))), .groups = 'drop')
    } else {
      tibble(STUDYID = character())
    }
    
    master_lb1 <- if (nrow(master_lb) > 0 && "STUDYID" %in% names(master_lb)) {
      master_lb %>% 
        group_by(STUDYID) %>% 
        summarise(across(where(is.numeric), ~ first(na.omit(.x))), .groups = 'drop')
    } else {
      tibble(STUDYID = character())
    }
    
    master_mi1 <- if (nrow(master_mi) > 0 && "STUDYID" %in% names(master_mi)) {
      master_mi %>% 
        group_by(STUDYID) %>% 
        summarise(across(where(is.numeric), ~ first(na.omit(.x))), .groups = 'drop')
    } else {
      tibble(STUDYID = character())
    }
    
    # Join the results, handling empty data frames gracefully
    result_individual <- master_bw1
    
    if (nrow(master_lb1) > 0 && "STUDYID" %in% names(master_lb1)) {
      result_individual <- result_individual %>% full_join(master_lb1, by = "STUDYID")
    }
    
    if (nrow(master_mi1) > 0 && "STUDYID" %in% names(master_mi1)) {
      result_individual <- result_individual %>% full_join(master_mi1, by = "STUDYID")
    }
    
    if (nrow(all_errors) > 0) {
      attr(result_individual, "errors") <- all_errors
    }
    
    return(result_individual)
    
  } else if (output_zscore_by_USUBJID) {
    combined_df <- master_usubjid_bw %>%
      full_join(master_usubjid_lb, by = c("STUDYID", "USUBJID")) %>%
      full_join(master_usubjid_mi, by = c("STUDYID", "USUBJID"))
    
    if (nrow(all_errors) > 0) {
      attr(combined_df, "errors") <- all_errors
    }
    
    return(combined_df)
    
  } else {
    FOUR_Kidney_Score_avg <- FOUR_Kidney_Score_avg %>%
      mutate(across(-STUDYID, ~ round(.x, 2)))
    return(FOUR_Kidney_Score_avg)
  }
}
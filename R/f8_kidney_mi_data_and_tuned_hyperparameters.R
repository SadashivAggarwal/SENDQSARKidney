#' @title Get Random Forest Data and Tuned Hyperparameters
#'
#' @description
#' The `get_ml_data_and_tuned_hyperparameters` function processes input data and metadata to prepare data for
#' random forest analysis. It includes steps for data preprocessing, optional imputation, rounding,
#' error correction, and hyperparameter tuning.
#'
#' @param Data data.frame. Input data frame containing scores, typically named `scores_df`.
#' First column is "STUDYID", followed by columns with score values.
#' @param studyid_metadata data.frame. Metadata containing `STUDYID` and `Target_Organ`.
#' @param Impute logical. Indicates whether to impute missing values in the dataset using random forest imputation. Default is `FALSE`.
#' @param Round logical. Specifies whether to round specific numerical columns according to predefined rules. Default is `FALSE`.
#' @param reps integer. Number of repetitions for cross-validation. A value of `0` skips repetition.
#' @param holdback numeric. Fraction of data to hold back for testing. A value of `1` performs leave-one-out cross-validation.
#' @param Undersample logical. Indicates whether to undersample the training data to balance the target classes. Default is `FALSE`.
#' @param hyperparameter_tuning logical. Specifies whether to perform hyperparameter tuning for the random forest model. Default is `FALSE`.
#' @param error_correction_method character. Specifies the method for error correction. Can be `"Flip"`, `"Prune"`, or `None`. Default is `None`.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{rfData}{The final processed data after preprocessing and error correction.}
#'   \item{best.m}{The best `mtry` hyperparameter determined for the random forest model.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' Data <- scores_df
#' studyid_metadata <- read.csv("path/to/study_metadata.csv")
#' result <- get_ml_data_and_tuned_hyperparameters(
#'   Data = f7_Data_result,
#'   studyid_metadata = studyid_metadata,
#'   Impute = TRUE,
#'   Round = TRUE,
#'   reps = 10,
#'   holdback = 0.75,
#'   Undersample = TRUE,
#'   hyperparameter_tuning = TRUE,
#'   error_correction_method = "Flip"
#' )
#' rfData <- result$rfData
#' best_mtry <- result$best.m
#' }

revised_f8_ml_data_and_tuned_hyperparameters <- function(Data,
                                                               studyid_metadata,
                                                               Impute = FALSE,
                                                               Round = FALSE,
                                                               reps,
                                                               holdback,
                                                               Undersample = FALSE,
                                                               hyperparameter_tuning = FALSE,
                                                               error_correction_method = "None") {
  # 1. Merge data and metadata
  merged_Data <- dplyr::left_join(
    studyid_metadata %>% dplyr::mutate(STUDYID = as.character(STUDYID)),
    Data,
    by = "STUDYID"
  )
  
  # 2. Prepare rfData
  rfData <- merged_Data[, -1]
  rfData$Target_Organ <- ifelse(rfData$Target_Organ == "Liver", 1, 0)
  rfData$Target_Organ <- factor(rfData$Target_Organ, levels = c(1, 0))
  
  # Diagnostic checks
  message("NAs in merged_Data$Target_Organ: ", sum(is.na(merged_Data$Target_Organ)))
  message("Unique values in merged_Data$Target_Organ: ", paste(unique(merged_Data$Target_Organ), collapse = ", "))
  message("NAs in rfData$Target_Organ after binarization: ", sum(is.na(rfData$Target_Organ)))
  
  # Remove NAs
  rfData <- stats::na.omit(rfData)
  
  # 3. Impute missing values if requested
  if (Impute) {
    rfData <- randomForest::rfImpute(Target_Organ ~ ., rfData)
  }
  
  # 4. Round and cap values if requested
  if (Round) {
    # Floor and cap z-score and liver columns
    zscoreIndex <- c(
      grep('avg_', colnames(rfData)),
      grep('liver', colnames(rfData))
    )
    for (i in zscoreIndex) {
      rfData[, i] <- floor(rfData[, i])
      rfData[rfData[, i] > 5, i] <- 5
    }
    # Ceiling for histopathology columns (excluding Target_Organ)
    histoIndex <- setdiff(
      which(substr(colnames(rfData), 1, 1) %in% LETTERS),
      which(colnames(rfData) == "Target_Organ")
    )
    for (i in histoIndex) {
      rfData[, i] <- ceiling(rfData[, i])
    }
  }
  
  # 5. Cross-validation and processing loop
  best.m <- NULL
  for (rep in seq_len(reps)) {
    message(rep / reps * 100, "% Complete...")
    
    # Split into training and test sets
    n <- nrow(rfData)
    if (holdback == 1) {
      ind <- sample(n, 1)
      train_indices <- setdiff(seq_len(n), ind)
      test_indices  <- ind
    } else {
      prob <- c(1 - holdback, holdback)
      grp <- sample(c(1, 2), n, replace = TRUE, prob = prob)
      train_indices <- which(grp == 1)
      test_indices  <- which(grp == 2)
    }
    train <- rfData[train_indices, , drop = FALSE]
    test  <- rfData[test_indices,  , drop = FALSE]
    train <- stats::na.omit(train)
    test  <- stats::na.omit(test)
    
    # 5c. Optional undersampling
    if (Undersample) {
      pos <- which(train$Target_Organ == 1)
      neg <- which(train$Target_Organ == 0)
      nPos <- length(pos)
      negSample <- if (length(neg) >= nPos) sample(neg, nPos) else neg
      train <- train[c(pos, negSample), , drop = FALSE]
    }
    
    # 5d. Hyperparameter tuning or default
    if (hyperparameter_tuning) {
      ctrl <- caret::trainControl(method = "repeatedcv", number = 10, repeats = 3)
      rfFit <- caret::train(
        Target_Organ ~ ., data = train,
        method = "rf", trControl = ctrl, tuneLength = 15
      )
      best.m <- rfFit$bestTune$mtry
    } else {
      best.m <- 4
    }
    
    # 5e. Train Random Forest model
    rf_model <- randomForest::randomForest(
      Target_Organ ~ ., data = train, mtry = best.m, ntree = 500
    )
    
    # 5f. Optional error correction stub
    if (error_correction_method != "None") {
      message("Error correction (", error_correction_method, ") not implemented.")
      rfData <- test  # placeholder
    }
  }
  
  # 6. Return processed data and best mtry
  list(rfData = rfData, best.m = best.m)
}

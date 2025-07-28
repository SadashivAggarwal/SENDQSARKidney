#' @title Random Forest Model with Cross-validation and Exclusion
#' @description Implements a Random Forest classification model with cross-validation,
#' optionally undersampling, handling indeterminate predictions, and calculating
#' sensitivity, specificity, PPV, NPV, prevalence, balanced accuracy, and indeterminate proportions across multiple repetitions.
#' @details
#' * **Sensitivity** (True Positive Rate): proportion of actual positive cases (Target_Organ = 1) correctly classified.
#' * **Specificity** (True Negative Rate): proportion of actual negative cases (Target_Organ = 0) correctly classified.
#' * **PPV** (Positive Predictive Value): among all cases predicted positive, the proportion that are truly positive.
#' * **NPV** (Negative Predictive Value): among all cases predicted negative, the proportion that are truly negative.
#' * **Prevalence**: proportion of positive cases in the test dataset, indicating baseline class frequency.
#' * **Balanced Accuracy**: average of sensitivity and specificity, mitigating effects of class imbalance.
#' * **Indeterminate Proportion** (`nRemoved`): fraction of predictions flagged as indeterminate (NA) when probabilities fall within the specified bounds.
#'
#' @param scores_data_df Data frame containing predictor features and the target column `Target_Organ`.
#' @param Undersample Logical; if TRUE, undersamples the majority class in each training fold to match the minority class.
#' @param best.m Numeric; number of variables to randomly sample at each split (`mtry`).
#' @param testReps Integer; number of cross-validation repetitions (must be >= 2).
#' @param indeterminateUpper Numeric; upper threshold for marking a prediction as indeterminate.
#' @param indeterminateLower Numeric; lower threshold for marking a prediction as indeterminate.
#' @param Type Integer; importance type: 1 for Mean Decrease Accuracy, 2 for Mean Decrease Gini (impurity).
#' @import randomForest caret
#' @return Invisibly returns a list containing:
#'   \describe{
#'     \item{performance_summary}{Named numeric vector of averaged metrics across repetitions.}
#'     \item{raw_results}{List of metric vectors and indeterminate counts per fold.}
#'   }
revised_f10_zone_exclusioned_rf_model_with_cv <- function(scores_data_df,
                                                  Undersample = FALSE,
                                                  best.m = NULL,
                                                  testReps,
                                                  indeterminateUpper,
                                                  indeterminateLower,
                                                  Type) {
  # Here, we make a local copy of the input data frame and name it `Data`. This lets us refer to `Data` throughout the function without altering the original `scores_data_df` object.
  Data <- scores_data_df
  
  # Initialize storage for metrics and indeterminate proportions
  # We initialize empty numeric vectors of length `testReps` for each performance metric.
  # Using `numeric(testReps)` ensures each vector starts filled with zeros and can store one value per CV repetition.
  Sensitivity   <- numeric(testReps)
  Specificity   <- numeric(testReps)
  PPV           <- numeric(testReps)
  NPV           <- numeric(testReps)
  Prevalence    <- numeric(testReps)
  Accuracy      <- numeric(testReps)
  nRemoved      <- numeric(testReps)
  
  # This creates an index vector `sampleIndicies` from 1 to the number of rows in `Data`, which we'll use to split observations into training and test sets across CV folds.
  sampleIndicies <- seq_len(nrow(Data))
  
  # Cross-validation loop
  # Begins the cross‑validation loop. We repeat the following code block `testReps` times, computing train/test splits and evaluating model performance each iteration.
  for (i in seq_len(testReps)) {
    # For all but the final iteration, randomly sample approximately `nrow(Data) / testReps - 1` rows for testing,
    # remove them from `sampleIndicies`, and use the remainder for training.
    # The final iteration simply uses whatever indices remain for testing.
    if (i < testReps) {
      ind <- sample(sampleIndicies, floor((nrow(Data) / testReps) - 1))
      sampleIndicies <- setdiff(sampleIndicies, ind)
    } else {
      ind <- sampleIndicies
    }
    
    # Split train and test sets
    # Creates the training set by excluding the indices in `ind`. `setdiff` ensures no overlap between train and test indices.
    train <- Data[setdiff(seq_len(nrow(Data)), ind), , drop = FALSE]
    test  <- Data[ind, , drop = FALSE]
    
    # If `Undersample` is `TRUE`, this block balances the training data by sampling the majority class (negative)
    # down to the size of the minority class (positive), preventing model bias toward the more frequent class.
    if (Undersample) {
      pos <- which(train$Target_Organ == 1)
      neg <- which(train$Target_Organ == 0)
      train <- train[c(pos, sample(neg, length(pos), replace = TRUE)), , drop = FALSE]
    }
    
    # This line trains a Random Forest classifier using your specified `mtry` (number of variables tried at each split),
    # number of trees (`ntree=500`), and returns variable importance and proximity measures for further analysis if needed.
    rf <- randomForest::randomForest(Target_Organ ~ ., data = train,
                                     mtry = best.m, importance = TRUE,
                                     ntree = 500, proximity = TRUE)
    
    # Generates predicted class probabilities for the positive class on the test set.
    # We extract the first column (positive class) to work with continuous scores before rounding.
    probs <- predict(rf, test, type = "prob")[, 1]
    
    # Identifies any predictions whose probabilities fall within the indeterminate range.
    # These are flagged for removal (set to `NA`) so they don’t contribute to performance metrics.
    idx <- which(probs > indeterminateLower & probs < indeterminateUpper)
    nRemoved[i] <- length(idx) / length(probs)
    probs[idx]  <- NA
    
    # Rounds the remaining probabilities to the nearest integer (0 or 1) to form discrete class predictions for the confusion matrix.
    preds <- round(probs)
    # Uses `caret` to compute the confusion matrix and extract performance metrics (sensitivity, specificity, etc.) comparing `preds` against the true `Target_Organ` labels.
    cm <- caret::confusionMatrix(
      factor(preds, levels = c(1, 0)),
      factor(test$Target_Organ, levels = c(1, 0))
    )
    
    # Stores each metric from the confusion matrix into the corresponding vector slot for the current repetition `i`.
    # After the loop, each vector holds one metric per fold.
    Sensitivity[i] <- cm$byClass["Sensitivity"]
    Specificity[i] <- cm$byClass["Specificity"]
    PPV[i]         <- cm$byClass["Pos Pred Value"]
    NPV[i]         <- cm$byClass["Neg Pred Value"]
    Prevalence[i]  <- cm$byClass["Prevalence"]
    Accuracy[i]    <- cm$byClass["Balanced Accuracy"]
  }
  
  # Aggregate performance across all repetitions
  # Aggregates the results across all folds by computing the column means of the combined metric vectors, excluding NAs from indeterminate predictions.
  Perf <- colMeans(
    cbind(Sensitivity, Specificity, PPV, NPV, Prevalence, Accuracy, nRemoved),
    na.rm = TRUE
  )
  # Prints the averaged performance metrics to the console for quick inspection.
  print(Perf)
  
  # Return results invisibly as a list
  # Wraps the performance summary and raw per-iteration results into a list and returns it invisibly.
  invisible(list(
    performance_summary = Perf,
    raw_results = list(
      Sensitivity   = Sensitivity,
      Specificity   = Specificity,
      PPV           = PPV,
      NPV           = NPV,
      Prevalence    = Prevalence,
      Accuracy      = Accuracy,
      indeterminate = nRemoved
    )
  ))
}

#' @title Random Forest with Cross-Validation
#'
#' @description
#' Builds a random forest model with cross-validation, computes metrics (sensitivity, specificity, accuracy),
#' and optionally balances classes via undersampling.
#'
#' @param scores_data_df Data frame. Must include a column named `Target_Organ` as the response.
#' @param Undersample Logical. If TRUE, undersamples majority class. Default FALSE.
#' @param best.m Numeric or NULL. Number of predictors sampled at each split. NULL uses default.
#' @param testReps Integer. Number of cross-validation repetitions (>=2).
#' @param Type Integer. Importance metric type: 1 = Mean Decrease Accuracy, 2 = Gini.
#'
#' @return A list with:
#' \describe{
#'   \item{performance_metrics}{Named vector of aggregated metrics: Sensitivity, Specificity, PPV, NPV, Prevalence, Accuracy.}
#'   \item{raw_results}{List of raw `sensitivity`, `specificity`, and `accuracy` vectors for each fold.}
#' }
#'
#' @export
revised_f9_get_rf_model_with_cv <- function(scores_data_df,
                                 Undersample = FALSE,
                                 best.m = NULL,
                                 testReps,
                                 Type) {
  # Input validation
  if (!is.data.frame(scores_data_df)) stop("scores_data_df must be a data.frame")
  if (testReps < 2) stop("testReps must be at least 2")
  
  n <- nrow(scores_data_df)
  # Preallocate metric vectors
  Sensitivity <- Specificity <- PPV <- NPV <- Prevalence <- Accuracy <- numeric(testReps)
  
  remaining <- NULL
  for (i in seq_len(testReps)) {
    if (i == 1L) {
      remaining <- seq_len(n)
    }
    if (i < testReps) {
      ind <- sample(seq_len(n), floor(n / testReps) - 1)
      remaining <- setdiff(remaining, ind)
    } else {
      ind <- remaining
    }
    
    train_idx <- setdiff(seq_len(n), ind)
    train <- scores_data_df[train_idx, , drop = FALSE]
    test  <- scores_data_df[ind, , drop = FALSE]
    
    if (isTRUE(Undersample)) {
      pos <- which(train[[1]] == 1)
      neg <- which(train[[1]] == 0)
      neg_samp <- sample(neg, length(pos), replace = TRUE)
      train <- train[c(pos, neg_samp), , drop = FALSE]
      test  <- rbind(train[-c(pos, neg_samp), , drop = FALSE], test)
    }
    
    rf <- randomForest::randomForest(
      Target_Organ ~ ., data = train,
      mytry = best.m,
      importance = TRUE,
      ntree = 500,
      proximity = TRUE
    )
    
    probs <- predict(rf, test, type = "prob")[, 1]
    preds <- factor(round(probs), levels = c(1, 0))
    truth <- factor(test$Target_Organ, levels = c(1, 0))
    cm <- caret::confusionMatrix(preds, truth)
    
    Sensitivity[i] <- cm$byClass["Sensitivity"]
    Specificity[i] <- cm$byClass["Specificity"]
    PPV[i]         <- cm$byClass[["Pos Pred Value"]]
    NPV[i]         <- cm$byClass[["Neg Pred Value"]]
    Prevalence[i]  <- cm$byClass[["Prevalence"]]
    Accuracy[i]    <- cm$byClass[["Balanced Accuracy"]]
  }
  
  PerformanceSummary <- colMeans(
    cbind(Sensitivity, Specificity, PPV, NPV, Prevalence, Accuracy),
    na.rm = TRUE
  )
  
  list(
    performance_metrics = PerformanceSummary,
    raw_results = list(
      sensitivity  = Sensitivity,
      specificity  = Specificity,
      accuracy     = Accuracy
    )
  )
}

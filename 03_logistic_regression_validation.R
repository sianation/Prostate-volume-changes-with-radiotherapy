################################################################################
# Script: 03_logistic_regression_validation.R
# Purpose: Logistic regression for binary outcomes and validation analysis
# Author: MOMENTUM Study Team
# Date: December 2024
################################################################################

# Load required packages --------------------------------------------------
library(tidyverse)
library(broom)        # For tidy model outputs
library(pROC)         # For ROC analysis
library(caret)        # For confusion matrix
library(boot)         # For bootstrapping

# Configuration -----------------------------------------------------------

# Significance level
ALPHA <- 0.05

# Bootstrap parameters
N_BOOTSTRAP <- 1000
BOOT_SEED <- 123

# Helper functions --------------------------------------------------------

#' Fit logistic regression model
#'
#' @param data Data frame with outcome and predictors
#' @param formula Model formula
#' @param family Family for glm (default: binomial)
#' @return Fitted glm object
fit_logistic <- function(data, formula, family = binomial(link = "logit")) {
  
  if (is.character(formula)) {
    formula <- as.formula(formula)
  }
  
  message("Fitting model: ", deparse(formula))
  
  model <- glm(formula, data = data, family = family)
  
  return(model)
}

#' Extract model coefficients with odds ratios
#'
#' @param model Fitted glm model
#' @param conf.level Confidence level (default: 0.95)
#' @return Data frame with coefficients, OR, CI, and p-values
get_odds_ratios <- function(model, conf.level = 1 - ALPHA) {
  
  # Tidy coefficients
  coefs <- tidy(model, conf.int = TRUE, conf.level = conf.level, exponentiate = FALSE)
  
  # Calculate odds ratios and CIs
  coefs <- coefs %>%
    mutate(
      odds_ratio = exp(estimate),
      or_lower = exp(conf.low),
      or_upper = exp(conf.high),
      sig = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        p.value < 0.1 ~ ".",
        TRUE ~ ""
      )
    )
  
  return(coefs)
}

#' Calculate model performance metrics
#'
#' @param model Fitted logistic regression model
#' @param data Data frame with actual outcomes
#' @param threshold Probability threshold for classification (default: 0.5)
#' @return List with AUC, accuracy, sensitivity, specificity, etc.
calculate_performance <- function(model, data = NULL, threshold = 0.5) {
  
  if (is.null(data)) {
    data <- model$data
  }
  
  # Get predictions
  predicted_prob <- predict(model, newdata = data, type = "response")
  predicted_class <- ifelse(predicted_prob >= threshold, 1, 0)
  
  # Get actual outcomes
  outcome_var <- all.vars(formula(model))[1]
  actual <- data[[outcome_var]]
  
  # ROC analysis
  roc_obj <- roc(actual, predicted_prob, quiet = TRUE)
  auc_val <- auc(roc_obj)
  
  # Confusion matrix
  cm <- confusionMatrix(
    factor(predicted_class, levels = c(0, 1)),
    factor(actual, levels = c(0, 1)),
    positive = "1"
  )
  
  # Extract metrics
  metrics <- list(
    auc = as.numeric(auc_val),
    accuracy = as.numeric(cm$overall["Accuracy"]),
    sensitivity = as.numeric(cm$byClass["Sensitivity"]),
    specificity = as.numeric(cm$byClass["Specificity"]),
    ppv = as.numeric(cm$byClass["Pos Pred Value"]),
    npv = as.numeric(cm$byClass["Neg Pred Value"]),
    f1 = as.numeric(cm$byClass["F1"]),
    confusion_matrix = cm$table,
    roc = roc_obj
  )
  
  return(metrics)
}

#' Bootstrap confidence intervals for performance metrics
#'
#' @param model Fitted model
#' @param data Data for bootstrapping
#' @param n_boot Number of bootstrap samples
#' @param seed Random seed
#' @return Data frame with bootstrapped CIs
bootstrap_performance <- function(model, data, n_boot = N_BOOTSTRAP, seed = BOOT_SEED) {
  
  set.seed(seed)
  
  message("Bootstrapping with ", n_boot, " samples...")
  
  # Bootstrap function
  boot_fn <- function(data, indices) {
    boot_data <- data[indices, ]
    boot_model <- update(model, data = boot_data)
    metrics <- calculate_performance(boot_model, boot_data)
    return(c(metrics$auc, metrics$accuracy, metrics$sensitivity, metrics$specificity))
  }
  
  # Run bootstrap
  boot_results <- boot(data = data, statistic = boot_fn, R = n_boot)
  
  # Extract CIs
  metrics_names <- c("AUC", "Accuracy", "Sensitivity", "Specificity")
  ci_results <- map_df(1:4, function(i) {
    ci <- boot.ci(boot_results, type = "bca", index = i)
    data.frame(
      metric = metrics_names[i],
      estimate = boot_results$t0[i],
      lower_ci = ci$bca[4],
      upper_ci = ci$bca[5]
    )
  })
  
  return(ci_results)
}

#' Plot ROC curve
#'
#' @param roc_obj ROC object from pROC
#' @param title Plot title
#' @return ggplot object
plot_roc_curve <- function(roc_obj, title = "ROC Curve") {
  
  roc_df <- data.frame(
    sensitivity = roc_obj$sensitivities,
    specificity = roc_obj$specificities,
    thresholds = roc_obj$thresholds
  )
  
  p <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(size = 1.2, color = "#2E86AB") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    annotate(
      "text", 
      x = 0.7, y = 0.3, 
      label = sprintf("AUC = %.3f", auc(roc_obj)),
      size = 5
    ) +
    labs(
      title = title,
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12)
    ) +
    coord_equal()
  
  return(p)
}

#' Calibration plot
#'
#' @param model Fitted logistic model
#' @param data Data frame
#' @param n_bins Number of bins for calibration (default: 10)
#' @return ggplot object
plot_calibration <- function(model, data, n_bins = 10) {
  
  # Get predictions and outcomes
  predicted_prob <- predict(model, newdata = data, type = "response")
  outcome_var <- all.vars(formula(model))[1]
  actual <- data[[outcome_var]]
  
  # Create bins
  cal_df <- data.frame(
    predicted = predicted_prob,
    actual = actual
  ) %>%
    mutate(
      bin = cut(predicted, breaks = n_bins, include.lowest = TRUE)
    ) %>%
    group_by(bin) %>%
    summarise(
      predicted_mean = mean(predicted),
      observed_mean = mean(actual),
      n = n(),
      .groups = "drop"
    )
  
  # Plot
  p <- ggplot(cal_df, aes(x = predicted_mean, y = observed_mean)) +
    geom_point(aes(size = n), alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = TRUE, color = "#2E86AB") +
    scale_size_continuous(range = c(2, 10)) +
    labs(
      title = "Calibration Plot",
      x = "Predicted Probability",
      y = "Observed Proportion",
      size = "N"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1))
  
  return(p)
}

# Validation workflow -----------------------------------------------------

#' Validate model on external cohort
#'
#' @param model Fitted model from training data
#' @param validation_data External validation dataset
#' @param outcome_var Name of outcome variable
#' @param plot_dir Directory to save plots (default: NULL)
#' @return List with validation metrics and plots
validate_model <- function(model, validation_data, outcome_var, plot_dir = NULL) {
  
  message("\nValidating model on external cohort...")
  message(paste("Validation N:", nrow(validation_data)))
  message(paste("Outcome prevalence:", mean(validation_data[[outcome_var]])))
  
  # Calculate performance
  metrics <- calculate_performance(model, validation_data)
  
  message("\nValidation Performance:")
  message(sprintf("  AUC: %.3f", metrics$auc))
  message(sprintf("  Accuracy: %.3f", metrics$accuracy))
  message(sprintf("  Sensitivity: %.3f", metrics$sensitivity))
  message(sprintf("  Specificity: %.3f", metrics$specificity))
  
  # Bootstrap CIs
  boot_metrics <- bootstrap_performance(model, validation_data)
  message("\nBootstrapped 95% CIs:")
  print(boot_metrics)
  
  # Generate plots
  roc_plot <- plot_roc_curve(metrics$roc, "Validation ROC Curve")
  cal_plot <- plot_calibration(model, validation_data)
  
  # Save plots
  if (!is.null(plot_dir)) {
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(plot_dir, "validation_roc.png"), roc_plot, width = 8, height = 8)
    ggsave(file.path(plot_dir, "validation_calibration.png"), cal_plot, width = 8, height = 8)
  }
  
  return(list(
    metrics = metrics,
    bootstrap_ci = boot_metrics,
    plots = list(roc = roc_plot, calibration = cal_plot)
  ))
}

#' Compare model performance between training and validation
#'
#' @param model Fitted model
#' @param training_data Training dataset
#' @param validation_data Validation dataset
#' @return Data frame with comparison
compare_performance <- function(model, training_data, validation_data) {
  
  train_metrics <- calculate_performance(model, training_data)
  val_metrics <- calculate_performance(model, validation_data)
  
  comparison <- data.frame(
    cohort = c("Training", "Validation"),
    n = c(nrow(training_data), nrow(validation_data)),
    auc = c(train_metrics$auc, val_metrics$auc),
    accuracy = c(train_metrics$accuracy, val_metrics$accuracy),
    sensitivity = c(train_metrics$sensitivity, val_metrics$sensitivity),
    specificity = c(train_metrics$specificity, val_metrics$specificity)
  )
  
  message("\nPerformance Comparison:")
  print(comparison)
  
  return(comparison)
}

# Example usage (uncomment to run) ----------------------------------------
# # Load data
# training_data <- read_csv("data/training_cohort.csv")
# validation_data <- read_csv("data/validation_cohort.csv")
# 
# # Define binary outcome (e.g., large volume decrease)
# training_data <- training_data %>%
#   mutate(large_decrease = as.numeric(volume_change_percent < -30))
# 
# validation_data <- validation_data %>%
#   mutate(large_decrease = as.numeric(volume_change_percent < -30))
# 
# # Fit model on training data
# model <- fit_logistic(
#   data = training_data,
#   formula = "large_decrease ~ baseline_volume + age + treatment_type"
# )
# 
# # Extract odds ratios
# or_table <- get_odds_ratios(model)
# print(or_table)
# 
# # Validate on external cohort
# validation_results <- validate_model(
#   model = model,
#   validation_data = validation_data,
#   outcome_var = "large_decrease",
#   plot_dir = "results/validation_plots"
# )
# 
# # Compare performance
# comparison <- compare_performance(model, training_data, validation_data)

################################################################################
# End of script
################################################################################

################################################################################
# Script: 02_linear_mixed_effects_modeling.R
# Purpose: Linear mixed effects models for prostate volume trajectories
# Author: MOMENTUM Study Team
# Date: December 2024
################################################################################

# Load required packages --------------------------------------------------
library(tidyverse)
library(lme4)         # For linear mixed models
library(lmerTest)     # For p-values in lmer
library(broom.mixed)  # For tidy model outputs
library(emmeans)      # For estimated marginal means
library(performance)  # For model diagnostics

# Configuration -----------------------------------------------------------

# Model convergence settings
OPTIMIZER <- "bobyqa"
MAX_ITER <- 1e5

# Significance level
ALPHA <- 0.05

# Helper functions --------------------------------------------------------

#' Fit linear mixed effects model for volume change
#'
#' @param data Data frame with patient_id, weeks_from_rt_start, volume_change_percent, and covariates
#' @param formula Model formula (as string or formula object)
#' @param optimizer Optimizer for lmer (default: "bobyqa")
#' @param max_iter Maximum iterations (default: 1e5)
#' @return Fitted lmer model object
fit_lmm <- function(data, formula, optimizer = OPTIMIZER, max_iter = MAX_ITER) {
  
  # Convert formula if needed
  if (is.character(formula)) {
    formula <- as.formula(formula)
  }
  
  message("Fitting model: ", deparse(formula))
  
  # Fit model with specified optimizer
  model <- lmer(
    formula, 
    data = data,
    control = lmerControl(
      optimizer = optimizer,
      optCtrl = list(maxfun = max_iter)
    )
  )
  
  return(model)
}

#' Extract model coefficients in tidy format
#'
#' @param model Fitted lmer model
#' @param conf.level Confidence level (default: 0.95)
#' @return Data frame with coefficients, SE, CI, and p-values
get_model_coefficients <- function(model, conf.level = 1 - ALPHA) {
  
  coefs <- tidy(model, conf.int = TRUE, conf.level = conf.level, effects = "fixed")
  
  # Add significance stars
  coefs <- coefs %>%
    mutate(
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

#' Calculate model diagnostics and fit statistics
#'
#' @param model Fitted lmer model
#' @return List with AIC, BIC, R-squared, ICC, and convergence info
get_model_diagnostics <- function(model) {
  
  # Basic fit statistics
  diag <- list(
    aic = AIC(model),
    bic = BIC(model),
    loglik = logLik(model)[1],
    nobs = nobs(model),
    converged = model@optinfo$conv$opt == 0
  )
  
  # R-squared
  r2 <- performance::r2(model)
  diag$r2_conditional <- r2$R2_conditional
  diag$r2_marginal <- r2$R2_marginal
  
  # Intraclass correlation
  icc <- performance::icc(model)
  diag$icc <- icc$ICC_adjusted
  
  return(diag)
}

#' Generate predicted trajectories for each patient
#'
#' @param model Fitted lmer model
#' @param newdata Data frame with prediction timepoints
#' @return Data frame with predictions and confidence intervals
predict_trajectories <- function(model, newdata = NULL) {
  
  if (is.null(newdata)) {
    # Use original data if no new data provided
    newdata <- model@frame
  }
  
  # Predictions with confidence intervals
  preds <- predict(model, newdata = newdata, re.form = ~(1 | patient_id))
  
  # Get prediction intervals using bootMer (computationally expensive)
  # For faster computation, use predictInterval from merTools package
  
  pred_df <- newdata %>%
    mutate(predicted_value = preds)
  
  return(pred_df)
}

#' Test model assumptions with diagnostic plots
#'
#' @param model Fitted lmer model
#' @param plot_dir Directory to save plots (default: NULL = no save)
#' @return List of ggplot objects
check_model_assumptions <- function(model, plot_dir = NULL) {
  
  # Extract residuals and fitted values
  model_df <- data.frame(
    fitted = fitted(model),
    residuals = residuals(model),
    std_residuals = residuals(model) / sd(residuals(model))
  )
  
  # 1. Residuals vs fitted
  p1 <- ggplot(model_df, aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = TRUE) +
    labs(
      title = "Residuals vs Fitted Values",
      x = "Fitted values",
      y = "Residuals"
    ) +
    theme_bw()
  
  # 2. Q-Q plot
  p2 <- ggplot(model_df, aes(sample = std_residuals)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    labs(
      title = "Normal Q-Q Plot",
      x = "Theoretical quantiles",
      y = "Standardized residuals"
    ) +
    theme_bw()
  
  # 3. Scale-location plot
  p3 <- ggplot(model_df, aes(x = fitted, y = sqrt(abs(std_residuals)))) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "loess", se = TRUE) +
    labs(
      title = "Scale-Location Plot",
      x = "Fitted values",
      y = expression(sqrt("|Standardized residuals|"))
    ) +
    theme_bw()
  
  # Save plots if directory specified
  if (!is.null(plot_dir)) {
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(plot_dir, "residuals_vs_fitted.png"), p1, width = 8, height = 6)
    ggsave(file.path(plot_dir, "qq_plot.png"), p2, width = 8, height = 6)
    ggsave(file.path(plot_dir, "scale_location.png"), p3, width = 8, height = 6)
  }
  
  return(list(residuals_fitted = p1, qq_plot = p2, scale_location = p3))
}

# Main modeling workflow --------------------------------------------------

#' Build and compare linear mixed effects models
#'
#' @param data Cleaned data frame from 01_data_cleaning script
#' @param outcome Outcome variable name (default: "volume_change_percent")
#' @param time_var Time variable name (default: "weeks_from_rt_start")
#' @param covariates Vector of covariate names
#' @param output_file Path to save model results (default: NULL)
#' @return List with models and diagnostics
build_lmm_models <- function(data, 
                             outcome = "volume_change_percent",
                             time_var = "weeks_from_rt_start",
                             covariates = NULL,
                             output_file = NULL) {
  
  message("\nBuilding linear mixed effects models...")
  message(paste("Outcome:", outcome))
  message(paste("Time variable:", time_var))
  message(paste("N patients:", n_distinct(data$patient_id)))
  message(paste("N observations:", nrow(data)))
  
  # Filter to complete cases
  vars_needed <- c("patient_id", outcome, time_var, covariates)
  data_complete <- data %>%
    select(all_of(vars_needed)) %>%
    drop_na()
  
  message(paste("N complete observations:", nrow(data_complete)))
  
  # Model 1: Null model (random intercept only)
  formula_null <- paste0(outcome, " ~ 1 + (1 | patient_id)")
  model_null <- fit_lmm(data_complete, formula_null)
  
  # Model 2: Time only (random intercept)
  formula_time <- paste0(outcome, " ~ ", time_var, " + (1 | patient_id)")
  model_time <- fit_lmm(data_complete, formula_time)
  
  # Model 3: Time with random slope
  formula_time_slope <- paste0(outcome, " ~ ", time_var, " + (1 + ", time_var, " | patient_id)")
  model_time_slope <- tryCatch(
    fit_lmm(data_complete, formula_time_slope),
    error = function(e) {
      message("Random slope model failed to converge, using random intercept only")
      return(model_time)
    }
  )
  
  # Model 4: Time + covariates (if provided)
  if (!is.null(covariates) && length(covariates) > 0) {
    cov_string <- paste(covariates, collapse = " + ")
    formula_cov <- paste0(outcome, " ~ ", time_var, " + ", cov_string, " + (1 | patient_id)")
    model_cov <- fit_lmm(data_complete, formula_cov)
  } else {
    model_cov <- NULL
  }
  
  # Compare models
  message("\nModel comparison:")
  model_list <- list(
    null = model_null,
    time = model_time,
    time_slope = model_time_slope
  )
  
  if (!is.null(model_cov)) {
    model_list$covariates <- model_cov
  }
  
  # Extract diagnostics for all models
  comparison <- map_df(names(model_list), function(name) {
    diag <- get_model_diagnostics(model_list[[name]])
    data.frame(
      model = name,
      aic = diag$aic,
      bic = diag$bic,
      loglik = diag$loglik,
      r2_marginal = diag$r2_marginal,
      r2_conditional = diag$r2_conditional,
      icc = diag$icc,
      converged = diag$converged
    )
  })
  
  print(comparison)
  
  # Select best model by BIC
  best_model_name <- comparison$model[which.min(comparison$bic)]
  best_model <- model_list[[best_model_name]]
  
  message("\nBest model (by BIC): ", best_model_name)
  
  # Extract coefficients from best model
  coefs <- get_model_coefficients(best_model)
  message("\nFixed effects:")
  print(coefs)
  
  # Save results if output file specified
  if (!is.null(output_file)) {
    results <- list(
      model_comparison = comparison,
      best_model_name = best_model_name,
      coefficients = coefs
    )
    saveRDS(results, output_file)
    message("\nResults saved to: ", output_file)
  }
  
  return(list(
    models = model_list,
    best_model = best_model,
    comparison = comparison,
    coefficients = coefs
  ))
}

# Example usage (uncomment to run) ----------------------------------------
# # Load cleaned data from script 01
# data <- read_csv("data/cleaned_volume_data.csv")
# 
# # Filter to quality-passed observations
# data_qc <- data %>%
#   filter(qc_flag == "pass")
# 
# # Build models
# model_results <- build_lmm_models(
#   data = data_qc,
#   outcome = "volume_change_percent",
#   time_var = "weeks_from_rt_start",
#   covariates = c("baseline_volume", "age", "treatment_type"),
#   output_file = "results/lmm_results.rds"
# )
# 
# # Check assumptions
# diagnostic_plots <- check_model_assumptions(
#   model_results$best_model,
#   plot_dir = "results/diagnostic_plots"
# )

################################################################################
# End of script
################################################################################

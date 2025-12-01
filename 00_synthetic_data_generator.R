################################################################################
# Script: 00_synthetic_data_generator.R
# Purpose: Generate synthetic volume data for testing pipeline
# Author: MOMENTUM Study Team
# Date: December 2024
################################################################################

# Load required packages --------------------------------------------------
library(tidyverse)
library(MASS)  # For mvrnorm

# Configuration -----------------------------------------------------------

# Sample size parameters
N_PATIENTS <- 200
N_TIMEPOINTS <- 10

# Random seed for reproducibility
SEED <- 42

# Synthetic data generation -----------------------------------------------

#' Generate synthetic prostate volume data
#'
#' @param n_patients Number of patients to simulate
#' @param n_timepoints Number of scans per patient
#' @param seed Random seed
#' @return Data frame with synthetic volume measurements
generate_synthetic_data <- function(n_patients = N_PATIENTS, 
                                   n_timepoints = N_TIMEPOINTS,
                                   seed = SEED) {
  
  set.seed(seed)
  
  message("Generating synthetic data...")
  message(paste("  Patients:", n_patients))
  message(paste("  Timepoints per patient:", n_timepoints))
  
  # Generate patient-level baseline characteristics
  patient_data <- tibble(
    patient_id = sprintf("P%03d", 1:n_patients),
    age = rnorm(n_patients, mean = 68, sd = 8),
    baseline_psa = rlnorm(n_patients, meanlog = 1.5, sdlog = 0.8),
    treatment_type = sample(c("IMRT", "SBRT", "Brachytherapy"), n_patients, replace = TRUE),
    # Baseline prostate volume (cm^3)
    baseline_volume = rnorm(n_patients, mean = 40, sd = 15)
  )
  
  # Ensure positive values
  patient_data <- patient_data %>%
    mutate(
      age = pmax(age, 40),
      baseline_psa = pmax(baseline_psa, 0.1),
      baseline_volume = pmax(baseline_volume, 15)
    )
  
  # Define timepoints (weeks from RT start)
  timepoints <- c(0, 1, 2, 3, 4, 5, 6, 8, 12, 24)
  
  # Generate trajectory patterns (3 clusters)
  # Cluster 1: Stable (minimal change)
  # Cluster 2: Gradual decrease
  # Cluster 3: Rapid decrease then plateau
  
  cluster_assignment <- sample(1:3, n_patients, replace = TRUE, prob = c(0.3, 0.4, 0.3))
  
  # Generate longitudinal data
  longitudinal_data <- map_df(1:n_patients, function(i) {
    patient <- patient_data[i, ]
    cluster <- cluster_assignment[i]
    
    # Base trajectory by cluster
    if (cluster == 1) {
      # Stable: small random fluctuations
      volume_changes <- rnorm(n_timepoints, mean = -2, sd = 5)
    } else if (cluster == 2) {
      # Gradual decrease: linear trend
      volume_changes <- -1.5 * timepoints[1:n_timepoints] + rnorm(n_timepoints, 0, 4)
    } else {
      # Rapid decrease then plateau
      volume_changes <- -3 * sqrt(timepoints[1:n_timepoints]) + rnorm(n_timepoints, 0, 3)
    }
    
    # Add treatment effect
    treatment_effect <- case_when(
      patient$treatment_type == "SBRT" ~ -5,
      patient$treatment_type == "Brachytherapy" ~ -8,
      TRUE ~ 0
    )
    
    volume_changes <- volume_changes + treatment_effect
    
    # Add baseline volume effect (larger prostates shrink more)
    volume_effect <- (patient$baseline_volume - 40) * 0.1
    volume_changes <- volume_changes + volume_effect
    
    # Calculate actual volumes
    volumes <- patient$baseline_volume * (1 + volume_changes / 100)
    volumes <- pmax(volumes, 5)  # Minimum 5 cm^3
    
    # Generate scan dates
    rt_start_date <- as.Date("2023-01-01") + sample(0:365, 1)
    scan_dates <- rt_start_date + 7 * timepoints[1:n_timepoints]
    
    # Create data frame
    tibble(
      patient_id = patient$patient_id,
      scan_date = scan_dates,
      scan_timepoint = ifelse(timepoints[1:n_timepoints] == 0, "Fraction_1",
                             paste0("Week_", timepoints[1:n_timepoints])),
      weeks_from_rt_start = timepoints[1:n_timepoints],
      roi_name = "Prostate",
      volume = volumes,
      voxel_count = round(volumes / 0.001),
      age = patient$age,
      baseline_psa = patient$baseline_psa,
      treatment_type = patient$treatment_type,
      baseline_volume = patient$baseline_volume,
      true_cluster = cluster
    )
  })
  
  message(paste("\nGenerated", nrow(longitudinal_data), "observations"))
  message("Cluster distribution:")
  print(table(longitudinal_data %>% 
               group_by(patient_id) %>% 
               slice(1) %>% 
               pull(true_cluster)))
  
  return(longitudinal_data)
}

#' Add quality control issues to synthetic data
#'
#' @param data Synthetic data frame
#' @param outlier_rate Proportion of outliers to add
#' @return Data with some QC issues
add_qc_issues <- function(data, outlier_rate = 0.05) {
  
  set.seed(SEED + 1)
  
  n_outliers <- round(nrow(data) * outlier_rate)
  outlier_indices <- sample(1:nrow(data), n_outliers)
  
  data_with_issues <- data
  
  # Add volume outliers
  data_with_issues$volume[outlier_indices] <- data_with_issues$volume[outlier_indices] * 
    runif(n_outliers, 1.5, 3)
  
  # Add some missing values
  n_missing <- round(nrow(data) * 0.02)
  missing_indices <- sample(1:nrow(data), n_missing)
  data_with_issues$volume[missing_indices] <- NA
  
  message(paste("\nAdded", n_outliers, "outliers and", n_missing, "missing values"))
  
  return(data_with_issues)
}

# Main execution ----------------------------------------------------------

if (interactive()) {
  
  message("\n=== Generating Synthetic MOMENTUM Data ===")
  
  # Generate clean synthetic data
  synthetic_data <- generate_synthetic_data()
  
  # Add some QC issues for realism
  synthetic_data_with_issues <- add_qc_issues(synthetic_data)
  
  # Create output directory
  dir.create("data", showWarnings = FALSE)
  
  # Save both versions
  write_csv(synthetic_data, "data/synthetic_volume_data_clean.csv")
  write_csv(synthetic_data_with_issues, "data/synthetic_volume_data.csv")
  
  message("\nSynthetic data saved to:")
  message("  - data/synthetic_volume_data_clean.csv (no QC issues)")
  message("  - data/synthetic_volume_data.csv (with QC issues for testing)")
  
  # Preview
  message("\nData preview:")
  print(head(synthetic_data_with_issues, 10))
  
  message("\nSummary statistics:")
  print(summary(synthetic_data_with_issues %>% select(age, baseline_volume, volume)))
  
  message("\n=== Ready to test pipeline! ===")
  message("Next steps:")
  message("  1. Run 01_data_cleaning_roi_harmonization.R")
  message("  2. Run 02_linear_mixed_effects_modeling.R")
  message("  3. Run 03_logistic_regression_validation.R")
  message("  4. Run 04_trajectory_cluster_analysis.R")
}

################################################################################
# End of script
################################################################################

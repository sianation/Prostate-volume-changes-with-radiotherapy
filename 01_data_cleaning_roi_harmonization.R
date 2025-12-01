################################################################################
# Script: 01_data_cleaning_roi_harmonization.R
# Purpose: Data cleaning, ROI harmonization, and volume calculation
# Author: MOMENTUM Study Team
# Date: December 2024
################################################################################

# Load required packages --------------------------------------------------
library(tidyverse)
library(lubridate)

# Configuration -----------------------------------------------------------
# Adjust these parameters as needed for your dataset

# Volume calculation parameters
VOXEL_VOLUME <- 0.001  # cm^3 per voxel (adjust based on imaging protocol)

# Quality control thresholds
MIN_VOLUME <- 5        # minimum plausible prostate volume (cm^3)
MAX_VOLUME <- 200      # maximum plausible prostate volume (cm^3)
MAX_VOLUME_CHANGE_RATE <- 0.5  # maximum fractional change per day

# ROI name harmonization mapping
# Add your institution's ROI naming variants here
ROI_MAPPING <- tribble(
  ~original_name, ~standard_name,
  "Prostate", "Prostate",
  "prostate", "Prostate",
  "PROSTATE", "Prostate",
  "Pros", "Prostate",
  "CTV_Prostate", "Prostate",
  "CTV_PROSTATE", "Prostate",
  # Add more variants as needed
)

# Helper functions --------------------------------------------------------

#' Calculate volume from voxel count
#'
#' @param voxel_count Integer vector of voxel counts
#' @param voxel_volume Volume per voxel in cm^3 (default: 0.001)
#' @return Numeric vector of volumes in cm^3
calculate_volume <- function(voxel_count, voxel_volume = VOXEL_VOLUME) {
  return(voxel_count * voxel_volume)
}

#' Harmonize ROI names to standard nomenclature
#'
#' @param roi_names Character vector of ROI names
#' @param mapping Data frame with original_name and standard_name columns
#' @return Character vector of harmonized ROI names
harmonize_roi_names <- function(roi_names, mapping = ROI_MAPPING) {
  # Create a named vector for fast lookup
  lookup <- setNames(mapping$standard_name, mapping$original_name)
  
  # Map names, keep original if not in mapping
  harmonized <- roi_names
  matched <- roi_names %in% names(lookup)
  harmonized[matched] <- lookup[roi_names[matched]]
  
  # Warn about unmapped names
  unmapped <- unique(roi_names[!matched])
  if (length(unmapped) > 0) {
    warning("Unmapped ROI names found: ", paste(unmapped, collapse = ", "))
  }
  
  return(harmonized)
}

#' Flag outliers based on volume range
#'
#' @param volume Numeric vector of volumes
#' @param min_vol Minimum plausible volume
#' @param max_vol Maximum plausible volume
#' @return Logical vector indicating outliers (TRUE = outlier)
flag_volume_outliers <- function(volume, 
                                 min_vol = MIN_VOLUME, 
                                 max_vol = MAX_VOLUME) {
  return(volume < min_vol | volume > max_vol)
}

#' Flag outliers based on rate of volume change
#'
#' @param data Data frame with patient_id, scan_date, and volume columns
#' @param max_rate Maximum fractional change per day
#' @return Logical vector indicating outliers (TRUE = outlier)
flag_change_rate_outliers <- function(data, max_rate = MAX_VOLUME_CHANGE_RATE) {
  data <- data %>%
    arrange(patient_id, scan_date) %>%
    group_by(patient_id) %>%
    mutate(
      days_between = as.numeric(difftime(scan_date, lag(scan_date), units = "days")),
      volume_change = volume - lag(volume),
      fractional_change = abs(volume_change) / lag(volume),
      change_rate = fractional_change / days_between,
      outlier = !is.na(change_rate) & change_rate > max_rate
    ) %>%
    ungroup()
  
  return(data$outlier)
}

#' Clean and validate volume data
#'
#' @param data Data frame with raw volume measurements
#' @return Data frame with cleaned data and QC flags
clean_volume_data <- function(data) {
  
  # Step 1: Harmonize ROI names
  message("Harmonizing ROI names...")
  data <- data %>%
    mutate(roi_name_harmonized = harmonize_roi_names(roi_name))
  
  # Step 2: Calculate volumes if voxel counts provided
  if ("voxel_count" %in% names(data) && !"volume" %in% names(data)) {
    message("Calculating volumes from voxel counts...")
    data <- data %>%
      mutate(volume = calculate_volume(voxel_count))
  }
  
  # Step 3: Parse and validate dates
  message("Validating dates...")
  if (is.character(data$scan_date)) {
    data <- data %>%
      mutate(scan_date = ymd(scan_date))
  }
  
  # Step 4: Flag volume outliers
  message("Flagging volume outliers...")
  data <- data %>%
    mutate(volume_outlier = flag_volume_outliers(volume))
  
  # Step 5: Calculate time from radiotherapy start
  message("Calculating time from RT start...")
  data <- data %>%
    group_by(patient_id) %>%
    arrange(scan_date) %>%
    mutate(
      rt_start_date = scan_date[which(scan_timepoint == "Fraction_1")[1]],
      days_from_rt_start = as.numeric(difftime(scan_date, rt_start_date, units = "days")),
      weeks_from_rt_start = days_from_rt_start / 7
    ) %>%
    ungroup()
  
  # Step 6: Flag rate-of-change outliers
  message("Flagging rate-of-change outliers...")
  change_outliers <- flag_change_rate_outliers(data)
  data <- data %>%
    mutate(change_rate_outlier = change_outliers)
  
  # Step 7: Create composite QC flag
  data <- data %>%
    mutate(
      qc_flag = case_when(
        volume_outlier ~ "volume_outlier",
        change_rate_outlier ~ "change_rate_outlier",
        is.na(volume) ~ "missing_volume",
        is.na(scan_date) ~ "missing_date",
        TRUE ~ "pass"
      )
    )
  
  # Report QC summary
  qc_summary <- data %>%
    count(qc_flag) %>%
    mutate(percent = 100 * n / sum(n))
  
  message("\nQuality control summary:")
  print(qc_summary)
  
  return(data)
}

#' Calculate baseline and follow-up volumes
#'
#' @param data Cleaned data frame
#' @return Data frame with baseline volumes and percentage changes
calculate_volume_changes <- function(data) {
  
  data <- data %>%
    group_by(patient_id) %>%
    arrange(scan_date) %>%
    mutate(
      baseline_volume = first(volume[scan_timepoint == "Fraction_1"]),
      volume_change_absolute = volume - baseline_volume,
      volume_change_percent = 100 * (volume - baseline_volume) / baseline_volume,
      volume_ratio = volume / baseline_volume
    ) %>%
    ungroup()
  
  return(data)
}

# Main workflow -----------------------------------------------------------

#' Process raw volume data through complete cleaning pipeline
#'
#' @param input_file Path to input CSV file with raw data
#' @param output_file Path to output CSV file for cleaned data
#' @return Data frame with cleaned and processed data
process_volume_data <- function(input_file, output_file = NULL) {
  
  message("Starting data processing pipeline...")
  message(paste("Input file:", input_file))
  
  # Read data
  message("\nReading input data...")
  raw_data <- read_csv(input_file, show_col_types = FALSE)
  message(paste("Rows read:", nrow(raw_data)))
  message(paste("Columns:", paste(names(raw_data), collapse = ", ")))
  
  # Required columns check
  required_cols <- c("patient_id", "scan_date", "roi_name", "scan_timepoint")
  missing_cols <- setdiff(required_cols, names(raw_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Clean data
  cleaned_data <- clean_volume_data(raw_data)
  
  # Calculate volume changes
  message("\nCalculating volume changes...")
  final_data <- calculate_volume_changes(cleaned_data)
  
  # Write output
  if (!is.null(output_file)) {
    message("\nWriting cleaned data to:", output_file)
    write_csv(final_data, output_file)
  }
  
  message("\nProcessing complete!")
  return(final_data)
}

# Example usage (uncomment to run) ----------------------------------------
# cleaned_data <- process_volume_data(
#   input_file = "data/raw_volume_data.csv",
#   output_file = "data/cleaned_volume_data.csv"
# )

################################################################################
# End of script
################################################################################

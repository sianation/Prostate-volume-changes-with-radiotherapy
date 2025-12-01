# Prostate Volume Changes with Radiotherapy

## MOMENTUM Study - Analysis Pipeline

**Repository for statistical analysis of prostate volume trajectories during radiotherapy**

---

## Overview

This repository contains production-ready R scripts for analyzing longitudinal prostate volume measurements during radiotherapy. The pipeline includes data cleaning, mixed effects modeling, validation analysis, and trajectory clustering.

### Key Features

- **Data Cleaning & QC**: ROI harmonization, outlier detection, quality control
- **Statistical Modeling**: Linear mixed effects models with multiple specifications
- **External Validation**: Logistic regression with ROC analysis and bootstrapping
- **Pattern Recognition**: Trajectory clustering to identify distinct volume change patterns
- **Synthetic Data**: Test dataset generator for pipeline validation

---

## Repository Structure

```
.
├── 00_synthetic_data_generator.R      # Generate test data
├── 01_data_cleaning_roi_harmonization.R    # Data QC and preparation  
├── 02_linear_mixed_effects_modeling.R      # LMM analysis
├── 03_logistic_regression_validation.R     # External validation
├── 04_trajectory_cluster_analysis.R        # Clustering analysis
└── README.md                               # This file
```

---

## Quick Start

### 1. Install Required Packages

```r
install.packages(c(
  "tidyverse", "lubridate",           # Data manipulation
  "lme4", "lmerTest", "broom.mixed",  # Mixed models
  "emmeans", "performance",           # Model diagnostics
  "pROC", "caret", "boot",            # Validation
  "kml", "lcmm", "factoextra",        # Clustering
  "ggrepel", "MASS"                   # Visualization
))
```

### 2. Generate Synthetic Test Data

```r
source("00_synthetic_data_generator.R")
# Creates: data/synthetic_volume_data.csv
```

### 3. Run Analysis Pipeline

```r
# Step 1: Data cleaning
source("01_data_cleaning_roi_harmonization.R")
cleaned_data <- process_volume_data(
  input_file = "data/synthetic_volume_data.csv",
  output_file = "data/cleaned_volume_data.csv"
)

# Step 2: Mixed effects modeling  
source("02_linear_mixed_effects_modeling.R")
model_results <- build_lmm_models(
  data = cleaned_data %>% filter(qc_flag == "pass"),
  outcome = "volume_change_percent",
  time_var = "weeks_from_rt_start"
)

# Step 3: Validation (if external cohort available)
source("03_logistic_regression_validation.R")
# See script for validation workflow

# Step 4: Trajectory clustering
source("04_trajectory_cluster_analysis.R")
trajectory_results <- perform_trajectory_analysis(
  data = cleaned_data,
  k = NULL  # Auto-select optimal clusters
)
```

---

## Script Documentation

### `00_synthetic_data_generator.R`

**Purpose**: Generate realistic synthetic prostate volume data for testing

**Features**:
- 200 patients, 10 timepoints each
- Three trajectory patterns (stable, gradual decrease, rapid decrease)
- Baseline characteristics: age, PSA, treatment type
- Optional QC issues for testing cleaning pipeline

**Usage**:
```r
source("00_synthetic_data_generator.R")  # Run interactively
```

---

### `01_data_cleaning_roi_harmonization.R`

**Purpose**: Clean raw volume data and perform quality control

**Key Functions**:
- `harmonize_roi_names()`: Standardize ROI nomenclature
- `calculate_volume()`: Convert voxel counts to volumes
- `flag_volume_outliers()`: Identify implausible measurements
- `flag_change_rate_outliers()`: Detect unrealistic rate of change
- `process_volume_data()`: Complete cleaning pipeline

**Output**:
- Cleaned dataset with QC flags
- Baseline volumes and percentage changes
- Time variables (days/weeks from RT start)

**Configuration**:
```r
VOXEL_VOLUME <- 0.001     # cm^3 per voxel
MIN_VOLUME <- 5           # Minimum plausible volume
MAX_VOLUME <- 200         # Maximum plausible volume
MAX_VOLUME_CHANGE_RATE <- 0.5  # Max fractional change per day
```

---

### `02_linear_mixed_effects_modeling.R`

**Purpose**: Fit linear mixed effects models for volume trajectories

**Key Functions**:
- `fit_lmm()`: Fit lmer models with convergence control
- `get_model_coefficients()`: Extract tidy coefficients with p-values
- `get_model_diagnostics()`: Calculate AIC, BIC, R², ICC
- `check_model_assumptions()`: Diagnostic plots
- `build_lmm_models()`: Compare multiple model specifications

**Models Tested**:
1. Null model (random intercept only)
2. Time effect (random intercept)
3. Time with random slope
4. Time + covariates

**Output**:
- Model comparison table (AIC/BIC)
- Fixed effects coefficients
- Diagnostic plots (residuals, Q-Q, scale-location)

---

### `03_logistic_regression_validation.R`

**Purpose**: Binary outcome modeling and external validation

**Key Functions**:
- `fit_logistic()`: Fit logistic regression models
- `get_odds_ratios()`: Extract ORs with confidence intervals
- `calculate_performance()`: AUC, accuracy, sensitivity, specificity
- `bootstrap_performance()`: Bootstrapped confidence intervals
- `validate_model()`: External cohort validation with plots
- `plot_roc_curve()`: ROC visualization
- `plot_calibration()`: Calibration plots

**Metrics**:
- AUC-ROC
- Accuracy, sensitivity, specificity
- Positive/negative predictive value
- Calibration assessment

---

### `04_trajectory_cluster_analysis.R`

**Purpose**: Identify distinct patterns of volume change

**Key Functions**:
- `prepare_trajectory_data()`: Convert to wide format for clustering
- `perform_kmeans_clustering()`: K-means on trajectories
- `evaluate_cluster_number()`: Elbow method for optimal k
- `plot_trajectory_clusters()`: Visualization with confidence bands
- `characterize_clusters()`: Summarize by baseline features
- `perform_trajectory_analysis()`: Complete clustering workflow

**Output**:
- Cluster assignments
- Elbow curve plot
- Trajectory visualization
- Cluster characteristics

---

## Data Format

### Input Data Requirements

Your raw data should have the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `patient_id` | character | Unique patient identifier |
| `scan_date` | Date/character | Scan date (YYYY-MM-DD) |
| `roi_name` | character | ROI name (e.g., "Prostate", "CTV_Prostate") |
| `scan_timepoint` | character | Timepoint label (e.g., "Fraction_1", "Week_4") |
| `volume` or `voxel_count` | numeric | Volume in cm³ or voxel count |
| Additional covariates | numeric/character | Age, treatment type, etc. |

### Example Data

```csv
patient_id,scan_date,roi_name,scan_timepoint,voxel_count,age,treatment_type
P001,2023-01-15,Prostate,Fraction_1,45231,68,IMRT
P001,2023-02-12,Prostate,Week_4,42187,68,IMRT
P002,2023-01-20,CTV_Prostate,Fraction_1,38942,72,SBRT
```

---

## Key Outputs

### 1. Cleaned Data
- Quality-controlled volume measurements
- Harmonized ROI names
- Calculated time variables and percentage changes

### 2. Model Results
- Fixed and random effects estimates
- Model fit statistics (AIC, BIC, R²)
- Predicted trajectories

### 3. Validation Metrics
- ROC curves and AUC
- Calibration plots
- Performance comparison (training vs validation)

### 4. Trajectory Clusters
- Distinct volume change patterns
- Cluster characteristics
- Individual and mean trajectory plots

---

## Citation

If you use this code in your research, please cite:

```
MOMENTUM Study Team (2024). 
Prostate Volume Changes with Radiotherapy: Statistical Analysis Pipeline.
GitHub repository: github.com/sianation/Prostate-volume-changes-with-radiotherapy
```

---

## License

This code is provided for academic and research purposes.

---

## Contact

For questions or issues, please open a GitHub issue or contact the MOMENTUM study team.

---

## Notes

- All scripts are designed to work with de-identified data only
- No hardcoded file paths or PHI
- Scripts include extensive documentation and examples
- Tested on R >= 4.0.0

---

**Last Updated**: December 2024

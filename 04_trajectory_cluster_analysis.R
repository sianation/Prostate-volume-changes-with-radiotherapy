################################################################################
# Script: 04_trajectory_cluster_analysis.R
# Purpose: Trajectory clustering and pattern identification
# Author: MOMENTUM Study Team
# Date: December 2024
################################################################################

# Load required packages --------------------------------------------------
library(tidyverse)
library(kml)          # For k-means longitudinal clustering
library(lcmm)         # For latent class mixed models
library(factoextra)   # For cluster visualization
library(ggrepel)      # For non-overlapping labels

# Configuration -----------------------------------------------------------

# Clustering parameters
MIN_CLUSTERS <- 2
MAX_CLUSTERS <- 6
N_START <- 50  # Number of random starts for k-means
SEED <- 123

# Helper functions --------------------------------------------------------

#' Prepare data for trajectory clustering
#'
#' @param data Long-format data with patient_id, time, and outcome
#' @param time_var Name of time variable
#' @param outcome_var Name of outcome variable
#' @return Wide-format matrix suitable for clustering
prepare_trajectory_data <- function(data, time_var = "weeks_from_rt_start", 
                                   outcome_var = "volume_change_percent") {
  
  # Convert to wide format
  wide_data <- data %>%
    select(patient_id, all_of(c(time_var, outcome_var))) %>%
    pivot_wider(
      names_from = all_of(time_var),
      values_from = all_of(outcome_var),
      names_prefix = "week_"
    )
  
  # Extract patient IDs
  patient_ids <- wide_data$patient_id
  
  # Convert to matrix
  trajectory_matrix <- wide_data %>%
    select(-patient_id) %>%
    as.matrix()
  
  rownames(trajectory_matrix) <- patient_ids
  
  message(paste("Prepared trajectories for", nrow(trajectory_matrix), "patients"))
  message(paste("Number of timepoints:", ncol(trajectory_matrix)))
  
  return(trajectory_matrix)
}

#' Perform k-means clustering on trajectories
#'
#' @param trajectory_matrix Matrix of trajectories (rows = patients, cols = timepoints)
#' @param k_range Range of cluster numbers to test
#' @param n_start Number of random starts
#' @param seed Random seed
#' @return List with kmeans results for each k
perform_kmeans_clustering <- function(trajectory_matrix, 
                                     k_range = MIN_CLUSTERS:MAX_CLUSTERS,
                                     n_start = N_START,
                                     seed = SEED) {
  
  set.seed(seed)
  
  message("\nPerforming k-means clustering...")
  
  # Cluster for each k
  kmeans_results <- map(k_range, function(k) {
    message(paste("  Testing k =", k))
    kmeans(trajectory_matrix, centers = k, nstart = n_start)
  })
  
  names(kmeans_results) <- paste0("k", k_range)
  
  return(kmeans_results)
}

#' Calculate optimal number of clusters
#'
#' @param kmeans_results List of kmeans objects
#' @return Data frame with metrics for each k
evaluate_cluster_number <- function(kmeans_results) {
  
  k_values <- as.numeric(str_extract(names(kmeans_results), "[0-9]+"))
  
  metrics <- map_df(seq_along(kmeans_results), function(i) {
    km <- kmeans_results[[i]]
    k <- k_values[i]
    
    # Within-cluster sum of squares
    wss <- km$tot.withinss
    
    # Between-cluster sum of squares
    bss <- km$betweenss
    
    # Silhouette (approximation)
    # For exact silhouette, use cluster::silhouette()
    
    data.frame(
      k = k,
      wss = wss,
      bss = bss,
      ratio = bss / km$totss  # Between-cluster variance ratio
    )
  })
  
  message("\nCluster evaluation metrics:")
  print(metrics)
  
  return(metrics)
}

#' Plot elbow curve for cluster selection
#'
#' @param metrics Data frame from evaluate_cluster_number
#' @return ggplot object
plot_elbow_curve <- function(metrics) {
  
  p <- ggplot(metrics, aes(x = k, y = wss)) +
    geom_line(size = 1.2, color = "#2E86AB") +
    geom_point(size = 3, color = "#2E86AB") +
    scale_x_continuous(breaks = metrics$k) +
    labs(
      title = "Elbow Curve for Optimal Cluster Number",
      x = "Number of Clusters (k)",
      y = "Within-Cluster Sum of Squares"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
  
  return(p)
}

#' Extract cluster assignments and add to original data
#'
#' @param data Original long-format data
#' @param kmeans_result Selected kmeans object
#' @return Data with cluster assignments
assign_clusters <- function(data, kmeans_result) {
  
  # Extract assignments
  cluster_df <- data.frame(
    patient_id = names(kmeans_result$cluster),
    cluster = as.factor(kmeans_result$cluster)
  )
  
  # Join to original data
  data_clustered <- data %>%
    left_join(cluster_df, by = "patient_id")
  
  message(paste("\nCluster sizes:"))
  print(table(cluster_df$cluster))
  
  return(data_clustered)
}

#' Plot trajectory clusters
#'
#' @param data Data with cluster assignments
#' @param time_var Time variable name
#' @param outcome_var Outcome variable name
#' @param add_individual Show individual trajectories (default: TRUE)
#' @return ggplot object
plot_trajectory_clusters <- function(data, 
                                    time_var = "weeks_from_rt_start",
                                    outcome_var = "volume_change_percent",
                                    add_individual = TRUE) {
  
  # Calculate cluster means
  cluster_means <- data %>%
    group_by(cluster, .data[[time_var]]) %>%
    summarise(
      mean_value = mean(.data[[outcome_var]], na.rm = TRUE),
      se = sd(.data[[outcome_var]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Base plot
  p <- ggplot(data, aes(x = .data[[time_var]], y = .data[[outcome_var]], 
                        color = cluster, group = interaction(patient_id, cluster)))
  
  # Add individual trajectories if requested
  if (add_individual) {
    p <- p + geom_line(alpha = 0.1, size = 0.3)
  }
  
  # Add cluster means
  p <- p +
    geom_line(
      data = cluster_means,
      aes(x = .data[[time_var]], y = mean_value, group = cluster),
      size = 2,
      inherit.aes = FALSE,
      mapping = aes(color = cluster)
    ) +
    geom_point(
      data = cluster_means,
      aes(x = .data[[time_var]], y = mean_value),
      size = 3,
      inherit.aes = FALSE,
      mapping = aes(color = cluster)
    ) +
    geom_ribbon(
      data = cluster_means,
      aes(x = .data[[time_var]], 
          ymin = mean_value - 1.96*se, 
          ymax = mean_value + 1.96*se,
          fill = cluster,
          group = cluster),
      alpha = 0.2,
      inherit.aes = FALSE
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "Volume Change Trajectory Clusters",
      x = "Weeks from RT Start",
      y = "Volume Change (%)",
      color = "Cluster",
      fill = "Cluster"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "right"
    )
  
  return(p)
}

#' Characterize clusters by baseline features
#'
#' @param data Data with cluster assignments and baseline covariates
#' @param cluster_var Name of cluster variable
#' @param covariate_names Vector of covariate names
#' @return Data frame with cluster characteristics
characterize_clusters <- function(data, 
                                 cluster_var = "cluster",
                                 covariate_names) {
  
  # Get unique patient-level data
  patient_data <- data %>%
    group_by(patient_id) %>%
    slice(1) %>%
    ungroup()
  
  # Summarize by cluster
  cluster_chars <- patient_data %>%
    group_by(.data[[cluster_var]]) %>%
    summarise(
      n = n(),
      across(
        all_of(covariate_names),
        list(
          mean = ~mean(., na.rm = TRUE),
          sd = ~sd(., na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
  
  message("\nCluster characteristics:")
  print(cluster_chars)
  
  return(cluster_chars)
}

# Main workflow -----------------------------------------------------------

#' Complete trajectory clustering analysis
#'
#' @param data Long-format cleaned data
#' @param time_var Time variable name
#' @param outcome_var Outcome variable name
#' @param k Specific number of clusters (if NULL, will optimize)
#' @param covariate_names Covariates for cluster characterization
#' @param plot_dir Directory to save plots
#' @return List with results
perform_trajectory_analysis <- function(data,
                                       time_var = "weeks_from_rt_start",
                                       outcome_var = "volume_change_percent",
                                       k = NULL,
                                       covariate_names = NULL,
                                       plot_dir = NULL) {
  
  message("\nStarting trajectory clustering analysis...")
  
  # Prepare data
  trajectory_matrix <- prepare_trajectory_data(data, time_var, outcome_var)
  
  # Perform clustering
  kmeans_results <- perform_kmeans_clustering(trajectory_matrix)
  
  # Evaluate cluster numbers
  metrics <- evaluate_cluster_number(kmeans_results)
  
  # Plot elbow curve
  elbow_plot <- plot_elbow_curve(metrics)
  
  # Select k if not specified
  if (is.null(k)) {
    # Simple heuristic: find elbow using second derivative
    metrics <- metrics %>%
      arrange(k) %>%
      mutate(
        wss_diff = lag(wss) - wss,
        wss_diff2 = lag(wss_diff) - wss_diff
      )
    k <- metrics$k[which.max(metrics$wss_diff2)]
    message(paste("\nOptimal k selected:", k))
  }
  
  # Get selected clustering
  selected_clustering <- kmeans_results[[paste0("k", k)]]
  
  # Assign clusters
  data_with_clusters <- assign_clusters(data, selected_clustering)
  
  # Plot trajectories
  trajectory_plot <- plot_trajectory_clusters(data_with_clusters, time_var, outcome_var)
  
  # Characterize clusters
  cluster_characteristics <- NULL
  if (!is.null(covariate_names)) {
    cluster_characteristics <- characterize_clusters(
      data_with_clusters,
      covariate_names = covariate_names
    )
  }
  
  # Save plots
  if (!is.null(plot_dir)) {
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(plot_dir, "elbow_curve.png"), elbow_plot, width = 8, height = 6)
    ggsave(file.path(plot_dir, "trajectory_clusters.png"), trajectory_plot, width = 12, height = 8)
  }
  
  return(list(
    data = data_with_clusters,
    kmeans = selected_clustering,
    metrics = metrics,
    characteristics = cluster_characteristics,
    plots = list(elbow = elbow_plot, trajectories = trajectory_plot)
  ))
}

# Example usage (uncomment to run) ----------------------------------------
# # Load cleaned data
# data <- read_csv("data/cleaned_volume_data.csv") %>%
#   filter(qc_flag == "pass")
# 
# # Perform clustering
# results <- perform_trajectory_analysis(
#   data = data,
#   time_var = "weeks_from_rt_start",
#   outcome_var = "volume_change_percent",
#   k = NULL,  # Auto-select optimal k
#   covariate_names = c("baseline_volume", "age"),
#   plot_dir = "results/trajectory_analysis"
# )
# 
# # Access results
# head(results$data)
# results$characteristics

################################################################################
# End of script
################################################################################

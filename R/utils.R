#' Return a List of Named NA Metrics for Zero-Point Tiles
#'
#' Provides default NA values for graph or canopy metrics when input data is insufficient.
#'
#' @param type Character string, either `"graph"`, `"canopy"` or `"trees"`. Determines which metric structure to return.
#'
#' @return A named list of metrics filled with `NA_real_` values appropriate to the type.
#' @export
named_zero_metrics <- function(type = 'graph') {
  if (type == 'graph') {
    return(list(
      mean_degree = 0,
      mean_strength = 0,
      mean_betweenness = 0,
      mean_closeness = NA_real_,
      n_components = NA_real_,
      avg_path_length = NA_real_,
      eigen_ratio = NA_real_,
      graph_density = 0,
      mean_abs_sa = NA_real_,
      n_m2 = 0
    ))
  } else if (type == 'canopy') {
    return(list(
      fractional_canopy_cover = 0,
      rumple_index = 0,
      LAI = 0,
      LAD_max = 0,
      LAD_mean = 0,
      LAD_z_max = 0,
      n_m2 = 0
    ))
  } else if (type == 'trees') {
    return(list(
      n_trees = 0,
      trees_per_acre = 0,
      n_strata_low = 0,
      n_strata_low_mid = 0,
      n_strata_mid = 0,
      n_strata_mid_upper = 0,
      n_strata_upper = 0,
      n_gt_6_1 = 0,
      n_gt_12_1 = 0,
      n_gt_24_1 = 0
    ))
  } else {
    stop("Invalid type. Must be 'graph' or 'canopy'.")
  }
}

#' Load Required Graph and LiDAR Dependencies
#'
#' Loads necessary packages used in tree detection and voxel-based graph computation.
#' Suppresses loading messages.
#'
#' @return Invisibly returns `TRUE` if all packages are successfully loaded.
#' @export
load_graph_deps <- function() {
  pkgs <- c("lidR", "dbscan", "igraph", "dplyr", "tibble", "purrr", "stats", "geometry", "lidar.metrics")
  for (pkg in pkgs) {
    suppressMessages(require(pkg, character.only = TRUE))
  }
  invisible(TRUE)
}


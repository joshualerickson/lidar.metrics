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
      vertical_complexity_index = 0,
      LAI = 0,
      LAD_max = 0,
      LAD_mean = 0,
      LAD_z_max = 0
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
      n_gt_24_1 = 0,
      mean_canopy_area = 0,
      median_canopy_area = 0,
      min_canopy_area = 0,
      max_canopy_area = 0,
      sd_canopy_area = 0,
      cv_canopy_area = 0,
      topo_residual_sd = NA_real_,
      topo_entropy = NA_real_,
      smoothness_score = NA_real_
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
  pkgs <- c("lidR", "dbscan", "igraph", "dplyr", "tibble", "purrr", "stats", "geometry")
  for (pkg in pkgs) {
    suppressMessages(require(pkg, character.only = TRUE))
  }
  invisible(TRUE)
}


#' Compute Voxel-Level Structure Metrics
#'
#' Calculates structural metrics for a voxel such as convex hull volume
#' (normalized), spatial variance, and point count.
#'
#' @param z Z coordinates (height) within the voxel.
#' @param x X coordinates within the voxel.
#' @param y Y coordinates within the voxel.
#'
#' @return A named list with `convex_hull_vol`, `norm_var`, and `point_count`.
#' @export

voxel_structure_metrics <- function(z, x, y) {
  coords <- cbind(x, y, z)

  # If not enough points, return fallback values
  if (nrow(coords) < 4) {
    var_xyz <- mean(apply(coords, 2, var), na.rm = TRUE)
    return(list(
      convex_hull_vol = 0,
      norm_var = var_xyz,
      point_count = nrow(coords)
    ))
  }

  # Try computing 3D convex hull volume
  hull_volume <- tryCatch({
    ch <- geometry::convhulln(coords, options = "FA")
    as.numeric(ch$vol)
  }, error = function(e) 0)

  # Compute normalized variance of spatial coordinates
  var_xyz <- mean(apply(coords, 2, var), na.rm = TRUE)

  return(list(
    convex_hull_vol = log1p(hull_volume/nrow(coords)),
    norm_var = var_xyz,
    point_count = nrow(coords)
  ))
}


#' Return a List of Named NA Metrics for Zero-Point Tiles
#'
#' Provides default NA values for graph or canopy metrics when input data is insufficient.
#'
#' @param type Character string, either `"graph"` or `"canopy"`. Determines which metric structure to return.
#'
#' @return A named list of metrics filled with `NA_real_` values appropriate to the type.
#' @export
named_zero_metrics <- function(type = 'graph') {
  if (type == 'graph') {
    return(list(
      mean_weighted_degree = 0,
      mean_betweenness = 0,
      mean_closeness = 0,
      n_components = 0,
      avg_path_length = 0,
      eigen_ratio = 0,
      graph_density = 0
    ))
  } else if (type == 'canopy') {
    return(list(
      n_trees = 0,
      trees_per_acre = 0,
      n_low = 0,
      n_low_mid = 0,
      n_mid = 0,
      n_mid_upper = 0,
      n_upper = 0,
      topo_residual_sd = NA_real_,
      topo_entropy = NA_real_,
      fractional_canopy_cover = NA_real_,
      LAI = NA_real_,
      LAD_max = NA_real_,
      LAD_mean = NA_real_,
      LAD_z_max = NA_real_
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
  pkgs <- c("lidR", "dbscan", "igraph", "dplyr", "tibble", "purrr", "stats")
  for (pkg in pkgs) {
    suppressMessages(require(pkg, character.only = TRUE))
  }
  invisible(TRUE)
}


#' Remove Named Objects Safely from Environment and Trigger Garbage Collection
#'
#' Removes objects passed via `...` if they exist in the current environment.
#'
#' @param ... Character names of objects to remove.
#'
#' @return Called for side effects. Triggers `gc()`.
#' @export
safe_cleanup <- function(...) {
  for (obj in list(...)) {
    if (exists(obj, inherits = FALSE)) {
      rm(list = obj, inherits = FALSE)
    }
  }
  gc()
}

#' Save CRS (Coordinate Reference System) to a Folder as Text File
#'
#' Writes the string representation of a CRS object to a `crs.txt` file in a folder.
#'
#' @param crs An object of class `crs` or string. Usually obtained via `terra::crs()` or `sf::st_crs()`.
#' @param out_dir Output directory path where the `crs.txt` file should be saved.
#'
#' @return Called for side effects. Creates a file at `out_dir/crs.txt`.
#' @export
save_crs_to_folder <- function(crs, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  crs_txt <- file.path(out_dir, "crs.txt")
  writeLines(as.character(crs), crs_txt)
}

#' @importFrom magrittr %>%
NULL


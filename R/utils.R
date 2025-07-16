
named_zero_metrics <- function(type = 'graph') {
  if(type == 'graph') {
    list(
      
      mean_weighted_degree = NA_real_,
      mean_betweenness = NA_real_,
      mean_closeness = NA_real_,
      n_components = NA_real_,
      avg_path_length = NA_real_,
      eigen_ratio = NA_real_,
      graph_density = NA_real_
    )
    
  } else if(type == 'canopy') {
    list(
      n_trees = NA_real_,
      trees_per_acre = NA_real_,
      n_low = NA_real_,
      n_low_mid = NA_real_,
      n_mid = NA_real_,
      n_mid_upper = NA_real_,
      n_upper = NA_real_,
      topo_residual_sd = NA_real_,
      topo_entropy = NA_real_,
      #rumple = NA_real_,
      fractional_canopy_cover = NA_real_,
      LAI = NA_real_,
      LAD_max = NA_real_,
      LAD_mean = NA_real_,
      LAD_z_max = NA_real_
    )
  }
}
load_graph_deps <- function() {
  pkgs <- c("lidR", "dbscan", "igraph", "dplyr", "tibble", "purrr", "stats")
  for (pkg in pkgs) {
    suppressMessages(require(pkg, character.only = TRUE))
  }
}

safe_cleanup <- function(...) {
  for (obj in list(...)) {
    if (exists(obj, inherits = FALSE)) rm(list = obj, inherits = FALSE)
  }
  gc()
}



save_crs_to_folder <- function(crs, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  crs_txt <- file.path(out_dir, "crs.txt")
  writeLines(as.character(crs), crs_txt)
}
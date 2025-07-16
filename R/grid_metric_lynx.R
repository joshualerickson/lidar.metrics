#' Compute Graph-Based Connectivity Metrics in Vertical Bins
#'
#' Calculates graph-theoretic metrics (e.g., betweenness, closeness, path length)
#' for voxels within user-defined height bins. This is intended to assess forest
#' structure connectivity relevant to wildlife habitat (e.g., for lynx).
#'
#' @param x X coordinates of LiDAR points.
#' @param y Y coordinates of LiDAR points.
#' @param z Z (height) values of LiDAR points.
#' @param edge_thresh_values Numeric vector of edge thresholds (e.g., `[3, 3]`) for each vertical bin.
#' @param z_1 Minimum height threshold (default: `0.3`).
#' @param z_20 Threshold separating bin 1 from bin 2 (default: `6.1`).
#' @param z_40 Maximum height for bin 2 (default: `12.1`).
#' @param voxel_res Numeric, voxel resolution in the XY plane (default: `3`).
#'
#' @return A named list of flattened graph metrics from two vertical bins.
#' @export
connectivity_metrics_binned <- function(x, y, z,
                                        edge_thresh_values = c(3, 3),
                                        z_1 = 0.3,
                                        z_20 = 6.1,
                                        z_40 = 12.1,
                                        voxel_res = 3) {

  metric_names <- c(
    "mean_weighted_degree",
    "mean_betweenness",
    "mean_closeness",
    "n_components",
    "avg_path_length",
    "eigen_ratio",
    "graph_density"
  )

  bins <- list(
    list(zmin = z_1, zmax = z_20, edge_thresh = edge_thresh_values[1], voxel_res = voxel_res, prefix = "bin_1_20_"),
    list(zmin = z_20, zmax = z_40, edge_thresh = edge_thresh_values[2], voxel_res = voxel_res, prefix = "bin_20_40_")
  )

  las_all <- suppressMessages(LAS(data.frame(X = x, Y = y, Z = z)))

  results <- purrr::map(bins, function(b) {
    tryCatch({
      out <- compute_graph_metrics(
        las = las_all,
        z_min = b$zmin,
        z_max = b$zmax,
        edge_thresh = b$edge_thresh,
        voxel_res = b$voxel_res
      )
      setNames(out, paste0(b$prefix, names(out)))
    }, error = function(cond) {
      setNames(named_zero_metrics(), paste0(b$prefix, names(named_zero_metrics())))
    })
  })

  return(do.call(c, results))
}


#' Compute Graph Metrics for a Filtered LAS Within a Vertical Bin
#'
#' Creates a voxel-based graph using a fixed edge threshold, and calculates graph
#' metrics such as weighted degree, betweenness, closeness, and component structure.
#'
#' @param las A LAS object (point cloud) from the `lidR` package.
#' @param z_min Minimum Z height threshold.
#' @param z_max Maximum Z height threshold.
#' @param edge_thresh Numeric, distance threshold for creating edges.
#' @param voxel_res Numeric, XY voxel resolution.
#'
#' @return A named list of graph-theoretic metrics.
#' @export
compute_graph_metrics <- function(las, z_min, z_max, edge_thresh, voxel_res) {
  load_graph_deps()

  las_filtered <- filter_poi(las, Z > z_min & Z <= z_max)
  if (is.empty(las_filtered) || length(las_filtered@data$Z) < 2) {
    return(named_zero_metrics())
  }

  count_voxel <- voxel_metrics(las_filtered, ~list(point_count = length(Z)), res = voxel_res)
  voxel_df <- as.data.frame(count_voxel)

  required_cols <- c("X", "Y", "Z", "point_count")
  if (!all(required_cols %in% names(voxel_df))) {
    return(named_zero_metrics())
  }

  idx <- which(voxel_df$Z > z_min & voxel_df$Z <= z_max)
  if (length(idx) < 3) return(named_zero_metrics())

  coords <- voxel_df[idx, c("X", "Y", "Z")]
  coords <- as.matrix(data.frame(lapply(coords, as.numeric)))
  pt_counts <- voxel_df$point_count[idx]

  nn <- dbscan::frNN(coords, eps = edge_thresh)
  edges <- tibble(
    from = rep(seq_along(nn$id), lengths(nn$id)),
    to = unlist(nn$id)
  ) %>%
    filter(from != to, from < to) %>%
    mutate(weight = pt_counts[from] + pt_counts[to])

  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  E(g)$weight <- edges$weight
  E(g)$inv_weight <- 1 / pmax(E(g)$weight, 1e-6)

  comps <- components(g)
  results <- list(
    mean_weighted_degree = mean(strength(g, weights = E(g)$inv_weight), na.rm = TRUE),
    mean_betweenness = tryCatch(mean(betweenness(g, weights = E(g)$inv_weight), na.rm = TRUE), error = function(e) 0),
    mean_closeness = tryCatch(mean(closeness(g, weights = E(g)$inv_weight), na.rm = TRUE), error = function(e) 0),
    n_components = comps$no,
    avg_path_length = tryCatch({
      d <- distances(g, weights = E(g)$inv_weight)
      mean(d[is.finite(d) & d > 0], na.rm = TRUE)
    }, error = function(e) 0),
    eigen_ratio = tryCatch({
      pca <- prcomp(coords)
      pca$sdev[1] / pca$sdev[2]
    }, error = function(e) 0),
    graph_density = tryCatch(edge_density(g), error = function(e) 0)
  )

  return(results)
}


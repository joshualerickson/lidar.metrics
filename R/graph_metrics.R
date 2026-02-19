#' Compute Binned Graph Connectivity Metrics from LiDAR Point Cloud
#'
#' Calculates graph-theoretic metrics (e.g., betweenness, closeness, path length)
#' for voxelized LiDAR data within vertical height bins. This function is used to
#' assess 3D structural connectivity (e.g., for wildlife habitat like lynx).
#'
#' @param x X coordinates of LiDAR points.
#' @param y Y coordinates of LiDAR points.
#' @param z Z (height) values of LiDAR points.
#' @param edge_thresh_values Numeric vector of edge thresholds (e.g., c(3, 3)) for each vertical bin.
#' @param z_1 Lower height threshold (default: 1 m).
#' @param z_20 Midstory lower bound (default: 6 m).
#' @param z_40 Midstory upper bound (default: 12 m).
#' @param voxel_res Numeric, voxel resolution in XY plane (default: 3).
#' @param psid unquoted `PointSourceID` from the LAS file.
#' @param QL1 logical for whether to QL1 methods or not. Default (FALSE)
#' @param n Numeric to pass to `lidR::random_per_voxel()` density argument.
#' @param res Numeric to pass to `lidR::random_per_voxel()` res argument.
#'
#' @return A named list of graph-theoretic metrics for each bin (understory, midstory).
#' @note The bounds for `z_1, ,z_20, and z_40` can be whatever you want but with always be labelled `understory_` and `midstory_`.
#' @export

connectivity_metrics_binned <- function(x, y, z,
                                        edge_thresh_values = c(3, 3),
                                        z_1 = 1,
                                        z_20 = 6,
                                        z_40 = 12,
                                        voxel_res = 3,
                                        psid,
                                        QL1 = FALSE,
                                        n = 1,
                                        res = 3) {

  bins <- list(
    list(zmin = z_1, zmax = z_20, edge_thresh = edge_thresh_values[1],
         voxel_res = voxel_res, prefix = "understory_", QL1 = QL1, n = n, res = res),
    list(zmin = z_20, zmax = z_40, edge_thresh = edge_thresh_values[2],
         voxel_res = voxel_res, prefix = "midstory_", QL1 = QL1, n = n, res = res)
  )


  las_all <- suppressMessages(lidR::LAS(data.frame(X = x,
                                                   Y = y,
                                                   Z = z,
                                                   PointSourceID = psid)))

  results <- purrr::map(bins, function(b) {
    tryCatch({
      out <- compute_graph_metrics(
        las = las_all,
        z_min = b$zmin,
        z_max = b$zmax,
        edge_thresh = b$edge_thresh,
        voxel_res = b$voxel_res,
        b$QL1,
        b$n,
        b$res
      )
      setNames(out, paste0(b$prefix, names(out)))
    }, error = function(cond) {
      setNames(named_zero_metrics(), paste0(b$prefix, names(named_zero_metrics())))
    })
  })


  return(do.call(c, results))
}


#' Compute Graph Metrics for a Voxelized Vertical Bin
#'
#' Calculates key graph-theoretic metrics for a filtered LAS object, including
#' strength, degree, betweenness, closeness, and more. Uses convex hull volume
#' as edge weights.
#'
#' @param las A LAS object (LiDAR point cloud) from the `lidR` package.
#' @param z_min Minimum Z height to include in bin.
#' @param z_max Maximum Z height to include in bin.
#' @param edge_thresh Distance threshold (in meters) for connecting centroids.
#' @param voxel_res Numeric voxel resolution for binning.
#' @param QL1 logical for whether to QL1 methods or not.
#' @param n Numeric to pass to `lidR::random_per_voxel()` density argument.
#' @param res Numeric to pass to `lidR::random_per_voxel()` res argument.
#'
#' @return A named list of graph metrics.
#' @export

compute_graph_metrics <- function(las, z_min, z_max, edge_thresh, voxel_res, QL1, n, res) {

  load_graph_deps()

  las_filtered <- lidR::filter_poi(las, Z >= z_min & Z <= z_max)

  if(QL1) {


  las_filtered <- lidR::decimate_points(las_filtered, algorithm = lidR::random_per_voxel(n = n, res = res))


  } else {


  psid_vec <- las_filtered@data$PointSourceID

  keep_ids <- as.integer(names(which.max(table(psid_vec))))

  las_filtered <- lidR::filter_poi(las_filtered, PointSourceID %in% keep_ids)


  }



  if (is.empty(las_filtered) || length(las_filtered@data$Z) < 2) {
    return(named_zero_metrics())
  }

  voxel_df <- lidR::voxel_metrics(
    las_filtered,
    func = ~lidar.metrics::voxel_structure_metrics(Z, X, Y),
    res = voxel_res
  )


  voxel_df <- as.data.frame(voxel_df)

  # Reorder to avoid idx before defining it
  idx <- which(voxel_df$Z > z_min & voxel_df$Z <= z_max)
  if (length(idx) < 3) return(named_zero_metrics())

  coords <- voxel_df[idx, c("X", "Y", "Z")]
  coords <- as.matrix(data.frame(lapply(coords, as.numeric)))


  required_cols <- c("X", "Y", "Z")
  if (!all(required_cols %in% names(voxel_df))) {
    return(named_zero_metrics())
  }

  nn <- dbscan::frNN(coords, eps = edge_thresh)


  vol <- voxel_df$convex_hull_vol[idx]

  edges <- dplyr::tibble(
    from = rep(seq_along(nn$id), lengths(nn$id)),
    to = unlist(nn$id)
  ) %>%
    dplyr::filter(from != to, from < to) %>%
    dplyr::mutate(weight = vol[from] + vol[to])

  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  E(g)$weight <- edges$weight
  E(g)$inv_weight <- 1 /  pmax(log1p(E(g)$weight), 1e-3)

  comps <- igraph::components(g)

  results <- list(
    mean_degree = mean(igraph::degree(g), na.rm = TRUE),
    mean_strength = mean(igraph::strength(g, weights = E(g)$inv_weight), na.rm = TRUE),
    mean_betweenness = tryCatch({
      mean(igraph::betweenness(g, weights = E(g)$inv_weight, normalized = FALSE), na.rm = TRUE)
    }, error = function(e) 0),

    mean_closeness = tryCatch({
      closeness_vals <- igraph::closeness(g, weights = E(g)$inv_weight, normalized = FALSE)
      # Normalize manually to [0, 1]
      if (max(closeness_vals, na.rm = TRUE) > 0) {
        closeness_norm <- closeness_vals / max(closeness_vals, na.rm = TRUE)
        mean(closeness_norm, na.rm = TRUE)
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_),

    n_components = comps$no,

    avg_path_length = tryCatch({
      d <- igraph::distances(g, weights = E(g)$inv_weight)
      d <- d[is.finite(d) & d > 0]
      if (length(d) > 0) mean(d) else 0
    }, error = function(e) 0),

    eigen_ratio = tryCatch({
      pca <- prcomp(coords)
      ratio <- pca$sdev[1] / max(pca$sdev[2], 1e-6)
      pmin(ratio, 3)  # Cap to avoid wild outliers
    }, error = function(e) NA_real_),


    graph_density = tryCatch(igraph::edge_density(g), error = function(e) 0),
    n_m2 = length(las_filtered@data$Z)/(30*30)
  )

  return(results)
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


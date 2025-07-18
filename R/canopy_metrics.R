#' Compute Canopy Cover and Tree Metrics
#'
#' Calculates canopy cover, tree strata metrics, and leaf area density metrics
#' for a small area of LiDAR point cloud data.
#'
#' @param x Numeric vector of X coordinates.
#' @param y Numeric vector of Y coordinates.
#' @param z Numeric vector of Z (height) values.
#' @param return_number Integer vector of return numbers (e.g., from LAS file).
#' @param window_func Function that returns variable window size for tree detection (e.g., based on height).
#' @param out_dir Directory path for saving optional outputs (currently unused).
#'
#' @return A named list with metrics including:
#' \itemize{
#'   \item \code{fractional_canopy_cover}
#'   \item \code{n_trees}, \code{trees_per_acre}, and strata-specific counts
#'   \item \code{LAI}, \code{LAD_max}, \code{LAD_mean}, \code{LAD_z_max}
#'   \item \code{topo_residual_sd} and \code{topo_entropy}
#' }
#' @export
canopy_cover_metrics <- function(x, y, z, return_number,
                                 window_func, out_dir) {
  if (length(z) < 5) return(named_zero_metrics(type = 'canopy'))

  trees <- tree_detection(x, y, z, window_func = window_func, out_dir = out_dir)

  first_above <- sum(z > 2 & return_number == 1)
  total_first <- sum(return_number == 1)
  canopy_cover <- ifelse(total_first == 0, NA_real_, first_above / total_first)
  cc <- list(fractional_canopy_cover = canopy_cover)

  lad <- lad_metrics(z, dz = 1, z0 = 2)

  return(c(trees, lad, cc))
}


#' Detect Trees and Compute Strata Metrics
#'
#' Detects individual trees using a variable window local maxima filter
#' and classifies trees into vertical strata. Also fits a spline surface
#' to derive topographic residual SD and entropy.
#'
#' @param x Numeric vector of X coordinates.
#' @param y Numeric vector of Y coordinates.
#' @param z Numeric vector of Z (height) values.
#' @param window_func Function that returns a window size for tree detection.
#' @param out_dir Path for optional output directory (currently unused).
#'
#' @return A named list of tree metrics:
#' \itemize{
#'   \item \code{n_trees}, \code{trees_per_acre}
#'   \item \code{n_low}, \code{n_low_mid}, \code{n_mid}, \code{n_mid_upper}, \code{n_upper}
#'   \item \code{topo_residual_sd}, \code{topo_entropy}
#' }
#' @export
tree_detection <- function(x, y, z,
                           window_func = function(x) {
                             w <- pmax(x, 1)
                             w <- log1p(w) * 1.5 + 2
                             w <- pmin(w, 6)
                             w[is.na(w) | w <= 0] <- 3
                             return(w)
                           },
                           out_dir) {
  load_graph_deps()

  tryCatch({
    if (length(x) < 5 || length(z) < 5) {
      return(setNames(rep(0, 9), c("n_trees", "trees_per_acre",
                                   "n_low", "n_low_mid", "n_mid",
                                   "n_mid_upper", "n_upper",
                                   "topo_residual_sd", "topo_entropy")))
    }

    las_data <- suppressMessages(LAS(data.frame(X = x, Y = y, Z = z)))

    trees <- tryCatch({
      lidR::locate_trees(las_data, lidR::lmf(ws = window_func))
    }, error = function(e) {
      message("locate_trees error: ", e$message)
      return(NULL)
    })

    if (is.null(trees) || nrow(trees) == 0) {
      return(setNames(rep(0, 9), c("n_trees", "trees_per_acre",
                                   "n_low", "n_low_mid", "n_mid",
                                   "n_mid_upper", "n_upper",
                                   "topo_residual_sd", "topo_entropy")))
    }

    coords <- sf::st_coordinates(trees)

    if(!is.null(out_dir)){
    out_path <- file.path(out_dir, paste0("X", round(mean(range(x)), 3), "_Y", round(mean(range(y)), 3), '.fst'))
    fst::write_fst(as.data.frame(coords), out_path)
    }
    z_vals <- trees$Z
    zbin <- cut(z_vals,
                breaks = c(0.3, 6.1, 12.1, 24.1, 36.1, Inf),
                labels = c("low", "low_mid", "mid", "mid_upper", "upper"))

    strata <- c("low", "low_mid", "mid", "mid_upper", "upper")
    strata_table <- table(zbin)
    strata_counts <- setNames(rep(0, length(strata)), strata)
    strata_counts[names(strata_table)] <- as.integer(strata_table)

    df <- data.frame(x = coords[, 1], y = coords[, 2], z = z_vals)
    residual_sd <- NA_real_
    topo_entropy <- NA_real_

    if (nrow(df) >= 5) {
      tryCatch({
        fit <- lm(z ~ splines::bs(x, df = 4) + splines::bs(y, df = 4), data = df)
        residuals <- df$z - predict(fit, newdata = df)
        residual_sd <- sd(residuals)
        h <- hist(residuals, breaks = 10, plot = FALSE)
        p <- h$counts / sum(h$counts)
        topo_entropy <- -sum(p * log(p + 1e-10))
      }, error = function(e) {
        message("Spline fit error: ", e$message)
      })
    }

    n_trees <- nrow(trees)
    trees_per_acre <- n_trees / 0.222395

    return(list(
      n_trees = as.numeric(n_trees),
      trees_per_acre = as.numeric(trees_per_acre),
      n_low = strata_counts["low"],
      n_low_mid = strata_counts["low_mid"],
      n_mid = strata_counts["mid"],
      n_mid_upper = strata_counts["mid_upper"],
      n_upper = strata_counts["upper"],
      topo_residual_sd = residual_sd,
      topo_entropy = topo_entropy
    ))
  }, error = function(e) {
    message("tree_detection() general error: ", e$message)
    return(setNames(rep(NA_real_, 9), c("n_trees", "trees_per_acre",
                                        "n_low", "n_low_mid", "n_mid",
                                        "n_mid_upper", "n_upper",
                                        "topo_residual_sd", "topo_entropy")))
  })
}

#' Leaf Area Density Metrics
#'
#' Computes LAD and LAI metrics from vertical profile of LiDAR Z values.
#'
#' @param z Numeric vector of Z heights.
#' @param dz Vertical resolution of profile bins (default = 1).
#' @param k Extinction coefficient (default = 0.5).
#' @param z0 Minimum Z height considered vegetation (default = 2).
#'
#' @return A named list of metrics:
#' \itemize{
#'   \item \code{LAI} Leaf Area Index
#'   \item \code{LAD_max}, \code{LAD_mean}
#'   \item \code{LAD_z_max} Height of maximum LAD
#' }
#' @export
lad_metrics <- function(z, dz = 1, k = 0.5, z0 = 2) {
  tryCatch({
    profile <- lidR::LAD(z, dz = dz, k = k, z0 = z0)

    if (nrow(profile) == 0) {
      return(setNames(rep(NA_real_, 4), c("LAI", "LAD_max", "LAD_mean", "LAD_z_max")))
    }

    LAI <- sum(profile$lad * dz, na.rm = TRUE)
    LAD_max <- max(profile$lad, na.rm = TRUE)
    LAD_mean <- mean(profile$lad, na.rm = TRUE)
    LAD_z_max <- profile$z[which.max(profile$lad)]

    return(list(
      LAI = LAI,
      LAD_max = LAD_max,
      LAD_mean = LAD_mean,
      LAD_z_max = LAD_z_max
    ))
  }, error = function(e) {
    message("LAD error: ", e$message)
    return(setNames(rep(NA_real_, 4), c("LAI", "LAD_max", "LAD_mean", "LAD_z_max")))
  })
}

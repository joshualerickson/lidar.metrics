
#' Detect Trees and Compute Structural Metrics from a LiDAR Tile
#'
#' Applies a local maxima filter with an adaptive window function to detect trees,
#' classifies them into vertical strata. Designed to be used within
#' `lidR::pixel_metrics()` for tile-based LAScatalog processing.
#'
#' @param x Numeric vector of X coordinates of LiDAR points.
#' @param y Numeric vector of Y coordinates of LiDAR points.
#' @param z Numeric vector of Z (height) values of LiDAR points.
#' @param window_func A function that returns a scalar window size given a vector of Z values.
#'        Used in `lidR::lmf()` for local maxima filtering. The default function is based
#'        on log-scaling and height-based branching logic.
#'
#' @return A named list containing:
#' \describe{
#'   \item{n_trees}{Total number of detected trees.}
#'   \item{trees_per_acre}{Number of detected trees per 0.222395 ha (≈0.55 acres).}
#'   \item{n_low, n_low_mid, n_mid, n_mid_upper, n_upper}{Counts of trees in fixed height strata: 0.3–6.1 m, 6.1–12.1 m, 12.1–24.1 m, 24.1–36.1 m, >36.1 m.}
#'   \item{n_0_6_1, n_0_12_1, n_0_24_1, n_0_36_1, n_0_Inf}{Cumulative counts of trees taller than 0.3 m up to given height thresholds.}
#'   \item{n_gt_6_1, n_gt_12_1, n_gt_24_1, n_gt_36_1}{Counts of trees taller than each specified threshold.}
#' }
#'
#' @seealso \code{\link[lidR]{lmf}}, \code{\link[lidR]{locate_trees}}, \code{\link{named_zero_metrics}}
#' @export
tree_detection <- function(las,
                           window_func = function(x) {ifelse(x*0.25 < 1, 1, x*0.25)}) {
  load_graph_deps()

   tryCatch({
    if (is.empty(las)) {
      return(named_zero_metrics('trees'))
    }

    trees <- tryCatch({
      lidR::locate_trees(las, lidR::lmf(ws = window_func))
    }, error = function(e) {
      return(NULL)

    })

    if (is.null(trees) || nrow(trees) == 0) {
      return(named_zero_metrics('trees'))
    }

    coords <- sf::st_coordinates(trees)
    z_vals <- trees$Z

    # Reverse cumulative bins: height >= threshold
    n_gt_6_1    <- as.numeric(sum(z_vals >= 6.1))
    n_gt_12_1   <- as.numeric(sum(z_vals >= 12.1))
    n_gt_24_1   <- as.numeric(sum(z_vals >= 24.1))


    # Stratified bins (non-cumulative)
    zbin <- cut(z_vals,
                breaks = c(0.3, 6.1, 12.1, 24.1, 36.1, Inf),
                labels = c("low", "low_mid", "mid", "mid_upper", "upper"))

    strata <- c("low", "low_mid", "mid", "mid_upper", "upper")
    strata_table <- table(zbin)
    strata_counts <- setNames(rep(0, length(strata)), strata)
    strata_counts[names(strata_table)] <- as.numeric(strata_table)

    n_trees <- nrow(trees)
    trees_per_acre <- n_trees / 0.222395

    return(list(
      n_trees = as.numeric(n_trees),
      trees_per_acre = as.numeric(trees_per_acre),
      n_strata_low = strata_counts["low"],
      n_strata_low_mid = strata_counts["low_mid"],
      n_strata_mid = strata_counts["mid"],
      n_strata_mid_upper = strata_counts["mid_upper"],
      n_strata_upper = strata_counts["upper"],
      n_gt_6_1 = n_gt_6_1,
      n_gt_12_1 = n_gt_12_1,
      n_gt_24_1 = n_gt_24_1
    ))
  }, error = function(e) {
    return(named_zero_metrics('trees'))
  })

}


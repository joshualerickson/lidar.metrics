
#' Detect Trees and Compute Structural Metrics from a LiDAR Tile
#'
#' Applies a local maxima filter with an adaptive window function to detect trees,
#' classifies them into vertical strata, and calculates both vertical structural
#' complexity and spatial smoothness metrics. Designed to be used within
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
#'   \item{topo_residual_sd}{Standard deviation of residuals from a 2D spline surface fit to tree heights (spatial topographic variation).}
#'   \item{topo_entropy}{Shannon entropy of the residual distribution (complexity of vertical deviation).}
#'   \item{smoothness_score}{Composite structural score combining scaled `topo_residual_sd` and `topo_entropy`.}
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

    if (is.null(chm)) return(named_zero_metrics("trees"))

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

    # Residual surface fitting and entropy
    df <- data.frame(x = coords[, 1], y = coords[, 2], z = z_vals)
    residual_sd <- NA_real_
    topo_entropy <- NA_real_
    smoothness_score <- NA_real_

    if (nrow(df) >= 5) {
      tryCatch({
        pca <- prcomp(df, center = TRUE, scale. = FALSE)
        pca_complexity <- sum((pca$sdev^2 / sum(pca$sdev^2))[2:3])

        # Spline model of Z ~ x + y
        fit <- lm(z ~ splines::bs(x, df = 4) + splines::bs(y, df = 4), data = df)
        residuals <- df$z - predict(fit, newdata = df)

        # Residual SD
        residual_sd <- sd(residuals)

        # Residual entropy
        h <- hist(residuals, breaks = 10, plot = FALSE)
        p <- h$counts / sum(h$counts)
        topo_entropy <- -sum(p * log(p + 1e-10))

        # Residual kurtosis
        if (!requireNamespace("e1071", quietly = TRUE)) stop("e1071 is required for kurtosis")
        residual_kurtosis <- e1071::kurtosis(residuals)

        # Normalize and compute final score
        scale_0_1 <- function(x, min_val, max_val) max(0, min(1, (x - min_val) / (max_val - min_val)))

        scaled_scores <- c(
          scale_0_1(residual_sd, 0, 20),
          scale_0_1(topo_entropy, 0, 3),
          scale_0_1(pca_complexity, 0, 0.5),
          scale_0_1(residual_kurtosis, 0, 10)
        )

        smoothness_score <- mean(scaled_scores)
      }, error = function(e) {
        residual_sd <- NA_real_
        topo_entropy <- NA_real_
        smoothness_score <- NA_real_
      })
    }

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
      n_gt_24_1 = n_gt_24_1,
      topo_residual_sd = residual_sd,
      topo_entropy = topo_entropy,
      smoothness_score = smoothness_score
    ))
  }, error = function(e) {
    return(named_zero_metrics('trees'))
  })

}


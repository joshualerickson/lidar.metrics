
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
#'
#' @return A named list of tree metrics:
#' \itemize{
#'   \item \code{n_trees}, \code{trees_per_acre}
#'   \item \code{n_low}, \code{n_low_mid}, \code{n_mid}, \code{n_mid_upper}, \code{n_upper}
#'   \item \code{topo_residual_sd}, \code{topo_entropy}
#' }
#' @export
tree_detection <- function(x, y, z,
                           window_func = function(x) {x*0.17 + 3}
                           ) {
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



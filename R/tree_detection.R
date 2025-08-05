tree_detection <- function(x, y, z,
                           window_func = function(x) {
                             w <- pmax(x, 1)
                             w <- log1p(w) * 1.5 + 2
                             w <- pmin(w, 6)
                             w[is.na(w) | w <= 0] <- 3
                             return(w)
                           }) {
  load_graph_deps()

  tryCatch({
    if (length(x) < 5 || length(z) < 5) {
      return(named_zero_metrics('trees'))
    }

    las_data <- suppressMessages(LAS(data.frame(X = x, Y = y, Z = z)))

    trees <- tryCatch({
      lidR::locate_trees(las_data, lidR::lmf(ws = window_func))
    }, error = function(e) {
      return(named_zero_metrics('trees'))
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

    # Residual surface fitting and entropy
    df <- data.frame(x = coords[, 1], y = coords[, 2], z = z_vals)
    residual_sd <- NA_real_
    topo_entropy <- NA_real_
    smoothness_score <- NA_real_

    if (nrow(df) >= 5) {
      tryCatch({
        fit <- lm(z ~ splines::bs(x, df = 4) + splines::bs(y, df = 4), data = df)
        residuals <- df$z - predict(fit, newdata = df)
        residual_sd <- sd(residuals)

        h <- hist(residuals, breaks = 10, plot = FALSE)
        p <- h$counts / sum(h$counts)
        topo_entropy <- -sum(p * log(p + 1e-10))

        # Composite smoothness score: scaled SD and entropy
        scale_0_1_scalar <- function(x, min_val, max_val) {
          x_scaled <- (x - min_val) / (max_val - min_val)
          return(max(0, min(1, x_scaled)))  # constrain to [0, 1]
        }

        sd_scaled      <- scale_0_1_scalar(residual_sd, 0, 20)
        entropy_scaled <- scale_0_1_scalar(topo_entropy, 0, 3)
        smoothness_score <- 0.5 * sd_scaled + 0.5 * entropy_scaled
      }, error = function(e) {
        residual_sd <- NA_real_
        topo_entropy <- NA_real_
        smoothness_score <- NA_real_
      })
    }

    n_trees <- as.numeric(nrow(trees))
    trees_per_acre <- as.numeric(n_trees / 0.222395)

    return(list(
      n_trees = n_trees,
      trees_per_acre = trees_per_acre,
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


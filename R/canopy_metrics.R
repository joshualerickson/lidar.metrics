
canopy_cover_metrics <- function(x, y, z, return_number,
                                 window_func, out_dir) {
  
  if (length(z) < 5) return(named_zero_metrics(type = 'canopy'))  # Skip if too few points
  
  trees <- tree_detection(x, y, z, window_func = window_func, out_dir = out_dir)
  
  
  # Canopy Cover
  
  first_above <- sum(z > 2 & return_number == 1)
  total_first <- sum(return_number == 1)
  canopy_cover <- ifelse(total_first == 0, NA_real_, first_above / total_first)
  cc <- list(fractional_canopy_cover = canopy_cover)
  
  # rumple <- tryCatch(
  #   suppressMessages(suppressWarnings(rumple_index(x, y, z))),
  #   error = function(e) NA_real_
  # )
  # 
  lad <- lad_metrics(z, dz = 1, z0 = 2)
  return(c(
    trees, 
    lad,
    cc
    #list(rumple = rumple))
  ))
  
}


tree_detection <- function(x,
                           y,
                           z,
                           window_func = function(x) {
                             w <- pmax(x, 1)
                             w <- log1p(w) * 1.5 + 2
                             w <- pmin(w, 6)
                             if (any(w <= 0 | is.na(w))) {
                               w[is.na(w) | w <= 0] <- 3  # fallback default
                             }
                             return(w)
                           },
                           out_dir) {
  load_graph_deps()
  tryCatch({
    if (length(x) < 5 || length(z) < 5) {
      return(list(
        n_trees = 0, trees_per_acre = 0,
        n_low = 0, n_low_mid = 0, n_mid = 0, n_mid_upper = 0, n_upper = 0,
        topo_residual_sd = 0, topo_entropy = 0
      ))
    }
    
    # Build point cloud
    las_data <- suppressMessages(LAS(data.frame(X = x, Y = y, Z = z)))
    
    
    # Locate trees
    trees <- tryCatch({
      lidR::locate_trees(las_data, lidR::lmf(ws = window_func))
    }, error = function(e) {
      message("locate_trees error: ", e$message)
      return(NULL)
    })
    
    if (is.null(trees) || nrow(trees) == 0) {
      return(list(
        n_trees = 0, trees_per_acre = 0,
        n_low = 0, n_low_mid = 0, n_mid = 0, n_mid_upper = 0, n_upper = 0,
        topo_residual_sd = 0, topo_entropy = 0
      ))
    }
    
    # Save coords as .fst
    coords <- sf::st_coordinates(trees)
    # out_path <- file.path(out_dir, paste0("X", round(mean(range(x)), 3), "_Y", round(mean(range(y)), 3), '.fst'))
    # fst::write_fst(as.data.frame(coords), out_path)
    # 
    # Tree metrics
    z_vals <- trees$Z
    zbin <- cut(z_vals,
                breaks = c(0.3, 6.1, 12.1, 24.1, 36.1, Inf),
                labels = c("low", "low_mid", "mid", "mid_upper", "upper"))
    
    strata <- c("low", "low_mid", "mid", "mid_upper", "upper")
    strata_table <- table(zbin)
    strata_counts <- setNames(rep(0, length(strata)), strata)
    strata_counts[names(strata_table)] <- as.integer(strata_table)
    
    # Topo smooth residuals
    df <- data.frame(x = coords[, 1], y = coords[, 2], z = z_vals)
    residual_sd <- NA_real_
    topo_entropy <- NA_real_
    
    if (nrow(df) >= 5) {
      tryCatch({
        fit <- lm(z ~ splines::bs(x, df = 4) + splines::bs(y, df = 4), data = df)
        z_pred <- predict(fit, newdata = df)
        residuals <- df$z - z_pred
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
    # ✅ Return named flat numeric list
    return(list(
      n_trees = as.numeric(n_trees),
      trees_per_acre = as.numeric(trees_per_acre),
      n_low = unname(strata_counts["low"]),
      n_low_mid = unname(strata_counts["low_mid"]),
      n_mid = unname(strata_counts["mid"]),
      n_mid_upper = unname(strata_counts["mid_upper"]),
      n_upper = unname(strata_counts["upper"]),
      topo_residual_sd = as.numeric(residual_sd),
      topo_entropy = as.numeric(topo_entropy)
    ))}, error = function(e) {
      message("tree_detection() general error: ", e$message)
      return(list(
        n_trees = NA_real_, trees_per_acre = NA_real_,
        n_low = NA_real_, n_low_mid = NA_real_, n_mid = NA_real_,
        n_mid_upper = NA_real_, n_upper = NA_real_,
        topo_residual_sd = NA_real_, topo_entropy = NA_real_
      ))
    })
  
}


lad_metrics <- function(z, dz = 1, k = 0.5, z0 = 2) {
  tryCatch({
    profile <- lidR::LAD(z, dz = dz, k = k, z0 = z0)
    
    if (nrow(profile) == 0) {
      return(list(
        LAI = NA_real_,
        LAD_max = NA_real_,
        LAD_mean = NA_real_,
        LAD_z_max = NA_real_
      ))
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
    return(list(
      LAI = NA_real_,
      LAD_max = NA_real_,
      LAD_mean = NA_real_,
      LAD_z_max = NA_real_
    ))
  })
}

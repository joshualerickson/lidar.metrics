#' Compute Canopy Cover
#'
#' Calculates canopy cover, and leaf area density metrics
#' for a small area of LiDAR point cloud data.
#'
#' @param x Numeric vector of X coordinates.
#' @param y Numeric vector of Y coordinates.
#' @param z Numeric vector of Z (height) values.
#' @param return_number Integer vector of return numbers (e.g., from LAS file).
#'
#' @return A named list with metrics including:
#' \itemize{
#'   \item \code{fractional_canopy_cover}
#'   \item \code{LAI}, \code{LAD_max}, \code{LAD_mean}, \code{LAD_z_max}
#' }
#' @export
canopy_cover_metrics <- function(x, y, z, return_number, scan_angle,
                                 voxel_res = 3,
                                 angle_scale, density_scale) {
  # Safety check for low point count
  if (length(z) < 5) return(named_zero_metrics(type = 'canopy'))

  # Create LAS
  las_all <- suppressMessages(lidR::LAS(data.frame(
    X = x,
    Y = y,
    Z = z,
    ReturnNumber = as.integer(ifelse(return_number > 7, 7L, return_number)),
    ScanAngle = scan_angle
  )))

  if (is.empty(las_all)) return(named_zero_metrics(type = 'canopy'))

  # Decimate
  las_all <- las_decimate_by_scan_and_density(
    las_all,
    voxel_res = voxel_res,
    angle_scale = angle_scale,
    density_scale = density_scale
  )

  # Check after decimation
  if (is.empty(las_all)) return(named_zero_metrics(type = 'canopy'))

  # Canopy cover
  first_above <- lidR::filter_poi(las_all, Z > 2 & ReturnNumber == 1)
  total_first <- lidR::filter_poi(las_all, ReturnNumber == 1)

  canopy_cover <- tryCatch({
    if (length(total_first@data$Z) == 0) NA_real_
    else length(first_above@data$Z) / length(total_first@data$Z)
  }, error = function(e) NA_real_)

  # LAD
  lad <- tryCatch(
    lad_metrics(las_all@data$Z, dz = 1, z0 = 2),
    error = function(e) named_zero_metrics(type = 'canopy')
  )

  # Final output as valid named list
  out <- c(
    lad,
    list(fractional_canopy_cover = canopy_cover)
  )

  return(out)
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
  default <- list(
    LAI = as.numeric(NA),
    LAD_max = as.numeric(NA),
    LAD_mean = as.numeric(NA),
    LAD_z_max = as.numeric(NA)
  )

  tryCatch({
    profile <- lidR::LAD(z, dz = dz, k = k, z0 = z0)

    if (nrow(profile) == 0) return(default)

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
    return(default)
  })
}

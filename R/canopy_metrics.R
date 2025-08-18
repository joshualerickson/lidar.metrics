#' Compute Canopy Cover
#'
#' Calculates canopy cover, and leaf area density metrics
#' for a small area of LiDAR point cloud data.
#'
#' @param x Numeric vector of X coordinates.
#' @param y Numeric vector of Y coordinates.
#' @param z Numeric vector of Z (height) values.
#' @param return_number Integer vector of return numbers (e.g., from LAS file).
#' @param zmax Numeric vector for VIC.
#'
#' @return A named list with metrics including:
#' \itemize{
#'   \item \code{fractional_canopy_cover}
#'   \item \code{LAI}, \code{LAD_max}, \code{LAD_mean}, \code{LAD_z_max}
#' }
#' @export
canopy_cover_metrics <- function(x, y, z, return_number, zmax) {
  # Safety check for low point count
  if (length(z) < 5) return(named_zero_metrics(type = 'canopy'))

  # Create LAS
  las_all <- suppressMessages(lidR::LAS(data.frame(
    X = x,
    Y = y,
    Z = z,
    ReturnNumber = as.integer(ifelse(return_number > 7, 7L, return_number))
  )))

  if (is.empty(las_all)) return(named_zero_metrics(type = 'canopy'))

  las_all <- lidR::filter_poi(las_all, ReturnNumber == 1)



  chm <- tryCatch({
    lidR::rasterize_canopy(las_all, res = 1, algorithm = lidR::pitfree(thresholds = c(0, 10, 20), max_edge =  c(0, 1.5)))
  }, error = function(e) NULL)

  canopy_cover <- tryCatch({
    vals <- values(chm)
    if (length(vals) == 0 || all(is.na(vals))) {
      NA_real_
    } else {
      sum(vals > 2, na.rm = TRUE) / sum(!is.na(vals))  # proportion of canopy pixels
    }
  }, error = function(e) NA_real_)

  rumple <- tryCatch({
    vals <- values(chm)
    if (length(vals) == 0 || all(is.na(vals))) {
      NA_real_
    } else {
      rumple_index(chm)  # proportion of canopy pixels
    }
  }, error = function(e) NA_real_)

  vci <- tryCatch({
      VCI(na.omit(las_all@data$Z), by = 1, zmax = zmax)  # proportion of canopy pixels

  }, error = function(e) NA_real_)


  # LAD
  lad <- tryCatch(
    lad_metrics(las_all@data$Z, dz = 1, z0 = 2),
    error = function(e) named_zero_metrics(type = 'canopy')
  )

  # Final output as valid named list
  out <- c(
    lad,
    list(fractional_canopy_cover = canopy_cover,
         rumple_index = rumple,
         vertical_complexity_index = vci)
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

#' Compute Canopy Cover
#'
#' Calculates canopy cover, and leaf area density metrics
#' for a small area of LiDAR point cloud data.
#'
#' @param x Numeric vector of X coordinates.
#' @param y Numeric vector of Y coordinates.
#' @param z Numeric vector of Z (height) values.
#' @param psid unquoted `PointSourceID` from the LAS file.
#' @param QL1 logical for whether to QL1 methods or not. Default (FALSE)
#' @param density Numeric to pass to `lidR::homogenize()` density argument.
#' @param dec_res Numeric to pass to `lidR::homogenize()` res argument.
#' @param ... Arguments to pass to `lidar.metrics::tree_detection()`.
#'
#' @return A named list with metrics including:
#' \itemize{
#'   \item \code{fractional_canopy_cover}
#'   \item \code{rumple_index}
#'   \item \code{LAI}, \code{LAD_max}, \code{LAD_mean}, \code{LAD_z_max}
#' }
#' @export
canopy_cover_metrics <- function(x, y, z,
                                 psid,
                                 QL1 = FALSE,
                                 n = 15,
                                 res = 3, ...) {
  # Safety check for low point count
  if (length(z) < 5) return(c(named_zero_metrics('trees'), named_zero_metrics('canopy')))

  # Create LAS
  las <- suppressMessages(lidR::LAS(data.frame(
    X = x,
    Y = y,
    Z = z,
    PointSourceID = psid
  )))

  if (is.empty(las)) return(c(named_zero_metrics('trees'), named_zero_metrics('canopy')))

  if(QL1) {

    las_filtered <- lidR::decimate_points(las, algorithm = lidR::homogenize(density = n, res = res))


  } else {

    psid_vec <- las@data$PointSourceID

    keep_ids <- as.integer(names(which.max(table(psid_vec))))

    las_filtered <- lidR::filter_poi(las, PointSourceID %in% keep_ids)

  }

  chm <- tryCatch({
    lidR::rasterize_canopy(las_filtered, res = 0.5, algorithm = lidR::pitfree(thresholds = c(0, 10, 20), max_edge =  c(0, 1.5)))
  }, error = function(e) NA_real_)


    canopy_cover <- tryCatch({
    vals <- terra::values(chm)

      if (length(vals) == 0L) NA_real_ else {
        denom <- sum(!is.na(vals))
        if (denom == 0L) NA_real_ else {
          num <- sum(vals > 2, na.rm = TRUE)
          as.numeric(num / denom)
        }
      }
    }, error = function(e) NA_real_)

    # Replace NaN/Inf defensively
    if (!is.finite(canopy_cover)) canopy_cover <- NA_real_

  rumple <- tryCatch({
    vals <- terra::values(chm)
    if (length(vals) == 0L) NA_real_ else {
      denom <- sum(!is.na(vals))
      if (denom == 0L) NA_real_ else {
       lidR::rumple_index(chm)
      }
    }
  }, error = function(e) NA_real_)


  if (!is.finite(rumple)) rumple <- NA_real_


  lad <- tryCatch(
    lad_metrics(las_filtered@data$Z, dz = 1, z0 = 2),
    error = function(e) named_zero_metrics(type = 'canopy')
  )

  trees <- tryCatch(tree_detection(las_filtered, ...),
                    error = function(e) named_zero_metrics(type = 'trees'))

  # Final output as valid named list
  out <- c(
    lad,
    list(fractional_canopy_cover = canopy_cover,
    rumple_index = rumple),
    trees,
    n_m2 = length(las_filtered@data$Z)/900
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



#' Calculate Summary LiDAR Metrics for Vegetation Structure
#'
#' This function computes a set of descriptive statistics (e.g., mean, standard deviation,
#' quantiles, skewness, kurtosis, etc.) from airborne LiDAR point cloud data. It calculates
#' metrics for all returns > 2 m and for first returns > 2 m, using both elevation (`z`)
#' and intensity (`i`) values. It is designed to be memory-efficient and avoids repeated
#' subsetting or reliance on external packages like `moments`.
#'
#' @param z Numeric vector. Elevation (`Z`) values of LiDAR returns (typically in meters).
#' @param i Numeric vector. `Intensity` values of LiDAR returns.
#' @param r Integer or logical vector. Return numbers (`ReturnNumber`), where 1 typically indicates a first return.
#'
#' @return A named list of LiDAR metrics including:
#' \itemize{
#'   \item Metrics for all returns above 2 m (e.g., mean height, P95, CV, skewness, etc.)
#'   \item Metrics for first returns above 2 m
#'   \item Intensity-based statistics for both groups
#'   \item Canopy relief ratio, interquartile range, and higher-order moments
#'   \item Percentages of returns above 2 m and above the mean
#' }
#'
#' @details
#' Metrics are grouped and prefixed as follows:
#' \describe{
#'   \item{ALL_RETURNS_}{All returns with elevation > 2 m}
#'   \item{FIRST_RETURNS_}{First returns with elevation > 2 m}
#' }
#'
#' This function avoids dependency on packages like `moments` for skewness and kurtosis by using
#' internal helper functions. It also avoids repeated subsetting by computing indices once.
#'
#' @note This function is a memory-efficient rewrite of older versions of `Metrics()` functions commonly used
#' in LiDAR processing pipelines. It is well-suited for use in `*_metrics()` or similar functions
#' in the `lidR` package.
#'
#' @examples
#' \dontrun{
#' metrics <- Metrics(z = las$Z, i = las$Intensity, r = las$ReturnNumber)
#' }
#'
#' @export
Metrics = function(z, i, r) {
  # Pre-filter data once to avoid repeated subsetting
  twoPlus_idx = which(z > 2)
  first_idx = which(z > 2 & r == 1)

  # Extract subsets once
  z_twoPlus = z[twoPlus_idx]
  i_twoPlus = i[twoPlus_idx]
  z_first = z[first_idx]
  i_first = i[first_idx]

  # Helper functions to replace moments package
  safe_skewness = function(x) {
    if(length(x) < 3) return(NA)
    n = length(x)
    m = mean(x, na.rm = TRUE)
    s = sd(x, na.rm = TRUE)
    if(s == 0 || is.na(s)) return(0)
    skew = sum(((x - m)/s)^3, na.rm = TRUE) / n
    return(skew)
  }

  safe_kurtosis = function(x) {
    if(length(x) < 4) return(NA)
    n = length(x)
    m = mean(x, na.rm = TRUE)
    s = sd(x, na.rm = TRUE)
    if(s == 0 || is.na(s)) return(0)
    kurt = sum(((x - m)/s)^4, na.rm = TRUE) / n - 3
    return(kurt)
  }

  # Simple mode function (most frequent value)
  safe_mode = function(x) {
    if(length(x) == 0) return(NA)
    ux = unique(x)
    if(length(ux) == 1) return(ux[1])
    tab = tabulate(match(x, ux))
    return(as.numeric(ux[which.max(tab)]))
  }

  # Pre-calculate common values to avoid repeated computation
  z_twoPlus_mean = mean(z_twoPlus)
  i_twoPlus_mean = mean(i_twoPlus)
  z_first_mean = mean(z_first)
  i_first_mean = mean(i_first)

  metrics = list(
    # ALL RETURNS metrics
    ALL_RETURNS_all_cnt_2plus = length(z_twoPlus),
    ALL_RETURNS_elev_AAD_2plus = mean(abs(z_twoPlus - z_twoPlus_mean)),
    ALL_RETURNS_elev_ave_2plus = z_twoPlus_mean,
    ALL_RETURNS_elev_canopy_relief_ratio = {
      z_min = min(z_twoPlus)
      z_max = max(z_twoPlus)
      (z_twoPlus_mean - z_min) / (z_max - z_min)
    },
    ALL_RETURNS_elev_CV_2plus = sd(z_twoPlus) / z_twoPlus_mean,
    ALL_RETURNS_elev_IQ_2plus = IQR(z_twoPlus),
    ALL_RETURNS_elev_kurtosis_2plus = safe_kurtosis(z_twoPlus),
    ALL_RETURNS_elev_max_2plus = max(z_twoPlus),
    ALL_RETURNS_elev_P01_2plus = quantile(z_twoPlus, 0.01, na.rm = TRUE),
    ALL_RETURNS_elev_P05_2plus = quantile(z_twoPlus, 0.05, na.rm = TRUE),
    ALL_RETURNS_elev_P25_2plus = quantile(z_twoPlus, 0.25, na.rm = TRUE),
    ALL_RETURNS_elev_P50_2plus = quantile(z_twoPlus, 0.50, na.rm = TRUE),
    ALL_RETURNS_elev_P75_2plus = quantile(z_twoPlus, 0.75, na.rm = TRUE),
    ALL_RETURNS_elev_P95_2plus = quantile(z_twoPlus, 0.95, na.rm = TRUE),
    ALL_RETURNS_elev_P99_2plus = quantile(z_twoPlus, 0.99, na.rm = TRUE),
    ALL_RETURNS_elev_skewness_2plus = safe_skewness(z_twoPlus),
    ALL_RETURNS_elev_stddev_2plus = sd(z_twoPlus),
    ALL_RETURNS_elev_variance_2plus = var(z_twoPlus),

    # ALL RETURNS intensity metrics
    ALL_RETURNS_int_AAD_2plus = mean(abs(i_twoPlus - i_twoPlus_mean)),
    ALL_RETURNS_int_ave_2plus = i_twoPlus_mean,
    ALL_RETURNS_int_CV_2plus = sd(i_twoPlus) / i_twoPlus_mean,
    ALL_RETURNS_int_IQ_2plus = IQR(i_twoPlus),
    ALL_RETURNS_int_kurtosis_2plus = safe_kurtosis(i_twoPlus),
    ALL_RETURNS_int_max_2plus = as.numeric(max(i_twoPlus)),
    ALL_RETURNS_int_min_2plus = as.numeric(min(i_twoPlus)),
    ALL_RETURNS_int_P01_2plus = quantile(i_twoPlus, 0.01, na.rm = TRUE),
    ALL_RETURNS_int_P05_2plus = quantile(i_twoPlus, 0.05, na.rm = TRUE),
    ALL_RETURNS_int_P25_2plus = quantile(i_twoPlus, 0.25, na.rm = TRUE),
    ALL_RETURNS_int_P50_2plus = quantile(i_twoPlus, 0.50, na.rm = TRUE),
    ALL_RETURNS_int_P75_2plus = quantile(i_twoPlus, 0.75, na.rm = TRUE),
    ALL_RETURNS_int_P95_2plus = quantile(i_twoPlus, 0.95, na.rm = TRUE),
    ALL_RETURNS_int_P99_2plus = quantile(i_twoPlus, 0.99, na.rm = TRUE),
    ALL_RETURNS_int_skewness_2plus = safe_skewness(i_twoPlus),
    ALL_RETURNS_int_stddev_2plus = sd(i_twoPlus),
    ALL_RETURNS_int_variance_2plus = var(i_twoPlus),

    # FIRST RETURNS metrics
    FIRST_RETURNS_all_cnt_2plus = length(z_first),
    FIRST_RETURNS_elev_AAD_2plus = mean(abs(z_first - z_first_mean)),
    FIRST_RETURNS_elev_ave_2plus = z_first_mean,
    FIRST_RETURNS_elev_canopy_relief_ratio = {
      z_first_min = min(z_first)
      z_first_max = max(z_first)
      (z_first_mean - z_first_min) / (z_first_max - z_first_min)
    },
    FIRST_RETURNS_elev_CV_2plus = sd(z_first) / z_first_mean,
    FIRST_RETURNS_elev_IQ_2plus = IQR(z_first),
    FIRST_RETURNS_elev_kurtosis_2plus = safe_kurtosis(z_first),
    FIRST_RETURNS_elev_max_2plus = max(z_first),
    FIRST_RETURNS_elev_P01_2plus = quantile(z_first, 0.01, na.rm = TRUE),
    FIRST_RETURNS_elev_P05_2plus = quantile(z_first, 0.05, na.rm = TRUE),
    FIRST_RETURNS_elev_P25_2plus = quantile(z_first, 0.25, na.rm = TRUE),
    FIRST_RETURNS_elev_P50_2plus = quantile(z_first, 0.50, na.rm = TRUE),
    FIRST_RETURNS_elev_P75_2plus = quantile(z_first, 0.75, na.rm = TRUE),
    FIRST_RETURNS_elev_P95_2plus = quantile(z_first, 0.95, na.rm = TRUE),
    FIRST_RETURNS_elev_P99_2plus = quantile(z_first, 0.99, na.rm = TRUE),
    FIRST_RETURNS_elev_skewness_2plus = safe_skewness(z_first),
    FIRST_RETURNS_elev_stddev_2plus = sd(z_first),
    FIRST_RETURNS_elev_variance_2plus = var(z_first),

    # FIRST RETURNS intensity metrics
    FIRST_RETURNS_int_AAD_2plus = mean(abs(i_first - i_first_mean)),
    FIRST_RETURNS_int_ave_2plus = i_first_mean,
    FIRST_RETURNS_int_CV_2plus = sd(i_first) / i_first_mean,
    FIRST_RETURNS_int_IQ_2plus = IQR(i_first),
    FIRST_RETURNS_int_kurtosis_2plus = safe_kurtosis(i_first),
    FIRST_RETURNS_int_max_2plus = as.numeric(max(i_first)),
    FIRST_RETURNS_int_min_2plus = as.numeric(min(i_first)),
    FIRST_RETURNS_int_mode_2plus = safe_mode(i_first),
    FIRST_RETURNS_int_P01_2plus = quantile(i_first, 0.01, na.rm = TRUE),
    FIRST_RETURNS_int_P05_2plus = quantile(i_first, 0.05, na.rm = TRUE),
    FIRST_RETURNS_int_P25_2plus = quantile(i_first, 0.25, na.rm = TRUE),
    FIRST_RETURNS_int_P50_2plus = quantile(i_first, 0.50, na.rm = TRUE),
    FIRST_RETURNS_int_P75_2plus = quantile(i_first, 0.75, na.rm = TRUE),
    FIRST_RETURNS_int_P95_2plus = quantile(i_first, 0.95, na.rm = TRUE),
    FIRST_RETURNS_int_P99_2plus = quantile(i_first, 0.99, na.rm = TRUE),
    FIRST_RETURNS_int_skewness_2plus = safe_skewness(i_first),
    FIRST_RETURNS_int_stddev_2plus = sd(i_first),
    FIRST_RETURNS_int_variance_2plus = var(i_first),

    # Percentage metrics
    Pct1stRtns_above_2 = length(z_first) / length(z[r == 1]) * 100,
    Pct1stRtns_above_mean = length(z[r == 1 & z > z_first_mean]) / length(z[r == 1]) * 100,
    PctAllRtns_above_2 = length(z_twoPlus) / length(z) * 100,
    PctAllRtns_above_mean = length(z[z > z_twoPlus_mean]) / length(z) * 100
  )

  return(metrics)
}


#' Calculate Stratum-Based LiDAR Metrics
#'
#' Computes vertical distribution metrics for LiDAR return heights, stratified by predefined
#' elevation bins (e.g., 0.5-6.1, ..., 64+ m). Returns statistics for both all returns
#' and first returns in each vertical bin. This is useful for assessing canopy structure, layering,
#' and vertical complexity in forest stands.
#'
#' @param z Numeric vector. Elevation (`Z`) values of LiDAR returns (meters above ground or DEM).
#' @param r Integer or logical vector. Return number (`ReturnNumber`) of each point (1 typically indicates first return).
#'
#' @return A named list of stratum metrics. For each height bin (stratum), the following statistics
#' are returned for both all and first returns:
#' \itemize{
#'   \item Coefficient of variation (CV)
#'   \item Kurtosis
#'   \item Skewness
#'   \item Standard deviation
#' }
#'
#' @details
#' The vertical strata are defined as follows (in meters):
#' \preformatted{
#'   0.5-6.1, 6.1–16, 16–32, 32–48, 48–64, 64+
#' }
#' Prefixes:
#' \describe{
#'   \item{ALL_RETURNS_strata_}{All returns in each stratum}
#'   \item{FIRST_RETURNS_strata_}{First returns in each stratum}
#' }
#'
#' Skewness and kurtosis are calculated internally to avoid reliance on external packages.
#' This function assumes that `z` is already normalized to height above ground.
#'
#' @note This is a memory-optimized alternative to more verbose stratum metric pipelines,
#' particularly suitable for processing in tiling frameworks such as `catalog_apply()` and `*_metrics()` in `lidR`.
#'
#' @examples
#' \dontrun{
#' strata <- StrataMetrics(z = las$Z, r = las$ReturnNumber)
#' }
#'
#' @export
StrataMetrics = function(z, r) {
  # Helper functions
  safe_skewness = function(x) {
    if(length(x) < 3) return(NA)
    n = length(x)
    m = mean(x, na.rm = TRUE)
    s = sd(x, na.rm = TRUE)
    if(s == 0 || is.na(s)) return(0)
    sum(((x - m)/s)^3, na.rm = TRUE) / n
  }

  safe_kurtosis = function(x) {
    if(length(x) < 4) return(NA)
    n = length(x)
    m = mean(x, na.rm = TRUE)
    s = sd(x, na.rm = TRUE)
    if(s == 0 || is.na(s)) return(0)
    sum(((x - m)/s)^4, na.rm = TRUE) / n - 3
  }

  safe_cv = function(x) {
    if(length(x) == 0) return(NA)
    m = mean(x, na.rm = TRUE)
    if(m == 0 || is.na(m)) return(NA)
    sd(x, na.rm = TRUE) / m
  }

  # Define strata breakpoints (new bins)
  bins = list(
    s0 = which(z > 0.5 & z < 6.1),
    s1 = which(z >= 6.1 & z < 16),
    s2 = which(z >= 16 & z < 32),
    s3 = which(z >= 32 & z < 48),
    s4 = which(z >= 48 & z < 64),
    s5 = which(z >= 64)
  )

  bins_first = list(
    s0 = which(z > 0.5 & z < 6.1 & r == 1),
    s1 = which(z >= 6.1 & z < 16 & r == 1),
    s2 = which(z >= 16 & z < 32 & r == 1),
    s3 = which(z >= 32 & z < 48 & r == 1),
    s4 = which(z >= 48 & z < 64 & r == 1),
    s5 = which(z >= 64 & r == 1)
  )

  # Bin names
  bin_labels = c("0p5to6p1M", "6p1to16M", "16to32M", "32to48M", "48to64M", "64M_plus")

  metrics = list()

  # Loop through all strata
  for (i in seq_along(bins)) {
    bin_data = z[bins[[i]]]
    bin_first = z[bins_first[[i]]]
    label = bin_labels[i]

    metrics[[paste0("ALL_RETURNS_strata_", label, "_CV")]] = safe_cv(bin_data)
    metrics[[paste0("ALL_RETURNS_strata_", label, "_kurtosis")]] = safe_kurtosis(bin_data)
    metrics[[paste0("ALL_RETURNS_strata_", label, "_skewness")]] = safe_skewness(bin_data)
    metrics[[paste0("ALL_RETURNS_strata_", label, "_stddev")]] = sd(bin_data, na.rm = TRUE)

    metrics[[paste0("FIRST_RETURNS_strata_", label, "_CV")]] = safe_cv(bin_first)
    metrics[[paste0("FIRST_RETURNS_strata_", label, "_kurtosis")]] = safe_kurtosis(bin_first)
    metrics[[paste0("FIRST_RETURNS_strata_", label, "_skewness")]] = safe_skewness(bin_first)
    metrics[[paste0("FIRST_RETURNS_strata_", label, "_stddev")]] = sd(bin_first, na.rm = TRUE)
  }

  # Clean return (ensure all numeric)
  metrics <- lapply(metrics, function(x) {
    if(is.logical(x)) return(as.numeric(x))
    if(is.null(x) || length(x) == 0) return(NA_real_)
    if(!is.numeric(x)) return(as.numeric(x))
    if(length(x) > 1) return(x[1])
    return(as.numeric(x))
  })

  return(metrics)
}

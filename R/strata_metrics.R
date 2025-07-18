

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
#' elevation bins (e.g., 0–0.5 m, 0.5–1 m, ..., 64+ m). Returns statistics for both all returns
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
#'   \item Total number of returns
#'   \item Proportion of returns in that bin
#' }
#'
#' @details
#' The vertical strata are defined as follows (in meters):
#' \preformatted{
#'   0–0.5, 0.5–1, 1–2, 2–4, 4–8, 8–16, 16–32, 32–48, 48–64, 64+
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

  # Helper functions to replace moments package (same as before)
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

  # Safe CV calculation
  safe_cv = function(x) {
    if(length(x) == 0) return(NA)
    m = mean(x, na.rm = TRUE)
    if(m == 0 || is.na(m)) return(NA)
    return(sd(x, na.rm = TRUE) / m)
  }

  # Pre-calculate all strata indices once
  gt0_idx = which(z > 0)
  s0_idx = which(z > 0 & z < 0.5)
  s1_idx = which(z > 0.5 & z < 1)
  s2_idx = which(z > 1 & z < 2)
  s3_idx = which(z > 2 & z < 4)
  s4_idx = which(z > 4 & z < 8)
  s5_idx = which(z > 8 & z < 16)
  s6_idx = which(z > 16 & z < 32)
  s7_idx = which(z > 32 & z < 48)
  s8_idx = which(z > 48 & z < 64)
  s9_idx = which(z > 64)

  # First return strata indices
  gt0first_idx = which(z > 0 & r == 1)
  s0first_idx = which(z > 0 & z < 0.5 & r == 1)
  s1first_idx = which(z > 0.5 & z < 1 & r == 1)
  s2first_idx = which(z > 1 & z < 2 & r == 1)
  s3first_idx = which(z > 2 & z < 4 & r == 1)
  s4first_idx = which(z > 4 & z < 8 & r == 1)
  s5first_idx = which(z > 8 & z < 16 & r == 1)
  s6first_idx = which(z > 16 & z < 32 & r == 1)
  s7first_idx = which(z > 32 & z < 48 & r == 1)
  s8first_idx = which(z > 48 & z < 64 & r == 1)
  s9first_idx = which(z > 64 & r == 1)

  # Extract subsets once
  z_gt0 = z[gt0_idx]
  z_s0 = z[s0_idx]
  z_s1 = z[s1_idx]
  z_s2 = z[s2_idx]
  z_s3 = z[s3_idx]
  z_s4 = z[s4_idx]
  z_s5 = z[s5_idx]
  z_s6 = z[s6_idx]
  z_s7 = z[s7_idx]
  z_s8 = z[s8_idx]
  z_s9 = z[s9_idx]

  z_gt0first = z[gt0first_idx]
  z_s0first = z[s0first_idx]
  z_s1first = z[s1first_idx]
  z_s2first = z[s2first_idx]
  z_s3first = z[s3first_idx]
  z_s4first = z[s4first_idx]
  z_s5first = z[s5first_idx]
  z_s6first = z[s6first_idx]
  z_s7first = z[s7first_idx]
  z_s8first = z[s8first_idx]
  z_s9first = z[s9first_idx]

  # Pre-calculate counts for proportions
  n_gt0 = length(z_gt0)
  n_gt0first = length(z_gt0first)

  metrics = list(
    # ALL RETURNS - Strata 0.5 to 1M
    ALL_RETURNS_strata_0p5to1M_CV = safe_cv(z_s1),
    ALL_RETURNS_strata_0p5to1M_kurtosis = safe_kurtosis(z_s1),
    ALL_RETURNS_strata_0p5to1M_return_proportion = length(z_s1) / n_gt0,
    ALL_RETURNS_strata_0p5to1M_skewness = safe_skewness(z_s1),
    ALL_RETURNS_strata_0p5to1M_stddev = sd(z_s1, na.rm = TRUE),
    ALL_RETURNS_strata_0p5to1M_total_return_cnt = length(z_s1),

    # ALL RETURNS - Strata 0 to 0.5M
    ALL_RETURNS_strata_0to0p5M_CV = safe_cv(z_s0),
    ALL_RETURNS_strata_0to0p5M_kurtosis = safe_kurtosis(z_s0),
    ALL_RETURNS_strata_0to0p5M_return_proportion = length(z_s0) / n_gt0,
    ALL_RETURNS_strata_0to0p5M_skewness = safe_skewness(z_s0),
    ALL_RETURNS_strata_0to0p5M_stddev = sd(z_s0, na.rm = TRUE),
    ALL_RETURNS_strata_0to0p5M_total_return_cnt = length(z_s0),

    # ALL RETURNS - Strata 16 to 32M
    ALL_RETURNS_strata_16to32M_CV = safe_cv(z_s6),
    ALL_RETURNS_strata_16to32M_kurtosis = safe_kurtosis(z_s6),
    ALL_RETURNS_strata_16to32M_return_proportion = length(z_s6) / n_gt0,
    ALL_RETURNS_strata_16to32M_skewness = safe_skewness(z_s6),
    ALL_RETURNS_strata_16to32M_stddev = sd(z_s6, na.rm = TRUE),
    ALL_RETURNS_strata_16to32M_total_return_cnt = length(z_s6),

    # ALL RETURNS - Strata 1 to 2M
    ALL_RETURNS_strata_1to2M_CV = safe_cv(z_s2),
    ALL_RETURNS_strata_1to2M_kurtosis = safe_kurtosis(z_s2),
    ALL_RETURNS_strata_1to2M_return_proportion = length(z_s2) / n_gt0,
    ALL_RETURNS_strata_1to2M_skewness = safe_skewness(z_s2),
    ALL_RETURNS_strata_1to2M_stddev = sd(z_s2, na.rm = TRUE),
    ALL_RETURNS_strata_1to2M_total_return_cnt = length(z_s2),

    # ALL RETURNS - Strata 2 to 4M
    ALL_RETURNS_strata_2to4M_CV = safe_cv(z_s3),
    ALL_RETURNS_strata_2to4M_kurtosis = safe_kurtosis(z_s3),
    ALL_RETURNS_strata_2to4M_return_proportion = length(z_s3) / n_gt0,
    ALL_RETURNS_strata_2to4M_skewness = safe_skewness(z_s3),
    ALL_RETURNS_strata_2to4M_stddev = sd(z_s3, na.rm = TRUE),
    ALL_RETURNS_strata_2to4M_total_return_cnt = length(z_s3),

    # ALL RETURNS - Strata 32 to 48M
    ALL_RETURNS_strata_32to48M_CV = safe_cv(z_s7),
    ALL_RETURNS_strata_32to48M_kurtosis = safe_kurtosis(z_s7),
    ALL_RETURNS_strata_32to48M_return_proportion = length(z_s7) / n_gt0,
    ALL_RETURNS_strata_32to48M_skewness = safe_skewness(z_s7),
    ALL_RETURNS_strata_32to48M_stddev = sd(z_s7, na.rm = TRUE),
    ALL_RETURNS_strata_32to48M_total_return_cnt = length(z_s7),

    # ALL RETURNS - Strata 48 to 64M
    ALL_RETURNS_strata_48to64M_CV = safe_cv(z_s8),
    ALL_RETURNS_strata_48to64M_kurtosis = safe_kurtosis(z_s8),
    ALL_RETURNS_strata_48to64M_return_proportion = length(z_s8) / n_gt0,
    ALL_RETURNS_strata_48to64M_skewness = safe_skewness(z_s8),
    ALL_RETURNS_strata_48to64M_stddev = sd(z_s8, na.rm = TRUE),
    ALL_RETURNS_strata_48to64M_total_return_cnt = length(z_s8),

    # ALL RETURNS - Strata 4 to 8M
    ALL_RETURNS_strata_4to8M_CV = safe_cv(z_s4),
    ALL_RETURNS_strata_4to8M_kurtosis = safe_kurtosis(z_s4),
    ALL_RETURNS_strata_4to8M_return_proportion = length(z_s4) / n_gt0,
    ALL_RETURNS_strata_4to8M_skewness = safe_skewness(z_s4),
    ALL_RETURNS_strata_4to8M_stddev = sd(z_s4, na.rm = TRUE),
    ALL_RETURNS_strata_4to8M_total_return_cnt = length(z_s4),

    # ALL RETURNS - Strata 64M plus
    ALL_RETURNS_strata_64M_plus_CV = safe_cv(z_s9),
    ALL_RETURNS_strata_64M_plus_kurtosis = safe_kurtosis(z_s9),
    ALL_RETURNS_strata_64M_plus_return_proportion = length(z_s9) / n_gt0,
    ALL_RETURNS_strata_64M_plus_skewness = safe_skewness(z_s9),
    ALL_RETURNS_strata_64M_plus_stddev = sd(z_s9, na.rm = TRUE),
    ALL_RETURNS_strata_64M_plus_total_return_cnt = length(z_s9),

    # ALL RETURNS - Strata 8 to 16M
    ALL_RETURNS_strata_8to16M_CV = safe_cv(z_s5),
    ALL_RETURNS_strata_8to16M_kurtosis = safe_kurtosis(z_s5),
    ALL_RETURNS_strata_8to16M_return_proportion = length(z_s5) / n_gt0,
    ALL_RETURNS_strata_8to16M_skewness = safe_skewness(z_s5),
    ALL_RETURNS_strata_8to16M_stddev = sd(z_s5, na.rm = TRUE),
    ALL_RETURNS_strata_8to16M_total_return_cnt = length(z_s5),

    # FIRST RETURNS - Strata 0.5 to 1M
    FIRST_RETURNS_strata_0p5to1M_CV = safe_cv(z_s1first),
    FIRST_RETURNS_strata_0p5to1M_kurtosis = safe_kurtosis(z_s1first),
    FIRST_RETURNS_strata_0p5to1M_return_proportion = length(z_s1first) / n_gt0first,
    FIRST_RETURNS_strata_0p5to1M_skewness = safe_skewness(z_s1first),
    FIRST_RETURNS_strata_0p5to1M_stddev = sd(z_s1first, na.rm = TRUE),
    FIRST_RETURNS_strata_0p5to1M_total_return_cnt = length(z_s1first),

    # FIRST RETURNS - Strata 0 to 0.5M
    FIRST_RETURNS_strata_0to0p5M_CV = safe_cv(z_s0first),
    FIRST_RETURNS_strata_0to0p5M_kurtosis = safe_kurtosis(z_s0first),
    FIRST_RETURNS_strata_0to0p5M_return_proportion = length(z_s0first) / n_gt0first,
    FIRST_RETURNS_strata_0to0p5M_skewness = safe_skewness(z_s0first),
    FIRST_RETURNS_strata_0to0p5M_stddev = sd(z_s0first, na.rm = TRUE),
    FIRST_RETURNS_strata_0to0p5M_total_return_cnt = length(z_s0first),

    # FIRST RETURNS - Strata 16 to 32M
    FIRST_RETURNS_strata_16to32M_CV = safe_cv(z_s6first),
    FIRST_RETURNS_strata_16to32M_kurtosis = safe_kurtosis(z_s6first),
    FIRST_RETURNS_strata_16to32M_return_proportion = length(z_s6first) / n_gt0first,
    FIRST_RETURNS_strata_16to32M_skewness = safe_skewness(z_s6first),
    FIRST_RETURNS_strata_16to32M_stddev = sd(z_s6first, na.rm = TRUE),
    FIRST_RETURNS_strata_16to32M_total_return_cnt = length(z_s6first),

    # FIRST RETURNS - Strata 1 to 2M
    FIRST_RETURNS_strata_1to2M_CV = safe_cv(z_s2first),
    FIRST_RETURNS_strata_1to2M_kurtosis = safe_kurtosis(z_s2first),
    FIRST_RETURNS_strata_1to2M_return_proportion = length(z_s2first) / n_gt0first,
    FIRST_RETURNS_strata_1to2M_skewness = safe_skewness(z_s2first),
    FIRST_RETURNS_strata_1to2M_stddev = sd(z_s2first, na.rm = TRUE),
    FIRST_RETURNS_strata_1to2M_total_return_cnt = length(z_s2first),

    # FIRST RETURNS - Strata 2 to 4M
    FIRST_RETURNS_strata_2to4M_CV = safe_cv(z_s3first),
    FIRST_RETURNS_strata_2to4M_kurtosis = safe_kurtosis(z_s3first),
    FIRST_RETURNS_strata_2to4M_return_proportion = length(z_s3first) / n_gt0first,
    FIRST_RETURNS_strata_2to4M_skewness = safe_skewness(z_s3first),
    FIRST_RETURNS_strata_2to4M_stddev = sd(z_s3first, na.rm = TRUE),
    FIRST_RETURNS_strata_2to4M_total_return_cnt = length(z_s3first),

    # FIRST RETURNS - Strata 32 to 48M
    FIRST_RETURNS_strata_32to48M_CV = safe_cv(z_s7first),
    FIRST_RETURNS_strata_32to48M_kurtosis = safe_kurtosis(z_s7first),
    FIRST_RETURNS_strata_32to48M_return_proportion = length(z_s7first) / n_gt0first,
    FIRST_RETURNS_strata_32to48M_skewness = safe_skewness(z_s7first),
    FIRST_RETURNS_strata_32to48M_stddev = sd(z_s7first, na.rm = TRUE),
    FIRST_RETURNS_strata_32to48M_total_return_cnt = length(z_s7first),

    # FIRST RETURNS - Strata 48 to 64M
    FIRST_RETURNS_strata_48to64M_CV = safe_cv(z_s8first),
    FIRST_RETURNS_strata_48to64M_kurtosis = safe_kurtosis(z_s8first),
    FIRST_RETURNS_strata_48to64M_return_proportion = length(z_s8first) / n_gt0first,
    FIRST_RETURNS_strata_48to64M_skewness = safe_skewness(z_s8first),
    FIRST_RETURNS_strata_48to64M_stddev = sd(z_s8first, na.rm = TRUE),
    FIRST_RETURNS_strata_48to64M_total_return_cnt = length(z_s8first),

    # FIRST RETURNS - Strata 4 to 8M
    FIRST_RETURNS_strata_4to8M_CV = safe_cv(z_s4first),
    FIRST_RETURNS_strata_4to8M_kurtosis = safe_kurtosis(z_s4first),
    FIRST_RETURNS_strata_4to8M_return_proportion = length(z_s4first) / n_gt0first,
    FIRST_RETURNS_strata_4to8M_skewness = safe_skewness(z_s4first),
    FIRST_RETURNS_strata_4to8M_stddev = sd(z_s4first, na.rm = TRUE),
    FIRST_RETURNS_strata_4to8M_total_return_cnt = length(z_s4first),

    # FIRST RETURNS - Strata 64M plus
    FIRST_RETURNS_strata_64M_plus_CV = safe_cv(z_s9first),
    FIRST_RETURNS_strata_64M_plus_kurtosis = safe_kurtosis(z_s9first),
    FIRST_RETURNS_strata_64M_plus_return_proportion = length(z_s9first) / n_gt0first,
    FIRST_RETURNS_strata_64M_plus_skewness = safe_skewness(z_s9first),
    FIRST_RETURNS_strata_64M_plus_stddev = sd(z_s9first, na.rm = TRUE),
    FIRST_RETURNS_strata_64M_plus_total_return_cnt = length(z_s9first),

    # FIRST RETURNS - Strata 8 to 16M
    FIRST_RETURNS_strata_8to16M_CV = safe_cv(z_s5first),
    FIRST_RETURNS_strata_8to16M_kurtosis = safe_kurtosis(z_s5first),
    FIRST_RETURNS_strata_8to16M_return_proportion = length(z_s5first) / n_gt0first,
    FIRST_RETURNS_strata_8to16M_skewness = safe_skewness(z_s5first),
    FIRST_RETURNS_strata_8to16M_stddev = sd(z_s5first, na.rm = TRUE),
    FIRST_RETURNS_strata_8to16M_total_return_cnt = length(z_s5first)
  )

  # Force all values to be numeric
  metrics <- lapply(metrics, function(x) {
    if(is.logical(x)) {
      return(as.numeric(x))  # Convert TRUE/FALSE to 1/0
    } else if(is.null(x) || length(x) == 0) {
      return(NA_real_)
    } else if(!is.numeric(x)) {
      return(as.numeric(x))  # Force conversion to numeric
    } else if(length(x) > 1) {
      return(x[1])  # Take first value if somehow a vector
    } else {
      return(as.numeric(x))  # Ensure it's numeric type
    }
  })
  return(metrics)
}

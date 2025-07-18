% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strata_metrics.R
\name{StrataMetrics}
\alias{StrataMetrics}
\title{Calculate Stratum-Based LiDAR Metrics}
\usage{
StrataMetrics(z, r)
}
\arguments{
\item{z}{Numeric vector. Elevation (\code{Z}) values of LiDAR returns (meters above ground or DEM).}

\item{r}{Integer or logical vector. Return number (\code{ReturnNumber}) of each point (1 typically indicates first return).}
}
\value{
A named list of stratum metrics. For each height bin (stratum), the following statistics
are returned for both all and first returns:
\itemize{
\item Coefficient of variation (CV)
\item Kurtosis
\item Skewness
\item Standard deviation
\item Total number of returns
\item Proportion of returns in that bin
}
}
\description{
Computes vertical distribution metrics for LiDAR return heights, stratified by predefined
elevation bins (e.g., 0–0.5 m, 0.5–1 m, ..., 64+ m). Returns statistics for both all returns
and first returns in each vertical bin. This is useful for assessing canopy structure, layering,
and vertical complexity in forest stands.
}
\details{
The vertical strata are defined as follows (in meters):
\preformatted{
  0–0.5, 0.5–1, 1–2, 2–4, 4–8, 8–16, 16–32, 32–48, 48–64, 64+
}
Prefixes:
\describe{
\item{ALL_RETURNS_strata_}{All returns in each stratum}
\item{FIRST_RETURNS_strata_}{First returns in each stratum}
}

Skewness and kurtosis are calculated internally to avoid reliance on external packages.
This function assumes that \code{z} is already normalized to height above ground.
}
\note{
This is a memory-optimized alternative to more verbose stratum metric pipelines,
particularly suitable for processing in tiling frameworks such as \code{catalog_apply()} and \verb{*_metrics()} in \code{lidR}.
}
\examples{
\dontrun{
strata <- StrataMetrics(z = las$Z, r = las$ReturnNumber)
}

}

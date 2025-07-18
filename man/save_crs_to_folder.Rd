% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{save_crs_to_folder}
\alias{save_crs_to_folder}
\title{Save CRS (Coordinate Reference System) to a Folder as Text File}
\usage{
save_crs_to_folder(crs, out_dir)
}
\arguments{
\item{crs}{An object of class \code{crs} or string. Usually obtained via \code{terra::crs()} or \code{sf::st_crs()}.}

\item{out_dir}{Output directory path where the \code{crs.txt} file should be saved.}
}
\value{
Called for side effects. Creates a file at \code{out_dir/crs.txt}.
}
\description{
Writes the string representation of a CRS object to a \code{crs.txt} file in a folder.
}

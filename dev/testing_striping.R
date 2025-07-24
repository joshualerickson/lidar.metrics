library(future)
library(lidR)
library(terra)
library(dplyr)
library(sf)
library(lidar.metrics)
# plan(multisession, workers=15)
# set_lidr_threads(2)

ctg_norm_clip<- readLAScatalog("/home/josh.erickson/Documents/projects/vbflood/dev/flathead_testing/")
load_graph_deps()

opt_chunk_size(ctg_norm_clip) <- 0
opt_chunk_buffer(ctg_norm_clip) <- 30
opt_chunk_alignment(ctg_norm_clip) <- c(0, 0)
opt_stop_early(ctg_norm_clip) <- FALSE

table(las@data$Classification)

opt_filter(ctg_norm_clip) <- paste(
  "-drop_class 7 9",            # drop noise and water
  "-drop_withheld",            # drop withheld (often invalid points)
  "-drop_overlap",
  "-drop_below 0.3",
  "-drop_above 12.1"
)

opt_filter(ctg_norm_clip) <- paste(
  "-drop_class 7 9",            # drop noise and water
  "-drop_withheld",
  "-drop_intensity_below 5",  # low signal = noise
  "-drop_overlap",
  "-drop_scan_angle_above 25",  # limit to near-nadir points
  "-drop_scan_angle_below -25"
)

testing_canopy_cover <- pixel_metrics(ctg_norm_clip[1,],~canopy_cover_metrics(X,Y,Z, return_number = ReturnNumber),
                                      res=30)


testing_canopy_cover <- pixel_metrics(ctg_norm_clip,~lidar.metrics::canopy_cover_metrics(X,Y,Z, ReturnNumber),
                                      res=30)
)


testing_density <- pixel_metrics(ctg_norm_clip,~list(density = length(Z)/(30*30)),
                                 res=30)

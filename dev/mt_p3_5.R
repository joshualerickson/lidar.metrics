
Sys.setenv(TMPDIR = "/mnt/nvmeDrive/data/LIDAR/tmp/flathead_north")

library(future)
library(lidR)
library(terra)
library(dplyr)
terra::terraOptions(tempdir = "/mnt/nvmeDrive/data/LIDAR/tmp/flathead_north")
library(sf)
plan(multisession, workers=20)
set_lidr_threads(1)

ctg_norm_sub <- readLAScatalog(c(
#"/mnt/nvmeDrive/data/LIDAR/MT_P3_4_B21/normalized",
                                 "/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_5_B21/normalized/"
 #                                "/mnt/nvmeDrive/data/LIDAR/MT_Statewide_Ravalli/MT_Statewide_Ravalli_6/normalized"
 ))
#ctg_norm_sub <- ctg_norm_sub[st_coordinates(st_centroid(st_as_sf(ctg_norm_sub)))[,1] < 4e5 & st_coordinates(st_centroid(st_as_sf(ctg_norm_sub)))[,2] > 400000, ]

#install if you haven't already
#remotes::install_github('joshualerickson/lidar.metrics')

library(lidar.metrics)
load_graph_deps()

#for graph metrics
opt_filter(ctg_norm_sub) <- paste(
  "-drop_class 7 9",            # drop noise and water
  "-drop_withheld",            # drop withheld (often invalid points)
  "-drop_z_below 1",         # avoid ground clutter
  "-drop_z_above 12.1",        # truncate canopy top
  "-drop_overlap"
)

# for tree/canopy metrics
# opt_filter(ctg_norm_sub) <- paste(
#   "-drop_class 7 9",            # drop noise and water
#   "-drop_withheld",            # drop withheld (often invalid points)
#   "-drop_overlap"
# )


opt_chunk_size(ctg_norm_sub) <- 0
opt_chunk_buffer(ctg_norm_sub) <- 30
opt_chunk_alignment(ctg_norm_sub) <- c(0, 0)
opt_stop_early(ctg_norm_sub) <- FALSE

#### make sure to change suffix!!!!!!!!!!!!!!!!!!
opt_output_files(ctg_norm_sub) <- paste0('/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_5_B21/graph_metrics/',"{ORIGINALFILENAME}_graph_metrics")

pixel_metrics(ctg_norm_sub,~lidar.metrics::connectivity_metrics_binned(X,Y,Z,
                                                                      scan_angle = ScanAngle,
                                                                      psid = PointSourceID,
                                                                      return_number = ReturnNumber,
                                                                      z_1 =  1,
                                                                      z_20 = 6,
                                                                      z_40 = 12,
                                                                      voxel_res = 3,
                                                                      edge_thresh_values = c(3,3),
                                                                      QL1 = FALSE,
                                                                      density = 4,
                                                                      dec_res = 30
                                                        ),
                         res=30)

#
# pixel_metrics(ctg_norm_sub,~lidar.metrics::tree_detection(X,Y,Z, return_number = ReturnNumber),
#               res=30)


mos <- list.files('/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_5_B21/graph_metrics/')

mos_col <- terra::sprc(paste0('/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_5_B21/graph_metrics/', mos))

mos <- mosaic(mos_col)

writeRaster(mos, '/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_5_B21/mosaic/MT_Statewide_P3_5_B21_graph_metrics_ql2.tif')
#

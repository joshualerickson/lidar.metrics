
Sys.setenv(TMPDIR = "/mnt/nvmeDrive/data/LIDAR/tmp/flathead_south")

library(future)
library(lidR)
library(terra)
library(sf)
library(lidar.metrics)
terraOptions(tempdir = "/mnt/nvmeDrive/data/LIDAR/tmp/flathead_south")
plan(multisession, workers=16)
set_lidr_threads(1)

sfh <-c("/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_1_B21/normalized"
        #"/mnt/nvmeDrive/data/LIDAR/MT_P3_4_B21/normalized"
        )
ctg_norm <- readLAScatalog(sfh)
ctg_norm <- ctg_norm[st_coordinates(st_centroid(st_as_sf(ctg_norm)))[,1] < 4e5 & st_coordinates(st_centroid(st_as_sf(ctg_norm)))[,2] < 400000, ]

load_graph_deps()

plot(ctg_norm)
opt_chunk_size(ctg_norm) <- 0
opt_chunk_buffer(ctg_norm) <- 30
opt_chunk_alignment(ctg_norm) <- c(0, 0)
opt_stop_early(ctg_norm) <- FALSE

#### make sure to change suffix!!!!!!!!!!!!!!!!!!
opt_output_files(ctg_norm) <- paste0('/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_1_B21/graph_metrics/',"{ORIGINALFILENAME}_graph_metrics")
#opt_output_files(ctg_norm) <- paste0('/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_1_B21/tree_metrics/',"{ORIGINALFILENAME}_tree_metrics")


opt_filter(ctg_norm) <- paste(
  "-drop_class 7 9",            # drop noise and water
  "-drop_withheld",            # drop withheld (often invalid points)
  "-drop_overlap",
  "-drop_below 1",
  "-drop_above 12.1"
)

opt_filter(ctg_norm) <- paste(
  "-drop_class 7 9",            # drop noise and water
  "-drop_withheld",            # drop withheld (often invalid points)
  "-drop_overlap"
)
mi

pixel_metrics(ctg_norm,~lidar.metrics::connectivity_metrics_binned(X,Y,Z,
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

pixel_metrics(ctg_norm,~lidar.metrics::tree_detection(X,Y,Z, return_number = ReturnNumber),
              res=30)


mos <- list.files('/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_1_B21/graph_metrics/')

mos_col <- terra::sprc(paste0('/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_1_B21/graph_metrics/', mos))

mos <- mosaic(mos_col)

writeRaster(mos, '/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_1_B21/mosaic/MT_Statewide_P3_1_B21_graph_metrics.tif')
#

mos <- list.files('/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_5_B21/graph_metrics/')

mos_col <- terra::sprc(paste0('/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_5_B21/graph_metrics/', mos))

mos <- terra::mosaic(mos_col)

terra::writeRaster(mos, '/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_5_B21/mosaic/MT_Statewide_P3_5_B21_graph_metrics_checkup.tif')
#

mos <- list.files('/mnt/nvmeDrive/data/LIDAR/MT_Lincoln/graph_metrics/')

mos_col <- terra::sprc(paste0('/mnt/nvmeDrive/data/LIDAR/MT_Lincoln/graph_metrics/', mos))

mos <- terra::mosaic(mos_col)

terra::writeRaster(mos, '/mnt/nvmeDrive/data/LIDAR/MT_Lincoln/mosaic/MT_Lincoln_graph_metrics_checkup.tif', overwrite = TRUE)
#

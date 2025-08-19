
Sys.setenv(TMPDIR = "/mnt/nvmeDrive/data/LIDAR/tmp/lc")

library(future)
library(lidR)
library(terra)
library(sf)
library(lidar.metrics)

terra::terraOptions(tempdir = "/mnt/nvmeDrive/data/LIDAR/tmp/lc")
plan(multisession, workers=20)
set_lidr_threads(1)


lincoln_mt <- "/mnt/nvmeDrive/data/LIDAR/MT_Lincoln/Lincoln_normalized_transfer/lincoln_normalized_meters/"

ctg_norm <- readLAScatalog(lincoln_mt)

load_graph_deps()

opt_chunk_size(ctg_norm) <- 0
opt_chunk_buffer(ctg_norm) <- 30
opt_chunk_alignment(ctg_norm) <- c(0, 0)
opt_stop_early(ctg_norm) <- FALSE
opt_output_files(ctg_norm) <- paste0('/mnt/nvmeDrive/data/LIDAR/MT_Lincoln/tree_metrics/',"{ORIGINALFILENAME}_tree_metrics")
# 👇 Here's the key line
#opt_restart(ctg_norm) <- 17947

opt_filter(ctg_norm) <- paste(
  "-drop_class 7 9",            # drop noise and water
  "-drop_withheld",            # drop withheld (often invalid points)
  "-drop_overlap"
)


pixel_metrics(ctg_norm,~lidar.metrics::tree_detection(X,Y,Z, return_number = ReturnNumber),
              res=30)
              
pixel_metrics(ctg_norm,~lidar.metrics::connectivity_metrics_binned(X,Y,Z, scan_angle = ScanAngle,
                                                         z_1 =  1,
                                                         z_20 = 6,
                                                         z_40 = 12,
                                                         voxel_res = 3,
                                                         edge_thresh_values = c(3,3)),
              res=30)



test


mos <- list.files('/mnt/nvmeDrive/data/LIDAR/MT_Lincoln/graph_metrics/')

mos_col <- terra::sprc(paste0('/mnt/nvmeDrive/data/LIDAR/MT_Lincoln/graph_metrics/', mos))

mos <- mosaic(mos_col)

writeRaster(mos, '/mnt/nvmeDrive/data/LIDAR/MT_Lincoln/mosaic/MT_Lincoln_graph_metrics.tif')

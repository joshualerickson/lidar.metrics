
Sys.setenv(TMPDIR = "/mnt/nvmeDrive/data/LIDAR/tmp/flathead_south")

library(future)
library(lidR)
library(terra)
library(sf)
library(lidar.metrics)
terraOptions(tempdir = "/mnt/nvmeDrive/data/LIDAR/tmp/flathead_south")
plan(multisession, workers=20)
set_lidr_threads(1)

sfh <-c("/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_1_B21/normalized",
        "/mnt/nvmeDrive/data/LIDAR/MT_P3_4_B21/normalized")
ctg_norm <- readLAScatalog(sfh)
ctg_norm <- ctg_norm[st_coordinates(st_centroid(st_as_sf(ctg_norm)))[,1] < 4e5 & st_coordinates(st_centroid(st_as_sf(ctg_norm)))[,2] < 400000, ]

dlist=list.files('/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_1_B21/tree_metrics/', full=F, pattern=".tif")

if(length(dlist)>0) {
	dlist=gsub("_tree_metrics.tif", ".laz", dlist)
	flist=basename(ctg_norm@data$filename)	
	testdone=which(flist %in% dlist)
	ctg_norm$processed <-TRUE 
	ctg_norm$processed[testdone] <-FALSE 
}

load_graph_deps()

opt_chunk_size(ctg_norm) <- 0
opt_chunk_buffer(ctg_norm) <- 30
opt_chunk_alignment(ctg_norm) <- c(0, 0)
opt_stop_early(ctg_norm) <- FALSE
#opt_output_files(ctg_norm) <- paste0('/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_1_B21/graph_metrics/',"{ORIGINALFILENAME}_graph_metrics")
opt_output_files(ctg_norm) <- paste0('/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_1_B21/tree_metrics/',"{ORIGINALFILENAME}_tree_metrics")
# 👇 Here's the key line
opt_restart(ctg_norm) <- 4242

opt_filter(ctg_norm) <- paste(
  "-drop_class 7 9",            # drop noise and water
  "-drop_withheld",            # drop withheld (often invalid points)
  "-drop_intensity_below 5",  # low signal = noise
  "-drop_overlap",
  "-drop_below 0.3",
  "-drop_above 12.1"
)

opt_filter(ctg_norm) <- paste(
  "-drop_class 7 9",            # drop noise and water
  "-drop_withheld",            # drop withheld (often invalid points)
  "-drop_overlap"
)


pixel_metrics(ctg_norm,~lidar.metrics::connectivity_metrics_binned(X,Y,Z, scan_angle = ScanAngle,
                                                                   z_1 =  1,
                                                                   z_20 = 6,
                                                                   z_40 = 9,
                                                                   voxel_res = 3,
                                                                   edge_thresh_values = c(3,3)),
              res=30)

pixel_metrics(ctg_norm,~lidar.metrics::tree_detection(X,Y,Z, return_number = ReturnNumber),
              res=30)


mos <- list.files('/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_1_B21/tree_metrics/')

mos_col <- terra::sprc(paste0('/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_1_B21/tree_metrics/', mos))

mos <- mosaic(mos_col)

writeRaster(mos, '/mnt/nvmeDrive/data/LIDAR/MT_Statewide_P3_1_B21/mosaic/MT_Statewide_P3_1_B21_tree_metrics.tif')

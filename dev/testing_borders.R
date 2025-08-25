library(future)
library(lidR)
library(terra)
library(dplyr)
library(sf)
library(lidar.metrics)
plan(multisession, workers=10)
set_lidr_threads(1)
ctg_norm_clip_fc<- readLAScatalog("/mnt/DataDrive1/data/LIDAR/flathead_testing/flathead")
ctg_norm_clip_lc<- readLAScatalog("/mnt/DataDrive1/data/LIDAR/flathead_testing/lc_border")
plot(ctg_norm_clip_fc)
plot(ctg_norm_clip_lc)

borders <- st_intersection(ctg_norm_clip_fc@data, st_transform(ctg_norm_clip_lc@data, st_crs(ctg_norm_clip_fc@data)))
borders_lc <- st_intersection(ctg_norm_clip_lc@data, st_transform(ctg_norm_clip_fc@data, st_crs(ctg_norm_clip_lc@data)))

lc_borders <- ctg_norm_clip_lc[ctg_norm_clip_lc@data$filename %in% borders_lc$filename,]
fc_borders <- ctg_norm_clip_fc[ctg_norm_clip_fc@data$filename %in% borders$filename,]

# for tree/canopy metrics
opt_filter(lc_borders) <- paste(
  "-drop_class 7 9",            # drop noise and water
  "-drop_withheld",            # drop withheld (often invalid points)
  "-drop_overlap"
)


opt_chunk_size(lc_borders) <- 0
opt_chunk_buffer(lc_borders) <- 30
opt_chunk_alignment(lc_borders) <- c(0, 0)
opt_stop_early(lc_borders) <- FALSE

testing_lc <- pixel_metrics(lc_borders,~connectivity_metrics_binned(X,Y,Z,
                                                                    scan_angle = ScanAngle,
                                                                    psid = PointSourceID,
                                                                    return_number = ReturnNumber,
                                                                    z_1 =  1,
                                                                    z_20 = 6,
                                                                    z_40 = 12,
                                                                    voxel_res = 3,
                                                                    edge_thresh_values = c(3,3),
                                                                    QL1 = TRUE,
                                                                    density = 15,
                                                                    dec_res = 3
),
res=30)
testing_lc_cc <- pixel_metrics(lc_borders,~canopy_cover_metrics(X,Y,Z,
                                                                psid = PointSourceID,
                                                                return_number = ReturnNumber,
                                                                QL1 = T,
                                                                n = 5,
                                                                res = 10

),
res=30)

plot(testing_fc_cc[[20]])
plot(testing_lc_cc[[20]])

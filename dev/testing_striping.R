library(future)
library(lidR)
library(terra)
library(dplyr)
library(sf)
library(lidar.metrics)
remotes::install_github('joshualerickson/lidar.metrics')
 plan(multisession, workers=15)
 set_lidr_threads(1)

ctg_norm_clip<- readLAScatalog("/home/josh.erickson/Documents/projects/vbflood/dev/flathead_testing/")
ctg_norm_clip<- readLAScatalog("/mnt/boreas/lidar_download/LincolnCounty/lincoln_normalized_meters/")
ctg_norm_clip<- readLAScatalog("/mnt/DataDrive1/data/LIDAR/MT_P3_3_B21/normalized/")
plot(ctg_norm_clip)
ctg_norm_clip
plot(ctg_norm_clip@data$geometry)
m <- mapview::mapview(ctg_norm_clip@data)@map
bb <- mapedit::drawFeatures(map = ) %>% st_transform(st_crs(ctg_norm_clip))

ctg_norm_clip2 <- lidR::catalog_intersect(ctg_norm_clip, bb)
plot(ctg_norm_clip2)
load_graph_deps()

opt_chunk_size(ctg_norm_clip2) <- 0
opt_chunk_buffer(ctg_norm_clip2) <- 30
opt_chunk_alignment(ctg_norm_clip2) <- c(0, 0)
opt_stop_early(ctg_norm_clip2) <- FALSE

opt_filter(ctg_norm_clip2) <- paste(
  "-drop_class 7 9",            # drop noise and water
  "-drop_withheld",            # drop withheld (often invalid points)
  "-drop_overlap",
  "-drop_below 0.3",
  "-drop_above 12.1"
)

opt_filter(ctg_norm_clip) <- paste(
  "-drop_class 7 9",            # drop noise and water
  "-drop_withheld",
  "-drop_overlap"
)

testing_canopy_cover <- pixel_metrics(ctg_norm_clip[1,],~canopy_cover_metrics(X,Y,Z, return_number = ReturnNumber, zmax = 50),
                                      res=30)

plot(testing_canopy_cover)

testing_density <- pixel_metrics(ctg_norm_clip2,~list(density = length(Z)/(30*30),
                                                     mean_scan_angle = mean(ScanAngle, na.rm = T)),
                                 res=30)

plot(ctg_norm_clip2[c(70:73,82:86),])

testing_graph_chv_lc<- pixel_metrics(ctg_norm_clip2,~lidar.metrics::connectivity_metrics_binned(X, Y, Z,
                                                                              scan_angle = ScanAngle,
                                                                              psid = PointSourceID,

                                                                              QL1 = TRUE),
                                 res=30)

testing_graph_chv_lc <- pixel_metrics(ctg_norm_clip2,~connectivity_metrics_binned(X, Y, Z, return_number = ReturnNumber,
                                                                              scan_angle = ScanAngle, density = 4,
                                                                              psid = PointSourceID, QL1 = TRUE),
                                 res=30)

testing_graph_chv_lc<- pixel_metrics(ctg_norm_clip2[10:30,],~list(psid_n = length(unique(PointSourceID)),
                                                                   psid_dom = as.integer(names(which.max(table(PointSourceID)))),
                                                                   psid_p = max(table(PointSourceID))/sum(table(PointSourceID)),
                                                         n = length(Z),
                                                         psid_n = max(table(PointSourceID))),
                                 res=30)

mapview::mapview(testing_graph_chv_lc[[9]])
plot(testing_graph_chv_lc[[4:5]])
plot(testing_graph_chv_lc[[c(2:3, 9)]])
plot(testing_graph_chv_lc[[c(19, 9)]])
plot(testing_graph_chv_lc[[c(1:10)]])
plot(testing_graph_chv_lc[[c(11:20)]])
plot(testing_graph_chv_lc[[c(9)]])
plot(testing_graph_chv_fc[[3]]< 0.6)
plot(testing_graph_chv_lc[[c(3,7, 9:10)]])
plot(testing_graph_chv_lc[[c(10)]])

plot(testing_graph_chv_lc2[[c(9, 3)]])
plot(testing_graph_chv_lc2[[c(3,7, 9:10)]])
plot(testing_graph_chv_lc2[[c(10)]])
plot(c(testing_graph_chv_lc[[1:2]], (testing_graph_pc_lc[[1:2]])))
plot(testing_graph[[c(1, 9)]])

names(testing_graph_pc) <- paste0(names(testing_graph_pc), '_pc')
names(testing_graph_chv) <- paste0(names(testing_graph_chv), '_chv')
testing_df <- as.data.frame(c(testing_graph, testing_graph_pc,testing_graph_chv_voxel_norm))
testing_df <- as.data.frame(c(testing_graph_chv, testing_graph_pc))

library(dplyr)
library(ggplot2)
library(GGally)

testing_df %>%
  select(starts_with('understory')) %>%
  #glimpse() %>%
  ggpairs()
testing_df %>% slice_sample(n = 5000) %>%
  # filter(understory_mean_var_xyz < 1, understory_mean_var_xyz > 0.5,
  #                     understory_n_m2 > 2.5) %>%
  filter(understory_mean_degree_chv > 0) %>%
  mutate(cut_bt = cut(understory_mean_betweenness_chv, 4)) %>%
  ggplot(aes(understory_mean_abs_sa_pc, understory_n_m2_pc
             )) + geom_point(aes(
                            color = understory_mean_betweenness_chv
                            )) +
  scale_color_gradientn(colors = hcl.colors(11, 'Zissou1'))  +
  geom_smooth()

  facet_wrap(~cut_bt)


  choose_swaths_robust <- function(psid, sa,
                                   p_thresh = 0.60,      # dominance level
                                   delta_min = 0.15,     # min gap between top1 and top2
                                   n_min = 500,          # need enough points to trust 1-swath
                                   max_keep = 2L,
                                   h_thresh = 0.60) {    # entropy threshold (normalized)
    tab <- sort(table(psid), decreasing = TRUE)
    n   <- sum(tab)

    if (length(tab) == 0) return(list(keep_ids = integer(0), p1 = NA, p2 = NA, n = 0))
    if (length(tab) == 1L) return(list(keep_ids = as.integer(names(tab)[1L]),
                                       p1 = 1, p2 = 0, n = n))

    p <- as.numeric(tab) / n
    p1 <- p[1]; p2 <- p[2]

    # normalized entropy (0 = one swath dominates, 1 = all equal)
    H  <- -sum(p * log(p))
    Hn <- H / log(length(p))

    one_swath_ok <- (n >= n_min) && (p1 >= p_thresh) && ((p1 - p2) >= delta_min) && (Hn <= h_thresh)

    if (one_swath_ok) {
      keep_ids <- as.integer(names(tab)[1L])
    } else {
      k <- min(max_keep, length(tab))
      keep_ids <- as.integer(names(tab)[seq_len(k)])

      # tie-break last slot by lowest mean |SA|
      cutoff <- tab[k]
      tied   <- names(tab)[tab == cutoff]
      if (length(tied) > 1L) {
        m_sa <- tapply(abs(sa), psid, function(v) mean(v, na.rm = TRUE))
        fixed <- as.integer(names(tab)[seq_len(k - 1L)])
        best  <- as.integer(names(sort(m_sa[tied], decreasing = FALSE))[1L])
        keep_ids <- unique(c(fixed, best))
      }
    }
    list(keep_ids = keep_ids, p1 = p1, p2 = p2, n = n, Hn = Hn)
  }


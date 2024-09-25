#####Understory vegetation data: Pre-processing#####

##Goal: Load understory vegetation data & calculate summaries for area-based analysis 
##DATE: 25 July 2023
##Author: Dr. Marc Mayes, SIG-NAL/USU/UCSB

#R packages
library(tidyverse)
library(sf)
library(lubridate)

#####1. Load line-intercept data#####
#data input directory:
indir <- "./data/"

#output directory
outdir <- "./output/"

#line-intercept data file:
line_intercept_infile <- "ebrook_lineintercept_Rin_20230724.csv"

tsect_lin <- read.csv(paste0(indir, line_intercept_infile))
head(tsect_lin)


#####2. Line-intercept preprocess calculations#####

#create grouping variable "sample_month" for month, then recode to match sampling period
tsect_lin$sample_month <- as.character(substr(tsect_lin$date, 5,6))
tsect_lin$obs <- dplyr::recode(tsect_lin$sample_month,
                               "05" = "initial_may",
                               "06" = "pg01_june",
                               "07" = "pg02_july")

#create variable for transect and obs period
tsect_lin$tsect_obs <- paste0(tsect_lin$transect,"_",
                              tsect_lin$obs)
#calculate observation lengths: subtract dist_stop_m - dist_start_m.
tsect_lin$obs_length_m <- tsect_lin$dist_stop_m - tsect_lin$dist_start_m

#TEST: if you sum observation lengths by observation period (tsect_obs), you should get about 20 except for transects that are of known shorter distances.
tsect_obsdist_check <- tsect_lin %>%
  dplyr::group_by(tsect_obs) %>%
  summarize(obs_length_tot_m = sum(obs_length_m))
#write.csv(tsect_obsdist_check, paste0(outdir, "transect_tot_distances_check.csv"))

#not as elegant as mutate...but: left-join total transect obs distances.
tsect_obsdist_calc <- left_join(tsect_lin,
                                tsect_obsdist_check,
                                by="tsect_obs")

#calculate the percentage of transect length for each record
tsect_obsdist_calc$pct_tsect <- round(tsect_obsdist_calc$obs_length_m/tsect_obsdist_calc$obs_length_tot_m*100, 2)

#calculate weighted-average height of ALL materials (including bare ground) by transect and obs period:
tsect_hts_wtavg_all <- tsect_obsdist_calc %>%
  group_by(tsect_obs) %>%
  summarize(ht_wtavg_m = sum(height_m*(pct_tsect/100),na.rm=T))
#write this out as a .csv:
#write.csv(tsect_hts_wtavg_all, paste0(outdir, "transect_height_weightedavg_allrecords.csv"))

#calculate average height and cover of materials:
tsects_cover_pct_ht <- tsect_obsdist_calc %>%
  group_by(tsect_obs, veg_type, .drop=FALSE) %>%
  summarize(cover_pct = sum(pct_tsect),
            ht_avg_m = mean(height_m, na.rm=T))
#write this out as .csv:
#write.csv(tsects_cover_pct_ht, paste0(outdir, "transect_coverpct-htm_allrecords.csv"))



##DRAFT##
# tsect_lin_test <- tsect_lin %>%
#   dplyr::group_by(tsect_obs, veg_type) %>%
#   summarize(veg_obs_length_m = sum(obs_length_m)) %>%
#   ungroup() %>%
#   group_by()
#     obs_length_pct = round((veg_obs_length_m/obs_length_tot_m)*100, 2))

            

#Introduction#
#This R function removes the suggested flags for pixel data extracted from Sentinel 3 OLCI Chlorophyll-a Neural Network Product as described in "PAPER NAME".
# Carla L. Mejias-Rivera
# carla.mejias@upr.edu



process_pixEx <- function(file_path) {
  
  # Load necessary packages
  library(readr)
  library(dplyr)
  library(tidyr)
  
  # Define the flags to filter
  Flags_NN <- c("WQSF_lsb_LAND", 
                "WQSF_lsb_CLOUD", 
                "WQSF_lsb_CLOUD_AMBIGUOUS",
                "WQSF_lsb_CLOUD_MARGIN",
                "WQSF_lsb_INVALID",
                "WQSF_lsb_COSMETIC", 
                "WQSF_lsb_SATURATED", 
                "WQSF_lsb_SUSPECT",
                "WQSF_lsb_HISOLZEN", 
                "WQSF_lsb_HIGHGLINT",
                "WQSF_lsb_OCNN_FAIL", 
                "WQSF_lsb_MEGLINT",
                "WQSF_lsb_COASTLINE",
                "WQSF_lsb_WHITECAPS",
                "WQSF_lsb_ADJAC")
  
  # Read and process the Extracted Pixel file
  processed_data <- read_csv(file_path, col_types = cols(`Date` = col_date(format = "%m/%d/%Y"))) %>% 
    mutate(ALL_FLAGS = rowSums((.[, Flags_NN]), na.rm = TRUE)) %>% 
    filter(ALL_FLAGS < 1) %>% 
    rename(Site = "Name",
           Kd490 = "KD490_M07") %>% 
    group_by(Date, Site) %>% 
    subset(select = c(Site, Latitude, Longitude, Date, CHL_NN, CHL_NN_err, Kd490, ALL_FLAGS)) %>% 
    dplyr::mutate(
      Optical_Depth = 1 / Kd490) %>%
    replace_with_na_all(condition = ~.x == -999)
  
  return(processed_data)
}

DESIRED_FILE_NAME <- process_pixEx("INSERT YOUR PIXEL EXTRACTION FILE PATH HERE")

PRCRMP_Flags <- process_pixEx("pixEx_PRCRMP_OL_2_WFR_measurements.csv")
Virtual_Flags <-process_pixEx("pixEx_Virtual_OL_2_WFR_measurements.csv")

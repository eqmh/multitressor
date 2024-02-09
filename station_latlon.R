library(dplyr)
setwd("~/Google Drive/My Drive/GDrive/proposals/2022_02_MultiStressor_NOAA/module_1")

df <- read.csv("WS_Cruise-Stn-Coords(2015-2023).csv")

# # Make sure all longitudes are negative (some values are positive in the original table)
# positive_vals<- which(df$dec_lon > 0)

unique_stations <- df %>% 
  group_by(station_id) %>% 
  summarise(
    mean_lon = mean(dec_lon), 
    mean_lat = mean(dec_lat), 
    std_lon = sd(dec_lon), 
    std_lat = sd(dec_lat)) %>% 
      ungroup()

write.csv(unique_stations, "unique_stations.csv", row.names = FALSE)
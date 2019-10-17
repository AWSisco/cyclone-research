library(ggplot2)
library(magrittr)
library(lubridate)
library(cowplot)
library(plyr)
library(dplyr)
#==================================================================================================
#Read cyclones.dat and add some additional fields

cyclones <- read.delim(file = "cyclones.dat", header = FALSE, sep = "",
                       col.names = c("year", "month", "day", "hour", "lat", "lon", "vort"))

dim(cyclones) #Raw data has 26414 cases

cyclones <- cyclones[-c(which(cyclones$hour==12)),] #Remove all 1200 UTC records from the dataset
#Removes 6677 cases for a new total of 19737

#Function to generate a date and time string
genDateTime <- function(x) {
  r <- rownames(x)
  paste(as.character(x[r,1]),
        as.character(x[r,2]),
        as.character(x[r,3]),
        as.character(x[r,4]), sep = "-") %>% ymd_h()
}

dateTime <- genDateTime(cyclones) #Generate date and time strings in YYYY-MM-DD HH:MM:SS UTC format
cyclones$dateTime <- dateTime #Attach a dateTime column

rm(dateTime, genDateTime)

#Subset into winter cyclones
winter_cyclones <- subset(cyclones, !(month >= 8 & month <= 10)) #13847 cases

#==================================================================================================
#Count the number of cyclones at each 2.5 degree point
#Write a 2D (25 x 144) text file containing the number of cyclones at each grid point
winter_counts <- ddply(winter_cyclones, .(winter_cyclones$lat, winter_cyclones$lon), nrow)
names(winter_counts) <- c("lat", "lon", "count")

lat <- seq(20.0, 80.0, 2.5)
lon<- seq(0,357.5, 2.5)
lat_lon<- expand.grid(lat,lon)
colnames(lat_lon) <- c("lat", "lon")
merged_counts <- merge(winter_counts, lat_lon, all = TRUE)
merged_counts$count[is.na(merged_counts$count)] <- 0
count_data <- matrix(merged_counts$count, nrow = 25, byrow = TRUE)
write.table(count_data, file = "winter_cyclone_counts.txt", sep = "\t", col.names = FALSE, row.names = FALSE)
rm(winter_counts, lat, lon, lat_lon, merged_counts, count_data)
#==================================================================================================

#Subset by basin
wnp_cyclones <- subset(winter_cyclones, (lat >= 30.0 & lat <= 45.0 & lon >= 135 & lon <= 195)) #1208
atl_cyclones <- subset(winter_cyclones, (lat >= 30.0 & lat <= 45.0 & lon >= 285 & lon <= 345)) #825

#==================================================================================================
#The box wi11 start at the lowest latitude and shift north 5 times by 2.5 degrees. This means I need to
#construct coordinates for 5 boxes.

wnp_boxes <- data.frame(S = seq(30,40,2.5), N = seq(35,45,2.5), W = rep(135,5), E = rep(195,5))
atl_boxes <- data.frame(S = seq(30,40,2.5), N = seq(35,45,2.5), W = rep(285,5), E = rep(345,5))
colnames(wnp_boxes) <- c("S", "N", "W", "E")
colnames(atl_boxes) <- c("S", "N", "W", "E")

count_cyclones <- function(basin_boxes, basin_cyclones){
  basin_count <- c()
  
  for (i in 1:dim(basin_boxes)[1]){
    
    box_subset <- subset(basin_cyclones, 
                         lat >= basin_boxes[i,1] & #Lat min
                           lat <= basin_boxes[i,2] & #Lat max
                           lon >= basin_boxes[i,3] & #Lon min
                           lon <= basin_boxes[i,4])  #Lon max
    basin_count[i] <- dim(box_subset)[1]
  }
  
  return(basin_count)
}

wnp_boxes$count <- count_cyclones(wnp_boxes, wnp_cyclones)
atl_boxes$count <- count_cyclones(atl_boxes, atl_cyclones)

wnp_plot <- ggplot(wnp_boxes, aes(y = count, x = S)) + geom_col() + coord_flip() + 
  theme_minimal_vgrid() + scale_y_continuous(expand = c(0,0), limits = c(0,750), breaks = seq(0,700,100)) +
  labs(title = "WNP")
atl_plot <- ggplot(atl_boxes, aes(y = count, x = S)) + geom_col() + coord_flip() + 
  theme_minimal_vgrid() + scale_y_continuous(expand = c(0,0), limits = c(0,750), breaks = seq(0,700,100)) +
  labs(title = "ATL")

plot_grid(wnp_plot, atl_plot, ncol = 2)

#==================================================================================================
#Extract the identified cyclones and prepare output files w/ metadata required to obtain
#ERA5 data

wnp_domain <- head(wnp_boxes[order(wnp_boxes$count, decreasing = T),], n = 1)
atl_domain <- head(atl_boxes[order(atl_boxes$count, decreasing = T),], n = 1)

wnp_cases <- subset(wnp_cyclones,
                    lat >= wnp_domain[,1] &
                    lat <= wnp_domain[,2] &
                    lon >= wnp_domain[,3] &
                    lon <= wnp_domain[,4])

atl_cases <- subset(atl_cyclones,
                    lat >= atl_domain[,1] &
                    lat <= atl_domain[,2] &
                    lon >= atl_domain[,3] &
                    lon <= atl_domain[,4])

#Construct a cyclone-specific domain for each case
wnp_cases$domain_lonMin <- wnp_cases$lon - 30.0
wnp_cases$domain_lonMax <- wnp_cases$lon + 90.0
wnp_cases$domain_latMin <- 30.0
wnp_cases$domain_latMax <- 60.0
wnp_cases$wnp_id <- paste0("wnp_", sprintf("%03d", 001:dim(wnp_cases)[1]))

atl_cases$domain_lonMin <- atl_cases$lon - 30.0
atl_cases$domain_lonMax <- atl_cases$lon + 90.0
atl_cases$domain_latMin <- 30.0
atl_cases$domain_latMax <- 60.0
atl_cases$atl_id <- paste0("atl_", sprintf("%03d", 001:dim(atl_cases)[1]))

wnp_era5_metadata <- data.frame(wnp_id = wnp_cases$wnp_id,
                            date = date(wnp_cases$dateTime),
                            time = substr(wnp_cases$dateTime, 12, 16),
                            lat = wnp_cases$lat,
                            lon = wnp_cases$lon,
                            domain_latMin = wnp_cases$domain_latMin,
                            domain_latMax = wnp_cases$domain_latMax,
                            domain_lonMin = wnp_cases$domain_lonMin,
                            domain_lonMax = wnp_cases$domain_lonMax)

atl_era5_metadata <- data.frame(atl_id = atl_cases$atl_id,
                            date = date(atl_cases$dateTime),
                            time = substr(atl_cases$dateTime, 12, 16),
                            lat = atl_cases$lat,
                            lon = atl_cases$lon,
                            domain_latMin = atl_cases$domain_latMin,
                            domain_latMax = atl_cases$domain_latMax,
                            domain_lonMin = atl_cases$domain_lonMin,
                            domain_lonMax = atl_cases$domain_lonMax)

write.table(wnp_era5_metadata, "wnp_era5_metadata.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(atl_era5_metadata, "atl_era5_metadata.txt", quote = FALSE, sep = "\t", row.names = FALSE)

#==================================================================================================
#Prepare output files with metadata needed to obtain GEFS reforecast files

#Generate a reforecast lead time for each cyclone in YYYYMMDDHH format
genLeadTime <- function(basin_df, num_days) {
  
  leadDateTime <- if_else((hour(basin_df$dateTime)==0 | hour(basin_df$dateTime)==6),
                          floor_date(basin_df$dateTime - days(num_days), unit = "day"),
                          ceiling_date(basin_df$dateTime - days(num_days), unit = "day"))
  
  basin_df$leadDateTime <- leadDateTime
  
  return(basin_df)
}

wnp_cases <- genLeadTime(wnp_cases, 2)
atl_cases <- genLeadTime(atl_cases, 2)

wnp_gefs_metadata <- data.frame(dateTime = format(wnp_cases$dateTime, "%Y%m%d%H"),
                                leadDateTime = format(wnp_cases$leadDateTime, "%Y%m%d%H"),
                                lon = format(round(wnp_cases$lon, 1), nsmall = 1))

atl_gefs_metadata <- data.frame(dateTime = format(atl_cases$dateTime, "%Y%m%d%H"),
                                leadDateTime = format(atl_cases$leadDateTime, "%Y%m%d%H"),
                                lon = format(round(atl_cases$lon, 1), nsmall = 1))

write.table(wnp_gefs_metadata ,"wnp_gefs_metadata.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(atl_gefs_metadata,"atl_gefs_metadata.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

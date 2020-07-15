# ====================================================================================================#
# ====================================================================================================#
# Convert presence absence rasters to dataframe
# ====================================================================================================#
# 
library(raster)
library(tidyverse)
library(rgdal)
library(scales)
library(dplyr)

#_____________________________
# Presence Background (Bootstrap method - 10% Sampled) 
# 25 arc min resolution (~45km^2 at equator)

setwd("/Users/rubendouglas/OneDrive - The Royal Botanic Gardens, Kew/Nutridiv R/Agg_threshold_binary/")

#list the tifs, convert the tifs to rastr and stack the raster files 
list_tifs <- as.list(list.files())
tifftorast <- lapply(list_tifs,raster)


PA_mat_25min <- curSpatialdf

for(i in 2:length(tifftorast)){

  curSpatialdf <- as.data.frame(rasterToPoints(tifftorast[[i]]))

  # Merge the current crop production matrix to the preallocated dataframe by xy coordinates
  PA_mat_25min <- left_join(PA_mat_25min, curSpatialdf, by = c("x","y"))
  }

#_____________________________

# copy file
Presabs_mat_25min <- PA_mat_25min

# Remove colnames '_...'
colnames(Presabs_mat_25min) <- sub('\\_.*','',colnames(Presabs_mat_25min))

#_____________________________

# output matrix to use in dbFD calculations:
view(Presabs_mat_25min)




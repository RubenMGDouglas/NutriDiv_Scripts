# ====================================================================================================#
# ====================================================================================================#
# Upscale Production files to 25 arc min (~45km^2 at equator)
# ====================================================================================================#
# The production.tif files contain: 'total  crop  production  in  metric  tons  on  the  land-area  mass  of  a  gridcell.
# Harvested  area  in  hectares  was  multiplied  by  yield  per  hectare  to  create  this data product.'

library(raster)
library(tidyverse)
library(rgdal)
library(scales)
library(dplyr)

#_____________________________

#Load in CHELSA current bio1 map, which has the correct resolution & origin
setwd( "/Users/rubendouglas/OneDrive - The Royal Botanic Gardens, Kew/Nutridiv R/Ruben environmental data/BioClim_25min/")
Bio1_25min <- raster('Bio01_CHELSA_25min.tif')
plot(Bio1_25min)

#_____________________________

# read-in the file names of our crop tifs from the local
# directory as a list that we will iterate over later
setwd( "/Users/rubendouglas/OneDrive - The Royal Botanic Gardens, Kew/Nutridiv R/Production.tif files/")
croptiffs <- list.files()

# Remove file 167 ('vert.txt')
croptiffs <- croptiffs[-167]

#_____________________________

# Then loop to convert each tif file to raster, resample that raster to the same resolution as the Bio1 (25 arc min),
#and then stack them together 
#The bilinear method means that the combined cells have values averaged (mean).

Production_25min <- stack()

for(i in 1:length(croptiffs)){
  
  one <- raster(croptiffs[i])
  
  two <- resample(one, Bio1_25min, method="bilinear")
   
  Production_25min <- stack(Production_25min,two)

}

#_____________________________

# Loop which takes i crop rasters to create a gridcell x crop abundance matrix
 

# pre-allocate some variable space for our output data
preallocate <- curSpatialdf

for(i in 2:length(Production_25min)) {
  
  # read in tif as raster
  curRaster <- Production_25min[[i]]
  
  # convert raster to a dataframe (the individual crop production matrix)
  curSpatialdf <- as.data.frame(rasterToPoints(curRaster))
  
  # Merge the current crop production matrix to the preallocated dataframe by xy coordinates
  preallocate <- merge(preallocate, curSpatialdf, by = c("x","y"))
  
}

# Check for NAs in each column
colSums(is.na(preallocate))

#_____________________________

# Rename to 'Production_mat_25min' & Remove '_...' from column names
Production_mat_25min <- preallocate
colnames(Production_mat_25min) <- sub('\\_.*','',colnames(Production_mat_25min))

#_____________________________

# output matrix to use in dbFD calculations:
Production_mat_25min


# ====================================================================================================#
# ====================================================================================================#
# Environmental & Socioeconomic Drivers of Nutritional FD patterns
# ====================================================================================================#
library(pcaMethods)
library(tidyverse)
library(raster)
library(rgdal)
library(scales)
library(dplyr)
library(sf)
library(sp)
library(stars)
library(Hmisc)
library(data.table)

# Load in PNAS data
setwd("/Users/rubendouglas/OneDrive - The Royal Botanic Gardens, Kew/Nutridiv R/Ruben environmental data/Human factors")

# 1st sheet
PNAS_tibble <- readxl::read_xlsx("pnas.1912710116.sd01.xlsx")

# 2nd sheet
FoodSov_tibble <- readxl::read_xlsx("pnas.1912710116.FoodSov_Indicators.xlsx")


# Harmonize names with previous country datsets
PNAS_tibble$Country <- str_replace_all(PNAS_tibble$Country, c("Libya" = 'Libyan Arab Jamahiriya',"Russian Federation" = 'Russia',
                                                                      "Macedonia" = 'The former Yugoslav Republic of Macedonia', "Viet nam" = 'Viet Nam'))

FoodSov_tibble$Country <- str_replace_all(FoodSov_tibble$Country, c("Libya" = 'Libyan Arab Jamahiriya',"Russian Federation" = 'Russia',
                                                              "Macedonia" = 'The former Yugoslav Republic of Macedonia', "Viet nam" = 'Viet Nam'))

PNAS_tibble <- PNAS_tibble[which(PNAS_tibble$Country %in% SES.FDPD$Country),]

FoodSov_tibble <- FoodSov_tibble[which(FoodSov_tibble$Country %in% SES.FDPD$Country),]

# =========================================

# Run a PCA on all continous variables to identify a handful of key variables to use in downstream analyses
# nhttps://www.r-bloggers.com/principal-component-analysis-in-r/

# pcamethods is favorable here because the 'pca' function can impute missing values via nipals PCA & other algorithms,
# it allows in-built scaling and centering, we can perform cross validation to inform number of PCs to retain and it has intuitive plotting

X<-cbind(PNAS_tibble[,32:78], FoodSov_tibble[,3:47])
rownames(X) <- FoodSov_tibble$Country

PNAS_PCAmeth <- pca(X, scale = "uv", center = TRUE, nPcs = 20, method = "nipals",cv="q2")


slplot(PNAS_PCAmeth, scoresLoadings =c(FALSE,TRUE))

PNAS_PCAmeth@R2cum

PNAS_PCA_loadings <- PNAS_PCAmeth@loadings %>% as.data.frame()

# Which PCs Should be retained?

# Eigenvalues >1 

plot(PNAS_PCAmeth)
# If cross-validation was done for the PCA the plot will also show the CV based statistics.
# A common rule-of-thumb for determining the optimal number of PCs is the PC where the CV diagnostic is at its maximum but not very far from R^2.

scores(PNAS_PCAmeth)

# NB I Should also test for colinearity between variables / PCs when modelling with FD indices

# =========================================

# Still need to account for climate, topography & soil nutrition. 

# Climate

# 4 uncorrelated bioclim variables from Dannys analyses:

# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO12 = Annual Precipitation
# BIO15 = Precipitation Seasonality (Coefficient of Variation)

setwd('/Users/rubendouglas/OneDrive - The Royal Botanic Gardens, Kew/Nutridiv R/Ruben environmental data/BioClim_25min/Danny Uncorrelated/')

Bio1_25min <- raster('Bio01_CHELSA_25min.tif')
Bio2_25min <- raster('Bio02_CHELSA_25min.tif')
Bio12_25min <- raster('Bio12_CHELSA_25min.tif')
Bio15_25min <- raster('Bio15_CHELSA_25min.tif')

# BIO1
out <- raster::extract(Bio1_25min, SpatialPolygons(wrld_simpl@polygons)) # Extract raster data from each of the country polygons in wrld_simpl
Bio1_Country <- data.frame(Country=wrld_simpl$NAME, mean=unlist(lapply(out, mean))) # unlist the data for each country polygon in a dataframe with the mean value in each country
names(Bio1_Country)[2] <- 'BIO1'

# BIO2
out <- raster::extract(Bio2_25min, SpatialPolygons(wrld_simpl@polygons)) # Extract raster data from each of the country polygons in wrld_simpl
Bio2_Country <- data.frame(Country=wrld_simpl$NAME, mean=unlist(lapply(out, mean))) # unlist the data for each country polygon in a dataframe with the mean value in each country
names(Bio2_Country)[2] <- 'BIO2'

# BIO12
out <- raster::extract(Bio12_25min, SpatialPolygons(wrld_simpl@polygons)) # Extract raster data from each of the country polygons in wrld_simpl
Bio12_Country <- data.frame(Country=wrld_simpl$NAME, mean=unlist(lapply(out, mean))) # unlist the data for each country polygon in a dataframe with the mean value in each country
names(Bio12_Country)[2] <- 'BIO12'

# BIO15
out <- raster::extract(Bio15_25min, SpatialPolygons(wrld_simpl@polygons)) # Extract raster data from each of the country polygons in wrld_simpl
Bio15_Country <- data.frame(Country=wrld_simpl$NAME, mean=unlist(lapply(out, mean))) # unlist the data for each country polygon in a dataframe with the mean value in each country
names(Bio15_Country)[2] <- 'BIO15'
  
# =========================================
# =========================================
# Topography

# Global elevation
#GloElev_5min.asc

setwd('/Users/rubendouglas/OneDrive - The Royal Botanic Gardens, Kew/Nutridiv R/Ruben environmental data/Soil & energy nutrition/GloElev_5min')

GloElev_5min <- raster('GloElev_5min.asc')

# Country Extract
out <- raster::extract(GloElev_5min, SpatialPolygons(wrld_simpl@polygons)) # Extract raster data from each of the country polygons in wrld_simpl
GloElev_Country <- data.frame(Country=wrld_simpl$NAME, mean=unlist(lapply(out, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))) # unlist the data for each country polygon in a dataframe with the mean value in each country
names(GloElev_Country)[2] <- 'ELEV'


# =========================================

# Soil

# Soil depth
# REF_DEPTH1.nc

setwd('/Users/rubendouglas/OneDrive - The Royal Botanic Gardens, Kew/Nutridiv R/Ruben environmental data/Soil & energy nutrition/')

REF_DEPTH1 <- raster('REF_DEPTH1.nc')

# Far too high resolution, aggregate by factor of 10 (100x fewer cells)
REF_DEPTH1 <- aggregate(REF_DEPTH1, fact = 10)

# Country Extract
out <- raster::extract(REF_DEPTH1, SpatialPolygons(wrld_simpl@polygons)) # Extract raster data from each of the country polygons in wrld_simpl
REF_DEPTH1_Country <- data.frame(Country=wrld_simpl$NAME, mean=unlist(lapply(out, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))) # unlist the data for each country polygon in a dataframe with the mean value in each country
names(REF_DEPTH1_Country)[2] <- 'DEPTH'

# =========================================

# ISRIC-WISE Soil properties
# https://stackoverflow.com/questions/23568899/access-data-base-import-to-r-installation-of-mdb-tools-on-mac

# WISE soil property summary file
setwd('/Users/rubendouglas/OneDrive - The Royal Botanic Gardens, Kew/Nutridiv R/Ruben environmental data/Soil & energy nutrition/wise_05min_v12/')
WISE5by5min <- mdb.get('WISE5by5min.mdb')
WISEsummary <- WISE5by5min$WISEsummaryFile

# Raster file to link with soil properties
setwd('/Users/rubendouglas/OneDrive - The Royal Botanic Gardens, Kew/Nutridiv R/Ruben environmental data/Soil & energy nutrition/wise_05min_v12/Grid/smw5by5min/')
dblbnd <- raster("dblbnd.adf")

# The gridcell ID corresponding to the WISE summary dataset are in the raster attribute table (RAT)
dblbndattr <- dblbnd$dblbnd@data@attributes[[1]]$ID

# Create a single raster layer with just the gridcell ID
derat_dblbnd <- deratify(dblbnd, 'ID')

# and convert this to a dataframe with lat, long and gridcell ID
rast2pts_dblbnd <- rasterToPoints(derat_dblbnd)#, spatial = TRUE) 
rast2pts_dblbnd <- rast2pts_dblbnd %>% as.data.frame()
rast2pts_dblbnd$ID.1 <- as.factor(rast2pts_dblbnd$ID.1)

# Aggregate the focal numeric varibales as mean values per  gridcell ID
WISEsummary_agg <- aggregate(cbind(TOTC,TOTN,CNrt,PHAQ,CECS,CECc,ECEC,BSAT,TCEQ,GYPS,ELCO,CFRAG,BULK,TAWC) ~ 
                               SUID, WISEsummary, mean )
# Get rid of first row (think this is global)
WISEsummary_agg <- WISEsummary_agg[-1,]

# Aggregate the rast2pts_dblbnd lat and long as mean values per gridcell ID
# cbind(vars, you, want, to, aggregate) ~ vars to aggregate by, dataframe, aggregate function 
rast2pts_dblbnd_agg <- aggregate(cbind(x,y) ~ 
                               ID.1, rast2pts_dblbnd, mean)

rast2pts_dblbnd_agg <- rast2pts_dblbnd_agg[-1,]

soilsummary <- merge(WISEsummary_agg, rast2pts_dblbnd_agg, by.x = 'SUID', by.y = 'ID.1')

soilsummary_spdf <- SpatialPointsDataFrame(coords=soilsummary[16:17], data=soilsummary[1:15])

?extract
proj4string(soilsummary_spdf) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")   
 #spTransform(soilsummary_spdf, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "))

soilsummary <- soilsummary[c(16:17,1:15)]
soilsummary <- data.table(as.data.frame(soilsummary, xy = TRUE))

plot(soilsummary)

# set up an 'empty' raster, here via an extent object derived from your data
e <- extent(soilsummary[,1:2])
#e <- e + 1000 # add this as all y's are the same

r <- raster(e, ncol=864, nrow=346)
# or r <- raster(xmn=, xmx=,  ...

# you need to provide a function 'fun' for when there are multiple points per cell
soil_rastbrick <- rasterize(soilsummary[, 1:2], r, soilsummary[,3:17], fun=mean)
proj4string(soil_rastbrick) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ") 
plot(soil_rastbrick[[10]])

# for each soil property raster layer in the brick, extract mean values for each country polygon in wrld_simpl

lapply(1:nrow(pts), function(i){extract(b, cbind(pts$x[i],pts$y[i]), layer=pts$layerindex[i], nl=1)})

Soildf <- Soil

for (i in 3:nlayers(soil_rastbrick)) {

  out <- extract(soil_rastbrick[[i]], SpatialPolygons(wrld_simpl@polygons)) # Extract raster data from each of the country polygons in wrld_simpl
  Soil <- data.frame(Country=wrld_simpl$NAME, mean=unlist(lapply(out, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA )))
  names(Soil)[2] <- names(soil_rastbrick[[i]])

Soildf <- left_join(Soildf, Soil, by = 'Country')
}

# ============
# Evapotranspiration
# LandFluxEVAL.merged.89-05.yearly.all.nc

setwd('/Users/rubendouglas/OneDrive - The Royal Botanic Gardens, Kew/Nutridiv R/Ruben environmental data/')
LandFluxEVAL <- raster('LandFluxEVAL.merged.89-05.yearly.all.nc')

# Country Extract
out <- raster::extract(LandFluxEVAL, SpatialPolygons(wrld_simpl@polygons)) # Extract raster data from each of the country polygons in wrld_simpl
LandFluxEVAL_Country <- data.frame(Country=wrld_simpl$NAME, mean=unlist(lapply(out, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))) # unlist the data for each country polygon in a dataframe with the mean value in each country
names(LandFluxEVAL_Country)[2] <- 'LandFLUXEval'

# =========================================

# Combine topography, Soil & Evapotranspiration data and perform PCA to select key variables

dfs <- list(
  GloElev_Country,
  Soildf,
  LandFluxEVAL_Country
)

func <- function(...){
  df1 = list(...)[[1]]
  df2 = list(...)[[2]]
  col1 = colnames(df1)[1]
  col2 = colnames(df2)[1]
  xxx = left_join(..., by = setNames(col2,col1))
  return(xxx)
}

TopoSoilEva_Country <- Reduce( func, dfs)

# And subset to just focal countries

TopoSoilEva_Country <- TopoSoilEva_Country[which(TopoSoilEva_Country$Country %in% SES.FDPD$Country),]

TopoSoilEva_Country <- TopoSoilEva_Country[complete.cases(TopoSoilEva_Country),]

# And run the PCA, with nipas & cross validation 

TSE <- TopoSoilEva_Country[2:17] %>% as_data_frame()
rownames(TSE) <- TopoSoilEva_Country$Country

TopoSoilEva_PCAmeth <- pca(TSE, scale = "uv", center = TRUE, nPcs = 10, method = "nipals",cv="q2")

slplot(TopoSoilEva_PCAmeth, scoresLoadings =c(FALSE,TRUE))

plot(TopoSoilEva_PCAmeth)

TopoSoilEva_PCAmeth@R2cum

TopoSoilEva_PCA_loadings <- TopoSoilEva_PCAmeth@loadings %>% as.data.frame()

library(PerformanceAnalytics)
chart.Correlation(TSE, histogram=TRUE, pch=19)

# Variables with greatest positive & negative contribution to first 4 PCs are:

TSE_2 <- TSE[c(1,2,5,8,10,14,16)]
chart.Correlation(TSE_2, histogram=TRUE, pch=19)

# ====================================================================================================#
# Constrained analysis to assess correlations between Enironmental variables and Nutritional FD SESs

# Because directinality is not assured, eg. Nutritional FD indices may be both response and predictor variables for environmental conditions
# then I have choice between canonical correlation analysis or multivariate multiple regression. Other constrained analyses,
# eg. redundancy analysis depends on the assurance of directionality. 

# Canonical correlation analysis
# https://stats.idre.ucla.edu/r/dae/canonical-correlation-analysis/

# Canonical correlation analysis (CCorA) is suitable when you wish to examine linear relationships between two data sets where it is unclear what are response
# and what are explanatory variables. It attempts to find axes of maximum linear correlation between two corresponding data matrices.
# As it treats all variables equally, asserting no causal structure, it is a symmetrical canonical analysis.

# can put all raw variables in

# Multivatiate normal distribution assumptions are required for both sets of variables.
# Canonical correlation analysis is not recommended for small samples.

# ====================================================================================================#

# Multivariate Multiple Regression
# https://data.library.virginia.edu/getting-started-with-multivariate-multiple-regression/

# Needed to demonstrate relatinships of Nutritional FD variables with specific environmental variables


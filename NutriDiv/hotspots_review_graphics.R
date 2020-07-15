library(tidyverse)
library(tidyselect)
library(tidyr)
library(ggplot2)
library(ggdendro)
library(dendextend)

dend <- iris[1:30,-5] %>% scale %>% dist %>% 
  hclust %>% as.dendrogram %>%
  dendextend::set("branches_k_col", k=3) %>% dendextend::set("branches_lwd", 1.2) %>%
  dendextend::set("labels_colors") %>% dendextend::set("labels_cex", c(.9,1.2)) %>% 
  #dendextend::set("leaves_pch", 19) %>% dendextend::set("leaves_col", "blue") %>% 
  dendextend::set("hang_leaves", 0.2) 

# plot the dend in usual "base" plotting engine:
plot(dend)

# Create a radial plot and remove labels
ggplot(dend, labels = FALSE) + 
  scale_y_reverse(expand = c(0.2, 0)) +
  coord_polar(theta="x")

ggplot(dend, labels = FALSE)+ 
  scale_y_reverse(expand = c(0.2, 0))

# ================================================================================
# Overlay Crop species richness raster with coastline shapefile using levelplot
# ================================================================================
shape_path <- "/Users/rubendouglas/OneDrive - The Royal Botanic Gardens, Kew/Nutridiv R/"
coast_shapefile <- paste(shape_path, "ne_50m_coastline/ne_50m_coastline.shp", sep="")

# find out kind of shapefile (lines vs. polygons)
layer <- ogrListLayers(coast_shapefile)
ogrInfo(coast_shapefile, layer=layer)

# read the shape file
coast_lines <- readOGR(coast_shapefile, layer=layer)

# levelplot

# Just california

#lat 30 - 42
#long -125 -114

ext2 <- extent(-125,-114,30,42)
cali <- crop(`rast nbsp`, ext2)

# just Asia Pacific

# 115 - 166 
# 29 - 55

ext3 <- extent(115,150,29,55)
asiap <- crop(`rast nbsp`, ext3)


# levelplot
mapTheme2 <- rasterTheme(region=rev(brewer.pal(8,"RdBu")))

pltcali <- levelplot(cali, margin=TRUE, par.settings=mapTheme2, main = expression(paste("Crop Species Richness (per 100km"^"2)")))
pltcali + latticeExtra::layer(sp.lines(coast_lines, col="black", lwd=1))

pltasiap <- levelplot(asiap, margin=TRUE, par.settings=mapTheme2, main = expression(paste("Crop Species Richness (per 100km"^"2)")))
pltasiap + latticeExtra::layer(sp.lines(coast_lines, col="black", lwd=1))

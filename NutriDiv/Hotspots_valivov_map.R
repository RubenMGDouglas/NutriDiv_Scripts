# ====================================================================================================#
# ====================================================================================================#
# Hotspots & Valivov centres map
# ====================================================================================================#

library(sp)
library(sf)
library(rgdal)
library(oce)
library(ggplot2)

w <- spTransform(wrld_simpl, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))

hotspots<-readOGR("OneDrive - The Royal Botanic Gardens, Kew/Written Documents, Correspondences & Reports/Agrobiodiversity Hotspots Review/James_shapefiles_hotspots_2016_1/hotspots_2016_1.dbf")

hotspots <- spTransform(hotspots, CRS = prj)

plot(hotspots)

hotspots_land <- hotspots[hotspots$Type == "hotspot area",]

plot(hotspots_land)

# Still issue with horizontal lines

# ==================================================
# https://gis.stackexchange.com/questions/151601/lines-on-reprojected-sp-objects-with-mollweide-projection

# The horizontal lines occur because some of the shapes exceed the 180 degree limit

library(raster)
library(maps)
library(ggmap)
library(mapdata)
library(ggplot2)
library(tmap)
library(tmaptools)

hotspotmap <- SpatialPolygons2map(hotspots)

worldmapLines <- map2SpatialLines(hotspotmap,
                                  proj4string=CRS("+proj=longlat +datum=WGS84"))

w1 <- crop( worldmapLines, extent(-180, 180,-90,90))
x <- spTransform(w1, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))

hotspots<-readOGR("OneDrive - The Royal Botanic Gardens, Kew/Written Documents, Correspondences & Reports/Agrobiodiversity Hotspots Review/James_shapefiles_hotspots_2016_1/hotspots_2016_1.dbf")

w2 <- crop( hotspots, extent(-180, 180,-90,90))
x2 <- spTransform(w2, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))

w <- st_make_valid(w)
x2 <- st_make_valid(x2)

hotspots_land <- x2[x2$Type == "hotspot area",]
hotspots_land <- st_make_valid(hotspots_land)

hotspots_lim <- x2[x2$Type == "outer limit",]
hotspots_lim <- st_make_valid(hotspots_lim)


palette_explorer()

dev.off()

  # Remove Polynesia-Micronesia name
  hotspots_land <- hotspots_land[-26,]
  
  tm_shape(w) +
  tm_fill() + 
  #tm_layout(bg.color="lightblue") +
  tm_shape(x2) +  
  #tm_shape(hotspots_land) +
  tm_fill(col = "darkorange", alpha = 0.4) +
  tm_borders(col = "darkorange", lwd = 0.3)
  #tm_polygons(col = "NAME", palette = "Set3", legend.show = FALSE, lwd = 0.1) +
  #tm_text("NAME", size = 0.6, col = "black", auto.placement = TRUE) 
  
  # Grids & graticules in Tmap
  #https://geocompr.github.io/post/2019/tmap-grid/
  
  tm_shape(hotspots_land) +
    tm_graticules() 
    

# 1-8 (primary centers) could be drawn with plain lines and dashed lines for 2A 2B and 7A (secondary centers).

# Regarding colours. I would make it very simple. Biodiversity hotspots in one colour and Vavilov centers in another 
# (playing with transparency maybe so that we see well the overlap between the 2).


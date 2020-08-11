# ====================================================================================================#
# ====================================================================================================#
# Null modelling - Mapping Standardized values
# ====================================================================================================#
library(sp)
library(sf)
library(rgdal)
library(tmap)

# Map the Standardized PD & FD indices using tmap

# First add a dataframe containing the standardized indices to the country spatial polygons file

SESwrld <- merge(wrld_simpl, SES.FDPD, by.x="NAME", by.y="Country")


SESwrld_Moll <- spTransform(SESwrld, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))


# SF plotting
#SESwrld_sf <- st_as_sf(SESwrld)
SESwrld_Moll <- st_as_sf(SESwrld_Moll)

#SESwrld_sf <- st_make_valid(SESwrld_sf)
SESwrld_Moll <- st_make_valid(SESwrld_Moll)
st_is_valid(SESwrld_sf)

tmaptools::palette_explorer()

dev.off()

#change namess
#map.nbsp <- tm_shape(FDPDwrld_sf) + tm_polygons("nbsp", title="Species Richness", palette = 'Spectral')
map.SES.FRic <- tm_shape(SESwrld_Moll) + tm_polygons("FRic", title="FRic SES", palette = 'Spectral') +
  tm_layout(title = "(a)", title.size = 2)
map.SES.FDiv <- tm_shape(SESwrld_Moll) + tm_polygons("FDiv", title="FDiv SES", palette = 'Spectral') +
  tm_layout(title = "(b)", title.size = 2)
map.SES.FDis <- tm_shape(SESwrld_Moll) + tm_polygons("FDis", title="FDis SES", palette = 'Spectral') +
  tm_layout(title = "(c)", title.size = 2)
map.SES.FEve <- tm_shape(SESwrld_Moll) + tm_polygons("FEve", title="FEve SES", palette = 'Spectral') +
  tm_layout(title = "(d)", title.size = 2)

tmap_arrange(#map.nbsp,
  map.SES.FRic, map.SES.FDiv, map.SES.FDis, map.SES.FEve)

# And phylogenetic SES's

map.SES.Faith <- tm_shape(SESwrld_Moll) + tm_polygons("Phylogenetic_Faith", title="PD Faith SES", palette = 'Spectral')+
  tm_layout(title = "(e)", title.size = 2)
map.SES.MPD <- tm_shape(SESwrld_Moll) + tm_polygons("Phylogenetic_MPD", title="PD MPD SES", palette = 'Spectral')+
  tm_layout(title = "(f)", title.size = 2)
map.SES.PSE <- tm_shape(SESwrld_Moll) + tm_polygons("Phylogenetic_sp_evennes", title="PD PSE SES", palette = 'Spectral')+
  tm_layout(title = "(g)", title.size = 2)
map.SES.WeightFaith <- tm_shape(SESwrld_Moll) + tm_polygons("Phylogenetic_Weighted_Faith", title="PD WeightFaith SES", palette = 'Spectral')+
  tm_layout(title = "(h)", title.size = 2)

tmap_arrange(map.SES.Faith, map.SES.MPD, map.SES.WeightFaith, map.SES.PSE)

tmap_arrange(map.SES.FRic, map.SES.FDiv, map.SES.FDis, map.SES.FEve,
             map.SES.Faith, map.SES.MPD, map.SES.WeightFaith, map.SES.PSE)


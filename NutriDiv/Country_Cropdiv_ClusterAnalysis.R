# ====================================================================================================#
# ====================================================================================================#
# Cluster analysis of countries by taxonomic, phylogenetic & functional diversity
# ====================================================================================================#

# https://www-pnas-org.ezproxy.nhm.ac.uk/content/116/52/26465#sec-15
# 
# For the identification of groups of countries, we carried out a PCA [covariance (n − 1)] based on the matrix of countries (observations) 
# and indicators (variables). We followed the Kaiser criterion (eigenvalue >1) to determine the significant number of components. 
# To identify possible groups of countries with similar values of FSv indicators, we carried out an HCA based on the Euclidean distance
# (percentage of distance similarity at a 95% level of confidence) and Ward’s agglomerative method (71) using the standardized coordinates 
# of the most significant factors of the PCA. 


# =========================================================================================
# Prepare idiversity indices data
# =========================================================================================
# Copy country x diversity matrix
countryxdiv <- FDPDindices

#add country rownames
rownames(countryxdiv) <- countryxdiv$Country

#remove Djibouti
countryxdiv <- countryxdiv[-44,]

# CWMs
countryxCWM <- countryxdiv[,9:40]
countryxCWMScaled <- scale(countryxCWM)

# Remove unwanted index cols 
countryxdiv <- countryxdiv[,-c(2,4,9:41,43:48,50:51,53:54)]
sapply(countryxdiv, class)

# First it is neccessary to standardise the variables, as they are measured on different scales
countryxdivScaled <- scale(countryxdiv)

# combine & standardize diversity indices and CWMs
divandCWMScaled <- cbind(countryxdiv,countryxCWM) %>% scale() %>% as.data.frame()
sapply(divandCWMScaled, class)

divandCWMS <- cbind(countryxdiv,countryxCWM)

# =========================================================================================
# HCPC - Hierarchical Clustering on Principal Components
# =========================================================================================
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/117-hcpc-hierarchical-clustering-on-principal-components-essentials/#why-hcpc

library(FactoMineR)
library(factoextra)

# Diversity indices excluding CWMS
DivPCA <- PCA(countryxdiv, scale.unit = TRUE, graph = TRUE)

# All Diversity indices & CWMS
AllDivPCA <- PCA(divandCWMS, scale.unit = TRUE, graph = TRUE)

# Hierarchical clustering on principal components:
Div.hcpc <- HCPC(countryxdiv, graph = FALSE) # no CWMs
AllDiv.hcpc <- HCPC(AllDivPCA, graph = FALSE) # with CWMs


fviz_dend(AllDiv.hcpc, 
          cex = 0.4,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)

fviz_cluster(AllDiv.hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)

# Suggests 3 clusters

# ========================================

# Get dataframe with assigned cluster with $data.clust
AllDiv.hcpc$data.clust

# Descriptive statistics tables for each variable in the clusters
AllDiv.hcpc$desc.var$quanti

# representative observations of each cluster can be extracted as follow:
AllDiv.hcpc$desc.ind$para

# ========================================

# Plot clusters using tmap

library(sp)
library(sf)
library(tmap)

countryclust <- AllDiv.hcpc$data.clust %>% as.data.frame()
countryclust$Country <- rownames(countryclust)

clusters <- countryclust[,42:43]

# need to merge cluster results column to FDPDwrld_sf

# Adding dataframe data to Spatialpolygonsdataframe
clusterwrld_sf <- merge(FDPDwrld_sf,clusters, by.x="NAME", by.y="Country")

clusterwrld_sf <- st_make_valid(clusterwrld_sf)
st_is_valid(clusterwrld_sf)

# and plot
tmaptools::palette_explorer()
tm_shape(clusterwrld_sf) + tm_polygons("clust", title="Cluster analysis of crop TD, FD & PD")

# ========================================
# Boxplots showing the value of each variable for each of the three clusters
library(ggplot2)
library(tidyverse)
library(plyr)

# clust on x axis
# value on y axis
# facet by variable

dfmelt <- melt(countryclust)

ggplot(dfmelt, aes(x=clust, y=log10(as.numeric(value)),fill=clust))+
  geom_boxplot()+
  facet_wrap(.~variable, ncol=5)+
  labs(x="Clusters")+
  theme(axis.text.x=element_text(angle=-90, vjust=0.4,hjust=1))

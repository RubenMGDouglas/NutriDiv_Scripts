# ====================================================================================================#
# ====================================================================================================#
# Calculate Gridcell FD indices with dbfd() from FD package
# ====================================================================================================#

# Paper on FD package
# http://adn.biol.umontreal.ca/~numericalecology/Reprints/Laliberte_Legendre_Ecology_2010.pdf
library(tidyverse)
library(FD)

# Gower dissimilarity matrix mostly is able to tolerate missing values
# NA values may result in the PCoA returning negative eigenvalues, which need to be corrected to avoid biased FD values.
# This can arise from missing values & sometime when the distance matrix cannot be represented in Euclidean space e.g. Gower

# FRic & FDiv can only be calculated where there are more species than traits in an assemblage. 

# FRic cannot integrate information on relative abundances. Consequently, rare species with extreme trait values will greatly inflate FRic.
# This may or may not be a desirable property, depending on the application.

# FDis, like FRic, captures the spread, or dispersion of of the N species in T-th dimentional traitspace,
# but is abundance weighted and less sensitive to outliers than FRic (which is the multivariate equivilent of the range).

# FDis calculates the dispersion of Species in T-th dimentional space as the mean distance of species to the centroid,
# but weights the distance of the centroid by relative abundance of species, 
#such that the centroid is closer to more abundant species. 

# FDiv is almost perfectly correlated to Raos Q (r = 0.966)

# Another addition is functional group richness (FGR), which is computed from an a posteriori classification of species based on their
# functional traits (i.e., the ‘‘data-defined’’ approach of Gitay and Noble [1997]). 
# NB This requres user input to select the number of clusters/ distance threshold for the Hierarchical clustering algorithm,
# so makes sense to only look at FGR for the Country-level analysis, where its reasonable to manually select number of a postereori clusters manually. 

# FRic and FDiv both rely on finding the minimum convex hull that includes all species (Villéger et al. 2008). 
# This requires more species than traits. To circumvent this problem, dbFD takes only a subset of the PCoA axes as traits via argument m.
# This, however, comes at a cost of loss of information. The quality of the resulting reduced-space representation is returned by
# qual.FRic, which is computed as described by Legendre and Legendre (1998) and can be interpreted as a R^2-like ratio.

# CWM can be calcualted from qualitative traits under this framework

# ====================================================================================================#

# Paper on Fric, Feve & FDiv
# file:///Users/rubendouglas/Downloads/Villegeretal.2008Ecology.pdf

#Functional evenness describes the evenness of abundance distribution in a functional trait space (Mason et al. 2005).
# Basically, the new index quantifies the regularity with which the functional space is filled by species, weighted by their abundance.
# FEve decreases either when abundance is less evenly distributed among species or when functional distances among species are less regular

# ====================================================================================================#

# Branch length FD (Petchey & Gaston 2006)
library(vegan)

# Function treedive() finds the treeheight for each site (row) of a community matrix.
# The function uses a subset of dendrogram for those species that occur in each site,
# and excludes the tree root if that is not needed to connect the species (Petchey and Gaston 2006).
# The subset of the dendrogram is found by first calculating cophenetic distances from the input dendrogram,
# then reconstructing the dendrogram for the subset of the cophenetic distance matrix for species occurring in each site.
# Diversity is 0 for one species, and NA for empty communities.

# Needs a 

# ====================================================================================================#
# ====================================================================================================#

# input files

# Crop x Nutrient matrices:
Nutrients5 # With NAs
Nona_Nutrients5 # No NAs

# Gridcell x Crop productin matrix
Production_mat_25min

# Gridcell x Crop presence absence matrix
Presabs_mat_25min

common_3 <- intersect(intersect(colnames(Production_mat_25min),colnames(Presabs_mat_25min)),Nutrients5$`Crop name`)
nona_common_3 <- intersect(intersect(colnames(Production_mat_25min),colnames(Presabs_mat_25min)),Nona_Nutrients5$`Crop name`)
#109, meaning all unique names are harmonized

# Just ensure all are the same (for missing values included)
Nutrients5 <- Nutrients5[Nutrients5$`Crop name` %in% common_3,]
# Production matrix
Production_mat_25min <- cbind(Production_mat_25min[,c(1:2)],Production_mat_25min[,colnames(Production_mat_25min) %in% common_3])
# presence absence matrix
Presabs_mat_25min <- cbind(Presabs_mat_25min[,c(1:2)],Presabs_mat_25min[,colnames(Presabs_mat_25min) %in% common_3])

# And for no NAs
Nona_Nutrients5 <- Nona_Nutrients5[Nona_Nutrients5$`Crop name` %in% nona_common_3,]
# Production matrix
Nona_Production_mat_25min <- cbind(Production_mat_25min[,c(1:2)],Production_mat_25min[,colnames(Production_mat_25min) %in% nona_common_3])
# presence absence matrix
Nona_Presabs_mat_25min <- cbind(Presabs_mat_25min[,c(1:2)],Presabs_mat_25min[,colnames(Presabs_mat_25min) %in% nona_common_3])



# ================================================================================
# (1) For each row in Pres /abs matrix, create a subsetted nutrient matrix 
# containing only species present in the focal gridcell
# ================================================================================
# ================================================================================

#preallocate a dataframe with primed column names for functional diversity indices
fd_indices <- data.frame(matrix(nrow = 0, ncol = 37))
colnames(fd_indices) <- colnames(df_fdres)


for (i in 1:nrow(Presabs_mat_25min)){
  
  # create a table from the current grid cell
  curr.cell <- as.data.frame(Nona_Presabs_mat_25min[i,])
  
  #transpose current cell x species matrix
  t.curr.cell <- t(curr.cell)
  
  # set column name
  colnames(t.curr.cell) <- "val"
  
  # create a subset cell x species matrix containing only present species
  cell_pres <- as.data.frame(t.curr.cell[which(t.curr.cell == 1),])
  
  curr_xy <- curr.cell[,1:2]
  
  length(rownames(cell_pres))
  
  if (length(rownames(cell_pres)) == 0){
    curr_zeros <- cbind(curr_xy, data.frame(matrix(NA, nrow=1, ncol=37)))
    colnames(curr_zeros) <- fdcolnames2
    # compute columns which are present in both the master dataframe and the current gridcell dataframe
    common_cols <- intersect(colnames(fd_indices), colnames(curr_zeros))
    # rbind the current gridcell with the master dataframe by shared columns
    fd_indices <- rbind(fd_indices[, common_cols], curr_zeros[, common_cols])
    #fd_indices <- rbind.fill(fd_indices, curr_zeros)
    
    next
    
  } else if (length(rownames(cell_pres)) == 1){
    curr_zeros <- cbind(curr_xy, data.frame(matrix(NA, nrow=1, ncol=37)))
    colnames(curr_zeros) <- fdcolnames2
    # compute columns which are present in both the master dataframe and the current gridcell dataframe
    common_cols <- intersect(colnames(fd_indices), colnames(curr_zeros))
    
    # rbind the current gridcell with the master dataframe by shared columns
    fd_indices <- rbind(fd_indices[, common_cols], curr_zeros[, common_cols])
    #fd_indices <- rbind.fill(fd_indices, curr_zeros)
    
    next
  }
  
  # convert crop names to multi-level factors for compatability with nutrient crop name values:
  selected <- as.factor(rownames(cell_pres))
  selected
  
  # use %in% to create a subsetted species x nutrient table for our sample grid cell: 
  cell_nutri <- Nona_Nutrients5[which(Nona_Nutrients5$`Crop name` %in% selected),]
  
  # make species labels in cell_nutri and cell_prod identical
  #tryCatch({
  
  cell_nutri <- cell_nutri %>%
    arrange(cell_nutri$`Crop name`)
  
  rownames(cell_nutri) <- selected
  
  
  #}, error = function(e){})
  
  # and subset so only continous variables are present
  cell_nutri <- cell_nutri[,c(3:29)]
  
  # ================================================================================
  # (2) combine that with a subsetted production matrix with only the species in the
  # subset nutrient matrix - called 'cell_prod'
  # ================================================================================
  
  #TEST VARIABLES
  #curr.cell <- Presabs_mat_25min[937,]
  #cell_pres <- cbind(curr.cell[,1:2],curr.cell[,curr.cell == 1])
  #selected <- as.factor(colnames(cell_pres[,3:NCOL(cell_pres)]))
  #TEST VARIABLES
  
  curr.production <- merge(Nona_Production_mat_25min, curr.cell[c('x','y')])
  
  cell_prod <- curr.production[,colnames(curr.production) %in% selected]
  
  dim <- dim(cell_prod)
  dimsum <- sum(dim)
  dimsum
  
  if (dimsum == 1){
    curr_zeros <- cbind(curr_xy, data.frame(matrix(NA, nrow=1, ncol=37)))
    colnames(curr_zeros) <- fdcolnames2
    # compute columns which are present in both the master dataframe and the current gridcell dataframe
    common_cols <- intersect(colnames(fd_indices), colnames(curr_zeros))
    
    # rbind the current gridcell with the master dataframe by shared columns
    fd_indices <- rbind(fd_indices[, common_cols], curr_zeros[, common_cols])
    #fd_indices <- rbind.fill(fd_indices, curr_zeros)
    next
  }
  
  # ================================================================================
  # (3) return a dataframe with all FD indices calculated using dbFD, for the 
  # grid-cell specific species x nutrient table
  # ================================================================================
  
  # NB: for dbFD function, the species labels in a and x must be identical and in the same order.
  
  tryCatch({
    # because NAs are present in the trait matrix, a distance matrix cannot be represented in Euclidean space,
    # a correction must be aplied to avoid biased FD estimates. Here sqrt transformation didn't work, so the
    # approach described by Cailliez (1983) is used:
    
    #fdres <- dbFD(cell_nutri,cell_prod)
    
    fdres <- dbFD(cell_nutri, cell_prod, w.abun = TRUE, calc.CWM = TRUE, calc.FDiv = TRUE)
    fdres
  }, error = function(e){})
  
  # bind fd indices results with xy coordinates
  df_fdres <- cbind(curr_xy,as.data.frame(fdres))
  
  # compute columns which are present in both the master dataframe and the current gridcell dataframe
  common_cols <- intersect(colnames(fd_indices), colnames(df_fdres))
  
  # rbind the current gridcell with the master dataframe by shared columns
  fd_indices <- rbind(fd_indices[, common_cols], df_fdres[, common_cols])
  #fd_indices <- rbind.fill(fd_indices, df_fdres)
  
} #End of loop
# ================================================================================
# ================================================================================












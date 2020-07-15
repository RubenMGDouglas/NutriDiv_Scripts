# ====================================================================================================# 
# FD influence 
# ====================================================================================================# 
# This function is adapted from Woods (2017) FD.influence() function
# https://github.com/swood-ecology/kdg_nutrient_diversity/blob/master/fd_influence.R
# to first calculate and return the FD for a community (sensu Petchey & Gaston 2002, 2006)
# ie. total dendrogram branch length, and then an influence statistic for each species in the community,
# which is the change in branch length when that species is removed stepwise from the community dendrogram. 
# ====================================================================================================# 
FDinfluence <- function(traitfile, dis.metric = "euclidean", clus.method = "average"){
  
  # Create 'communities' where all species are present but one
  sp = rownames(traitfile)
  com <- melt(sapply(sp, function(x) sp[sp!=x]))[,c(3,2)]
  names(com) <- c("species", "community")
  
  # compute distance matrix for 'global community' in traitfile
  distances = dist(traitfile, method = dis.metric)
  
  # if there is more than one NA in the distance matrix, return the warning "incomparable species (check trait file)"
  if(length(distances[is.na(distances)]) > 1) { NA; warning("incomparable species (check trait file)")}
  
  # else, (1) compute a hierarchical cluster analysis (dendrogram) on the distance matrix 
  else { tree = hclust(distances, method = clus.method)
  
  # (2) calculate FD (total dendrogram branch length) for subset communities using vegan's treeheight function
  FD <- tapply(com$species, com$community,
               function(y) treeheight(hclust(dist(traitfile[y,], method = dis.metric),method = clus.method)
               )
  )
  
  # Calculate FD of 'global community' dendrogram
  full <- treeheight(tree)
  
  # calculate absolute or relative difference between full dendrogram and subset dendrogram
  absdiff <- full - FD
  reldiff <- 100*(full - FD ) / full
  }
  
  reldiff[rownames(traitfile)]
  
}

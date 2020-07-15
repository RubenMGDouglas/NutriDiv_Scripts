# ====================================================================================================#
# ====================================================================================================#
# Calculate Country level FD & PD indices 
# ====================================================================================================#
library(tidyverse)
library(ape)
library(FD)
library(adiv)
library(vegan)

# Read in plain text Phylogeny 
Aoife_tree <- read.tree('/Users/rubendouglas/OneDrive - The Royal Botanic Gardens, Kew/Nutridiv R/Aoifes Phylogeny/crop_tree_Ruben_common_v1')
plot(Aoife_tree, type = "fan")

# ==================================================================
# Harmonize crop names from phylogeny with the spatial crop matrices

tips <- Aoife_tree$tip.label %>% sort.default 
int <- intersect(tips,common_3) # only 88 intersect.. 

setdiff(tips,common_3) 
#[1] "aniseed"   "chille"    "cucumber"  "melon"     "peach"     "pumpkin"   "tangerine"

setdiff(common_3, tips) 
# [1] "aniseetc"       "banana"         "cauliflower"    "chilleetc"      "cucumberetc"    "fonio"          "garlic"         "grapefruitetc" 
# [9] "greenbroadbean" "groundnut"      "melonetc"       "mushroom"       "onion"          "peachetc"       "peppermint"     "plantain"      
# [17] "pumpkinetc"     "strawberry"     "stringbean"     "tangetc"        "triticale" 

# Rename the non-overlapping 7 crop tips so they correspond to those in the Crop PA, production & nutrient matrices
for (i in 1:length(Aoife_tree$tip.label)){
  if (Aoife_tree$tip.label[i]=="aniseed"){Aoife_tree$tip.label[i] <- "aniseetc"}
  else if (Aoife_tree$tip.label[i]=="chille"){Aoife_tree$tip.label[i] <- "chilleetc"}
  else if (Aoife_tree$tip.label[i]=="cucumber"){Aoife_tree$tip.label[i] <- "cucumberetc"}
  else if (Aoife_tree$tip.label[i]=="melon"){Aoife_tree$tip.label[i] <- "melonetc"}
  else if (Aoife_tree$tip.label[i]=="peach"){Aoife_tree$tip.label[i] <- "peachetc"}
  else if (Aoife_tree$tip.label[i]=="pumpkin"){Aoife_tree$tip.label[i] <- "pumpkinetc"}
  else if (Aoife_tree$tip.label[i]=="tangerine"){Aoife_tree$tip.label[i] <- "tangetc"}
}

tips <- Aoife_tree$tip.label %>% sort.default 
plot(Aoife_tree, type = "fan")
common_phy <- intersect(tips,common_3) # now all 95 tips intersect with common_3

# ========================================================================
# Subset the crop matrices to only have columns for crops in the phylogeny

prod_PDFD_country <- Prod_country 

prod_PDFD_country <- cbind(prod_PDFD_country$Country,prod_PDFD_country[,which(colnames(prod_PDFD_country) %in% common_phy)])
colnames(prod_PDFD_country)[1] <- "Country"

# =================================================================

# Remove outliers in production!

boxplot(prod_PDFD_country)
boxplot(prod_PDFD_country[,20:30])

# chicory has many extremely outlying values. 
prod_PDFD_country <- prod_PDFD_country[,-25]


# ================================================================================
# (1) For each row in Production matrix, create a subsetted nutrient matrix 
# containing only species present in the focal gridcell
# ================================================================================
# ================================================================================

#preallocate a dataframe with primed column names for functional diversity indices
fd_indices <- data.frame(matrix(nrow = 0, ncol = 36))
colnames(fd_indices) <- colnames(df_fdres)

for (i in 1:nrow(prod_PDFD_country)){
  
  # create a table from the current country
  curr.cell <- as.data.frame(prod_PDFD_country[3,])
  
  #country name
  curr_country <- as.data.frame(curr.cell[1,1])
  colnames(curr_country)[1] <- "Country"
  
  #subset all crops with production greater than 0
  cell_prod <- curr.cell[,2:ncol(curr.cell)]
  cell_prod <- cell_prod[,which(cell_prod > 0)]
  
  length(colnames(cell_prod))
  
  if (length(colnames(cell_prod)) == 0){
    curr_zeros <- cbind(curr_country, data.frame(matrix(NA, nrow=1, ncol=36)))
    colnames(curr_zeros) <- colnames(df_fdres)
    # compute columns which are present in both the master dataframe and the current gridcell dataframe
    common_cols <- intersect(colnames(fd_indices), colnames(curr_zeros))
    # rbind the current gridcell with the master dataframe by shared columns
    fd_indices <- rbind(fd_indices[, common_cols], curr_zeros[, common_cols])
    #fd_indices <- rbind.fill(fd_indices, curr_zeros)
    
    next
    
  } else if (length(colnames(cell_prod)) == 1){
    curr_zeros <- cbind(curr_country, data.frame(matrix(NA, nrow=1, ncol=36)))
    colnames(curr_zeros) <- colnames(df_fdres)
    # compute columns which are present in both the master dataframe and the current gridcell dataframe
    common_cols <- intersect(colnames(fd_indices), colnames(curr_zeros))
    
    # rbind the current gridcell with the master dataframe by shared columns
    fd_indices <- rbind(fd_indices[, common_cols], curr_zeros[, common_cols])
    #fd_indices <- rbind.fill(fd_indices, curr_zeros)
    
    next
  }
  
  # convert crop names to multi-level factors for compatability with nutrient crop name values:
  selected <- as.factor(colnames(cell_prod))
  selected
  
  # use %in% to create a subsetted species x nutrient table for our sample grid cell: 
  cell_nutri <- Nutrients5[which(Nutrients5$`Crop name` %in% selected),]
  
  # make species labels in cell_nutri and cell_prod identical
  #tryCatch({
  cell_nutri <- cell_nutri %>%
    arrange(cell_nutri$`Crop name`)
  
  rownames(cell_nutri) <- selected
  #}, error = function(e){})
  
  # and subset so only continous variables are present
  cell_nutri <- cell_nutri[,c(3:29)]
  
  # ================================================================================
  # (3) return a dataframe with all FD indices calculated using dbFD, for the 
  # grid-cell specific species x nutrient table
  # ================================================================================
  
  # NB: for dbFD function, the species labels in a and x must be identical and in the same order.
  
  tryCatch({
    # - FDis, Rao's Q, FEve, FDiv, and CWM are weighted by the relative abundances of species, so diversity is not inflated by superabundant crops
    # - All traits will be standardised to mean 0 and unit variance to allow relative comparison of FD between countries
    # - Because NAs are present in the trait matrix, a distance matrix cannot be represented in Euclidean space,
    # a correction must be aplied to avoid biased FD estimates. Here sqrt transformation didn't work, so the
    # approach described by Cailliez (1983) is used.
    # - FRic values were extremely low in previous iterations, so I will standardize FRic between 0 - 1 to allow comparison of relative values across countries
    
    fdres <- dbFD(cell_nutri, cell_prod, w.abun = TRUE, stand.x = TRUE, corr = "cailliez",stand.FRic = TRUE, calc.CWM = TRUE, calc.FDiv = TRUE)
  }, error = function(e){})
  
  # bind fd indices results with curr_country
  df_fdres <- cbind(curr_country,as.data.frame(fdres))
  
  # compute columns which are present in both the master dataframe and the current gridcell dataframe
  common_cols <- intersect(colnames(fd_indices), colnames(df_fdres))
  
  # rbind the current gridcell with the master dataframe by shared columns
  fd_indices <- rbind(fd_indices[, common_cols], df_fdres[, common_cols])
  #fd_indices <- rbind.fill(fd_indices, df_fdres)

  
} #End of loop
# ================================================================================
# ================================================================================

#set rownames as country names
rownames(fd_indices) <- prod_PDFD_country$Country



# =================================================================
# (3) integrate FD & PD index calculations from adiv & dbFD

# dbFD: FGR 
# treedive: FD (total branch length) 
# adiv: PD & Uniqueness 

# =================================================================

comm <- prod_PDFD_country
rownames(comm) <- prod_PDFD_country$Country
comm <- comm[,2:ncol(comm)]

# remove country rows with all zeros
comm_noze <- comm[rowSums(comm)!=0, ]

colnames(comm_noze)

#and NA row
comm_noze <- comm_noze[!rownames(comm_noze)=='NA',]

# ================================================================================

# Try dbFD with whole matrices 

# Need to prepare:
# (1) Subset of Nutrients5 with just rows for common_phy, rownames are crop names and just continous traits

#first remove chicory
common_phy <- common_phy[-24]

Wholenutrients <- Nutrients5[Nutrients5$`Crop name` %in% common_phy,]
rownames(Wholenutrients) <- Wholenutrients$`Crop name`
Wholenutrients <- Wholenutrients[,-c(1:2)]

# order rownames alphabetically
Wholenutrients <- Wholenutrients[order(rownames(Wholenutrients)),]


# (2)  Copy of prod_PDFD with country for rownames rather than a column for country
Wholeproduction <- prod_PDFD_country

rownames(Wholeproduction) <- Wholeproduction$Country
Wholeproduction <- Wholeproduction[,-1]

# order column names alphabetically
Wholeproduction <- Wholeproduction[ , order(names(Wholeproduction))]

# remove rows with all Nas
Wholeproduction <- Wholeproduction[complete.cases(Wholeproduction),]

# and rows with zero-sum production
Wholeproduction <- Wholeproduction[as.logical(rowSums(Wholeproduction != 0)), ]

WholedbFD <- dbFD(Wholenutrients, Wholeproduction, w.abun = TRUE, stand.x = TRUE, corr = "cailliez", calc.FRic = FALSE, calc.CWM = TRUE, m = "max", messages = TRUE) 
WholedbFD <- WholedbFD %>% as.data.frame()

WholedbFD2 <- dbFD(Wholenutrients, Wholeproduction, w.abun = TRUE, stand.x = TRUE, corr = "cailliez", calc.FRic = TRUE, calc.FDiv = FALSE, calc.CWM = TRUE, m = "min", messages = TRUE) 
WholedbFD2 <- WholedbFD2 %>% as.data.frame()

WholedbFD3 <- dbFD(Wholenutrients, Wholeproduction, w.abun = TRUE, stand.x = TRUE, corr = "cailliez", calc.FRic = TRUE, calc.FDiv = TRUE, calc.CWM = TRUE, m = "min", messages = TRUE) 
WholedbFD3 <- WholedbFD3 %>% as.data.frame()

# Copy as fd_indices for downstream leftjoins
fd_indices <- WholedbFD3
crops <- rownames(fd_indices)
fd_indices$Country <- crops


# ================================================================================
# ===============================
# Package "adiv"
#
# Listed below all of the relevant Phylogenetic & functional diversity functions in adiv which take as input:
#
# 'Comm':
# a data frame or a matrix typically with communities as rows, species as columns
# and an index of abundance as entries. Species should be labeled as in the phylogenetic tree where they are the tips.
#
# some functions additionally take a phylogeny or distance matrix 
#

# ===============================
# EH
# Faith Phylogenetic Diversity

# phyl	
# 
# select	= a vector containing the numbers of the leaves (species) which must be considered in the computation of Phylogenetic Diversity (PD)
# (or merely sum of branch lengths on the tree). This argument allows the calculation of PD for a subset of species (including the branch between the subtree and the most ancient node of the full tree).
# 
 adiv:::EH(Aoife_tree)# Full tree [1] 5908.207
 adiv:::EH(Aoife_tree, as.character(selected))# Community subset [1] 3843.444
#   
# ===============================
#
# evodiv 
# Indices of Phylogenetic Diversity
#
# The function evodiv calculates diversity indices that rely on the relative or absolute abundance of
# features on a phylogenetic tree, with the assumption that the number of features on a given branch
# of a phylogenetic tree is equal to the length of this branch (see Pavoine 2016).
#
# Faith, Gini, Simpson, Shannon... indices (7 in total)
#
# Arguments
#
#Phyl
#
#Comm
#
  evodiv(Aoife_tree, cell_prod, method = "full")
  # richness GiniSimpson  Simpson  Shannon Margalef Menhinick  McIntosh
  # 246 3843.444   0.9962287 265.1592 6.035701 355.0322  17.16243 0.9427989
  
  evodiv <- as.data.frame(evodiv(Aoife_tree, Wholeproduction, method = "full"))
  evodiv$Country <- rownames(evodiv)
  
  evodiv <- evodiv[,c(8,1,2,3,4,5,6,7)]
  
FDPDindices <- left_join(fd_indices, evodiv, by = "Country")
  
# ===============================
# evodivparam 
#Parametric Indices of Phylogenetic Diversity
#
# Function evodivparam calculates phylogenetic diversity in communities using parametric indices
# derived from Tsallis and Hill compositional indices. It can also be applied to functional trees rather
# than phylogenies, to calculate a functional diversity. The function plot.evodivparam plots the
# results of function evodivparam.
#
# only able to calculate one index at a time # method = c("hill", "tsallis", "renyi")
#
# Arguments
#
#Phyl
#
#Comm
#
#
#

  compcomm <- comm[,whi]
  
  # q	a vector with nonnegative value(s) for parameter q of functions qfeveHCDT, qfeveHill, and qfeveRenyi.
  # q is the parameter that increases with the importance given to abundant species compared to rare species in diversity.
  # If only one value of q is given, a vector with the phylogenetic diversity of each community is returned. If more than one value of q is given, a list of two objects is returned:
  # q	the vector of values for q;
  # div a data frame with the phylogenetic diversity of each community calculated for all values of q.
  
  hillparam <- evodivparam(Aoife_tree, comm_noze, method = "hill") 
  tsallisparam <- evodivparam(Aoife_tree, comm_noze, method = "tsallis") 
  renyiparam <- evodivparam(Aoife_tree, comm_noze, method = "renyi") 
  
  plot(x)
  
  # data(batcomm)
  # phy <- read.tree(text=batcomm$tre)
  # ab <- batcomm$ab[, phy$tip.label]
  # plot(evodivparam(phy, ab))
  # plot(evodivparam(phy, ab, q=seq(0, 10, length=20)))
  
# ===============================
#
# evoeveparam 
# Parametric Indices of Phylogenetic evenness
# 
# Function evoeveparam calculates phylogenetic evenness (evenness in features, which are branch units of a phylogenetic tree) in communities.
# It uses parametric indices derived from Tsallis and Hill compositional indices, and named qfeveHCDT, qfeveHill, and qfeveRenyi in Pavoine and Ricotta (2019).
# evoeveparam can also be applied to functional trees rather than phylogenies, to calculate a functional evenness.
# The function plot.evoeveparam plots the results of function evoeveparam.
#
# method = c("hill", "tsallis", "renyi")
#
# Arguments
#
#Phyl
#
#Comm
#
  evoeveparam(Aoife_tree, comm, method = "hill") # Error in .local(x, ...) : Only the root can have NA as edge length. 
  evoeveparam(Aoife_tree, cell_prod, method = "tsallis") # 
  evoeveparam(Aoife_tree, cell_prod, method = "renyi") #  
  
# ===============================
# evouniparam 
# Parametric Indices of Phylogenetic Uniqueness
# 
# Function evouniparam calculates phylogenetic uniqueness in communities using parametric indices derived from Tsallis and Hill compositional indices. evouniparam can also be applied to functional trees rather than phylogenies, to calculate a functional uniqueness. The function plot.evouniparam
# plots the results of function evouniparam.
# 
# 
  
  ab <- cell_prod[, Aoife_tree$tip.label]
  
  evouniparam(Aoife_tree, cell_prod, method = "hill") 
  # Error in .local(x, ...) : Only the root can have NA as edge length. 
  
  
# # ===============================
# rlqESLTP
# Linking Patterns in Phylogeny, Traits, Abiotic Variables and Space
# 
# An extention of the RLQ approach to identify potential environmental filters in species’ traits in a
# phylogenetic and spatial context.
# 
# speciesdiv Indices of Species Diversity
# The function speciesdiv calculates diversity indices that rely on relative or absolute species abundance.
# 
# Arguments
# comm 
# 
# # ===============================
# specieseve 
# Indices of Species Evenvess
# The function specieseve calculates evenness indices that rely on relative or absolute species abundance.
# 
# Arguments
# 
# comm 
# 
  
  specieseve(cell_prod)
  # GiniSimpson    Simpson   Shannon       Heip  McIntosh SmithWilson
  #   0.6472781 0.04660295 0.3992672 0.07058301 0.4563412  0.06297021
  
# # ===============================
# treeUniqueness 
# Community-level phylogenetic (or functional) redundancy
# 
# The function treeUniqueness calculates community-level phylogenetic (or tree-based) redundancy
# taking into account the branching pattern of the underlying phylogenetic tree (or any other tree, like a functional dendrogram).
# 
# Arguments
# 
# phyl = phylo or dendrogram
# 
# comm
# 
  
  treeUniqueness(Aoife_tree, cell_prod, index = "richness")
  #     Dk       Dp Phylogenetic Uniqueness Phylogenetic Redundancy
  # 246 59 28.27888               0.4793031               0.5206969
  
  treeUniqueness(Aoife_tree, cell_prod, index = "GiniSimpson")
  
  treeUniqueness(Aoife_tree, cell_prod, index = "Shannon")
  
  treeUniqueness(PDFDdend, cell_prod, index = "richness")
  #     Dk      Dp   Functional Uniqueness   Functional Redundancy
  # 246 59 5.86559              0.09941677               0.9005832
  
  treeUniqueness(PDFDdend, cell_prod, index = "GiniSimpson")
  
  treeUniqueness(PDFDdend, cell_prod, index = "Shannon")
 
  # # ===============================
  
  FDPDindices <- FDPDindices[,-c(49:54)]
  
  # Phylogenetic uniqueness (Richness) 
  PhyUniqRic <- as.data.frame(treeUniqueness(Aoife_tree, Wholeproduction, index = "richness"))
  PhyUniqRic$Country <- rownames(PhyUniqRic)
  PhyUniqRic <- PhyUniqRic[,3:5]
  colnames(PhyUniqRic)[1:2] <- c("Phylogenetic Uniqueness_Richness", "Phylogenetic Redundancy_Richness")
  FDPDindices <- left_join(FDPDindices, PhyUniqRic, by = "Country")
  
  # Phylogenetic uniqueness (GiniSimpson) 
  PhyUniqGisi <- as.data.frame(treeUniqueness(Aoife_tree, Wholeproduction, index = "GiniSimpson"))
  PhyUniqGisi$Country <- rownames(PhyUniqGisi)
  PhyUniqGisi <- PhyUniqGisi[,3:5]
  colnames(PhyUniqGisi)[1:2] <- c("Phylogenetic Uniqueness_GiniSimpson", "Phylogenetic Redundancy_GiniSimpson")
  FDPDindices <- left_join(FDPDindices, PhyUniqGisi, by = "Country")
  
  # Phylogenetic uniqueness (Shannon)
  PhyUniqSha <- as.data.frame(treeUniqueness(Aoife_tree, Wholeproduction, index = "Shannon"))
  PhyUniqSha$Country <- rownames(PhyUniqSha)
  PhyUniqSha <- PhyUniqSha[,3:5]
  colnames(PhyUniqSha)[1:2] <- c("Phylogenetic Uniqueness_Shannon", "Phylogenetic Redundancy_Shannon")
  FDPDindices <- left_join(FDPDindices, PhyUniqSha, by = "Country")
  
  # ==========
  
  # Functional uniqueness (Richness) 
  FunUniqRic <- as.data.frame(treeUniqueness(PDFDdend, Wholeproduction, index = "richness"))
  FunUniqRic$Country <- rownames(FunUniqRic)
  FunUniqRic <- FunUniqRic[,3:5]
  colnames(FunUniqRic)[1:2] <- c("Functional Uniqueness_Richness", "Functional Redundancy_Richness")
  FDPDindices <- left_join(FDPDindices, FunUniqRic, by = "Country")
  
  # Functional uniqueness (GiniSimpson) 
  FunUniqGisi <- as.data.frame(treeUniqueness(PDFDdend, Wholeproduction, index = "GiniSimpson"))
  FunUniqGisi$Country <- rownames(FunUniqGisi)
  FunUniqGisi <- FunUniqGisi[,3:5]
  colnames(FunUniqGisi)[1:2] <- c("Functional Uniqueness_GiniSimpson", "Functional Redundancy_GiniSimpson")
  FDPDindices <- left_join(FDPDindices, FunUniqGisi, by = "Country")
  
  # Functional uniqueness (Shannon)
  FunUniqSha <- as.data.frame(treeUniqueness(PDFDdend, Wholeproduction, index = "Shannon"))
  FunUniqSha$Country <- rownames(FunUniqSha)
  FunUniqSha <- FunUniqSha[,3:5]
  colnames(FunUniqSha)[1:2] <- c("Functional Uniqueness_Shannon", "Functional Redundancy_Shannon")
  FDPDindices <- left_join(FDPDindices, FunUniqSha, by = "Country")
  
  # Remove all redundancy indices
  FDPDindices1 <- FDPDindices[,-c(50,52,54,56,58,60)]
  
  
  # # ===============================
  
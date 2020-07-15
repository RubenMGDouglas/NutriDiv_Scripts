# ====================================================================================================#
# ====================================================================================================#
# Null modelling - Calculating observed values
# ====================================================================================================#

library(vegan)
library(picante)
library(phytools)
library(tidyverse)
library(geiger)

# (1) FRic captures the total amount of variation  in trait values,
# making it conceptually analogous to PDFaith (Tucker et al. 2017).

my.sample <- Wholeproduction

observed.pd <- pd(my.sample, Aoife_tree, include.root = FALSE)
observed.pd <- observed.pd[,-2] %>% as.data.frame()
observed.pd$Country <- rownames(Wholeproduction)
names(observed.pd)[1] <- 'Phylogenetic_Faith'


# ====================================================================================================#

# (2) Functional Divergence is analogous to phylogenetic MPD

# picante package has a mpd() function to calculate mpd for all sites from a site by species matrix and phylo distance matrix.

pd.matrix <- cophenetic(Aoife_tree)

observed.mpd <- mpd(samp = my.sample, dis = pd.matrix, abundance.weighted = TRUE) %>% as.data.frame() # you can weight by sp abundance also.
observed.mpd$Country <- rownames(Wholeproduction)
names(observed.mpd)[1] <- 'Phylogenetic_MPD'

#join them into a dataframe for all observed & SES columns
FDPDses <- left_join(observed.pd, observed.mpd, by = "Country")
FDPDses <- FDPDses[c(2,1,3)]

# ====================================================================================================#

# (3) Functional Eveness is analogous to VPD (Tucker et al. 2017), however there is very little information
# on how to implement this index in R. Another distance based measure of phylogenetic evenness is Helmus et.
# als (2007) Phylogenetic sp. evennes (PSE) and this blog has info about its implementation in R:
# https://daijiang.name/en/2014/05/04/notes-func-phylo-book-1/

# PSE is a measure of how evenly spread the individuals in a sample are, in terms of their evolutionary history. 

Aoife_Ultramet <- force.ultrametric(Aoife_tree)

observed.pse <- pse(my.sample, Aoife_Ultramet) %>% as.data.frame()
observed.pse$Country <- rownames(observed.pse)
observed.pse <- observed.pse[,-2]
names(observed.pse)[1] <- 'Phylogenetic_sp_evennes'

FDPDses <- left_join(FDPDses, observed.pse, by = "Country")

# ====================================================================================================#

# (4) Functional Dispersion is analygous to phylogenetic Rao's Q, which can be calculated using picante's raoD()

Phy_raoq <- raoD(my.sample, Aoife_Ultramet)  %>% as.data.frame()
Phy_raoq <- Phy_raoq[,1]
Phy_raoq <- Phy_raoq %>% as.data.frame()
Phy_raoq$Country <- rownames(Wholeproduction)
names(Phy_raoq)[1] <- 'Phylogenetic_Raos'

FDPDindices <- left_join(FDPDindices, Phy_raoq, by = "Country") # and join to original country x div index matrix

# HOWEVER this function is computationally exhaustive as it performs inter-species comparisons for each site. 

# The Abundance weighted Faith index is a more practical index for measuring phylogenetic dispersion for this application:

# Swenson https://github.com/NGSwenson/lefse_0.5/blob/master/R/weighted.faith.R 
# has written a function which calculated weighted faith index for all communities in  sample:

weighted.faith <-
  function(my.phylo, my.sample){
    
    weighted.faith.function = function(my.sub.sample){
      
      ## extract the names of species in a community with an abundance greater than zero and use that information to make a pruned phylogeny for that community.
      tmp.tree = treedata(my.phylo, my.sub.sample[my.sub.sample > 0],warnings=F)$phy
      
      
      ## Create empty branches matrix
      branches = matrix(NA, nrow = nrow(tmp.tree$edge), ncol = 4)
      
      ## Fill first two columns of the matrix with node numbers defining each edge.
      branches[,1:2] = tmp.tree$edge
      
      ## Fill the third column with the length of each branch
      branches[,3] = tmp.tree$edge.length
      
      get.leaves<-function(x){
        leaves.node<-tips(tmp.tree,x[2])
      }
      
      ## Apply the get.leaves() function to each row in the branches object. This will retrieve species names subtended by each branch (i.e. ## row) in the branches matrix
      leaves = apply(branches, MARGIN = 1, get.leaves)
      
      ## Now quickly loop through each set of leaves to ## calculate the mean abundance (Ai) for those species.
      for(i in 1:length(leaves)){
        branches[i,4] = mean(my.sub.sample[leaves[[i]]], na.rm = T) 
        
      }
      
      ## Lastly calculated the Weighted Faith’s Index
      nrow(tmp.tree$edge) * ((sum(branches[,3] * branches[,4])) / sum(branches[,4]))
      
      
    }
    outt = apply(my.sample, MARGIN = 1, weighted.faith.function)
    outt
    
  }

observed.Faith.weighted <- weighted.faith(Aoife_tree, my.sample) %>% as.data.frame()

observed.Faith.weighted$Country <- rownames(observed.Faith.weighted)
names(observed.Faith.weighted)[1] <- 'Phylogenetic_Weighted_Faith'

FDPDses <- left_join(FDPDses, observed.Faith.weighted, by = "Country")

# ====================================================================================================#
# ====================================================================================================#
# ====================================================================================================#
# Null models for standardized diversity measures - Tip Shuffle randomisation
# ====================================================================================================#
# Null models for phylogenetic and functional analysis:
# https://daijiang.name/en/2014/05/29/null-model-for-phylogenetic-functional-diversity/

library(vegan)
library(picante)
library(phytools)
library(tidyverse)
library(geiger)
library(parallel)
library(FD)

# ====================================================================================================#

# Copy original trait & community matrices
my.sample <- Wholeproduction
traits <- Wholenutrients

# ===================================================
# Phylogenetic diversity (FRic equivalent)

# Faith PD
observed.pd

rand.Faith.fun <- function(x) {
  
  tmp.phylo <- tipShuffle(x)
  
  # and caclulate Faiths PD:
  pd(my.sample, tmp.phylo)[,1]
  
}

# Next we can replicate this function ten times to produce ten random MPD values
# (columns of the output) for each community (rows of the output).

null.output <- replicate(1000, rand.Faith.fun(Aoife_tree))
hist(null.output[1, ]) # make sure the values are changing each time #place a vertical line where your observed value lands in the histogram.

# Although plotting the null distribution and observed value is useful, we also wish to calculate summary statistics that can be reported 
# and used for quantitatively based inferences. Specifically, we are interested in computing the standardize effect size (S.E.S .)
# and quantile score. The quantile score, or where the observed value ranks in the null distribution, can be used to calculate a P -value.
# The S.E.S. is calculated as the observed minus the mean of the null distribution and this value divided by the standard deviation of
# the null distribution. 

# the apply() function can be used to calculate the S.E.S values for all communities simultaneously
ses.Faith <- (observed.pd[1] - 
                apply(null.output, MARGIN = 1, mean)) /
  apply(null.output, MARGIN = 1, sd) %>% as.data.frame()

# Here we use the apply() function to calculate the mean and standard deviation of the values in each row of the null model output matrix.

# ===================================================
# Phylogenetic Mean Pairwise Distances (FDiv equivalent)

rand.mpd.fun <- function(x) {
  
  tmp.phylo <- tipShuffle(x)
  
  # and caclulate Faiths PD:
  mpd(my.sample, cophenetic(tmp.phylo))
}

# Next we can replicate this function ten times to produce ten random MPD values
# (columns of the output) for each community (rows of the output).

null.output <- replicate(1000, rand.mpd.fun(Aoife_tree))

hist(null.output[1, ]) 
abline(v = mpd(my.sample, cophenetic(Aoife_tree))[1],
       col = "red", lwd = 2)

ses.mpd <- (observed.mpd[1] -
              apply(null.output, MARGIN = 1, mean)) /
  apply(null.output, MARGIN = 1, sd) %>% as.data.frame()

# ===================================================
# Phylogenetic Species Evenness (FEve Equivalent)

rand.pse.fun <- function(x) {
  
  tmp.phylo <- tipShuffle(x) %>% force.ultrametric() # need to make tree ultrametric for pse 
  
  # and calculate pse:
  pse(my.sample, tmp.phylo)[,1]
  
}

null.output <- replicate(1000, rand.pse.fun(Aoife_tree)) 
hist(null.output[1, ])

ses.pse <- (observed.pse[1] - 
              apply(null.output, MARGIN = 1, mean)) /
  apply(null.output, MARGIN = 1, sd) %>% as.data.frame()

# ===================================================
# Abundance Weighted Faith (FDis Equivalent)

rand.weightFaith.fun <- function(x) {
  
  tmp.phylo <- tipShuffle(x)
  
  # and calculate Weighted Faith index
  weighted.faith(tmp.phylo, my.sample)
  
}

null.output <- replicate(1000, rand.weightFaith.fun(Aoife_tree))
hist(null.output[1, ]) 

ses.weightFaith <- (observed.Faith.weighted[1] -
                      apply(null.output, MARGIN = 1, mean)) /
  apply(null.output, MARGIN = 1, sd) %>% as.data.frame()

# ====================================================================================================#
# ====================================================================================================#

# Parrelelize functional diversity randomisations 
# (much faster than replicate)


# Calculate the number of cores
no_cores <- detectCores() - 1

cl <- makeCluster(no_cores, type="FORK")
clusterSetRNGStream(cl)

#... then parallel replicate...
FRic.null.1000 <- parSapply(cl, 1:1000, function(i,...){
  
  x <- traits
  
  rownames(x) <- sample(rownames(x),length(rownames(x)),replace = FALSE)
  
  x <- x[match(colnames(my.sample), rownames(x)),] # This ensures all species and trait labels match those in the community matrix
  
  dbFD(x, my.sample, w.abun = TRUE, stand.x = TRUE, corr = "cailliez",
       calc.FRic = TRUE, calc.FDiv = FALSE, calc.CWM = FALSE, m = "min", messages = TRUE)$FRic
} )


#stop the cluster
stopCluster(cl)

# And calculate SES

SES.FRic <- (WholedbFD3$FRic -
               apply(FRic.null.1000, MARGIN = 1, mean)) /
  apply(FRic.null.1000, MARGIN = 1, sd) %>% as.data.frame()



# ===================================================
# FEve

start_time <- Sys.time()

cl <- makeCluster(no_cores, type="FORK")
clusterSetRNGStream(cl)

FEve.null.1000 <- parSapply(cl, 1:1000, function(i,...){
  
  x <- traits
  
  rownames(x) <- sample(rownames(x),length(rownames(x)),replace = FALSE)
  
  x <- x[match(colnames(my.sample), rownames(x)),] # This ensures all species and trait labels match those in the community matrix
  
  dbFD(x, my.sample, w.abun = TRUE, stand.x = TRUE, corr = "cailliez",
       calc.FRic = FALSE, calc.FDiv = FALSE, calc.CWM = FALSE, m = "min", messages = TRUE)$FEve
} )

end_time <- Sys.time()
end_time - start_time

# Time difference of 6.192757 hours

SES.FEve <- (WholedbFD3$FEve -
               apply(FEve.null.1000, MARGIN = 1, mean)) /
  apply(FEve.null.1000, MARGIN = 1, sd) %>% as.data.frame()

# ===================================================
# FDis


start_time <- Sys.time()

cl <- makeCluster(no_cores, type="FORK")
clusterSetRNGStream(cl)

FDis.null.1000 <- parSapply(cl, 1:1000, function(i,...){
  
  x <- traits
  
  rownames(x) <- sample(rownames(x),length(rownames(x)),replace = FALSE)
  
  x <- x[match(colnames(my.sample), rownames(x)),] # This ensures all species and trait labels match those in the community matrix
  
  dbFD(x, my.sample, w.abun = TRUE, stand.x = TRUE, corr = "cailliez",
       calc.FRic = FALSE, calc.FDiv = FALSE, calc.CWM = FALSE, m = "min", messages = TRUE)$FDis
} )

end_time <- Sys.time()
end_time - start_time

# Time difference of 6.192757 hours

SES.FDis <- (WholedbFD3$FDis -
               apply(FDis.null.1000, MARGIN = 1, mean)) /
  apply(FDis.null.1000, MARGIN = 1, sd) %>% as.data.frame()

# ===================================================
# FDiv

start_time <- Sys.time()

cl <- makeCluster(no_cores, type="FORK")
clusterSetRNGStream(cl)

FDiv.null.1000 <- parSapply(cl, 1:1000, function(i,...){
  
  x <- traits
  
  rownames(x) <- sample(rownames(x),length(rownames(x)),replace = FALSE)
  
  x <- x[match(colnames(my.sample), rownames(x)),] # This ensures all species and trait labels match those in the community matrix
  
  dbFD(x, my.sample, w.abun = TRUE, stand.x = TRUE, corr = "cailliez",
       calc.FRic = TRUE, calc.FDiv = TRUE, calc.CWM = FALSE, m = "min", messages = TRUE)$FDiv
} )

end_time <- Sys.time()
end_time - start_time

# Time difference ~10 hours

# ===================================================

FDiv.null.1000copy <- FDiv.null.1000[!rownames(FDiv.null.1000)=='Djibouti',]

WholedbFD4 <- WholedbFD3[!rownames(WholedbFD3)=='Djibouti',]

SES.FDiv <- (WholedbFD3$FDiv -
               apply(FDiv.null.1000, MARGIN = 1, function(x) mean(x, na.rm=TRUE))) /
  apply(FDiv.null.1000, MARGIN = 1, function(x) sd(x, na.rm=TRUE)) %>% as.data.frame()

# Results in all Nas, because theres at least one Na in every country across the 1000 randomizations.

# Need to convert these Nas to numbers in an unbiased way



# ===================================================

# Create a dataframe with all Standardized indices

rownames(ses.Faith) <- rownames(SES.FRic)
rownames(ses.mpd) <- rownames(SES.FRic)
rownames(ses.pse) <- rownames(SES.FRic)
rownames(ses.weightFaith) <- rownames(SES.FRic)

names(SES.FRic)[1] <- 'FRic'
names(SES.FDiv)[1] <- 'FDiv'
names(SES.FEve)[1] <- 'FEve'
names(SES.FDis)[1] <- 'FDis'

l <- list(ses.Faith,ses.mpd,ses.pse,ses.weightFaith,SES.FRic,SES.FDiv,SES.FEve,SES.FDis)    

SES.FDPD <-Reduce(merge, lapply(l, function(x) data.frame(x, rn = row.names(x))))

names(SES.FDPD)[1] <- 'Country'

# Remove Djibouti 
SES.FDPD <- SES.FDPD[!SES.FDPD$Country=='Djibouti',]
  


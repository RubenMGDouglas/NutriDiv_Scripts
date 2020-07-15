# ====================================================================================================#
# ====================================================================================================#
# Calculate FD indices with dbfd() from FD package
# ====================================================================================================#

# Paper on FD package
# http://adn.biol.umontreal.ca/~numericalecology/Reprints/Laliberte_Legendre_Ecology_2010.pdf

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

# FRic and FDiv both rely on finding the minimum convex hull that includes all species (Villéger et al. 2008). 
# This requires more species than traits. To circumvent this problem, dbFD takes only a subset of the PCoA axes as traits via argument m.
# This, however, comes at a cost of loss of information. The quality of the resulting reduced-space representation is returned by
# qual.FRic, which is computed as described by Legendre and Legendre (1998) and can be interpreted as a R^2-like ratio.

# ====================================================================================================#

# Paper on Fric, Feve & FDiv
# file:///Users/rubendouglas/Downloads/Villegeretal.2008Ecology.pdf

#Functional evenness describes the evenness of abundance distribution in a functional trait space (Mason et al. 2005).
# Basically, the new index quantifies the regularity with which the functional space is filled by species, weighted by their abundance.
# FEve decreases either when abundance is less evenly distributed among species or when functional distances among species are less regular

# ====================================================================================================#
# Branch length FD (Petchey & Gaston)
# Functional uniqueness?







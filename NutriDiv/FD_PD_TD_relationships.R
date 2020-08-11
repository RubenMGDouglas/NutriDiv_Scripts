# ====================================================================================================#
# ====================================================================================================#
# Relationships of Nutritional Functional Diversity with Phylogenetic & Taxonomic Diversity
# ====================================================================================================#

library(ggplot2)
library(maptools)
library(PerformanceAnalytics)
library(corrplot)

# Chart corr vis
chart.Correlation(SES.FDPD[,2:10], histogram=TRUE, pch=19)

# corrplot vis
res1 <- cor.mtest(SES.FDPD[,2:10], conf.level = 0.95)
M <- cor(SES.FDPD[,2:10])
corrplot(M, p.mat = res1$p, method = "color",
         sig.level = c(.001, .01, .05), pch.cex = .9,
         insig = "label_sig", pch.col = "black", order = "AOE")


Bivar_nbspFEve <- ggplot(SESwrld_sf, aes(x=nbsp, y=FEve)) + 
  theme_bw() +
  geom_smooth(method="lm") + labs(x="Sp. Richness", y = "FEve SES") +
  geom_point(aes(size = AREA)) 

Bivar_nbspFDiv <- ggplot(SESwrld_sf, aes(x=nbsp, y=FDiv)) + 
  theme_bw() +
  geom_smooth(method="lm") + labs(x="Sp. Richness", y = "FDiv SES") +
  geom_point(aes(size = AREA)) 

  # 
  # FDPDEvenness <- ggplot(SESwrld_sf, aes(x=Phylogenetic_sp_evennes, y=FEve, color = nbsp)) + 
  #   ggtitle("Evenness") + theme_bw() +
  #   geom_smooth(method="lm") + labs(x="PSE SES", y = "FEve SES") +
  #   geom_point(aes(size = AREA)) 
  # 
  # FDPDDivergence <- ggplot(SESwrld_sf, aes(x=Phylogenetic_MPD, y=FDiv, color = nbsp)) + 
  #   ggtitle("Divergence") + theme_bw() +
  #   geom_smooth(method="lm") + labs(x="PMPD SES", y = "FDiv SES") +
  #   geom_point(aes(size = AREA))
  # 
  # FDPDDispersion <- ggplot(SESwrld_sf, aes(x=Phylogenetic_Weighted_Faith, y=FDis, color = nbsp)) + 
  #   ggtitle("Dispersion") + theme_bw() +
  #   geom_smooth(method="lm") + labs(x="Weighted Faith SES", y = "FDis SES") +
  #   geom_point(aes(size = AREA))

# Faith PD
# SES.FDPD[,c(2:6,10:11)] %>%
#   gather(-FRic, -nbsp, -AREA, key = "var", value = "value") %>% 
#   ggplot(aes(x = value, y = FRic, color = nbsp)) +
#   geom_point(aes(size = AREA)) +
#   facet_wrap(~ var, scales = "free") +
#   theme_bw() +
#   geom_smooth(method="lm")
  
 # Plot FDis ~ MPD
  Bivar_FDisFaith <- ggplot(SESwrld_sf, aes(x=Phylogenetic_MPD, y=FDis, color = nbsp)) + 
    theme_bw() +
    geom_smooth(method="lm") + labs(x="Phylogenetic MPD SES", y = "FDis SES") +
    geom_point(aes(size = AREA))  
  
# ====================================================================================================#

# Resampling methods to quantify uncertainty & adequacy of fit with associated
# https://afit-r.github.io/resampling_methods

library(purrr)

# Use Leave-one-out cross-validation (LOOCV) to asess which polynomial models 1 through 3 fits the data best
# by choosing the model with the lowest LOOCV MSE.

# ===========================================================#

# FEve ~ Sp. Richness
  
# create function that computes LOOCV MSE based on specified polynomial degree
loocv_error <- function(x) {
  glm.fit <- glm(FEve ~ poly(nbsp, x), data = SES.FDPD)
  cv.glm(SES.FDPD, glm.fit)$delta[1]
}

# compute LOOCV prediction error (MSE) for polynomial degrees 1-5
1:5 %>% map_dbl(loocv_error) 
# [1] 1.0386095 1.0342967 0.9744319 0.9733113 0.9831380
# 4th order polynomial results in lowest MSE (0.9733113), 
# However no significant reduction from 3rd polynoial on, so will just use 3rd

# Run the Cubic model on the whole dataset and plot with R2 and P value. 
FEvenbspSES_Cu <- lm(FEve ~ poly(nbsp, 3), data = SES.FDPD)

hist(FEvenbspSES_Cu$residuals)
# residuals approx. normally distributed

summary(FEvenbspSES_Quar)

# labels for R^2 & P values
lb1 <- paste("R^2 == 0.09")
lb2 <- paste("P < 0.001")

Bivar_nbspFEve <- ggplot(SESwrld_sf, aes(x=nbsp, y=FEve)) + theme_bw() +
  stat_smooth(aes(y = FEve),method = "lm", formula = y ~ poly(x,3), size = 1) +
  labs(x="Sp. Richness", y = "FEve SES") + geom_point(aes(size = AREA)) +
  annotate(geom = "text", x = 23, y = 3, label=lb1, parse=TRUE) +
  annotate(geom = "text", x = 23, y = 2.5, label=lb2, parse=TRUE) +
  theme(legend.position="none")

# ===========================================================#

# FDiv ~ Sp. Richness

# create function that computes LOOCV MSE based on specified polynomial degree
loocv_error <- function(x) {
  glm.fit <- glm(FDiv ~ poly(nbsp, x), data = SES.FDPD)
  cv.glm(SES.FDPD, glm.fit)$delta[1]
}

# compute LOOCV prediction error (MSE) for polynomial degrees 1-5
1:5 %>% map_dbl(loocv_error) 
# [1] 0.9682789 0.9708273 0.9744601 0.9390882 0.9482262
# Quartic model has lowest test error

# Run the Quartic model on the whole dataset and plot with R2 and P value. 
FDivnbspSES_Quar <- lm(FDiv ~ poly(nbsp, 4), data = SES.FDPD)

hist(FDivnbspSES_Quar$residuals)
# residuals normally distributed

summary(FDivnbspSES_Quar)

# labels for R^2 & P values
lb3 <- paste("R^2 == 0.09")
lb4 <- paste("P < 0.001")

Bivar_nbspFDiv <- ggplot(SESwrld_sf, aes(x=nbsp, y=FDiv)) + theme_bw() +
  stat_smooth(aes(y = FDiv),method = "lm", formula = y ~ poly(x,4), size = 1) +
  labs(x="Sp. Richness", y = "FDiv SES") + geom_point(aes(size = AREA)) +
  annotate(geom = "text", x = 23, y = 3, label=lb3, parse=TRUE) +
  annotate(geom = "text", x = 23, y = 2.5, label=lb4, parse=TRUE) +
  theme(legend.position="none")

# ===========================================================#

# FRic ~ Faith PD 

# create function that computes LOOCV MSE based on specified polynomial degree
loocv_error <- function(x) {
  glm.fit <- glm(FRic ~ poly(Phylogenetic_Faith, x), data = SES.FDPD)
  cv.glm(SES.FDPD, glm.fit)$delta[1]
}

# compute LOOCV MSE for polynomial degrees 1-5
1:5 %>% map_dbl(loocv_error)
#[1] 1.057426 1.065150 1.057697 1.100477 1.092929
# Test error (MSE) is lowest in linear model (1.057426)

# Run the linear model on the whole dataset and plot with R2 and P value. 
FRicFaithSES_OLS <- lm(FRic ~ poly(Phylogenetic_Faith, 1), data = SES.FDPD)

hist(FRicFaithSES_OLS$residuals)
# residuals normally distributed

summary(FRicFaithSES_OLS)

# labels for R^2 & P values
lb5 <- paste("R^2 == 0.19")
lb6 <- paste("P < 0.001")

Bivar_FRicFaith <- ggplot(SESwrld_Moll, aes(x=Phylogenetic_Faith, y=FRic, color = nbsp)) + 
  theme_classic() + 
  geom_smooth(method="lm") + labs(x="Faith SES", y = "FRic SES") +
  geom_point(aes(size = AREA)) + 
  annotate(geom = "text", x = -4.5, y = 1.5, label=lb5, parse=TRUE) +
  annotate(geom = "text", x = -4.5, y = 1.1, label=lb6, parse=TRUE) +
  scale_size(range=c(1,5)) +
  theme(legend.position =  c(0.91, 0.5))


# ===========================================================#

# FDis ~ Phylogenetic MPD

# create function that computes LOOCV MSE based on specified polynomial degree
loocv_error <- function(x) {
  glm.fit <- glm(FDis ~ poly(Phylogenetic_MPD, x), data = SES.FDPD)
  cv.glm(SES.FDPD, glm.fit)$delta[1]
}

theme_classic()+# compute LOOCV MSE for polynomial degrees 1-5
1:5 %>% map_dbl(loocv_error)
# [1] 0.4767913 0.4756030 0.4949305 0.4768153 0.5325689
# linear model

FDisMPDSES_OLS <- lm(FDis ~ poly(Phylogenetic_MPD, 1), data = SES.FDPD)

hist(FDisMPDSES_OLS$residuals)
# residuals normally distributed

summary(FDisMPDSES_OLS)

# labels for R^2 & P values
lb7 <- paste("R^2 == 0.36")
lb8 <- paste("P < 0.0001")


Bivar_FDisMPD <- ggplot(SESwrld_Moll, aes(x=Phylogenetic_MPD, y=FDis, color = nbsp)) + 
  theme_bw() +
  geom_smooth(method="lm") + labs(x="Phylogenetic MPD SES", y = "FDis SES") +
  geom_point(aes(size = AREA)) + 
  annotate(geom = "text", x = -7, y = 1, label=lb7, parse=TRUE) +
  annotate(geom = "text", x = -7, y = 0.5, label=lb8, parse=TRUE) +
  theme_classic()+
  scale_size(range=c(1,5)) 
# ===========================================================#

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(Bivar_nbspFEve, Bivar_nbspFDiv, Bivar_FRicFaith, Bivar_FDisMPD, cols=2)

# ===========================================================#
# Bootstrapping to test model adequacy
# statistic <- function(data, index) {
#   lm.fit <- lm(FRic ~ Phylogenetic_Faith, data = SES.FDPD, subset = index)
#   coef(lm.fit)
# }
# 
# boot(SES.FDPD, statistic, 1000)
# summary(lm(mpg ~ poly(horsepower, 2), data = auto))
# 
# quad.statistic <- function(data, index) {
#   lm.fit <- lm(FRic ~ poly(Phylogenetic_Faith, 2), data = SES.FDPD, subset = index)
#   coef(lm.fit)
# }
# 
# boot(SES.FDPD, quad.statistic, 1000)
# summary(lm(mpg ~ poly(horsepower, 2), data = auto))
# ===========================================================#

# To ensure there is no heterogenetity in country size, can model residuals against country size
  

  
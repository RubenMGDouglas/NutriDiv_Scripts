# ====================================================================================================#
# ====================================================================================================#
# Canonical Correlation Analysis of Environmental drivers of Distance-based Nutritional FD (Country Level)
# ====================================================================================================#
library(tidyverse)
library(ggplot2)
library(GGally)
library(lme4)
library(CCA) #facilitates canonical correlation analysis
library(CCP) #facilitates checking the significance of the canonical 
library(kableExtra)
library(vegan)
library(candisc)
library(heplots)
library(yacca)

# Subset just FD SES
SES.FD <- SES.FDPD[c(1,6:10)]

# Combine all Environmental variables into one dataframe

All_Socioeconomic <- cbind(FoodSov_tibble[c(1,6,7,15,16,22,25,26,29:33,38:40,42,45:47)],
                           PNAS_tibble[c(34,38,40,53,54,56,57,61,65,69,74,77,78)])


dfs <- list(
TopoSoilEva_Country,
Bio1_Country,
Bio2_Country,
Bio12_Country,
Bio15_Country,
All_Socioeconomic
)

All_Env_vars <- Reduce(func, dfs)

# Now use matcor to look at correlations within and between the sets of variables

# annoyingly need to have country rownames so only numeric variables present
rownames(SES.FD) <- SES.FD$Country
SES.FD <- SES.FD[-1]

SES.FD <- SES.FD[which(rownames(SES.FD) %in% rownames(All_Env_vars)),]
All_Env_vars <- All_Env_vars[ order(row.names(All_Env_vars)), ]

# remove patchiest socioeconomic variables
Comp_Env_vars <- All_Env_vars[-c(27,33,39,51,52)]

Comp_Env_vars <- Comp_Env_vars[complete.cases(Comp_Env_vars),]
Comp.SES.FD <- SES.FD[which(rownames(SES.FD) %in% rownames(Comp_Env_vars)),]
Comp.SES.FD <- Comp.SES.FD[-c(5,6)]

# ========================================================#

# CCA in R
# https://medium.com/analytics-vidhya/canonical-correlation-analysis-cca-in-r-a-non-technical-primer-b67d9bdeb9dd

#checking the between and within set associations
cormat <- matcor(Comp.SES.FD, Comp_Env_vars)

#extracting the within study correlations for set 1 and set 2 and #between set cor
round(cormat$Ycor, 4) # Set 1

round(cormat$Xcor, 4) # Set 2

#between set associations
between <- round(cormat$XYcor, 4)
between %>% mutate(cell_spec(FRic, "html", color = ifelse(FRic > 0.2, "red", "blue")))

between %>%
  kable() %>%
  kable_styling()

#obtaining the canonical correlations
can_cor1 <- cc(Comp.SES.FD, Comp_Env_vars)
can_cor1$cor

#raw canonical coefficients
raw_coeff_y <- can_cor1[4] %>% as.data.frame() # Shows for example, that with a 1 unit increase in population density 
# results in an increase of 17.38 units of canonical coefficent 1 for the environmental set. 

#compute the canonical loadings
can_cor2 <- comput(Comp.SES.FD, Comp_Env_vars, can_cor1)
can_cor2[3:6]  #the canonical loadings

corr.Y.xscores <- D %>% as.data.frame()
corr.X.yscores <- can_cor2[5] %>% as.data.frame()

plt.cc(can_cor1, d1 = 1, d2 = 2, type = "b", var.label = T)

can_cor1$cor
#===========================

# make a table with the highest variable cross loadings

# cross loadings of Environmental variables on Nutritional FD canonical dimensions
corr.Y.xscores %>%
  mutate(
    Variable = row.names(.),
    Dim.1 = cell_spec(corr.Y.xscores.1, "html", color = ifelse(corr.Y.xscores.1 > 0.2 | corr.Y.xscores.1 < -0.2, "red", "black")),
    Dim.2 = cell_spec(corr.Y.xscores.2, "html", color = ifelse(corr.Y.xscores.2 > 0.2 | corr.Y.xscores.2 < -0.2, "red", "black"))
      ) %>%
  select(Variable, Dim.1, Dim.2) %>%
  kable(format = "html", escape = FALSE) %>%
  kable_styling()

# ====================================================================================================#

# test of canonical dimensions
rho <- can_cor1$cor
##defining the number of observations, no of variables in first set,
#and number of variables in second set
n=dim(Comp.SES.FD)[1]
p=length(Comp.SES.FD)
q=length(Comp_Env_vars)
##Calculating the F approximations using different test statistics
#using wilks test statistic

CCA_Wilks <- p.asym(rho,n,p,q,tstat="Wilks") %>% 
  as.data.frame() 

CCA_Wilks[1:6] %>%
  kable() %>%
  kable_styling()

# only the first two canonical dimensions are significant at the 0.05 alpha level.

#standardizing the first set of canonical coefficients (FD indices)
std_coef1<-diag(sqrt(diag(cov(Comp.SES.FD))))
std_coef1%*%can_cor1$xcoef %>%#[,1:2] %>%
  kable() %>%
  kable_styling()

##standardizing the coeficents of the second set (Environmental vars)
std_coef2<-diag(sqrt(diag(cov(Comp_Env_vars))))

std_coef2%*%can_cor1$ycoef[,1:2] %>%
  kable() %>%
  kable_styling()


?cancor

out <- CCorA(Comp.SES.FD, Comp_Env_vars)

biplot(out, "b", cex = c(0.6,0.6), expand = 1.25) 

out

cc$structure$Y.xscores[,1:2] %>%
  kable() %>%
  kable_styling()

# ====================================================================================================#

cc <- cancor(Comp.SES.FD[1:4], Comp_Env_vars, set.names=c("NFD SES", "Environmental"),
             xcenter = TRUE, ycenter = TRUE, xscale = TRUE, yscale = TRUE
             )

# Cross loadings of Env variables on NFD variables
corr.Y.xscores <- cc$structure$Y.xscores[,1:2] %>% as.data.frame()

# cross loadings of Environmental variables on Nutritional FD canonical dimensions
corr.Y.xscores %>%
  mutate(
    Variable = row.names(.),
    Dim.1 = cell_spec(Xcan1, "html", color = ifelse(Xcan1 > 0.2 | Xcan1 < -0.2, "red", "black")),
    Dim.2 = cell_spec(Xcan2, "html", color = ifelse(Xcan2 > 0.2 | Xcan2 < -0.2, "red", "black"))
  ) %>%
  select(Variable, Dim.1, Dim.2) %>%
  kable(format = "html", escape = FALSE) %>%
  kable_styling()

# Loadings
corr.X.xscores <- cc$structure$X.xscores[,1:2] %>% as.data.frame()
corr.X.xscores %>%
  mutate(
    Variable = row.names(.),
    Dim.1 = cell_spec(Xcan1, "html", color = ifelse(Xcan1 > 0.4 | Xcan1 < -0.4, "red", "black")),
    Dim.2 = cell_spec(Xcan2, "html", color = ifelse(Xcan2 > 0.2 | Xcan2 < -0.4, "red", "black"))
  ) %>%
  select(Variable, Dim.1, Dim.2) %>%
  kable(format = "html", escape = FALSE) %>%
  kable_styling()

heplot(cc, fill.alpha=0.2,
       var.cex=0.8, var.col="orange", var.lwd=1
)


# ==============================================#

# Followed by Redundancy analysis (RDA) to determine the total variance in the first two canonical
# functions of the FD set explained by the environmental set (R^2)

cc_RDA <- redundancy(cc)
cc_RDA$Xcan.redun
# 37%  of variance in NFD values accounted for by the Environmental variables
# if we just take into account the first two (significant) canonical dimensions

cc_RDA$Ycan.redun
# 9% of variance in Environmental values accounted for by NFD values in the
# two significant canonical dimensions. 

# ====================================================================================================#
# 
# # Yacca package for CCA and helio plots
# library(yacca)
#  
# 
# sapply(SES.FDPD, class)
# 
# yacca.cca <- cca(SES.FDPD[2:9], Comp_Env_vars[1:6])#, xcenter = TRUE, 
#                  ycenter = TRUE, xscale = TRUE, yscale = TRUE, use = "complete.obs", na.rm = TRUE)
# 
# helio.plot(yacca.cca)

# ====================================================================================================#
# Second iteration of CCA

# Eliminate variables which are redundant or showed little contribution to the canonical dimensions.

Env_vars_2 <- All_Env_vars[,c(1:4,7:10,13,15:20,23,34,37,40,43,45,50)]
Env_vars_2 <- Env_vars_2[complete.cases(Env_vars_2),]
sapply(Env_vars_2, class)

SES.PDFD_2 <- SES.FDPD
rownames(SES.PDFD_2) <- SES.PDFD_2$Country
SES.PDFD_2 <- SES.PDFD_2[2:10]
sapply(SES.PDFD_2, class)
SES.PDFD_2$nbsp <- as.numeric(SES.PDFD_2$nbsp)

#subset to only countriesw with complete environmental data
SES.PDFD_2 <- SES.PDFD_2[which(rownames(SES.PDFD_2) %in% rownames(Env_vars_2)),]

CCA_2 <- cancor(SES.PDFD_2, Env_vars_2, set.names=c("NFD SES", "Environmental"),
                xcenter = TRUE, ycenter = TRUE, xscale = TRUE, yscale = TRUE
)


# test of canonical dimensions
rho <- CCA_2$cancor
##defining the number of observations, no of variables in first set,
#and number of variables in second set
n=dim(SES.PDFD_2)[1]
p=length(SES.PDFD_2)
q=length(Env_vars_2)
##Calculating the F approximations using different test statistics
#using wilks test statistic

CCA_Wilks <- p.asym(rho,n,p,q,tstat="Wilks") %>% 
  as.data.frame() 

CCA_Wilks[1:6] %>%
  kable() %>%
  kable_styling()


# cross loadings of Environmental variables on Diversity canonical dimensions
corr.Y.xscores %>%
  mutate(
    Variable = row.names(.),
    Dim.1 = cell_spec(Xcan1, "html", color = ifelse(Xcan1 > 0.2 | Xcan1 < -0.2, "red", "black")),
    Dim.2 = cell_spec(Xcan2, "html", color = ifelse(Xcan2 > 0.2 | Xcan2 < -0.2, "red", "black"))
  ) %>%
  select(Variable, Dim.1, Dim.2) %>%
  kable(format = "html", escape = FALSE) %>%
  kable_styling()

# Loadings
corr.X.xscores <- CCA_2$structure$X.xscores[,1:2] %>% as.data.frame()
corr.X.xscores %>%
  mutate(
    Variable = row.names(.),
    Dim.1 = cell_spec(Xcan1, "html", color = ifelse(Xcan1 > 0.4 | Xcan1 < -0.4, "red", "black")),
    Dim.2 = cell_spec(Xcan2, "html", color = ifelse(Xcan2 > 0.2 | Xcan2 < -0.4, "red", "black"))
  ) %>%
  select(Variable, Dim.1, Dim.2) %>%
  kable(format = "html", escape = FALSE) %>%
  kable_styling()

heplot(cc, fill.alpha=0.2,
       var.cex=0.8, var.col="orange", var.lwd=1
)

# Redundancy analysis
CCA_2_RDA <- redundancy(CCA_2)
CCA_2_RDA$Xcan.redun
# 37%  of variance in NFD values aCCA_2ounted for by the Environmental variables
# if we just take into aCCA_2ount the first two (significant) canonical dimensions

CCA_2_RDA$Ycan.redun
# 9% of variance in Environmental values aCCA_2ounted for by NFD values in the
# two significant canonical dimensions. 



# ====================================================================================================#
# Tidying Jenny's Crop x Nutrient matrix 
# ====================================================================================================#

library(tidyverse)
library(dplyr)

setwd("/Users/rubendouglas/OneDrive - The Royal Botanic Gardens, Kew/Nutridiv R/Crop x Nutrition matrices/")

Tidy_CN_mat <- readxl::read_xlsx("Ruben Update.xlsx")


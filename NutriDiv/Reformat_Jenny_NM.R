# ====================================================================================================#
# ====================================================================================================#
# Reformat Nutrition matrix 
# ====================================================================================================#
library(tidyverse)
library(tidyr)
library(dplyr)

# We prioritize retaining crop species over nutrient traits, so need to set a reasonable threshold of NAs
# However for dbFD() some missing values can be tolerated. 
#_____________________________
# Import Jenny's nutrient data

# NB: Nutrient data retrieved via Langual: https://www.langual.org/langual_linkcategory.asp?CategoryID=4&Category=Food+Composition
setwd("/Users/rubendouglas/OneDrive - The Royal Botanic Gardens, Kew/Nutridiv R/Crop x Nutrition matrices/")
Tidy_CN_mat <- readxl::read_xlsx("Ruben Update.xlsx")
#_____________________________

#Convert columns 12:77 to numeric
num <- sapply(Tidy_CN_mat[,12:77],as.numeric)
Tidy_CN_mat <- cbind(Tidy_CN_mat[,1:11],num)

#Remove columns with all or mostly zero values
Tidy_CN_mat <- Tidy_CN_mat[,-c(44,46,49,50,51,57,63,64,66,76,77)]

# Remove taxonomic columns
Tidy_CN_mat <- Tidy_CN_mat[,-c(4:9)]

# Edible part row 98 need to be renamed from "NA" to "seed"
Tidy_CN_mat[98,5] <- 'seed'
#_____________________________

# Harmonize names with those in presence absence matrix & production matrix 

# get the set of crop names from the list of Presabs .tif files
setwd("/Users/rubendouglas/OneDrive - The Royal Botanic Gardens, Kew/Nutridiv R/Agg_threshold_binary/")
#list the tifs, convert the tifs to rastr and stack the raster files 
list_tifs <- as.list(list.files())
PAcropnames  <- sub('\\_.*','',list_tifs)

# Create a vector of intersecting crop names 
Crop_names <- intersect(colnames(Production_mat_25min),PAcropnames)

# Create a column with just the Crop name in Tidy_CN_mat:
Tidy_CN_mat$`Crop name` <- sub('\\,.*','',Tidy_CN_mat$`Food name`)

Tidy_CN_mat$`Crop name` <- str_replace_all(Tidy_CN_mat$`Crop name`, c("green broad beans (fava beans)" = "greenbroadbean", 
                                                                      "green chili pepper" = "chilleetc", "green pea" = "pea"))

# Take away all strings after and including a space 
Tidy_CN_mat$`Crop name` <- sub('\\ .*','',Tidy_CN_mat$`Crop name`)

# de-capitalize A,E,L,P,Q
Tidy_CN_mat$`Crop name` <- str_replace_all(Tidy_CN_mat$`Crop name`, c("A" = "a", "E" = "e", "L" = "l", "P" = "p", "Q" = 'q'))

#_____________________________

Tidy_CN_mat$`Crop name` <- str_replace_all(Tidy_CN_mat$`Crop name`, c("anise" = 'aniseetc',"chickpeas" = 'chickpea',
                                                                       "blueberries" = 'blueberry', "cherries" = 'cherry',
                                                                       "clementines" = "tangetc", "colanut" = 'kolanut',
                                                                       "corn" = 'maize',"cranberries" = 'cranberry',
                                                                       "cucumber" = 'cucumberetc', "dates" = 'date',
                                                                       "flaxseed" = 'flax',"grapefruit" = "grapefruitetc",
                                                                       "onions" = "onion", "green" = "greenbroadbean",
                                                                       "oranges" = "orange", "peaches" = "peachetc",
                                                                       "pears" = "pear", "persimmons" = "persimmon",
                                                                       "pigeon" = "pigeonpea", "plantains" = "plantain",
                                                                       "plums" = "plum", "potatoes" = "potato", 
                                                                       "pumpkin" = "pumpkinetc", "squash" = "pumpkinetc",
                                                                       "gourd" = "pumpkinetc", "quinces" = "quince",
                                                                       "raspberries" = "rasberry","sour" = "sourcherry",
                                                                       "strawberries" = "strawberry", "string" = "stringbean",
                                                                       "sweet" = "sweetpotato", "tangerines" = "tangetc",
                                                                       "walnuts" = "walnut", "lemons" = "citrusnes", 
                                                                       "limes" = "citrusnes", "lentils" = "lentil", 
                                                                       "lupins" = "lupin", "mangos" = "mango", "guava" = "mango",
                                                                       "melons" = "melonetc", "oil" = "oilpalm", "olives" = "olive",
                                                                       "kiwifruit" = "kiwi"))


intersect(Crop_names, Tidy_CN_mat$`Crop name`)
# 109 intersecting names (length of unique crop names in new column)

#_____________________________

# add FAO grop group using left_join with Old Allnutrients matrix
FAO <- Allnutrients[,1:2]

Tidy_CN_mat <- left_join(Tidy_CN_mat, FAO, by = "Crop name")

# assign FAO group to the few remaining NAs

for (i in 1:nrow(Tidy_CN_mat)){
  
if (Tidy_CN_mat$`Crop name`[i] == "aniseetc") {
  Tidy_CN_mat$`Group (FAO)`[i] <- "other crops"
}
else if (Tidy_CN_mat$`Crop name`[i] == "cucumberetc") {
  Tidy_CN_mat$`Group (FAO)`[i] <- "vegetables & melons"
}
  else if (Tidy_CN_mat$`Crop name`[i] == "grapefruitetc") {
    Tidy_CN_mat$`Group (FAO)`[i] <- "fruit"
  }
  else if (Tidy_CN_mat$`Crop name`[i] == "citrusnes") {
    Tidy_CN_mat$`Group (FAO)`[i] <- "fruit"
  }
  else if (Tidy_CN_mat$`Crop name`[i] == "melonetc") {
    Tidy_CN_mat$`Group (FAO)`[i] <- "vegetables & melons"
  }
  else if (Tidy_CN_mat$`Crop name`[i] == "peachetc") {
    Tidy_CN_mat$`Group (FAO)`[i] <- "fruit"
  }
  else if (Tidy_CN_mat$`Crop name`[i] == "pumpkinetc") {
    Tidy_CN_mat$`Group (FAO)`[i] <- "vegetables & melons"
  }
  else if (Tidy_CN_mat$`Crop name`[i] == "tangetc") {
    Tidy_CN_mat$`Group (FAO)`[i] <- "fruit"
  }
  }


#_____________________________

# get rid of maize row with NAs
Tidy_CN_mat <- Tidy_CN_mat[-61,]

#_____________________________

# omit columns with ≥5% NAs (6 missing values or more)
Tidy_CN_mat5 <- Tidy_CN_mat[, colSums(is.na(Tidy_CN_mat)) <= 5]
# leaves 35 columns

# omit columns with ≥10% NAs (12 missing values or more)
Tidy_CN_mat10 <- Tidy_CN_mat[, colSums(is.na(Tidy_CN_mat)) <= 12]
# leaves 46 columns (Optimal)

# omit columns with ≥20% NAs
Tidy_CN_mat20 <- Tidy_CN_mat[, colSums(is.na(Tidy_CN_mat)) <= 24]
# leaves 49 columns

#_____________________________
# Remove Rows with NA values

Nona_Nutrients5 <- Tidy_CN_mat5[complete.cases(Tidy_CN_mat5),]
# 100 crops, 30 continous nutrient variables

Nona_Nutrients10 <- Tidy_CN_mat10[complete.cases(Tidy_CN_mat10),]
# 94 crops, 46 continous nutrient variables

#_____________________________

# Aggregate rows with the same Crop file name, so mean nutrient values are given
# (e.g. mango & guava only have a single class 'mango' in the PA & production maps)

Nutrients5 <- aggregate(Tidy_CN_mat5, by = list(Tidy_CN_mat5$`Crop name`, Tidy_CN_mat5$`Group (FAO)`), FUN = mean)

# Remove cols 3-7 + 40-41

Nutrients5 <- Nutrients5[,-c(3:7,40:41)]

# Rename first two cols
names(Nutrients5)[names(Nutrients5)=="Group.1"] <- "Crop name"
names(Nutrients5)[names(Nutrients5)=="Group.2"] <- "Group (FAO)"

sapply(Nutrients5, class)

# make a list of species which were aggregated for reference:

# "Tangetc" : Citrus reticulata, tangerines, raw | Citrus clementina hort. ex Tanaka clementines, raw
# "Pumpkinetc" : Pumpkin, raw | Cucurbita spp, squash, summer, all varieties, raw |  Lagenaria siceraria, gourd, white-flowered, raw (bottle gourd, elongate)
# "Pea" : Pisum sativum, pea, edible, podded, raw | Pisum sativum, green pea
# "Onion": Allium cepa, onions, raw | Allium cepa, onions, green, spring (including tops and bulbs)
# "Mangos": Psidium guajava, guava, common, raw | Mangifera indica, mangos, raw
# "Citrusnes": Citrus x limon, lemons, raw, without peel | Citrus x latifolia, limes, raw
#_____________________________

#And do the same for the no NA matrix too:
Nona_Nutrients5 <- aggregate(Nona_Nutrients5, by = list(Nona_Nutrients5$`Crop name`, Nona_Nutrients5$`Group (FAO)`), FUN = mean)

# Remove cols 3-7 + 40-41
Nona_Nutrients5 <- Nona_Nutrients5[,-c(3:7,40:41)]

# Rename first two cols
names(Nona_Nutrients5)[names(Nona_Nutrients5)=="Group.1"] <- "Crop name"
names(Nona_Nutrients5)[names(Nona_Nutrients5)=="Group.2"] <- "Group (FAO)"

#_____________________________

# Two output matrices to use in dbFD calculations:
# (1) 
Nutrients5
# (2) 
Nona_Nutrients5

# (all nutrient variables are numeric)



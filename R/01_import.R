#Load TransPlantNetwork data
load("data/TransPlantNetwork_060824.RData")

# Get the Koppen-Gieger climate region maps. Retrieved on 16th of July 2024: https://www.gloh2o.org/koppen/
kgmap = rast("data/koppen_geiger_0p1.tif")
kglegend = read.csv2("data/koppen_geiger_0p1_legend.txt", header=FALSE)%>%
  dplyr::select(-V4)%>%
  rename(Climate = V1, Climate_code = V2, Climate_explanation = V3)


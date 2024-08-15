library(tidyverse)
library(emmeans)
library(DescTools) #AUC
library(lme4)
library(sf)
library(terra)

# set colors
dataset.colors = c("#D55E00" , "#0072B2")
names(dataset.colors) = c("TP", "TRY")

odt.colors = c("#648FFF", "#FE6100", "#DC267F")
names(odt.colors) = c("originControls","destControls","warmed")

contrast.colors = c("grey30", "#FE6100", "#648FFF")
names(contrast.colors) = c("between_controls","to_destination","to_origin")

source("R/01_import.R")
source("R/02_organize.R")
source("R/03_calculate_winners_losers.R")
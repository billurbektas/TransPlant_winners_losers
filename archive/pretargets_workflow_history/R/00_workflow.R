library(tidyverse)
library(emmeans)
library(DescTools) #AUC
library(lme4)
library(sf)
library(terra)
library(metafor)
library(ggpubr)
library(patchwork)
library(factoextra)
library(FactoMineR)
library(lmerTest)
library(VGAM)
library(geomtextpath)

# set colors
dataset.colors = c("#D55E00" , "#0072B2")
names(dataset.colors) = c("TP", "TRY")

odt.colors = c("#648FFF", "#FE6100", "#DC267F")
names(odt.colors) = c("originControls","destControls","warmed")

contrast.colors = c("grey30", "#FE6100", "#648FFF")
names(contrast.colors) = c("between_controls","to_destination","to_origin")

cat.colors = c("#118AB2","#06D6A0","#EF476F", "grey80")
names(cat.colors) = c("colonizing_winner", "winner", "loser", "nc")

source("R/01_import.R")
source("R/02_organize.R")
#source("R/03_calculate_winners_losers.R")
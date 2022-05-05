#### 0. Load R packages ####
library(rgl)
library(geomorph)
library(devtools)
install_github("marta-vidalgarcia/morpho.tools.GM")
library(morpho.tools.GM)
library(Morpho)
library(Rvcg)
library(magick)
library(ggplot2)
library(vegan)

#### 1. QUESTIONS ASYMMETRY ####

# 1. Determine how much of the shape variance is explained by the asymmetric component
# 2. Look at fluctuating asymmetry vs directional asymmetry
# 3. Mirror each side of the face, so double up on specs, run GMM on these new datasets, by treatment
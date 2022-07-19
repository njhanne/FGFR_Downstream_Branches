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
# Use Morpho's GPA instead of geomorph's

# 2. Look at fluctuating asymmetry vs directional asymmetry

# Questions 1 & 2 are going to be very straight-forward, 3 more difficult

# Make plots with vector arrows (perhaps even 3*SD, to stress out differences) to see direction of shape changes on each side


# 3. Mirror each side of the face, so double up on specs, run GMM on these new datasets, by treatment
# We will end up with these funny-looking specimens, 4 treatments (CTRL-left, CTRL-right, MUT-left, MUT-right)

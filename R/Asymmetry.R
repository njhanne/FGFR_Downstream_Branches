#### 0. Load R packages ####
library(rgl)
library(geomorph)
library(devtools)
install_github("marta-vidalgarcia/morpho.tools.GM", force = TRUE)
library(morpho.tools.GM)
install_github("marta-vidalgarcia/symmetry")
library(symmetry)
library(mesh_process)
library(Morpho)
library(Rvcg)
library(magick)
library(Evomorph)
library(ggplot2)
library(vegan)
install_github("vqv/ggbiplot")
library(ggbiplot)
library(factoextra)
library(gt)
library(abind)

#### 1. QUESTIONS ASYMMETRY ####

# 1. Determine how much of the shape variance is explained by the asymmetric component
# Use Morpho's GPA instead of geomorph's


# 2. Look at fluctuating asymmetry vs directional asymmetry

# Questions 1 & 2 are going to be very straight-forward, 3 more difficult

# Make plots with vector arrows (perhaps even 3*SD, to stress out differences) to see direction of shape changes on each side


# 3. Mirror each side of the face, so double up on specs, run GMM on these new datasets, by treatment
# We will end up with these funny-looking specimens, 4 treatments (CTRL-left, CTRL-right, MUT-left, MUT-right)

#### 2. LOAD DATA ####
# Classifiers
classifiers_unord  <- read.csv("./data/classifiers.csv", header = TRUE)
head(classifiers_unord)
tail(classifiers_unord)

str(classifiers_unord)
classifiers_unord$treatment <- as.factor(classifiers_unord$treatment)
row.names(classifiers_unord) <- classifiers_unord$id
classifiers_unord

# Landmark array & GPA
head_array <- readRDS("./data/Head_LM_array_FGF_embryos.rds")
GPA_geomorph <- readRDS("./data/GPA_FGF_embryos.rds")

classifiers <- classifiers_unord[match(dimnames(head_array)[[3]], row.names(classifiers_unord)),]

# Curveslide
curveslide_all <- read.csv("./data/curveslide.csv")

# Head surface
head_surface.lm <- c(34:51)

# Atlases
head_mesh <- geomorph::read.ply("./data/ATLAS_chick_ctr_23_smooth_ascii_no_back.ply") # Not sure what is wrong with this mesh
head_lowres <- vcgQEdecim(head_mesh, percent = 0.15)

atlas_head_lm <- head_array[,, which(dimnames(head_array)[[3]] == "chick_ctr_23")]



# FIND LANDMARK PAIRS
?detect.symmetry
detect.symmetry(head_array[1:12,,], sym.plane = "yz", plot = TRUE)
detect.symmetry(head_array[1:33,,], sym.plane = "yz", plot = TRUE)
detect.symmetry(head_array[,,], sym.plane = "yz", plot = TRUE)


non.sym <- c(9, 10, 13:15)
side.1 <- c(2,4, 6, 8, 12, 16:24, 34:42)
side.2 <- c(1, 3, 5, 7, 11, 25:33, 43:51)

pairedLM <- cbind(side.1, side.2)

pairedLM

#### 3. ALL - ANALYSIS OF BILATERAL SYMMETRY ####

SYM_FGF <- bilat.symmetry(head_array, side = NULL, replicate = NULL, object.sym = TRUE, 
                          ind = dimnames(head_array )[[3]], land.pairs = pairedLM, iter = 999, seed = NULL, RRPP = TRUE)

summary(SYM_FGF)
cat("SYM_FGF", capture.output(summary(SYM_FGF)), 
    file="./output/SYM_FGF.txt", sep="\n", append=TRUE)

SYM_FGF$symm.shape

str(SYM_FGF$FA.component)
str(SYM_FGF$DA.component)


#### 3.1. SYMMETRIC COMPONENT ####

PCA_SYM_FGF <- gm.prcomp(SYM_FGF$symm.shape)
classifiers$treatment
palette(c("navy", "darkorange"))

png("./figs/PCA_symmetric_component_FGF.png", width = 750, height = 600)
plot(PCA_SYM_FGF, pch = 19, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_SYM_FGF, classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("bottomleft", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA FGF SYMMETRIC COMPONENT")
dev.off()

pdf("./figs/PCA_symmetric_component_FGF.pdf", width = 7.5, height = 6)
plot(PCA_SYM_FGF, pch = 19, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_SYM_FGF, classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("bottomleft", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA FGF SYMMETRIC COMPONENT")
dev.off()

svg("./figs/PCA_symmetric_component_FGF.svg", width = 7.5, height = 6)
plot(PCA_SYM_FGF, pch = 19, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_SYM_FGF, classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("bottomleft", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA FGF SYMMETRIC COMPONENT")
dev.off()

ANOVA_ALL_sym <- procD.lm(SYM_FGF$symm.shape ~ classifiers$treatment, 
                          iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_ALL_sym)
cat("ANOVA_SYM_FGF", capture.output(summary(ANOVA_ALL_sym)), 
    file="./output/ANOVA_symmetric_component_FGF.txt", sep="\n", append=TRUE)

### 3.2. ASYMMETRIC COMPONENT ####

PCA_ASYM_FGF <- gm.prcomp(SYM_FGF$asymm.shape)
classifiers$treatment
palette()
png("./figs/PCA_asymmetric_component_FGF.png", width = 750, height = 600)
plot(PCA_ASYM_FGF, pch = 19, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_ASYM_FGF, classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA FGF ASYMMETRIC COMPONENT")
dev.off()

pdf("./figs/PCA_asymmetric_component_FGF.pdf", width = 7.5, height = 6)
plot(PCA_ASYM_FGF, pch = 19, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_ASYM_FGF, classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA FGF ASYMMETRIC COMPONENT")
dev.off()

svg("./figs/PCA_asymmetric_component_FGF.svg", width = 7.5, height = 6)
plot(PCA_ASYM_FGF, pch = 19, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_ASYM_FGF, classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA FGF ASYMMETRIC COMPONENT")
dev.off()

identify(x = PCA_ASYM_FGF$x[,1], y = PCA_ASYM_FGF$x[,2], labels=row.names(PCA_ASYM_FGF$x))

ANOVA_ASYM_FGF  <- procD.lm(SYM_FGF$asymm.shape ~ classifiers$treatment, 
                            iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_ASYM_FGF)
cat("ANOVA_ASYM_FGF", capture.output(summary(ANOVA_ASYM_FGF)), 
    file="./output/ANOVA_asymmetric_component_FGF.txt", sep="\n", append=TRUE)


#### 3.3. FLUCTUATING ASYMMETRY COMPONENT ####
PCA_FA_FGF <- gm.prcomp(gpagen(SYM_FGF$FA.component)$coords)
classifiers$treatment
palette()
png("./figs/PCA_fluctuating_asymmetry_component_FGF.png", width = 750, height = 600)
plot(PCA_FA_FGF, pch = 19, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_FA_FGF, classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA FGF FLUCTUATING ASYMMETRY COMPONENT")
dev.off()

pdf("./figs/PCA_fluctuating_asymmetry_component_FGF.pdf", width = 7.5, height = 6)
plot(PCA_FA_FGF, pch = 19, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_FA_FGF, classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA FGF FLUCTUATING ASYMMETRY COMPONENT")
dev.off()

svg("./figs/PCA_fluctuating_asymmetry_component_FGF.svg", width = 7.5, height = 6)
plot(PCA_FA_FGF, pch = 19, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_FA_FGF, classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA FGF FLUCTUATING ASYMMETRY COMPONENT")
dev.off()

# identify(x = PCA_FA_FGF$x[,1], y = PCA_FA_FGF$x[,2], labels=row.names(PCA_FA_FGF$x))


ANOVA_FA_FGF  <- procD.lm(SYM_FGF$FA.component ~ classifiers$treatment, 
                          iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_FA_FGF)
cat("ANOVA_FA_FGF", capture.output(summary(ANOVA_ASYM_FGF)), 
    file="./output/ANOVA_Fluctuating_Asymmetry_FGF.txt", sep="\n", append=TRUE)

#### 3.4. DIRECTIONAL ASYMMETRY COMPONENT ####
summary(SYM_FGF)
SYM_FGF$DA.component # array with side 1 & side 2

ANOVA_DA_FGF  <- procD.lm(SYM_FGF$DA.component ~ classifiers$treatment, 
                          iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_FA_FGF)

#### 3.5. TO DO NEXT ####
#### 3.6. PCA ####
# Make morphs for each component (PCA)
# Make panels PC1-PC4
# Make heatmap for DA left vs right
# Make table with proportion of variance explained by each component
# Remake PCA plots like Costello comparisons

#### 3.7. CVA ####
# Same as above




#### 4. MIRRORING ####

# Mirror each side of the face, so double up on specs, run GMM on these new datasets, by treatment
# We will end up with these funny-looking specimens, 4 treatments (CTRL-left, CTRL-right, MUT-left, MUT-right)


## QUESTIONS ASYMMETRY #

# Mirror each side of the face, so double up on specs, run GMM on these new datasets, by treatment
# We will end up with these funny-looking specimens, 4 treatments (CTRL-left, CTRL-right, MUT-left, MUT-right)

#### 4.1 LOAD DATA ####
# Classifiers
classifiers_unord  <- read.csv("./data/classifiers.csv", header = TRUE)
head(classifiers_unord)
tail(classifiers_unord)

str(classifiers_unord)
classifiers_unord$treatment <- as.factor(classifiers_unord$treatment)
row.names(classifiers_unord) <- classifiers_unord$id
classifiers_unord

# Landmark array & GPA
head_array <- readRDS("./data/Head_LM_array_FGF_embryos.rds")
GPA_geomorph <- readRDS("./data/GPA_FGF_embryos.rds")

classifiers <- classifiers_unord[match(dimnames(head_array)[[3]], row.names(classifiers_unord)),]

# Curveslide
curveslide_all <- read.csv("./data/curveslide.csv")

# Head surface
head_surface.lm <- c(34:51)

# Atlases
head_mesh <- geomorph::read.ply("./data/ATLAS_chick_ctr_23_smooth_ascii_no_back.ply") # Not sure what is wrong with this mesh
head_lowres <- vcgQEdecim(head_mesh, percent = 0.15)

atlas_head_lm <- head_array[,, which(dimnames(head_array)[[3]] == "chick_ctr_23")]



# 4.2 FIND LANDMARK PAIRS ####
# ?detect.symmetry
# detect.symmetry(GPA_geomorph$coords, sym.plane = "yz", plot = TRUE)
# detect.symmetry(head_array[1:33,,], sym.plane = "yz", plot = TRUE)
# detect.symmetry(head_array[,,], sym.plane = "yz", plot = TRUE)


non.sym <- c(9, 10, 13:15)
side.1 <- c(2,4, 6, 8, 12, 16:24, 34:42) # LEFT
side.2 <- c(1, 3, 5, 7, 11, 25:33, 43:51) # RIGHT



open3d(zoom = 0.75, windowRect = c(0, 0, 700, 700)) 
rgl::shade3d(head_lowres, color = "gray", alpha =0.9)
rgl::plot3d(atlas_head_lm, aspect = "iso", type = "s", size=1.2, col = "darkblue", add = T)
rgl::text3d(x = atlas_head_lm[,1],
            y = atlas_head_lm[,2],
            z = atlas_head_lm[,3],
            texts = c(1:dim(atlas_head_lm)[1]),
            cex = 1.5, offset = 0.5, pos = 1)


# 4.3. MIRRORING LANDMARKS ####
# Mirror landmarks on both sides and generate two new arrays
# the non-symmetrical landmarks also remain the same

array_side1_mirrored <- GPA_geomorph$coords
array_side2_mirrored <- GPA_geomorph$coords

half_array_side1 <- sweep(GPA_geomorph$coords[side.1,,1], MARGIN = 2, c(-1,1,1), `*`) 

open3d()
plot3d(GPA_geomorph$coords[non.sym, , 1], col = "black", type = "s", aspect = "iso", 
       size = 1, add = TRUE, xlab = "x", ylab = "y", zlab = "z")
plot3d(GPA_geomorph$coords[side.1, , 1], col = "green", type = "s", aspect = "iso", 
       size = 1, add = TRUE, xlab = "x", ylab = "y", zlab = "z")
plot3d(half_array_side1, col = "blue", type = "s", aspect = "iso", 
       size = 1, add = TRUE, xlab = "x", ylab = "y", zlab = "z")
plot3d(half_array_side1, col = "red", type = "s", aspect = "iso", 
       size = 1, add = TRUE, xlab = "x", ylab = "y", zlab = "z")



for (i in 1:dim(head_array)[3]){
  array_side1_mirrored[side.2,,i] <- sweep(GPA_geomorph$coords[side.1,,1], 
                                           MARGIN = 2, c(-1,1,1), `*`)
  
  array_side2_mirrored[side.1,,i] <- sweep(GPA_geomorph$coords[side.2,,1], 
                                           MARGIN = 2, c(-1,1,1), `*`)
}

dorsal <- par3d()$userMatrix
open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = dorsal) 
plot3d(GPA_geomorph$coords[, , 1], col = "black", type = "s", aspect = "iso", 
       size = 1, add = TRUE, xlab = "x", ylab = "y", zlab = "z")
plot3d(array_side1_mirrored[,,1], col = "chartreuse", type = "s", aspect = "iso", 
       size = 0.75, add = TRUE, xlab = "x", ylab = "y", zlab = "z")

which(classifiers_mirrored$treatment == "triple")

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = dorsal) 
plot3d(GPA_geomorph$coords[, , 27], col = "black", type = "s", aspect = "iso", 
       size = 1, add = TRUE, xlab = "x", ylab = "y", zlab = "z")
plot3d(array_side1_mirrored[,,27], col = "chartreuse", type = "s", aspect = "iso", 
       size = 0.75, add = TRUE, xlab = "x", ylab = "y", zlab = "z")

# Change dimnames so we know which way they were mirrored
dimnames(array_side1_mirrored)[[3]] <- paste0("LEFT_", dimnames(array_side1_mirrored)[[3]])
dimnames(array_side2_mirrored)[[3]] <- paste0("RIGHT_", dimnames(array_side2_mirrored)[[3]])

# And join both arrays

mirrored_array <- abind(array_side1_mirrored, array_side2_mirrored, along = 3)

dimnames(mirrored_array)[[3]]

# Finally make new classifiers matrix

classifiers_mirrored <- rbind(classifiers, classifiers)
row.names(classifiers_mirrored) <- dimnames(mirrored_array)[[3]]
classifiers_mirrored$id <- dimnames(mirrored_array)[[3]]
classifiers_mirrored$treatment_mirror <- vector(mode = "character", length = dim(classifiers_mirrored)[1])
classifiers_mirrored$treatment_mirror[1:dim(classifiers)[1]] <- paste0("LEFT_", classifiers_mirrored$treatment[1:dim(classifiers)[1]])
classifiers_mirrored$treatment_mirror[(1+dim(classifiers)[1]):dim(classifiers_mirrored)[1]] <- paste0("RIGHT_", classifiers_mirrored$treatment[(1+dim(classifiers)[1]):dim(classifiers_mirrored)[1]])

classifiers_mirrored$treatment_mirror <- as.factor(classifiers_mirrored$treatment_mirror)


# 4.4. GPA & PCA - both sides mirrored ####
surface_semis <- c(34:51)
GPA_mirrored_double <- geomorph::gpagen(A = mirrored_array*c(GPA_geomorph$Csize,GPA_geomorph$Csize), curves = as.matrix(curveslide_all), 
                                        surfaces = surface_semis)

GPA_mirrored_LEFT <- geomorph::gpagen(A = mirrored_array[,,1:dim(array_side1_mirrored)[3]]*GPA_geomorph$Csize, curves = as.matrix(curveslide_all), 
                                      surfaces = surface_semis)
GPA_mirrored_RIGHT <- geomorph::gpagen(A = mirrored_array[,,(1+dim(array_side1_mirrored)[3]):dim(mirrored_array)[3]]*GPA_geomorph$Csize, curves = as.matrix(curveslide_all), 
                                       surfaces = surface_semis)

saveRDS(mirrored_array, "./data/mirrored_array_both_sides_May2023.rds")
saveRDS(GPA_mirrored_double, "./data/mirrored_array_both_sides_May2023.rds")
saveRDS(GPA_mirrored_LEFT, "./data/mirrored_array_LEFT_May2023.rds")
saveRDS(GPA_mirrored_RIGHT, "./data/mirrored_array_RIGHT_May2023.rds")
write.csv(classifiers_mirrored, "./data/classifiers_mirrored_both_sides_May2023.csv")

# Both sides
PCA_both_sides_May2023 <- gm.prcomp(GPA_mirrored_double$coords)
summary(PCA_both_sides_May2023)
str(PCA_both_sides_May2023)

# Delete file if it exists
if (file.exists("./output/PCA_both_sides_May2023.txt")) {
  file.remove("./output/PCA_both_sides_May2023.txt")
}
cat("PCA shape variables raw - both sides mirrored", capture.output(summary(PCA_both_sides_May2023)), 
    file="./output/PCA_both_sides_May2023.txt", sep="\n", append=TRUE)


levels(classifiers_mirrored$treatment_mirror)
palette(c("navy", "darkorange", "cornflowerblue", "goldenrod1"))

pdf("./figs/PCA_head_shape_treatment_sides_May2023.pdf", width = 8.25, height = 6)
plot(PCA_both_sides_May2023, pch = 19, col = classifiers_mirrored$treatment_mirror, cex = 1.25)
ordiellipse(PCA_both_sides_May2023, classifiers_mirrored$treatment_mirror, kind="ehull",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("bottomright", pch = 19, col = palette(), legend = levels(classifiers_mirrored$treatment_mirror))
title("PCA of shape coordinates - mirrored side and treatment")
dev.off()


png("./figs/PCA_head_shape_treatment_sides_May2023.png", width = 750, height = 600)
plot(PCA_both_sides_May2023, pch = 19, col = classifiers_mirrored$treatment_mirror, cex = 1.25)
ordiellipse(PCA_both_sides_May2023, classifiers_mirrored$treatment_mirror, kind="ehull",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("bottomright", pch = 19, col = palette(), legend = levels(classifiers_mirrored$treatment_mirror))
title("PCA of shape coordinates - mirrored side and treatment")
dev.off()

t<- geomorph.data.frame(GPA_mirrored_double, treatment = classifiers_mirrored$treatment_mirror, Pdist = Pdist)

ANOVA_both_mirrored  <- procD.lm(coords ~ treatment, data=t, 
                          iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_both_mirrored)
ANOVA_both_mirrored_pw <- pairwise(ANOVA_both_mirrored, groups = t$treatment)
summary(ANOVA_both_mirrored_pw)

# 4.5. PCA Left side ####

PCA_LEFT_May2023 <- gm.prcomp(GPA_mirrored_LEFT$coords)
summary(PCA_LEFT_May2023)
str(PCA_LEFT_May2023)

# Delete file if it exists
if (file.exists("./output/PCA_LEFT_May2023.txt")) {
  file.remove("./output/PCA_LEFT_May2023.txt")
}
cat("PCA shape variables raw - LEFT side mirrored", capture.output(summary(PCA_LEFT_May2023)), 
    file="./output/PCA_LEFT_May2023.txt", sep="\n", append=TRUE)


levels(classifiers$treatment)
palette(c("navy", "darkorange"))

pdf("./figs/PCA_head_shape_treatment_mirrored_LEFT_May2023.pdf", width = 8.25, height = 6)
plot(PCA_LEFT_May2023, pch = 19, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_LEFT_May2023, classifiers$treatment, kind="ehull",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("bottomright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA of shape coordinates - mirrored LEFT and treatment")
dev.off()


png("./figs/PCA_head_shape_treatment_LEFT_May2023.png", width = 750, height = 600)
plot(PCA_LEFT_May2023, pch = 19, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_LEFT_May2023, classifiers$treatment, kind="ehull",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("bottomright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA of shape coordinates - mirrored LEFT and treatment")
dev.off()

t_LEFT <- geomorph.data.frame(GPA_mirrored_LEFT, treatment = classifiers$treatment)

ANOVA_LEFT_mirrored  <- procD.lm(coords ~ treatment, data=t_LEFT, 
                                 iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_LEFT_mirrored)


# 4.6. Right side ####

PCA_RIGHT_May2023 <- gm.prcomp(GPA_mirrored_RIGHT$coords)
summary(PCA_RIGHT_May2023)
str(PCA_RIGHT_May2023)

# Delete file if it exists
if (file.exists("./output/PCA_RIGHT_May2023.txt")) {
  file.remove("./output/PCA_RIGHT_May2023.txt")
}
cat("PCA shape variables raw - RIGHT side mirrored", capture.output(summary(PCA_RIGHT_May2023)), 
    file="./output/PCA_RIGHT_May2023.txt", sep="\n", append=TRUE)


levels(classifiers$treatment)
palette(c("cornflowerblue", "goldenrod1"))

pdf("./figs/PCA_head_shape_treatment_mirrored_RIGHT_May2023.pdf", width = 8.25, height = 6)
plot(PCA_RIGHT_May2023, pch = 19, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_RIGHT_May2023, classifiers$treatment, kind="ehull",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("bottomright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA of shape coordinates - mirrored RIGHT and treatment")
dev.off()


png("./figs/PCA_head_shape_treatment_RIGHT_May2023.png", width = 750, height = 600)
plot(PCA_RIGHT_May2023, pch = 19, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_RIGHT_May2023, classifiers$treatment, kind="ehull",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("bottomright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA of shape coordinates - mirrored RIGHT and treatment")
dev.off()

t_RIGHT <- geomorph.data.frame(GPA_mirrored_RIGHT, treatment = classifiers$treatment)

ANOVA_RIGHT_mirrored  <- procD.lm(coords ~ treatment, data=t_RIGHT, 
                                 iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_RIGHT_mirrored)

# 4.7. Procrustes distance ####


Pdist <- ShapeDist(GPA_mirrored_double$coords, GPA_mirrored_double$consensus)

gdf_mirrored <- geomorph.data.frame(GPA_mirrored_double, 
                                treatment = classifiers_mirrored$treatment, 
                                treatment_mirror = classifiers_mirrored$treatment_mirror, 
                                Pdist = Pdist)

ggplot_df <- as.data.frame(cbind(as.character(gdf_mirrored$treatment_mirror), 
                                 as.character(gdf_mirrored$treatment), 
                                 as.character(gdf_mirrored$Pdist)))
colnames(ggplot_df) <- c("treatment_mirror", "treatment", "Pdist")

row.names(ggplot_df) <- dimnames(gdf_mirrored$coords)[[3]]

head(ggplot_df)

ggplot_df$treatment <- as.factor(ggplot_df$treatment)
ggplot_df$treatment_mirror <- as.factor(ggplot_df$treatment_mirror)
ggplot_df$Pdist <- as.numeric(as.character(ggplot_df$Pdist))

str(ggplot_df)


pdf("./figs/Pdist_treatment_mirrored.pdf", width = 6.5, height = 6.5)
ggplot(ggplot_df, aes(Pdist, fill = treatment_mirror)) +
  scale_fill_manual(values = c("navy", "darkorange", "cornflowerblue", "goldenrod1")) + geom_density(alpha = 0.75) + 
  ggtitle("Procrustes distances head shape - mirrored") + xlab("Procrustes distance") + ylab("Relative density") +
  theme(plot.title = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"), legend.text=element_text(size=10), 
        legend.title=element_text(size=12, face = "bold"))
dev.off()

png("./figs/Pdist_treatment_mirrored.png", width = 650, height = 650)
ggplot(ggplot_df, aes(Pdist, fill = treatment_mirror)) +
  scale_fill_manual(values = c("navy", "darkorange", "cornflowerblue", "goldenrod1")) + geom_density(alpha = 0.75) + 
  ggtitle("Procrustes distances head shape - mirrored") + xlab("Procrustes distance") + ylab("Relative density") +
  theme(plot.title = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"), legend.text=element_text(size=10), 
        legend.title=element_text(size=12, face = "bold"))
dev.off()


# 5. MORPHS MIRRORING ####
load("./figs/RGL_head_heatmaps_pos.rdata")
# frontal <- par3d()$userMatrix
# 
# rgl.close()
# 
# save(dorsal, frontal, file = "./figs/RGL_head_heatmaps_pos.rdata")

# load mesh of the specimen closest to the mean shape
head_mesh <- geomorph::read.ply("./data/ATLAS_chick_ctr_23_smooth_ascii_only_face.ply") 
# The mesh should only be surface, and be just have the frontal side here
head_lowres <- vcgQEdecim(head_mesh, percent = 0.15)
atlas_head_lm <- head_array[,, which(dimnames(head_array)[[3]] == "chick_ctr_23")]

levels(gdf_mirrored$treatment_mirror)

# Mean shape by treatment and mirrored side
# LEFT_control_mean_shape <- mshape(mirrored_array[,,which(gdf_mirrored$treatment_mirror == "LEFT_control")])
# RIGHT_control_mean_shape <- mshape(mirrored_array[,,which(gdf_mirrored$treatment_mirror == "RIGHT_control")])
# LEFT_triple_mean_shape <- mshape(mirrored_array[,,which(gdf_mirrored$treatment_mirror == "LEFT_triple")])
# RIGHT_triple_mean_shape <- mshape(mirrored_array[,,which(gdf_mirrored$treatment_mirror == "RIGHT_triple")])

LEFT_control_mean_shape <- mshape(gdf_mirrored$coords[,,which(gdf_mirrored$treatment_mirror == "LEFT_control")])
RIGHT_control_mean_shape <- mshape(gdf_mirrored$coords[,,which(gdf_mirrored$treatment_mirror == "RIGHT_control")])
LEFT_triple_mean_shape <- mshape(gdf_mirrored$coords[,,which(gdf_mirrored$treatment_mirror == "LEFT_triple")])
RIGHT_triple_mean_shape <- mshape(gdf_mirrored$coords[,,which(gdf_mirrored$treatment_mirror == "RIGHT_triple")])



# Create morphed meshes
LEFT_control_mesh <- tps3d(head_lowres, as.matrix(atlas_head_lm), LEFT_control_mean_shape, threads = 1)
RIGHT_control_mesh <- tps3d(head_lowres, as.matrix(atlas_head_lm), RIGHT_control_mean_shape, threads = 1)
LEFT_triple_mesh <- tps3d(head_lowres, as.matrix(atlas_head_lm), LEFT_triple_mean_shape, threads = 1)
RIGHT_triple_mesh <- tps3d(head_lowres, as.matrix(atlas_head_lm), RIGHT_triple_mean_shape, threads = 1)


# Plot morphs
open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = frontal)
rgl.pop("lights")
light3d(specular="black")
shade3d(LEFT_control_mesh, color="gray", alpha=0.9)
rgl.snapshot("./figs/Morph_LEFT_control_frontal_head.png", top = TRUE)
writePLY("./output/LEFT_control_mesh.ply")
rgl.close()

open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = frontal)
rgl.pop("lights")
light3d(specular="black")
shade3d(RIGHT_control_mesh, color="gray", alpha=0.9)
rgl.snapshot("./figs/Morph_RIGHT_control_frontal_head.png", top = TRUE)
writePLY("./output/RIGHT_control_mesh.ply")
rgl.close()

open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = frontal)
rgl.pop("lights")
light3d(specular="black")
shade3d(LEFT_triple_mesh, color="gray", alpha=0.9)
rgl.snapshot("./figs/Morph_LEFT_triple_frontal_head.png", top = TRUE)
writePLY("./output/LEFT_triple_mesh.ply")
rgl.close()

open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = frontal)
rgl.pop("lights")
light3d(specular="black")
shade3d(RIGHT_triple_mesh, color="gray", alpha=0.9)
rgl.snapshot("./figs/Morph_RIGHT_triple_frontal_head.png", top = TRUE)
writePLY("./output/RIGHT_triple_mesh.ply")
rgl.close()


# HEATMAPS
open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = frontal)
x11(width=1.7, height=8)
meshDist(LEFT_control_mesh, LEFT_triple_mesh, rampcolors = c("blue", "white", "red"), sign = TRUE)
rgl.snapshot("./figs/Heatmap_LEFT_ctrl_triple_mirrored_frontal.png", top = TRUE)
dev.print(pdf, "./figs/Heatmap_LEFT_ctrl_triple_mirrored_scale.pdf", width=2, height=9.5)
dev.off()
clear3d()
rgl.close()

open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = frontal)
x11(width=1.7, height=8)
meshDist(RIGHT_control_mesh, RIGHT_triple_mesh, rampcolors = c("blue", "white", "red"), sign = TRUE)
rgl.snapshot("./figs/Heatmap_RIGHT_ctrl_triple_mirrored_frontal.png", top = TRUE)
dev.print(pdf, "./figs/Heatmap_RIGHT_ctrl_triple_mirrored_scale.pdf", width=2, height=9.5)
dev.off()
clear3d()
rgl.close()

open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = frontal)
x11(width=1.7, height=8)
meshDist(LEFT_control_mesh, RIGHT_control_mesh, rampcolors = c("blue", "white", "red"), sign = TRUE)
rgl.snapshot("./figs/Heatmap_LEFT_RIGHT_ctrlmirrored_frontal.png", top = TRUE)
dev.print(pdf, "./figs/Heatmap_LEFT__RIGHT_ctrl_mirrored_scale.pdf", width=2, height=9.5)
dev.off()
clear3d()
rgl.close()


open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = frontal)
x11(width=1.7, height=8)
meshDist(RIGHT_triple_mesh, LEFT_triple_mesh, rampcolors = c("blue", "white", "red"), sign = TRUE)
rgl.snapshot("./figs/Heatmap_LEFT_RIGHT_triple_mirrored_frontal.png", top = TRUE)
dev.print(pdf, "./figs/Heatmap_LEFT__RIGHT_triple_mirrored_scale.pdf", width=2, height=9.5)
dev.off()
clear3d()
rgl.close()



#### 5.2. PC morphs mirroring ####
# Have a look at the heatmap, and change colour if needed
open3d(zoom=0.75, userMatrix = frontal, windowRect= c(0,0,1000,700))
x11(width=1.7, height=8)
PC_min <- tps3d(head_lowres, as.matrix(atlas_head_lm), PCA_both_sides_May2023$shapes[[1]]$min, threads = 1)
PC_max <- tps3d(head_lowres, as.matrix(atlas_head_lm), PCA_both_sides_May2023$shapes[[1]]$max, threads = 1)

meshDist(PC_min, PC_max, rampcolors = c("darkblue", "blue", "white", "red", "darkred"), sign = TRUE)

summary(PCA_both_sides_May2023)
n_dimensions <- 10 # number of PCs to include in the figure, decide depending on % of variance explained in PCs

# This loop automatically positions faces in frontal only for the 10 PCs
for (i in 1:n_dimensions){
  PC_min <- tps3d(head_lowres, as.matrix(atlas_head_lm), PCA_both_sides_May2023$shapes[[i]]$min, threads = 1)
  PC_max <- tps3d(head_lowres, as.matrix(atlas_head_lm), PCA_both_sides_May2023$shapes[[i]]$max, threads = 1)
  PC <- paste0("PC", i)
  
  #frontal views
  open3d(zoom=0.75, userMatrix = frontal, windowRect= c(0,0,1000,700))
  shade3d(PC_max, color="grey")
  rgl.snapshot(paste0("./figs/pc_morphs/mirrored/",PC,"_","max_frontal.png"), top = TRUE )
  clear3d()
  rgl.close()
  
  
  #frontal views
  open3d(zoom=0.75, userMatrix = frontal, windowRect= c(0,0,1000,700))
  shade3d(PC_min, color="grey")
  rgl.snapshot(paste0("./figs/pc_morphs/mirrored/",PC,"_","min_frontal.png"), top = TRUE )
  clear3d()
  rgl.close()
  
  
  #Heatmaps
  open3d(zoom=0.75, userMatrix = frontal, windowRect= c(0,0,1000,700))
  x11(width=1.7, height=8)
  meshDist(PC_min, PC_max, rampcolors = c("darkblue", "blue", "white", "red", "darkred"), sign = TRUE)
  rgl.snapshot(paste0("./figs/pc_morphs/mirrored/heat_",PC,"_","frontal.png"), top = TRUE)
  # dev.copy(pdf, paste0("./figs/pc_morphs/mirrored/heat_",PC,"_scale.pdf"), width=2, height=9.5)
  dev.print(pdf, paste0("./figs/pc_morphs/mirrored/heat_",PC,"_scale.pdf"), width=2, height=9.5)
  dev.off()
  clear3d()
  rgl.close()
  
  rm(PC_max, PC_min, PC) 
}
#END LOOP





#Create page with PCs
i <- 1
PC_min_img_sup <- image_read(paste0("./figs/pc_morphs/mirrored/PC",i,"_min_frontal.png"))
PC_max_img_sup <- image_read(paste0("./figs/pc_morphs/mirrored/PC",i,"_max_frontal.png"))
PC_HT_img_sup <- image_read(paste0("./figs/pc_morphs/mirrored/heat_PC",i,"_frontal.png"))


PC_col <- c(PC_min_img_sup,PC_max_img_sup,PC_HT_img_sup)

imgs <- image_append(image_scale(PC_col, "x200"), stack = TRUE)
# imgs <- image_border(imgs, "white", "0x15")
# imgs <- image_annotate(imgs, paste0("PC",i), font = "Times",
#                        location = "+100+0", size = 60)
imgs <- image_border(imgs, "white", "85x15")
imgs <- image_annotate(imgs, paste0("PC",i), font = "Times",
                       location = "+195+0", size = 60)
imgs <- image_annotate(imgs, paste0("PCmin"), font = "Times",
                       location = "+0+95", size = 32)
imgs <- image_annotate(imgs, paste0("PCmax"), font = "Times",
                       location = "+0+295", size = 32)
imgs <- image_annotate(imgs, paste0("Heatmap"), font = "Times",
                       location = "+0+495", size = 32)


image_browse(imgs)

for (i in 2:n_dimensions){
  PC_min_img_sup <- image_read(paste0("./figs/pc_morphs/mirrored/PC",i,"_min_frontal.png"))
  PC_max_img_sup <- image_read(paste0("./figs/pc_morphs/mirrored/PC",i,"_max_frontal.png"))
  PC_HT_img_sup <- image_read(paste0("./figs/pc_morphs/mirrored/heat_PC",i,"_frontal.png"))
  
  
  PC_row <- c(PC_min_img_sup,PC_max_img_sup,PC_HT_img_sup)
  img2 <-image_append(image_scale(PC_row, "x200"), stack = TRUE)
  img2 <-image_border(img2, "white", "0x15")
  img2 <-image_annotate(img2, paste0("PC",i), font = "Times", 
                        location = "+100+0", size = 60)
  
  imgs <-c(imgs,img2)
  
  imgs <-image_append(image_scale(imgs))
}

i

image_browse(imgs)
image_write(imgs, path = paste0("./figs/pc_morphs/mirrored_SHAPE_heatmap_morphs_PC1-", i, ".png"), format = "png")

# We are done with these analyses



# 6. Face integration CTRL vs Treatment ####

?integration.test

non.sym <- c(9, 10, 13:15)
side.1 <- c(2,4, 6, 8, 12, 16:24, 34:42) # LEFT
side.2 <- c(1, 3, 5, 7, 11, 25:33, 43:51) # RIGHT

side <- vector(mode = "character", length = 51)
side[side.1] <- "left"
side[side.2] <- "right"
side <- side[-non.sym]

# Integration face all
face_integration <- integration.test(GPA_geomorph$coords[-non.sym,,], partition.gp = side, iter = 999)
summary(face_integration) # Test summary
plot(face_integration) # PLS plot

# Compare integration of the face between treatment and control
face_integration_CTRL <- integration.test(GPA_geomorph$coords[-non.sym, , which(classifiers$treatment == "control")], 
                                          partition.gp = side, iter = 999)

face_integration_TREATMENT <- integration.test(GPA_geomorph$coords[-non.sym, , which(classifiers$treatment == "triple")], 
                                               partition.gp = side, iter = 999)

summary(face_integration_CTRL) # Test summary
summary(face_integration_TREATMENT) # Test summary
plot(face_integration_CTRL) # PLS plot
plot(face_integration_TREATMENT)

PLS_comparison <- compare.pls(CTRL = face_integration_CTRL, TRIPLE = face_integration_TREATMENT)

# Delete file if it exists
if (file.exists("./output/Integration_face_comparisons_CTRL_TRIPLE_May2023.txt")) {
  file.remove("./output/Integration_face_comparisons_CTRL_TRIPLE_May2023.txt")
}
cat("INTEGRATION on the face comparison - CTRL vs TREATMENT", capture.output(summary(PLS_comparison)), 
    file="./output/Integration_face_comparisons_CTRL_TRIPLE_May2023.txt", sep="\n", append=TRUE)



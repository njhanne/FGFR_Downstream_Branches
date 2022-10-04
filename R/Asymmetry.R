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


#### 3.1 SYMMETRIC COMPONENT ####

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

### 3.2 ASYMMETRIC COMPONENT ####

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


#### 3.3 FLUCTUATING ASYMMETRY COMPONENT ####
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

#### 3.4 DIRECTIONAL ASYMMETRY COMPONENT ####
summary(SYM_FGF)
SYM_FGF$DA.component # array with side 1 & side 2

#### 3.5 TO DO NEXT ####
# Make morphs for each component (PCA)
# Make panels PC1-PC4
# Make heatmap for DA left vs right
# Make table with proportion of variance explained by each component
# Remake PCA plots like Costello comparisons

# Mirror each side of the face, so double up on specs, run GMM on these new datasets, by treatment
# We will end up with these funny-looking specimens, 4 treatments (CTRL-left, CTRL-right, MUT-left, MUT-right)


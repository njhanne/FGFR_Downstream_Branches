#### Load R packages ####
# install.packages('devtools')
# library(devtools)
# install_github("marta-vidalgarcia/morpho.tools.GM", force = TRUE)
# install_github("marta-vidalgarcia/symmetry")
# install.packages(c('rgl', 'geomorph' 'devtools', 'Morpho', 'Rvcg', 'magick', 'Evomorph', 'ggplot2', 'vegan', 'factoextra', 'gt'))

library(rgl)
library(geomorph)
library(morpho.tools.GM)
library(symmetry)
# install_github("marta-vidalgarcia/mesh_process")
# library(mesh_process)
library(Morpho)
library(Rvcg)
library(magick)
library(Evomorph)
library(ggplot2)
library(vegan)
# install_github("vqv/ggbiplot")
# library(ggbiplot)
library(factoextra)
library(gt)

library(abind)
library(stringr)
library(dplyr)


#### 0 Helpers ####

#### 0.2 Mesh and 3D plot helpers ####
get_dec_mesh <- function(mesh='face') {
  if (mesh == 'head'){
    if (file.exists('chick_ctr_1_decimated.ply')) {
      decimated_mesh <- vcgImport("chick_ctr_1_decimated.ply")
    } else {
      full_mesh <- vcgImport("chick_ctr_1.ply") # load 3d mesh of contol embryo
      decimated_mesh <- vcgQEdecim(full_mesh, percent = 0.15) # decimate mesh to reduce complexity
      decimated_mesh <- vcgQEdecim(decimated_mesh, percent = 0.25)
      vcgPlyWrite(decimated_mesh, 'chick_ctr_1_decimated.ply')
    }
  } else {
    if (file.exists('ATLAS_chick_ctr_face_decimated.ply')) {
      decimated_mesh <- vcgImport("ATLAS_chick_ctr_face_decimated.ply")
    } else {
      full_mesh <- vcgImport("ATLAS_chick_ctr_face.ply")
      decimated_mesh <- vcgQEdecim(full_mesh, percent = 0.15)
      vcgPlyWrite(decimated_mesh, 'ATLAS_chick_ctr_face_decimated.ply')
    }
  }
  return(decimated_mesh)
}


get_palette <- function(lights=FALSE) {
  if (lights) {
    return(c('#dddddd', '#89cced', '#ccddaa', '#ccbb43', '#994556', '#bbbbbb', '#0177bb', '#13783d','#989936', '#882256'))
  } else {
    return(c('#bbbbbb', '#0177bb', '#13783d','#989936', '#882256'))
  }
}


#### Main ####
#### 1 Load data ####
# R doesn't have a good way to get the path of where this file is located, unfortunately
# if you are running this code in Rstudio, try this:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd() #check our working directory
setwd("../../../data/Morphology/3D")
getwd() #check our working directory
# if you aren't using rstudio, use the command setwd() 
# and point it to the data/Morphology/3D directory

# Classifiers
#### 1.2 Load and match classifiers ####
setwd("./lm_data/Landmarks/")
sample_names <- list.files(path = getwd(), pattern="fcsv$") # get all fiducial files
sample_names <- str_remove(sample_names, "_Fiducials.fcsv$") # remove file extension
classifiers_unord <- data.frame(sample_names) # make dataframe
colnames(classifiers_unord) <- "id" # rename first column
# pull treatment from filename
classifiers_unord$treatment <- str_match(classifiers_unord$id, "_(.+)_")[,2] 
# rename treatments
classifiers_unord <- classifiers_unord %>% mutate(treatment = case_when(treatment =='ctr' ~ 'DMSO',
                                                                        treatment =='expredo' ~ 'Triple',
                                                                        treatment =='LY' ~ 'LY294002',
                                                                        treatment =='U0' ~ 'U0126',
                                                                        treatment =='U73' ~ 'U73122'))
# make treatment a factor and set level order for graphing
classifiers_unord$treatment <- factor(classifiers_unord$treatment, levels = c('DMSO', 'U0126', 'LY294002', 'U73122', 'Triple'))
row.names(classifiers_unord) <- classifiers_unord$id


# match the landmark sample names to the classifiers csv
# classifiers <- classifiers_unord[match(dimnames(LMs)[[3]], row.names(classifiers_unord)),]

# Landmark array & GPA
setwd("../..")

head_array <- readRDS("./lm_data/Head_LM_array_FGF_embryos.rds")
GPA_geomorph <- readRDS("./lm_data/GPA_FGF_embryos.rds") #'gpa_head' from gmm_analysis.R
og_lm <- head_array[,, which(dimnames(head_array)[[3]] == "chick_ctr_23")]
atlas_head_lm <- head_array[,, which(dimnames(head_array)[[3]] == "chick_ctr_23")]

# remove PC1 from gpa and redo PLS
PCA_head <- gm.prcomp(GPA_geomorph$coords)
PC1_regression <- procD.lm(GPA_geomorph$coords~PCA_head$x[,1])
new_shapes <- colMeans(two.d.array(GPA_geomorph$coords)) + t(PC1_regression$residuals)
head_array <- arrayspecs(t(new_shapes), 51,3)
GPA_geomorph$coords <- head_array
# end PC1 code

classifiers <- classifiers_unord[match(dimnames(head_array)[[3]], row.names(classifiers_unord)),]

# Curveslide
curveslide_all <- read.csv("./lm_data/curveslide.csv")

# Head surface
head_surface.lm <- c(34:51)

# Find landmark pairs across midline
non.sym <- c(9, 10, 13:15)
side.1 <- c(2,4, 6, 8, 12, 16:24, 34:42) # contralateral
side.2 <- c(1, 3, 5, 7, 11, 25:33, 43:51) # treated

pairedLM <- cbind(side.1, side.2)

# Atlases
setwd('./lm_data/Meshes/')
head_lowres <- head_mesh_spec1_dec <- get_dec_mesh('face')
setwd("../../")

# plot landmarks 
open3d(zoom = 0.75, windowRect = c(0, 0, 700, 700)) 
rgl::shade3d(head_lowres, color = "gray", alpha =0.9)
rgl::plot3d(atlas_head_lm, aspect = "iso", type = "s", size=1.2, col = "darkblue", add = T)
rgl::text3d(x = atlas_head_lm[,1],
            y = atlas_head_lm[,2],
            z = atlas_head_lm[,3],
            texts = c(1:dim(atlas_head_lm)[1]),
            cex = 1.5, offset = 0.5, pos = 1)
rgl::close3d()


#### 3 Analysis of bilateral symmetry ####
SYM_FGF <- bilat.symmetry(head_array, side = NULL, replicate = NULL, object.sym = TRUE, 
                          ind = dimnames(head_array)[[3]], land.pairs = pairedLM, iter = 999, seed = NULL, RRPP = TRUE)

summary(SYM_FGF)
cat("SYM_FGF", capture.output(summary(SYM_FGF)), 
    file="./output/SYM_FGF.txt", sep="\n", append=TRUE)


#### 3.1. Symmetric component ####
PCA_SYM_FGF <- gm.prcomp(SYM_FGF$symm.shape)

pal <- get_palette()

pdf("./figs/PCA_symmetric_component_noPC1.pdf", width = 7.5, height = 6)
plot(PCA_SYM_FGF, pch = 16, col = pal[as.numeric(classifiers$treatment)], cex = 0.8)
# text(PCA_SYM_FGF[["x"]][,1], PCA_SYM_FGF[["x"]][,2], dimnames(head_array)[[3]])
ordiellipse(PCA_SYM_FGF, classifiers$treatment, kind="sd",conf=0.95, border = pal,
            draw = "polygon", alpha = 0, lty = 1)
ordiellipse(PCA_SYM_FGF, classifiers$treatment, kind="se",conf=0.95, border = pal,
            draw = "polygon", alpha = 0, lty = 1, lwd = 4)
legend("bottomleft", pch = 16, col = pal, legend = levels(classifiers$treatment))
dev.off()


# this gives pvalue used in manuscript for symmetric shape change
ANOVA_ALL_sym <- procD.lm(SYM_FGF$symm.shape ~ classifiers$treatment, 
                          iter=999, RRPP=TRUE, print.progress = FALSE)

summary(ANOVA_ALL_sym)
summary(ANOVA_ALL_sym, test.type = "var")

treatment_ph <- pairwise(ANOVA_ALL_sym, groups = classifiers$treatment)
summary(treatment_ph)
summary(treatment_ph, test.type = "var", confidence = 0.95, stat.table = TRUE)

cat("ANOVA_SYM_FGF", capture.output(summary(ANOVA_ALL_sym)),
    file="./output/ANOVA_symmetric_component_FGF.txt", sep="\n", append=TRUE)


### 3.2. Asymmetric component ####
PCA_ASYM_FGF <- gm.prcomp(SYM_FGF$asymm.shape)

# this figure is used in manuscript
pdf("./figs/PCA_asymmetric_component_noPC1.pdf", width = 7.5, height = 6)
plot(PCA_ASYM_FGF, pch = 16, col = pal[as.numeric(classifiers$treatment)], cex = 0.8)
# text(PCA_ASYM_FGF[["x"]][,1], PCA_ASYM_FGF[["x"]][,2], dimnames(head_array)[[3]], cex=0.5)
ordiellipse(PCA_ASYM_FGF, classifiers$treatment, kind="sd",conf=0.95, border = pal,
            draw = "polygon", alpha = 0, lty = 1)
ordiellipse(PCA_ASYM_FGF, classifiers$treatment, kind="se",conf=0.95, border = pal,
            draw = "polygon", alpha = 0, lty = 1, lwd = 4)
legend("topright", pch = 16, col = pal, legend = levels(classifiers$treatment))
dev.off()

# this is pvalue used in manuscript or asymmetry
ANOVA_ASYM_FGF  <- procD.lm(SYM_FGF$asymm.shape ~ classifiers$treatment, 
                            iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_ASYM_FGF)
summary(ANOVA_ASYM_FGF, test.type = "var")

treatment_ph <- pairwise(ANOVA_ASYM_FGF, groups = classifiers$treatment)
summary(treatment_ph)
summary(treatment_ph, test.type = "var", confidence = 0.95, stat.table = TRUE)


cat("ANOVA_ASYM_FGF", capture.output(summary(ANOVA_ASYM_FGF)),
    file="./output/ANOVA_asymmetric_component_FGF.txt", sep="\n", append=TRUE)


#### 3.3. Fluctuating asymmetry component ####
PCA_FA_FGF <- gm.prcomp(gpagen(SYM_FGF$FA.component)$coords)

pdf("./figs/PCA_fluctuating_asymmetry_component.pdf", width = 7.5, height = 6)
plot(PCA_FA_FGF, pch = 16, col = pal[as.numeric(classifiers$treatment)], cex = 1.25)
ordiellipse(PCA_FA_FGF, classifiers$treatment, kind="sd",conf=0.95, border = pal,
            draw = "polygon", alpha = 0, lty = 1)
ordiellipse(PCA_FA_FGF, classifiers$treatment, kind="se",conf=0.95, border = pal,
            draw = "polygon", alpha = 0, lty = 1, lwd = 4)
legend("topright", pch = 16, col = pal, legend = levels(classifiers$treatment))
dev.off()


# variance covariance matrix?
vcm <- cov(t(two.d.array(SYM_FGF$asymm.shape)))
total_variance <- sum(diag(vcm))
DMSO_var <- sum(diag(vcm[1:25,1:25]))

U0_var <- sum(diag(vcm[78:101,78:101]))


ANOVA_FA_FGF  <- procD.lm(SYM_FGF$FA.component ~ classifiers$treatment, 
                          iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_FA_FGF)
summary(ANOVA_FA_FGF, test.type = "var")

treatment_ph <- pairwise(ANOVA_FA_FGF, groups = classifiers$treatment)
summary(treatment_ph)
summary(treatment_ph, test.type = "var", confidence = 0.95, stat.table = TRUE)

cat("ANOVA_FA_FGF", capture.output(summary(ANOVA_ASYM_FGF)),
    file="./output/ANOVA_Fluctuating_Asymmetry_FGF.txt", sep="\n", append=TRUE)


#### 3.4. DIRECTIONAL ASYMMETRY COMPONENT ####
summary(SYM_FGF)
SYM_FGF$DA.component # array with side 1 & side 2

# this isn't a component, it is the mean direction
# not needed for our analysis
ANOVA_DA_FGF  <- procD.lm(SYM_FGF$DA.component ~ classifiers$treatment, 
                          iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_DA_FGF)


#### 3.5. TO DO NEXT ####
#### 3.6. PCA ####
# Make morphs for each component (PCA)
# load mesh of the specimen closest to the mean shape
setwd('./lm_data/Meshes/')
face_mesh <- get_dec_mesh('face')
setwd("../../")

# get landmarks for each treatment group
Pdist <- ShapeDist(GPA_geomorph$coords, GPA_geomorph$consensus)
# make the dataframe for better ggplot plotting
gdf_head <- geomorph.data.frame(GPA_geomorph, treatment = classifiers$treatment, Pdist = Pdist)

new_ctrl_coords <- gdf_head$coords[,,which(classifiers$treatment == "DMSO")]
new_ctrl_mean_shape <- mshape(new_ctrl_coords)
new_ctrl_mesh <- tps3d(face_mesh, as.matrix(atlas_head_lm), new_ctrl_mean_shape, threads = 1)

open3d(zoom = 0.75, windowRect = c(0, 0, 700, 700)) 
rgl::shade3d(new_ctrl_mesh, color = "gray", alpha =0.9)
rgl::plot3d(new_ctrl_mean_shape, aspect = "iso", type = "s", size=1.2, col = "darkblue", add = T)
rgl::text3d(x = atlas_head_lm[,1],
            y = atlas_head_lm[,2],
            z = atlas_head_lm[,3],
            texts = c(1:dim(atlas_head_lm)[1]),
            cex = 1.5, offset = 0.5, pos = 1)
rgl::close3d()

# get warp for certain PC values rather than mean shapes
# symmetric
sym_PC1 = PCA_SYM_FGF$x[,1] # get pc1
sym_PC2 = PCA_SYM_FGF$x[,2] # get pc2

sym_PC1_ctrl <- tps3d(new_ctrl_mesh, as.matrix(new_ctrl_mean_shape), shape.predictor(GPA_geomorph$coords, x= sym_PC1, Intercept = FALSE, pred1 = .01)[[1]], threads=1)
sym_PC1_trt <- tps3d(new_ctrl_mesh, as.matrix(new_ctrl_mean_shape), shape.predictor(GPA_geomorph$coords, x= sym_PC1, Intercept = FALSE, pred1 = -.02)[[1]], threads=1)


sym_PC2_ctrl <- tps3d(new_ctrl_mesh, as.matrix(new_ctrl_mean_shape), shape.predictor(GPA_geomorph$coords, x= sym_PC2, Intercept = FALSE, pred1 = .025)[[1]], threads=1)
sym_PC2_trt <- tps3d(new_ctrl_mesh, as.matrix(new_ctrl_mean_shape), shape.predictor(GPA_geomorph$coords, x= sym_PC2, Intercept = FALSE, pred1 = -.025)[[1]], threads=1)

open3d(zoom = 0.75,  windowRect = c(0, 0, 1000, 700)) 
meshDist(sym_PC1_ctrl, sym_PC1_trt, from= -0.015, to= 0.015, rampcolors = c("blue", "white", "red"), sign = TRUE)
meshDist(sym_PC2_ctrl, sym_PC2_trt, from= -0.015, to= 0.015, rampcolors = c("blue", "white", "red"), sign = TRUE)
close3d()


#asymmetry
asym_PC1 = PCA_ASYM_FGF$x[,1] # get pc1
asym_PC2 = PCA_ASYM_FGF$x[,2] # get pc2

asym_PC1_ctrl <- tps3d(face_mesh, as.matrix(atlas_head_lm), shape.predictor(GPA_geomorph$coords, x= asym_PC1, Intercept = FALSE, pred1 = .05)[[1]], threads=1)
asym_PC1_trt <- tps3d(face_mesh, as.matrix(atlas_head_lm), shape.predictor(GPA_geomorph$coords, x= asym_PC1, Intercept = FALSE, pred1 = -.05)[[1]], threads=1)

asym_PC2_ctrl <- tps3d(face_mesh, as.matrix(atlas_head_lm), shape.predictor(GPA_geomorph$coords, x= asym_PC2, Intercept = FALSE, pred1 = .025)[[1]], threads=1)
asym_PC2_trt <- tps3d(face_mesh, as.matrix(atlas_head_lm), shape.predictor(GPA_geomorph$coords, x= asym_PC2, Intercept = FALSE, pred1 = -.025)[[1]], threads=1)

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700)) 
meshDist(asym_PC1_ctrl, asym_PC1_trt, from= -0.015, to= 0.015, rampcolors = c("blue", "white", "red"), sign = TRUE)
meshDist(asym_PC2_ctrl, asym_PC2_trt, from= -0.015, to= 0.015, rampcolors = c("blue", "white", "red"), sign = TRUE)
close3d()

# Make heatmap for DA left vs right
# Make table with proportion of variance explained by each component
# Remake PCA plots like Costello comparisons

#### 3.7. CVA ####
# Same as above




#### 4 Mirrored faces ####
# Mirror each side of the face, so double up on specs, run GMM on these new datasets, by treatment
# We will end up with these funny-looking specimens, 4 treat-side combos:
# (CTRL-contralateral, CTRL-treated, triple-contralateral, triple-treated)

# 4.1 Mirroring landmarks ####
# Mirror landmarks on both sides and generate two new arrays
# the non-symmetrical landmarks also remain the same

array_side1_mirrored <- GPA_geomorph$coords # contralateral
array_side2_mirrored <- GPA_geomorph$coords # treated

half_array_side1 <- sweep(GPA_geomorph$coords[side.1,,1], MARGIN = 2, c(-1,1,1), `*`) # mirror side 1 (right)

open3d()
plot3d(GPA_geomorph$coords[non.sym, , 1], col = "black", type = "s", aspect = "iso", 
       size = 1, add = TRUE, xlab = "x", ylab = "y", zlab = "z")
plot3d(GPA_geomorph$coords[side.1, , 1], col = "green", type = "s", aspect = "iso", 
       size = 1, add = TRUE, xlab = "x", ylab = "y", zlab = "z")
plot3d(half_array_side1, col = "blue", type = "s", aspect = "iso", 
       size = 1, add = TRUE, xlab = "x", ylab = "y", zlab = "z")
plot3d(GPA_geomorph$coords[side.2, , 1], col = "red", type = "s", aspect = "iso",
       size = 1, add = TRUE, xlab = "x", ylab = "y", zlab = "z")
rgl::close3d()


for (i in 1:dim(head_array)[3]){
  array_side1_mirrored[side.2,,i] <- sweep(GPA_geomorph$coords[side.1,,i],
                                           MARGIN = 2, c(-1,1,1), `*`)
  
  array_side2_mirrored[side.1,,i] <- sweep(GPA_geomorph$coords[side.2,,i],
                                           MARGIN = 2, c(-1,1,1), `*`)
}

setwd('./lm_data/Meshes/')
head_lowres <- head_mesh_spec1_dec <- get_dec_mesh('head')
setwd("../../")

# load("./lm_data/RGL_head_pos.rdata")
open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700))
plot3d(GPA_geomorph$coords[, , 1], col = "black", type = "s", aspect = "iso",
       size = 1, add = TRUE, xlab = "x", ylab = "y", zlab = "z")
plot3d(array_side1_mirrored[,,1], col = "chartreuse", type = "s", aspect = "iso",
       size = 0.75, add = TRUE, xlab = "x", ylab = "y", zlab = "z")
rgl::close3d()
# 
# # which(classifiers_mirrored$treatment == "triple")
# 
open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700))
# rgl::shade3d(head_lowres, color = "gray", alpha =1,specular = "black", add=T)
# plot3d(GPA_geomorph$coords[, , 14], col = "black", type = "s", aspect = "iso",
#       size = 1, add = TRUE, xlab = "x", ylab = "y", zlab = "z")
plot3d(array_side1_mirrored[,,14], col = "chartreuse", type = "s", aspect = "iso",
       size = 0.75, add = TRUE, xlab = "x", ylab = "y", zlab = "z")
rgl::text3d(x = array_side1_mirrored[,1,14],
            y = array_side1_mirrored[,2,14],
            z = array_side1_mirrored[,3,14],
            texts = c(1:dim(array_side1_mirrored)[1]),
            cex = 1.5, offset = 0.5, pos = 1)

# Change dimnames so we know which way they were mirrored
dimnames(array_side2_mirrored)[[3]] <- paste0("treated_", dimnames(array_side2_mirrored)[[3]]) # this one has the treated side mirrored
dimnames(array_side1_mirrored)[[3]] <- paste0("contralateral_", dimnames(array_side1_mirrored)[[3]]) # this one has contralateral side mirrored

# And join both arrays

mirrored_array <- abind(array_side1_mirrored, array_side2_mirrored, along = 3)

# Finally make new classifiers matrix

classifiers_mirrored <- rbind(classifiers, classifiers)
row.names(classifiers_mirrored) <- dimnames(mirrored_array)[[3]]
classifiers_mirrored$id <- dimnames(mirrored_array)[[3]]
classifiers_mirrored$treatment_mirror <- vector(mode = "character", length = dim(classifiers_mirrored)[1])
classifiers_mirrored$treatment_mirror[1:dim(classifiers)[1]] <- paste0("contra_", classifiers_mirrored$treatment[1:dim(classifiers)[1]])
classifiers_mirrored$treatment_mirror[(1+dim(classifiers)[1]):dim(classifiers_mirrored)[1]] <- paste0("treat_", classifiers_mirrored$treatment[(1+dim(classifiers)[1]):dim(classifiers_mirrored)[1]])

classifiers_mirrored$treatment_mirror <- as.factor(classifiers_mirrored$treatment_mirror)

classifiers_mirrored$treatment_mirror <- factor(classifiers_mirrored$treatment_mirror, levels = c('contra_DMSO', 'contra_U0126', 'contra_LY294002', 'contra_U73122', 'contra_Triple',
                                                                                                  'treat_DMSO', 'treat_U0126', 'treat_LY294002', 'treat_U73122', 'treat_Triple'))



# 4.4. Mirrored GPA & PCA ####
surface_semis <- c(34:51)
# GPA_mirrored_double <- geomorph::gpagen(A = mirrored_array*c(GPA_geomorph$Csize,GPA_geomorph$Csize), curves = as.matrix(curveslide_all), 
#                                         surfaces = surface_semis)
GPA_mirrored_double <- geomorph::gpagen(A = mirrored_array, curves = as.matrix(curveslide_all), 
                                        surfaces = surface_semis)

# GPA_mirrored_contra <- geomorph::gpagen(A = mirrored_array[,,1:dim(array_side1_mirrored)[3]]*GPA_geomorph$Csize, curves = as.matrix(curveslide_all), 
#                                       surfaces = surface_semis)
GPA_mirrored_contra <- geomorph::gpagen(A = mirrored_array[,,1:dim(array_side1_mirrored)[3]], curves = as.matrix(curveslide_all),
                                        surfaces = surface_semis)

# GPA_mirrored_treat <- geomorph::gpagen(A = mirrored_array[,,(1+dim(array_side1_mirrored)[3]):dim(mirrored_array)[3]]*GPA_geomorph$Csize, curves = as.matrix(curveslide_all), 
#                                        surfaces = surface_semis)
GPA_mirrored_treat <- geomorph::gpagen(A = mirrored_array[,,(1+dim(array_side1_mirrored)[3]):dim(mirrored_array)[3]], curves = as.matrix(curveslide_all),
                                       surfaces = surface_semis)

saveRDS(mirrored_array, "./lm_data/mirrored_array_both_sides.rds")
saveRDS(GPA_mirrored_double, "./lm_data/mirrored_gpa_both_sides.rds")
saveRDS(GPA_mirrored_contra, "./lm_data/mirrored_gpa_contralateral.rds")
saveRDS(GPA_mirrored_treat, "./lm_data/mirrored_gpa_treated.rds")
write.csv(classifiers_mirrored, "./lm_data/classifiers_mirrored_both_sides.csv")

# Both sides
PCA_both_sides <- gm.prcomp(GPA_mirrored_double$coords)

# Delete file if it exists
if (file.exists("./output/PCA_both_sides.txt")) {
  file.remove("./output/PCA_both_sides.txt")
}
cat("PCA shape variables raw - both sides mirrored", capture.output(summary(PCA_both_sides)), 
    file="./output/PCA_both_sides.txt", sep="\n", append=TRUE)

# plot PC1vPC2
levels(classifiers_mirrored$treatment_mirror)
pal_light <- get_palette(lights=TRUE)
# this plot is in manuscript
pdf("./figs/mirrored_PCA_treatment_sides_noPC1.pdf", width = 8.25, height = 6)
plot(PCA_both_sides, pch = 16, col = pal_light[as.numeric(classifiers_mirrored$treatment_mirror)], cex = 0.6)
ordiellipse(PCA_both_sides, classifiers_mirrored$treatment_mirror, kind="sd",conf=0.95, border = pal_light,
            draw = "polygon", alpha = 0, lty = 1)
ordiellipse(PCA_both_sides, classifiers_mirrored$treatment_mirror, kind="se",conf=0.95, border = pal_light,
            draw = "polygon", alpha = 50, lty = 1, lwd=4)
legend("bottomleft", pch = 16, col = pal_light, legend = levels(classifiers_mirrored$treatment_mirror))
dev.off()


Pdist <- ShapeDist(GPA_mirrored_double$coords, GPA_mirrored_double$consensus)
t <- geomorph.data.frame(GPA_mirrored_double, treatment = classifiers_mirrored$treatment_mirror, Pdist = Pdist)

ANOVA_both_mirrored  <- procD.lm(coords ~ treatment, data=t, 
                          iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_both_mirrored)
ANOVA_both_mirrored_pw <- pairwise(ANOVA_both_mirrored, groups = t$treatment)
# these are pvalues used in manuscript
summary(ANOVA_both_mirrored_pw)


# perform CVA
# need Morpho package for the cva
cva_head_mirrored <- CVA(GPA_mirrored_double$coords, classifiers_mirrored$treatment_mirror)

pdf("./figs/cva_mirrored_mean_shape.pdf", width = 8.25, height = 6)
plot(cva_head_mirrored$CVscores[,1:2], col=pal_light[as.numeric(classifiers_mirrored$treatment_mirror)], pch=16, typ="p",asp=1)
# text(cva_head$CVscores, as.character(classifiers$treatment), col=as.numeric(classifiers$treatment), cex=.7)
# plot(cva_head$CVscores[,c(3,4)], bg=classifiers$treatment, pch=21, typ="p",asp=1)


# https://rdrr.io/cran/Morpho/man/CVA.html
# add chull (merge groups)
# for(jj in 1:length(levels(classifiers$treatment))){
#       ii=levels(classifiers$treatment)[jj]
#   kk=chull(cva_head$CVscores[classifiers$treatment==ii,1:2])
#   lines(cva_head$CVscores[classifiers$treatment==ii,1][c(kk, kk[1])],
#   cva_head$CVscores[classifiers$treatment==ii,2][c(kk, kk[1])], col=pal[jj])
#   }

# add 95% ellipses
ordiellipse(cva_head_mirrored$CVscores, classifiers_mirrored$treatment_mirror, kind="se",conf=0.95, border = pal_light,
            draw = "polygon", alpha = 0, lty = 1, lwd = 4)
ordiellipse(cva_head_mirrored$CVscores, classifiers_mirrored$treatment_mirror, kind="sd",conf=0.95, border = pal_light,
            draw = "polygon", alpha = .4, lty = 1, lwd = 1)
dev.off()




# 4.5. PCA contralateral side ####
PCA_contra <- gm.prcomp(GPA_mirrored_contra$coords)
# Delete file if it exists
if (file.exists("./output/PCA_contra.txt")) {
  file.remove("./output/PCA_contra.txt")
}
cat("PCA shape variables raw - contralateral side mirrored", capture.output(summary(PCA_contra)), 
    file="./output/PCA_contra.txt", sep="\n", append=TRUE)

pal_light_only <- pal_light[1:5]
pdf("./figs/mirrored_PCA_contralateral_treatment_noPC1.pdf", width = 8.25, height = 6)
plot(PCA_contra, pch = 16, col = pal_light_only[as.numeric(classifiers$treatment)], cex = .6)
ordiellipse(PCA_contra, classifiers$treatment, kind="sd",conf=0.95, border = pal_light_only,
            draw = "polygon", alpha = 0, lty = 1)
ordiellipse(PCA_contra, classifiers$treatment, kind="se",conf=0.95, border = pal_light_only,
            draw = "polygon", alpha = 50, lty = 1, lwd=4)
legend("bottomright", pch = 16, col = pal_light_only, legend = levels(classifiers$treatment))
dev.off()

t_contra <- geomorph.data.frame(GPA_mirrored_contra, treatment = classifiers$treatment)
ANOVA_contra_mirrored  <- procD.lm(coords ~ treatment, data=t_contra, 
                                 iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_contra_mirrored)
ANOVA_contra_mirrored_pw <- pairwise(ANOVA_contra_mirrored, groups = t_contra$treatment)
summary(ANOVA_contra_mirrored_pw)


#### 4.6. PCA treated side #####
PCA_treat <- gm.prcomp(GPA_mirrored_treat$coords)
# Delete file if it exists
if (file.exists("./output/PCA_treat.txt")) {
  file.remove("./output/PCA_treat.txt")
}
cat("PCA shape variables raw - treated side mirrored", capture.output(summary(PCA_treat)), 
    file="./output/PCA_treat.txt", sep="\n", append=TRUE)

pal_dark_only <- pal_light[6:10]
pdf("./figs/mirrored_PCA_treated_treatment_noPC1.pdf", width = 8.25, height = 6)
plot(PCA_treat, pch = 16, col = pal_dark_only[as.numeric(classifiers$treatment)], cex = .6)
ordiellipse(PCA_contra, classifiers$treatment, kind="sd",conf=0.95, border = pal_dark_only,
            draw = "polygon", alpha = 0, lty = 1)
ordiellipse(PCA_treat, classifiers$treatment, kind="se",conf=0.95, border = pal_dark_only,
            draw = "polygon", alpha = 50, lty = 1, lwd=4)
legend("bottomright", pch = 16, col = pal, legend = levels(classifiers$treatment))
dev.off()

t_treat <- geomorph.data.frame(GPA_mirrored_treat, treatment = classifiers$treatment)

ANOVA_treat_mirrored  <- procD.lm(coords ~ treatment, data=t_treat, 
                                 iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_treat_mirrored)
ANOVA_treat_mirrored_pw <- pairwise(ANOVA_treat_mirrored, groups = t_treat$treatment)
summary(ANOVA_treat_mirrored_pw)

#### 4.7. Procrustes distance ####
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
ggplot_df$treatment <- as.factor(ggplot_df$treatment)
# levels(ggplot_df$treatment) <- c('DMSO', 'U0126', 'LY294002', 'U73122', 'Triple')
ggplot_df$treatment_mirror <- as.factor(ggplot_df$treatment_mirror)
# levels(ggplot_df$treatment_mirror) <- c( 'treat_DMSO', 'treat_U0126', 'treat_LY294002', 'treat_U73122', 'treat_Triple',
                                       #  'contra_DMSO', 'contra_U0126',  'contra_LY294002', 'contra_U73122', 'contra_Triple')
ggplot_df$Pdist <- as.numeric(as.character(ggplot_df$Pdist))

pdf("./figs/mirrored_Pdist_treatment.pdf", width = 6.5, height = 6.5)
ggplot(ggplot_df, aes(Pdist, fill = treatment_mirror)) +
  scale_fill_manual(values = pal_light) + geom_density(alpha = 0.5) + 
  ggtitle("Procrustes distances head shape - mirrored") + xlab("Procrustes distance") + ylab("Relative density") +
  theme(plot.title = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"), legend.text=element_text(size=10), 
        legend.title=element_text(size=12, face = "bold"))
dev.off()


#### 5 MORPHS MIRRORING ####
load("./lm_data/RGL_heat_head_pos.rdata")
load("./lm_data/RGL_head_pos.rdata")
# frontal <- par3d()$userMatrix
# 
# rgl.close()
# 
# save(dorsal, frontal, file = "./figs/RGL_head_heatmaps_pos.rdata")

# load mesh of the specimen closest to the mean shape
# head_mesh <- vcgImport("./lm_data/Meshes/ATLAS_chick_ctr_face_decimated.ply") 
# The mesh should only be surface, and be just have the frontal side here
# head_lowres <- vcgQEdecim(head_mesh, percent = 0.15)

new_ctrl_mirrored_coords <- gdf_mirrored$coords[,,which(classifiers$treatment == "DMSO")]
new_ctrl_mirrored_mean_shape <- mshape(new_ctrl_mirrored_coords)
new_ctrl_mirrored_mesh <- tps3d(face_mesh, as.matrix(atlas_head_lm), new_ctrl_mirrored_mean_shape, threads = 1)

rgl::shade3d(new_ctrl_mirrored_mesh, color="gray", alpha=1, specular="black")
rgl::close3d()


levels(gdf_mirrored$treatment_mirror)

contra_DMSO_mean_shape <- mshape(gdf_mirrored$coords[,,which(gdf_mirrored$treatment_mirror == "contra_DMSO")])
treat_DMSO_mean_shape <- mshape(gdf_mirrored$coords[,,which(gdf_mirrored$treatment_mirror == "treat_DMSO")])

contra_LY_mean_shape <- mshape(gdf_mirrored$coords[,,which(gdf_mirrored$treatment_mirror == "contra_LY294002")])
treat_LY_mean_shape <- mshape(gdf_mirrored$coords[,,which(gdf_mirrored$treatment_mirror == "treat_LY294002")])

contra_Triple_mean_shape <- mshape(gdf_mirrored$coords[,,which(gdf_mirrored$treatment_mirror == "contra_Triple")])
treat_Triple_mean_shape <- mshape(gdf_mirrored$coords[,,which(gdf_mirrored$treatment_mirror == "treat_Triple")])


# contra_triple_mean_shape <- mshape(gdf_mirrored$coords[,,which(gdf_mirrored$treatment_mirror == "contra_triple")])
# treat_triple_mean_shape <- mshape(gdf_mirrored$coords[,,which(gdf_mirrored$treatment_mirror == "treat_triple")])

model_head_lm <- GPA_mirrored_contra$coords[,, which(dimnames(GPA_mirrored_contra$coords)[[3]] == "treated_chick_ctr_23")] # used for warping mesh
alt_lm <- gdf_mirrored$coords[,, which(dimnames(gdf_mirrored$coords)[[3]] == "treated_chick_ctr_23")]
# alt_lm <- gdf_mirrored$coords[,, which(dimnames(gdf_mirrored$coords)[[3]] == "treated_treated_chick_ctr_23")]

# Create morphed meshes
# model_DMSO_mesh <- tps3d(head_mesh, as.matrix(og_lm), atlas_head_lm, threads = 1)
# model_DMSO_mesh <- tps3d(head_mesh, as.matrix(og_lm), alt_lm, threads = 1)

contra_DMSO_mesh <- tps3d(new_ctrl_mirrored_mesh, as.matrix(new_ctrl_mirrored_mean_shape), contra_DMSO_mean_shape, threads = 1)
treat_DMSO_mesh <- tps3d(new_ctrl_mirrored_mesh, as.matrix(new_ctrl_mirrored_mean_shape), treat_DMSO_mean_shape, threads = 1)

contra_LY_mesh <- tps3d(new_ctrl_mirrored_mesh, as.matrix(new_ctrl_mirrored_mean_shape), contra_LY_mean_shape, threads = 1)
treat_LY_mesh <- tps3d(new_ctrl_mirrored_mesh, as.matrix(new_ctrl_mirrored_mean_shape), treat_LY_mean_shape, threads = 1)

contra_Triple_mesh <- tps3d(new_ctrl_mirrored_mesh, as.matrix(new_ctrl_mirrored_mean_shape), contra_Triple_mean_shape, threads = 1)
treat_Triple_mesh <- tps3d(new_ctrl_mirrored_mesh, as.matrix(new_ctrl_mirrored_mean_shape), treat_Triple_mean_shape, threads = 1)

PC1 = PCA_both_sides$x[,1] # get pc1
PC2 = PCA_both_sides$x[,2] # get pc2

PC1_ctrl <- tps3d(new_ctrl_mirrored_mesh, as.matrix(new_ctrl_mirrored_mean_shape), shape.predictor(GPA_mirrored_double$coords, x= PC1, Intercept = FALSE, pred1 = -.03)[[1]], threads=1)
PC1_trt <-  tps3d(new_ctrl_mirrored_mesh, as.matrix(new_ctrl_mirrored_mean_shape), shape.predictor(GPA_mirrored_double$coords, x= PC1, Intercept = FALSE, pred1 = .05)[[1]], threads=1)

PC2_ctrl <- tps3d(new_ctrl_mirrored_mesh, as.matrix(new_ctrl_mirrored_mean_shape), shape.predictor(GPA_mirrored_double$coords, x= PC2, Intercept = FALSE, pred1 = -.02)[[1]], threads=1)
PC2_trt <-  tps3d(new_ctrl_mirrored_mesh, as.matrix(new_ctrl_mirrored_mean_shape), shape.predictor(GPA_mirrored_double$coords, x= PC2, Intercept = FALSE, pred1 = .03)[[1]], threads=1)

# PC for contra only
PC2 = PCA_contra$x[,2] # get pc2

PC2_contra_ctrl <- tps3d(new_ctrl_mirrored_mesh, as.matrix(new_ctrl_mirrored_mean_shape), shape.predictor(GPA_mirrored_contra$coords, x= PC2, Intercept = FALSE, pred1 = -.015)[[1]], threads=1)
PC2_contra_trt <-  tps3d(new_ctrl_mirrored_mesh, as.matrix(new_ctrl_mirrored_mean_shape), shape.predictor(GPA_mirrored_contra$coords, x= PC2, Intercept = FALSE, pred1 = .02)[[1]], threads=1)



# create CVA mesh
# cvvis 1-4 is the cv number (like PC num)
# groupmean is the mean shape for that group in all cvs (very similar to mean shape above)
ctrl_mean_shape_CV1 <- 5*matrix(cva_head_mirrored$CVvis[,1], nrow(cva_head_mirrored$groupmeans), 3) + cva_head_mirrored$Grandm
# ctrl_mean_shape_CV1 <- matrix(cva_head$CVvis[,1], nrow(cva_head$groupmeans), 3) + cva_head$groupmeans[,,1]
ctrl_mesh_CV1 <- tps3d(new_ctrl_mirrored_mesh, as.matrix(new_ctrl_mirrored_mean_shape), ctrl_mean_shape_CV1, threads = 1)

triple_mean_shape_CV1 <- -3*matrix(cva_head_mirrored$CVvis[,1], nrow(cva_head_mirrored$groupmeans), 3) + cva_head_mirrored$Grandm
# triple_mean_shape_CV1 <- matrix(cva_head$CVvis[,1], nrow(cva_head$groupmeans), 3) + cva_head$groupmeans[,,5]
triple_mesh_CV1 <- tps3d(new_ctrl_mirrored_mesh, as.matrix(new_ctrl_mirrored_mean_shape), triple_mean_shape_CV1, threads = 1)

trt_mean_shape_CV2 <- 4*matrix(cva_head_mirrored$CVvis[,2], nrow(cva_head_mirrored$groupmeans), 3) + cva_head_mirrored$Grandm
# ctrl_mean_shape_CV1 <- matrix(cva_head$CVvis[,1], nrow(cva_head$groupmeans), 3) + cva_head$groupmeans[,,1]
trt_mesh_CV2 <- tps3d(new_ctrl_mirrored_mesh, as.matrix(new_ctrl_mirrored_mean_shape), trt_mean_shape_CV2, threads = 1)

contra_mean_shape_CV2 <- -2*matrix(cva_head_mirrored$CVvis[,2], nrow(cva_head_mirrored$groupmeans), 3) + cva_head_mirrored$Grandm
# triple_mean_shape_CV1 <- matrix(cva_head$CVvis[,1], nrow(cva_head$groupmeans), 3) + cva_head$groupmeans[,,5]
contra_mesh_CV2 <- tps3d(new_ctrl_mirrored_mesh, as.matrix(new_ctrl_mirrored_mean_shape), contra_mean_shape_CV2, threads = 1)




# Plot morphs
open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = frontal)
rgl::shade3d(treat_DMSO_mesh, color="gray", alpha=1, specular="black")
# rgl::plot3d(rotonto(as.matrix(og_lm), as.matrix(model_head_lm), scale=TRUE)$yrot, col = "black", type = "s", aspect = "iso", size = 1, add = TRUE)
rgl::plot3d(alt_lm[side.1,], col = "black", type = "s", aspect = "iso", size = 1, add = TRUE)
rgl::plot3d(alt_lm[non.sym,], col = "blue", type = "s", aspect = "iso", size = 1, add = TRUE)
rgl::plot3d(alt_lm[side.2,], col = "red", type = "s", aspect = "iso", size = 1, add = TRUE)
rgl::plot3d(alt_lm, col = "green", type = "s", aspect = "iso", size = 1, add = TRUE)
axes3d(tick=T)

rgl.snapshot("./figs/Morph_treated_DMSO_frontal_head.png", top = TRUE)
writePLY("./output/contralateral_DMSO_mesh.ply")
rgl::close3d()

open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = frontal)
shade3d(contra_DMSO_mesh, color="gray", alpha=0.5, specular='black')
rgl.snapshot("./figs/Morph_contralateral_DMSO_frontal_head.png", top = TRUE)
writePLY("./output/contralateral_DMSO_mesh.ply")
rgl::close3d()

open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = frontal)
shade3d(treat_DMSO_mesh, color="red", alpha=0.5, specular='black')
rgl.snapshot("./figs/Morph_treated_DMSO_frontal_head.png", top = TRUE)
writePLY("./output/treated_DMSO_mesh.ply")
rgl.close()

open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = frontal)
shade3d(contra_LY_mesh, color="gray", alpha=0.5, specular='black')
rgl.snapshot("./figs/Morph_contralateral_LY_frontal_head.png", top = TRUE)
writePLY("./output/contralateral_triple_mesh.ply")
rgl.close()

open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = frontal)
shade3d(treat_LY_mesh, color="gray", alpha=0.5, specular='black')
rgl::plot3d(test[,,1], col = "green", type = "s", aspect = "iso", size = 1, add = TRUE)
rgl.snapshot("./figs/Morph_treated_LY_frontal_head.png", top = TRUE)
writePLY("./output/treated_triple_mesh.ply")
rgl::close3d()


# HEATMAPS
# plot the two first PC heatmaps for mean shape - these are used in manuscript
open3d(zoom = 0.75,  windowRect = c(0, 0, 1000, 700)) 
pdf("./figs/heatmap_mirrored_PC1_-pt03_to_pt05_legend.pdf", width = 2.5, height = 6.5)
meshDist(contra_DMSO_mesh, treat_DMSO_mesh, from= -.015, to= 0.015, rampcolors = c("blue", "white", "red"), sign = TRUE)
meshDist(contra_DMSO_mesh, contra_LY_mesh, from= -.015, to= 0.015, rampcolors = c("blue", "white", "red"), sign = TRUE)
meshDist(contra_DMSO_mesh, contra_Triple_mesh, from= -.015, to= 0.02, rampcolors = c("blue", "white", "red"), sign = TRUE)

meshDist(treat_DMSO_mesh, treat_LY_mesh, from= -.015, to= 0.015, rampcolors = c("blue", "white", "red"), sign = TRUE)

meshDist(ctrl_mesh_CV1, triple_mesh_CV1, rampcolors = c("blue", "white", "red"), sign = TRUE)
meshDist(contra_mesh_CV2, trt_mesh_CV2, rampcolors = c("blue", "white", "red"), sign = TRUE)

meshDist(PC1_ctrl, PC1_trt, from= -.015, to= 0.02, rampcolors = c("blue", "white", "red"), sign = TRUE)
meshDist(PC2_contra_ctrl, PC2_contra_trt, rampcolors = c("blue", "white", "red"), sign = TRUE)
meshDist(PC2_ctrl, PC2_trt, from= -.015, to= 0.015,  rampcolors = c("blue", "white", "red"), sign = TRUE)

meshDist(contra_DMSO_mesh, PC2_trt,  rampcolors = c("blue", "white", "red"), sign = TRUE)
rgl.snapshot("./figs/heatmap_DMSO-contra_to_Triple-contra.png", top = TRUE) # this one captures 3d output
rgl::close3d() # this one captures the heatmap legend as pdf
dev.off()

open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = frontal)
x11(width=1.7, height=8)
meshDist(treat_DMSO_mesh, treat_triple_mesh, rampcolors = c("blue", "white", "red"), sign = TRUE)
rgl.snapshot("./figs/Heatmap_treated_ctrl-vs-triple_mirrored_frontal.png", top = TRUE)
dev.print(pdf, "./figs/Heatmap_treated_ctrl-vs-triple_mirrored_scale.pdf", width=2, height=9.5)
dev.off()
clear3d()
rgl.close()

open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = frontal)
x11(width=1.7, height=8)
meshDist(contra_DMSO_mesh, treat_DMSO_mesh, rampcolors = c("blue", "white", "red"), sign = TRUE)
rgl.snapshot("./figs/Heatmap_contra-vs-treat_DMSO_mirrored_frontal.png", top = TRUE)
dev.print(pdf, "./figs/Heatmap_contra-vs-treat_DMSO_mirrored_scale.pdf", width=2, height=9.5)
dev.off()
clear3d()
rgl.close()


open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = frontal)
x11(width=1.7, height=8)
meshDist(treat_triple_mesh, contra_triple_mesh, rampcolors = c("blue", "white", "red"), sign = TRUE)
rgl.snapshot("./figs/Heatmap_contra-vs-treat_triple_mirrored_frontal.png", top = TRUE)
dev.print(pdf, "./figs/Heatmap_contra-vs-treat_triple_mirrored_scale.pdf", width=2, height=9.5)
dev.off()
clear3d()
rgl.close()



#### 5.2. PC morphs mirroring ####
# Have a look at the heatmap, and change colour if needed
open3d(zoom=0.75, userMatrix = frontal, windowRect= c(0,0,1000,700))
x11(width=1.7, height=8)
PC_min <- tps3d(head_lowres, as.matrix(atlas_head_lm), PCA_both_sides$shapes[[1]]$min, threads = 1)
PC_max <- tps3d(head_lowres, as.matrix(atlas_head_lm), PCA_both_sides$shapes[[1]]$max, threads = 1)

meshDist(PC_min, PC_max, rampcolors = c("darkblue", "blue", "white", "red", "darkred"), sign = TRUE)

summary(PCA_both_sides)
n_dimensions <- 10 # number of PCs to include in the figure, decide depending on % of variance explained in PCs

# This loop automatically positions faces in frontal only for the 10 PCs
for (i in 1:n_dimensions){
  PC_min <- tps3d(head_lowres, as.matrix(atlas_head_lm), PCA_both_sides$shapes[[i]]$min, threads = 1)
  PC_max <- tps3d(head_lowres, as.matrix(atlas_head_lm), PCA_both_sides$shapes[[i]]$max, threads = 1)
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



# 6. Face integration CTRL vs Treatments ####

?integration.test

non.sym <- c(9, 10, 13:15)
side.2 <- c(2,4, 6, 8, 12, 16:24, 34:42) # LEFT
side.1 <- c(1, 3, 5, 7, 11, 25:33, 43:51) # RIGHT

side <- vector(mode = "character", length = 51)
side[side.2] <- "contralateral"
side[side.1] <- "atreated" # they will go in alphabetical order, we want treated to be 'x'
side <- side[-non.sym]

# Integration face all
face_integration <- integration.test(head_array[-non.sym,,], partition.gp = side, iter = 999)
summary(face_integration) # Test summary
p<- plot(face_integration) # PLS plot
make_ggplot(p)
plotx <- p$plot.args$x
ploty <- p$plot.args$y
plot_df <- stack(plotx)
plot_df <- plot_df %>% rename(x = values)
plot_df <- plot_df %>% mutate(y = ploty)
plot_df <- plot_df %>% mutate(treatment = classifiers$treatment)

pdf("./figs/integration_pls_noPC1.pdf", width = 7.5, height = 6)
ggplot(plot_df, aes(x=x, y=y, color = treatment)) +
  geom_point(shape=16) +
  scale_fill_manual(values=get_palette(FALSE)) +
  geom_smooth(method=lm, se=FALSE) +
  geom_smooth(method=lm, se=FALSE,aes(group=1), color='black') +
  coord_fixed()
dev.off()

# Compare integration of the face between treatment and control
face_integration_CTRL <- integration.test(GPA_geomorph$coords[-non.sym, , which(classifiers$treatment == "DMSO")], 
                                          partition.gp = side, iter = 999)

face_integration_Triple <- integration.test(GPA_geomorph$coords[-non.sym, , which(classifiers$treatment == "Triple")], 
                                               partition.gp = side, iter = 999)
face_integration_U0 <- integration.test(GPA_geomorph$coords[-non.sym, , which(classifiers$treatment == "U0126")], 
                                            partition.gp = side, iter = 999)
face_integration_LY <- integration.test(GPA_geomorph$coords[-non.sym, , which(classifiers$treatment == "LY294002")], 
                                            partition.gp = side, iter = 999)
face_integration_U73 <- integration.test(GPA_geomorph$coords[-non.sym, , which(classifiers$treatment == "U73122")], 
                                            partition.gp = side, iter = 999)

# these p-values are used in manuscript
summary(face_integration_CTRL) # Test summary
summary(face_integration_Triple) # Test summary
plot(face_integration_CTRL) # PLS plot
plot(face_integration_Triple)
summary(face_integration_U0) # Test summary
summary(face_integration_LY) # Test summary
summary(face_integration_U73) # Test summary

PLS_comparison_Triple <- compare.pls(CTRL = face_integration_CTRL, TRIPLE = face_integration_Triple)
summary(PLS_comparison_Triple)
PLS_comparison_U0 <- compare.pls(CTRL = face_integration_CTRL, U0126 = face_integration_U0)
summary(PLS_comparison_U0)
PLS_comparison_LY <- compare.pls(CTRL = face_integration_CTRL, LY294002 = face_integration_LY)
summary(PLS_comparison_LY)
PLS_comparison_U73 <- compare.pls(CTRL = face_integration_CTRL, U73122 = face_integration_U73)
summary(PLS_comparison_U73)

# Delete file if it exists
if (file.exists("./output/Integration_face_comparisons_CTRL_TRIPLE_May2023.txt")) {
  file.remove("./output/Integration_face_comparisons_CTRL_TRIPLE_May2023.txt")
}
cat("INTEGRATION on the face comparison - CTRL vs TREATMENT", capture.output(summary(PLS_comparison)), 
    file="./output/Integration_face_comparisons_CTRL_TRIPLE_May2023.txt", sep="\n", append=TRUE)

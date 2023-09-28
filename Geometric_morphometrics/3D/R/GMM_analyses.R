#### 0. Load R packages ####
# install.packages('devtools')
# library(devtools)
# install_github("marta-vidalgarcia/morpho.tools.GM", force = TRUE)
# install.packages(c('rgl', 'geomorph' 'devtools', 'Morpho', 'Rvcg', 'magick', 'Evomorph', 'ggplot2', 'vegan', 'factoextra', 'gt'))

library(rgl)
library(geomorph)
library(morpho.tools.GM)
# install_github("marta-vidalgarcia/mesh_process")
# library(mesh_process)
library(Morpho)
library(Rvcg)
library(magick)
library(Evomorph)
library(ggplot2)
library(vegan)
# install_github("vqv/ggbiplot")
library(ggbiplot)
library(factoextra)
library(gt)
library(abind)

#### 0 Helpers ####
#### 0.1 Landmark loading helpers ####
fcsv_number_check <- function(dir, pattern) {
  fcsv_list <- dir(dir, pattern = pattern)
  n_land <- vector("numeric", length = length(fcsv_list))
  for (i in 1:length(fcsv_list)) {
    n_land[i] <- length(count.fields(fcsv_list[i]))
  }
  checker_df <- do.call(rbind.data.frame, Map('c', fcsv_list, n_land))
  return(checker_df)
}


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


plot_3d_LMs <- function(LMs, color) {
  # this will plot the landmarks in the specifid color
  # need to already have a rgl window open
  rgl::plot3d(LMs[,,1], aspect = "iso", type = "s", size=.5, col = color, add = T)
  rgl::text3d(x = LMs[,1,1],
              y = LMs[,2,1],
              z = LMs[,3,1],
              texts = c(1:dim(LMs)[1]),
              cex = 1.5, offset = 0.5, pos = 1)
}


#### Main ####
#### 1.0 Setting up the DF ####
# R doesn't have a good way to get the path of where this file is located, unfortunately
# if you are running this code in Rstudio, try this:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd() #check our working directory
setwd("../../../data/Morphology/3D")
getwd() #check our working directory
# if you aren't using rstudio, use the command setwd() 
# and point it to the data/Morphology/3D directory


# Then save this script inside the R folder
# Copy landmark data to the Prop_LMs folder
# Copy atlases inside the data/atlas folder. We will need:
# PLY of head, endocast & mandible. 
# CURVESLIDE FILES for endocast, head, mandible
# ATLAS TAG file for head, endocast, mandible
# Copy classifiers file inside data and make sure it is a csv and there are no issues


#### 1.1 Load LM data and classifiers ####
classifiers_unord  <- read.csv("./lm_data/classifiers.csv", header = TRUE)
classifiers_unord$treatment <- as.factor(classifiers_unord$treatment)
row.names(classifiers_unord) <- classifiers_unord$id


#### 1.2 Import landmark data ####
# 1. Import all fcsv files into separate arrays (all specimens for east LM set)
# Traditional landmarks
setwd("./lm_data/Landmarks/") # directory with all the fcsv files in it
LMs <- fcsv2array(string_del = "_Fiducials") # cleanup filenames, from morphotools library

# match the landmark sample names to the classifiers csv
classifiers <- classifiers_unord[match(dimnames(LMs)[[3]], row.names(classifiers_unord)),]


# Curve semilandmarks
setwd("../Semi_Curves/")
# finds all files with the pattern in the name
curve_semis_center <- fcsv2array(pattern = "*center*", string_del = "_center_semi-curve")
curve_semis_L <- fcsv2array(pattern = "*_L_semi-curve*", string_del = "_L_semi-curve")
curve_semis_R <- fcsv2array(pattern = "*_R_semi-curve*", string_del = "_R_semi-curve")


# Surface semilandmarks grid
setwd("../Semi_Points/")
surf_semis_L <- fcsv2array(pattern = "*_L_semi-lm*", string_del = "_L_semi-lm")
surf_semis_R <- fcsv2array(pattern = "*_R_semi-lm*", string_del = "_R_semi-lm")

setwd("../../") # back to '3D' dir

#### 1.3 Combine LM data ####
# Combine all arrays into a single one, in a particular order that makes sense for our GMM
# Centre curve
# Delete positions #1 & #5 in curve_semis_center as they are traditional landmarks
curve_semis_center <- curve_semis_center[-c(1,5),,]
# combine LMs, should have an rray of 51 landmarks and 54 samples
head_array <- abind(LMs, curve_semis_center, curve_semis_L, curve_semis_R, surf_semis_L, 
                    surf_semis_R, along = 1)


#### 1.4 Curveslide ####
# Generate a curveslide file for the GPA that tells us how the curve semis have to slide (constrained)
# this is the same as what we had to do in the 2D code
# curve semis positions, based on the abind for the head_array, above
semis_center <- (dim(LMs)[1]+1):(dim(LMs)[1]+dim(curve_semis_center)[1])
semis_L <- (dim(LMs)[1]+dim(curve_semis_center)[1]+1):(dim(LMs)[1]+dim(curve_semis_center)[1]+dim(curve_semis_L)[1])
semis_R <- (dim(LMs)[1]+dim(curve_semis_center)[1]+dim(curve_semis_L)[1]+1):(dim(LMs)[1]+dim(curve_semis_center)[1]+dim(curve_semis_L)[1]+dim(curve_semis_R)[1])

# Plot the landmarks to check placement of semis
# this loads in a simplified mesh of a control embryo
# it can take a long time if it's the first time you've run it
setwd("./lm_data/Meshes/")
head_mesh_spec1_dec <- get_dec_mesh('head')
setwd("../../")

# make a 3D plot
open3d(zoom = 0.75, windowRect = c(0, 0, 700, 700)) 

# plot the decimated head mesh
rgl::shade3d(head_mesh_spec1_dec, color = "gray", alpha =0.9)
# plot the landmarks in blue
plot_3d_LMs(LMs, 'darkblue')
# plot the center curve in orange
plot_3d_LMs(curve_semis_center, 'orange')
# plot the right curve in orange
plot_3d_LMs(curve_semis_R, 'orange')
# plot the left curve in orange
plot_3d_LMs(curve_semis_L, 'orange')

rgl::close3d() # close plot


# Center curve - sliding landmarks
# We need to create a curveslide matrix to know where from to where to the semis slide
# landmarks 9 & 10 will be treated as fixed landmarks!
# semi1 slides between LM9 & semi2
# semi2 slides between semi1 & semi3
# semi3 slides between semi2 & LM10
curve_c_left <- c(9, semis_center[c(1,2)])h
curve_c_right <- c(semis_center[c(2,3)], 10)
curveslide_c <- cbind(curve_c_left, semis_center, curve_c_right)

# Right and left curves - positions from top to bottoms
# Similar to the center sliders, the LNP landmarks 16 & 24 will be treated as fixed landmarks!
curve_L_left <- semis_L[c(1:(length(semis_L)-2))] 
curve_L_right <- semis_L[c(3:length(semis_L))]
curve_L_center <- semis_L[c(2:(length(semis_L)-1))]
curveslide_L <- cbind(curve_L_left, curve_L_center, curve_L_right)

# Similar to the center sliders, the LNP landmarks 25 & 33 will be treated as fixed landmarks!
curve_R_left <- semis_R[c(1:(length(semis_R)-2))] 
curve_R_right <- semis_R[c(3:length(semis_R))]
curve_R_center <- semis_R[c(2:(length(semis_R)-1))]
curveslide_R <- cbind(curve_R_left, curve_R_center, curve_R_right)


# all our curveslide matrices
curveslide_list <- lapply(ls(pattern = "curveslide*"), get)
curveslide_all <- do.call(rbind, curveslide_list)
curveslide_all <- as.data.frame(curveslide_all)
colnames(curveslide_all) <- c("left", "sliding", "right")
write.csv(curveslide_all, "./lm_data/curveslide.csv", row.names = FALSE)


#### 1.5 Surface semi-landmarks ####
head_surface.lm <- (dim(LMs)[1]+dim(curve_semis_center)[1]+dim(curve_semis_L)[1]+dim(curve_semis_R)[1]+1):dim(head_array)[1]


#### 2 GPA and plots ####
#### 2.1 GPA ####
# Use this to find out who is the closest to the shape sphere centroid
GPA_head <- geomorph::gpagen(A = head_array, curves = as.matrix(curveslide_all), 
                             surfaces = head_surface.lm)
# check for outliers
# looks fine, outlier is just heavily affected
outlier <- plotOutliers_percentile(A = GPA_head$coords, percentile = 0.99, save.plot = FALSE)
min(outlier$All_Proc_d$`Proc. d. from mean`) # actually I got lazy & did it this easier way, same concept as above though
row.names(outlier$All_Proc_d[which(outlier$All_Proc_d$`Proc. d. from mean` == min(outlier$All_Proc_d$`Proc. d. from mean`)),])

# save progress
saveRDS(head_array, "./lm_data/Head_LM_array_FGF_embryos.rds")
saveRDS(GPA_head, "./lm_data/GPA_FGF_embryos.rds")


#### 2.2 Atlas head and plots ####
setwd('./lm_data/Meshes/')
# in this case I have already made a cleaned up mesh of just the face
# I used chick_ctr_23 as the mesh, so it will be out 'atlas'
face_mesh <- get_dec_mesh('face')
setwd("../../")

# load up the landmarks for our atlas mesh
atlas_head_lm <- head_array[,, which(dimnames(head_array)[[3]] == "chick_ctr_23")]
  
# Divide the data into type of landmark
head_fixed.lm <- c(1:dim(atlas_head_lm)[1])[-c(as.matrix(curveslide_all)[,2], head_surface.lm)]
head_curves.lm <- as.matrix(curveslide_all)[,2]
head_surface.lm <- head_surface.lm


# Plot the mesh with the landmarks, curve semi-landmarks, and surface semi-landmarks
open3d(zoom = 0.75, windowRect = c(0, 0, 700, 700)) # bigger rgl window
shade3d(face_mesh, color = "gray", alpha = 0.8) # alpha for transparency

# you need to click and move mouse around in the figure until it is en-face with
# the window, then run the line below to save the orientation
frontal <- par3d()$userMatrix
rgl.close3d()

save(frontal, file = "./lm_data/RGL_head_pos.rdata")
load("./lm_data/RGL_head_pos.rdata")

# Frontal view
# this figure is in the manuscript
open3d(zoom = 0.75, userMatrix = frontal, windowRect = c(0, 0, 1000, 700)) 
rgl::shade3d(face_mesh, color = "gray", alpha =1, specular = "black" )
rgl::plot3d(atlas_head_lm[head_fixed.lm,], aspect = "iso", type = "s", size=.8, col = "black", add = T)
rgl::plot3d(atlas_head_lm[head_curves.lm,], aspect = "iso", type = "s", size=0.8, col = "orange", add = T)
rgl::plot3d(atlas_head_lm[head_surface.lm,], aspect = "iso", type = "s", size=0.8, col = "turquoise2", add = T)
rgl::snapshot("./figs/head_LM_frontal.png", top = TRUE)
writeASY(title = "head_LM_frontal", prc = FALSE)
rgl::close3d()


#### 4. PCA plots HEAD ####
PCA_head <- gm.prcomp(GPA_head$coords)
summary(PCA_head)

# Delete file if it exists
if (file.exists("./output/PCA_head_shape_coords.txt")) {
  file.remove("./output/PCA_head_shape_coords.txt")
}
if (!dir.exists('./output/')) dir.create('./output/')
cat("PCA shape variables raw", capture.output(summary(PCA_head)), 
    file="./output/PCA_head_shape_coords.txt", sep="\n", append=TRUE)


palette(c("navy", "darkorange")) # this is fine, we'll change it in illustrator
classifiers$treatment <- as.character(classifiers$treatment)
classifiers$treatment <- as.factor(classifiers$treatment)

# this figure is in the manuscript
pdf("./figs/PCA_head_shape_treatment_raw.pdf", width = 8.25, height = 6)
plot(PCA_head, pch = 19, col = classifiers$treatment, cex = 1.25, xlim=c(-.1,.15))
ordiellipse(PCA_head, classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA of shape coordinates - CTRL vs treatment")
dev.off()


# PC1 - PC4
png("./figs/PCA_head_treatment_PC1-4.png", width = 750, height = 1200)
pdf("./figs/PCA_head_treatment_PC1-4.pdf", width = 8, height = 12)
par(mfrow=c(3,2))
plot(PCA_head, pch = 19, axis1 = 1, axis2 = 2, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_head, classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PC1 ~ PC2")

plot(PCA_head, pch = 19, axis1 = 1, axis2 = 3, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_head$x[,c(1,3)], classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PC1 ~ PC3")

plot(PCA_head, pch = 19, axis1 = 1, axis2 = 4, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_head$x[,c(1, 4)], classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PC1 ~ PC4")

plot(PCA_head, pch = 19, axis1 = 2, axis2 = 3, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_head$x[,c(2,3)], classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("bottomleft", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PC2 ~ PC3")

plot(PCA_head, pch = 19, axis1 = 2, axis2 = 4, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_head$x[,c(2,4)], classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PC2 ~ PC4")

plot(PCA_head, pch = 19, axis1 = 3, axis2 = 4, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_head$x[,3:4], classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0.2, lty = 0)
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PC3 ~ PC4")

dev.off()

par(mfrow=c(1,1))

PCA_comp <- PCA_head
class(PCA_comp) <- "princomp"

png("./figs/PCA_head_shape_scree_plot.png", width = 300, height = 300)
pdf("./figs/PCA_head_shape_scree_plot.pdf", height = 5, width = 5)
pca_scree <- fviz_eig(PCA_comp, addlabels=TRUE, hjust = -0.3,
                      barfill="darkgrey", barcolor ="black",
                      linecolor ="blue") + ylim(0, 85) + 
  theme_classic()

print(pca_scree)
dev.off()

#### 5. BOXPLOTS HEAD ####

classifiers

Pdist <- ShapeDist(GPA_head$coords, GPA_head$consensus)

gdf_head <- geomorph.data.frame(GPA_head, treatment = classifiers$treatment, Pdist = Pdist)


ggplot_df <- as.data.frame(cbind(as.character(gdf_head$Csize), 
                                 as.character(gdf_head$treatment), 
                                 as.character(gdf_head$Pdist)))
colnames(ggplot_df) <- c("Csize", "treatment", "Pdist")

row.names(ggplot_df) <- dimnames(gdf_head$coords)[[3]]

head(ggplot_df)

ggplot_df$treatment <- as.factor(ggplot_df$treatment)
ggplot_df$Csize <- as.numeric(as.character(ggplot_df$Csize))
ggplot_df$Pdist <- as.numeric(as.character(ggplot_df$Pdist))

str(ggplot_df)
log(ggplot_df$Csize)

pdf("./figs/Csize_boxplot_HEAD.pdf", width = 6.5, height = 6.5)
ggplot(ggplot_df, aes(x=treatment, y=log(Csize), fill=treatment)) + 
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("navy", "darkorange")) +
  geom_jitter(width = 0.1, size = 1.25) +
  theme(plot.title = element_text(size = 20, face = "bold"), axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"), legend.text=element_text(size=12), 
        legend.title=element_text(size=15, face = "bold"))
dev.off()


pdf("./figs/Pdist_treatment.pdf", width = 6.5, height = 6.5)
ggplot(ggplot_df, aes(Pdist, fill = treatment)) +
  scale_fill_manual(values=c("navy", "darkorange")) + geom_density(alpha = 0.65) + 
  ggtitle("Procrustes distances head shape") + xlab("Procrustes distance") + ylab("Relative density") +
  theme(plot.title = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"), legend.text=element_text(size=10), 
        legend.title=element_text(size=12, face = "bold"))
dev.off()

#### 6. ANOVAs & ALLOMETRY - HEAD ####

allometry_all <- procD.lm(coords ~ Csize, data = gdf_head, iter = 999, RRPP = TRUE)
summary(allometry_all)
# Delete file if it exists
if (file.exists("./output/HEAD_allometry_all.txt")) {
  file.remove("./output/HEAD_allometry_all.txt")
}

cat("Allometry (coords ~ Csize)", capture.output(summary(allometry_all)), 
    file="./output/HEAD_allometry_all.txt", sep="\n", append=TRUE)

shape_residuals <- allometry_all$residuals # let's save this for later. We will run a PCA again on the residuals
write.csv(shape_residuals, "./output/HEAD_shape_residuals.csv")

allometry_treatment <- procD.lm(coords ~ Csize * treatment, data = gdf_head, iter = 999, RRPP = TRUE)
summary(allometry_treatment)

# Delete file if it exists
if (file.exists("./output/HEAD_allometry_treatment.txt")) {
  file.remove("./output/HEAD_allometry_treatment.txt")
}

cat("Allometry * treatment (coords ~ Csize * treatment)", capture.output(summary(allometry_treatment)), 
    file="./output/HEAD_allometry_treatment.txt", sep="\n", append=TRUE)


# ANOVA treatment
treatment <- procD.lm(coords ~ treatment, data = gdf_head, RRPP = TRUE)
summary(treatment)

# Delete file if it exists
if (file.exists("./output/ANOVA_HEAD_shape_treatment.txt")) {
  file.remove("./output/ANOVA_HEAD_shape_treatment.txt")
}

cat("ANOVA shape - treatment (coords ~ treatment)", capture.output(summary(treatment)), 
    file="./output/ANOVA_HEAD_shape_treatment.txt", sep="\n", append=TRUE)


one_way <- aov(log(gdf_head$Csize) ~ gdf_head$treatment)
summary(one_way)

# Delete file if it exists
if (file.exists("./output/ANOVA_HEAD_Csize_treatment.txt")) {
  file.remove("./output/ANOVA_HEAD_Csize_treatment.txt")
}

cat("ANOVA Csize - treatment (Csize ~ treatment)", capture.output(summary(one_way)), 
    file="./output/ANOVA_HEAD_Csize_treatment.txt", sep="\n", append=TRUE)

#### 7. HEATMAPS - HEAD ####
#### 7.1. HEAD HEATMAPS group means ####

# load mesh of the specimen closest to the mean shape
face_mesh <- geomorph::read.ply("./data/ATLAS_chick_ctr_23_smooth_ascii_only_face.ply") 
# The mesh should only be surface, and be just have the frontal side here
head_lowres <- vcgQEdecim(face_mesh, percent = 0.15)
atlas_head_lm <- head_array[,, which(dimnames(head_array)[[3]] == "chick_ctr_23")]

levels(gdf_head$treatment)

who_is_MUT <- which(gdf_head$treat == "triple")
MUT_coords <- gdf_head$coords[,,who_is_MUT]
dim(MUT_coords)

CTRL_coords <- gdf_head$coords[,,-who_is_MUT]
dim(CTRL_coords)

dim(gdf_head$coords)
# we are going to calculate the mean shape on the shape coordinates. Everything will be in Procrustes dist.

MUT_mean_shape <- mshape(MUT_coords)

CTRL_mean_shape <- mshape(CTRL_coords)
str(CTRL_mean_shape)

# Create morphed meshes
MUT_mesh <- tps3d(head_lowres, as.matrix(atlas_head_lm), MUT_mean_shape, threads = 1)
CTRL_mesh <- tps3d(head_lowres, as.matrix(atlas_head_lm), CTRL_mean_shape, threads = 1)

open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = dorsal)
# meshDist(CTRL_mesh, MUT_mesh, rampcolors = diverge_hsv(n = 3), sign = TRUE)
meshDist(CTRL_mesh, MUT_mesh, rampcolors = c("blue", "white", "red"), sign = TRUE)
rgl.snapshot("./figs/Heatmap_CTRL_MUT_head_dorsal.png", top = TRUE)

writeWebGL(dir = "webGL", filename = file.path("./figs/Heatmap_CTRL_MUT_head.html"),
           template = system.file(file.path("WebGL", "template.html"), package = "rgl"),
           snapshot = TRUE, width = 1000, height = 700)

# set positions again
dorsal <- par3d()$userMatrix

rgl.close()

save(dorsal, file = "./figs/RGL_head_heatmaps_pos.rdata")

load("./figs/RGL_head_heatmaps_pos.rdata")


open3d(zoom=0.75, userMatrix = dorsal, windowRect = c(0,0,1000,700)) 
shade3d(MUT_mesh, color="gray", alpha=0.9)
rgl.snapshot("./figs/Morph_MUT_dorsal_head.png", top = TRUE)
writePLY("./output/MUT_mesh.ply")
rgl.close()


open3d(zoom=0.75, userMatrix = dorsal, windowRect = c(0,0,1000,700)) 
shade3d(CTRL_mesh, color="gray", alpha=0.9)
rgl.snapshot("./figs/Morph_CTRL_dorsal_head.png", top = TRUE)
writePLY("./output/CTRL_mesh.ply")
rgl.close()



#Create page with MUT vs CTRL morphs & heatmaps
# dorsal
Morph_CTRL_dorsal <- image_read(paste0("./figs/Morph_CTRL_dorsal_head.png"))
Morph_MUT_dorsal <- image_read(paste0("./figs/Morph_MUT_dorsal_head.png"))
Heatmap_dorsal <- image_read(paste0("./figs/Heatmap_CTRL_MUT_head_dorsal.png"))

Morph_CTRL_dorsal <- image_annotate(Morph_CTRL_dorsal, "CTRL", font = "times", location = "+400+20", size = 100)
Morph_MUT_dorsal <- image_annotate(Morph_MUT_dorsal, "MUT", font = "times", location = "+425+20", size = 100)
Heatmap_dorsal <- image_annotate(Heatmap_dorsal, "Heatmap", font = "times", location = "+350+20", size = 100)




stack_img <- image_append(image_scale(c(Morph_CTRL_dorsal, Morph_MUT_dorsal, Heatmap_dorsal), "x200"), stack = TRUE)

image_browse(stack_img)

image_write(stack_img, path = paste0("./figs/CTRL_MUT_morphs_heatmaps_head.png"), format = "png")


####7.2. HEATMAPS - PC1-10 ####

PC_min <- tps3d(face_mesh, as.matrix(atlas_head_lm), PCA_head$shapes[[1]]$min, threads = 1)
open3d(zoom=0.75, userMatrix = dorsal, windowRect= c(0,0,1000,700))
shade3d(PC_min, color="grey")

PC_max <- tps3d(face_mesh, as.matrix(atlas_head_lm), PCA_head$shapes[[1]]$max, threads = 1)
open3d(zoom=0.75, userMatrix = dorsal, windowRect= c(0,0,1000,700))
shade3d(PC_max, color="grey")


# Have a look at the heatmap, and change colour if needed
open3d(zoom=0.75, userMatrix = dorsal, windowRect= c(0,0,1000,700))
x11(width=1.7, height=8)
meshDist(PC_max, PC_min, rampcolors = c("blue", "white", "red"), sign = FALSE)

summary(PCA_head)
n_dimensions <- 10 # number of PCs to include in the figure, decide depending on % of variance explained in PCs

# This loop automatically positions faces in dorsal only for the 10 PCs
for (i in 1:n_dimensions){
  PC_min <- tps3d(head_lowres, as.matrix(atlas_head_lm), PCA_head$shapes[[i]]$min, threads = 1)
  PC_max <- tps3d(head_lowres, as.matrix(atlas_head_lm), PCA_head$shapes[[i]]$max, threads = 1)
  PC <- paste0("PC", i)
  
  #dorsal views
  open3d(zoom=0.75, userMatrix = dorsal, windowRect= c(0,0,1000,700))
  shade3d(PC_max, color="grey")
  rgl.snapshot(paste0("./figs/pc_morphs/head/",PC,"_","max_dorsal.png"), top = TRUE )
  clear3d()
  rgl.close()
  
  
  #dorsal views
  open3d(zoom=0.75, userMatrix = dorsal, windowRect= c(0,0,1000,700))
  shade3d(PC_min, color="grey")
  rgl.snapshot(paste0("./figs/pc_morphs/head/",PC,"_","min_dorsal.png"), top = TRUE )
  clear3d()
  rgl.close()
  
  
  #Heatmaps
  open3d(zoom=0.75, userMatrix = dorsal, windowRect= c(0,0,1000,700))
  x11(width=1.7, height=8)
  meshDist(PC_min, PC_max, rampcolors = c("darkblue", "blue", "white", "red", "darkred"), sign = TRUE)
  rgl.snapshot(paste0("./figs/pc_morphs/head/heat_",PC,"_","dorsal.png"), top = TRUE)
  # dev.copy(pdf, paste0("./figs/pc_morphs/head/heat_",PC,"_scale.pdf"), width=2, height=9.5)
  dev.print(pdf, paste0("./figs/pc_morphs/head/heat_",PC,"_scale.pdf"), width=2, height=9.5)
  dev.off()
  clear3d()
  rgl.close()
  
  rm(PC_max, PC_min, PC) 
}
#END LOOP





#Create page with PCs
i <- 1
PC_min_img_sup <- image_read(paste0("./figs/pc_morphs/head/PC",i,"_min_dorsal.png"))
PC_max_img_sup <- image_read(paste0("./figs/pc_morphs/head/PC",i,"_max_dorsal.png"))
PC_HT_img_sup <- image_read(paste0("./figs/pc_morphs/head/heat_PC",i,"_dorsal.png"))


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
  PC_min_img_sup <- image_read(paste0("./figs/pc_morphs/head/PC",i,"_min_dorsal.png"))
  PC_max_img_sup <- image_read(paste0("./figs/pc_morphs/head/PC",i,"_max_dorsal.png"))
  PC_HT_img_sup <- image_read(paste0("./figs/pc_morphs/head/heat_PC",i,"_dorsal.png"))
  
  
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
image_write(imgs, path = paste0("./figs/pc_morphs/head_SHAPE_heatmap_morphs_PC1-", i, ".png"), format = "png")

# We are done with these analyses

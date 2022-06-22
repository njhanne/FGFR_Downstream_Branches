#### 0. Load R packages ####
# install.packages(c('rgl', 'geomorph' 'devtools', 'Morpho', 'Rvcg', 'magick', 'Evomorph', 'ggplot2', 'vegan', 'factoextra', 'gt'))
library(rgl)
library(geomorph)
library(devtools)
install_github("marta-vidalgarcia/morpho.tools.GM", force = TRUE)
library(morpho.tools.GM)
# install_github("marta-vidalgarcia/mesh_process")
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

# setwd("~/Documents/GITHUB_repos/FGFR-Branches-GM/")
# setwd("C:/Users/nhanne/Box/FGF_inhibitor_paper_5-26-2020/Morphology/3D_data")

# # Create directories.
### WILL IT WORK ON WINDOWS???
# getwd() # where are we? We should always be in the project main directory
# folder_structure <- c("./figs", "./figs/pc_morphs", "./figs/pc_morphs/head", 
#                       "./data", "./data/atlas", "./data/Prop_LMs", "./data/Prop_LMs/LM_head", "./output", "./R")
# 
# dir.exists(folder_structure) # it should all be false unless you created folders manually. I recommend doing it this way
# 
# for (i in 1:length(folder_structure)){
#   dir.create(folder_structure[i], mode = "0777") # create these folders
# }


# Then save this script inside the R folder
# Copy landmark data to the Prop_LMs folder
# Copy atlases inside the data/atlas folder. We will need:
# PLY of head, endocast & mandible. 
# CURVESLIDE FILES for endocast, head, mandible
# ATLAS TAG file for head, endocast, mandible
# Copy classifiers file inside data and make sure it is a csv and there are no issues

#### 1. LOAD LM DATA & CLASSIFIERS ####
classifiers_unord  <- read.csv("./data/classifiers.csv", header = TRUE)
head(classifiers_unord)
tail(classifiers_unord)

str(classifiers_unord)
classifiers_unord$treatment <- as.factor(classifiers_unord$treatment)
row.names(classifiers_unord) <- classifiers_unord$id
classifiers_unord

#### 2. ANALYSES PREP EMBRYO FACES ####
#### 2.1 IMPORT LANDMARK DATA ####
# 1. Import all fcsv files into separate arrays (all specimens for east LM set)
?morpho.tools.GM::fcsv2array

# LANDMARKS
setwd("./data/Landmarks/")
dir()
LMs <- fcsv2array(string_del = "_Fiducials")
str(LMs)
dimnames(LMs)[[3]]
row.names(classifiers_unord)

classifiers <- classifiers_unord[match(dimnames(LMs)[[3]], row.names(classifiers_unord)),]

setwd("../../")

# CURVE SEMILANDMARKS
fcsv_number_check <- function(dir, pattern) {
  fcsv_list <- dir(dir, pattern = pattern)
  n_land <- vector("numeric", length = length(fcsv_list))
  for (i in 1:length(fcsv_list)) {
    n_land[i] <- length(count.fields(fcsv_list[i]))
  }
  checker_df <- do.call(rbind.data.frame, Map('c', fcsv_list, n_land))
  return(checker_df)
}

setwd("./data/Semi_Curves/")
dir(pattern = "center")
dir(pattern = "L")
dir(pattern = "R")

curve_semis_center <- fcsv2array(pattern = "*center*", string_del = "_center_semi-curve")
curve_semis_L <- fcsv2array(pattern = "*_L_semi-curve*", string_del = "_L_semi-curve")
curve_semis_R <- fcsv2array(pattern = "*_R_semi-curve*", string_del = "_R_semi-curve")

str(curve_semis_center)
str(curve_semis_L)
str(curve_semis_R)

dimnames(curve_semis_center)[[3]]
dimnames(curve_semis_L)[[3]]
dimnames(curve_semis_R)[[3]]


# SURFACE SEMILANDMARKS

setwd("../Semi_Points/")

dir()
surf_semis_L <- fcsv2array(pattern = "*_L_semi-lm*", string_del = "_L_semi-lm")
surf_semis_R <- fcsv2array(pattern = "*_R_semi-lm*", string_del = "_R_semi-lm")

setwd("../../")

#### 2.2. NEW ARRAY ####
# Combine all arrays into a single one, in a particular order that makes sense for our GMM
# CENTRE CURVE
# Delete positions #1 & #5 in curve_semis_center
curve_semis_center <- curve_semis_center[-c(1,5),,]

head_array <- abind(LMs, curve_semis_center, curve_semis_L, curve_semis_R, surf_semis_L, 
                    curve_semis_R, along = 1)

head_array_no_side_curve <- abind(LMs, curve_semis_center, surf_semis_L, curve_semis_R, along = 1)

dimnames(head_array)[[3]]

dim(LMs)[1]
dim(curve_semis_center)[1]
dim(curve_semis_L)[1]
dim(curve_semis_R)[1]
dim(surf_semis_L)[1]
dim(surf_semis_R)[1]


#### 2.3. CURVESLIDE ####
# Generate a curveslide file for the GPA that tells us how the curve semis have to slide (constrained)

# curve semis positions
semis_center <- (dim(LMs)[1]+1):(dim(LMs)[1]+dim(curve_semis_center)[1])
semis_L <- (dim(LMs)[1]+dim(curve_semis_center)[1]+1):(dim(LMs)[1]+dim(curve_semis_center)[1]+dim(curve_semis_L)[1])
semis_R <- (dim(LMs)[1]+dim(curve_semis_center)[1]+dim(curve_semis_L)[1]+1):(dim(LMs)[1]+dim(curve_semis_center)[1]+dim(curve_semis_L)[1]+dim(curve_semis_R)[1])

# PLOTTING THE LANDMARKS TO CHECK THE SEMIS
head_mesh_spec1 <- vcgImport("./data/Meshes/chick_ctr_1.ply") # Not sure what is wrong with this mesh
head_mesh_spec1_dec <- vcgQEdecim(head_mesh_spec1, percent = 0.15)
head_mesh_spec1_dec <- vcgQEdecim(head_mesh_spec1_dec, percent = 0.25)

open3d(zoom = 0.75, windowRect = c(0, 0, 700, 700)) 
rgl::shade3d(head_mesh_spec1_dec, color = "gray", alpha =0.9)
rgl::plot3d(LMs[,,1], aspect = "iso", type = "s", size=1.2, col = "darkblue", add = T)
rgl::text3d(x = LMs[,1,1],
            y = LMs[,2,1],
            z = LMs[,3,1],
            texts = c(1:dim(LMs)[1]),
            cex = 1.5, offset = 0.5, pos = 1)

rgl::plot3d(curve_semis_center[,,1], aspect = "iso", type = "s", size=0.5, col = "orange", add = T)
rgl::text3d(x = curve_semis_center[,1,1],
            y = curve_semis_center[,2,1],
            z = curve_semis_center[,3,1],
            texts = c(1:dim(curve_semis_center)[1]),
            cex = 0.75, offset = 0.5, pos = 2)

rgl::plot3d(curve_semis_R[,,1], aspect = "iso", type = "s", size=0.75, col = "orange", add = T)
rgl::text3d(x = curve_semis_R[,1,1],
            y = curve_semis_R[,2,1],
            z = curve_semis_R[,3,1],
            texts = c(1:dim(curve_semis_R)[1]),
            cex = 0.75, offset = 0.5, pos = 2)

rgl::plot3d(curve_semis_L[,,1], aspect = "iso", type = "s", size=0.75, col = "orange", add = T)
rgl::text3d(x = curve_semis_L[,1,1],
            y = curve_semis_L[,2,1],
            z = curve_semis_L[,3,1],
            texts = c(1:dim(curve_semis_L)[1]),
            cex = 0.75, offset = 0.5, pos = 2)
rgl.close()


# CENTER CURVE - CURVESLIDE
# We need to create a curveslide matrix to know where from to where to the semis slide
# semi1 slides between LM9 & semi2
# semi2 slides between semi1 & semi3
# semi3 slides between semi2 & LM10
curve_c_left <- c(9, semis_center[c(1,2)])
curve_c_right <- c(semis_center[c(2,3)], 10)
curveslide_c <- cbind(curve_c_left, semis_center, curve_c_right)

# RIGHT & LEFT CURVES - positions from top to bottoms
curve_L_left <- semis_L[c(1:(length(semis_L)-2))] # remember that landmarks 16 & 24 will be treated as landmarks!
# Change them later for the final analyses to slide between the nasal pits and the next semis
curve_L_right <- semis_L[c(3:length(semis_L))]
curve_L_center <- semis_L[c(2:(length(semis_L)-1))]
curveslide_L <- cbind(curve_L_left, curve_L_center, curve_L_right)

curve_R_left <- semis_R[c(1:(length(semis_R)-2))] 
curve_R_right <- semis_R[c(3:length(semis_R))]
curve_R_center <- semis_R[c(2:(length(semis_R)-1))]
curveslide_R <- cbind(curve_R_left, curve_R_center, curve_R_right)


# all our curveslide matrices
ls(pattern = "curveslide*")
curveslide_list <- lapply(ls(pattern = "curveslide*"), get)
str(curveslide_list)
curveslide_all <- do.call(rbind, curveslide_list)
curveslide_all <- as.data.frame(curveslide_all)
colnames(curveslide_all) <- c("left", "sliding", "right")
write.csv(curveslide_all, "./data/curveslide.csv")

# Surface landmarks
head_surface.lm <- (dim(LMs)[1]+dim(curve_semis_center)[1]+dim(curve_semis_L)[1]+dim(curve_semis_R)[1]+1):dim(head_array)[1]
head_surface.lm2 <- (dim(LMs)[1]+dim(curve_semis_center)[1]+1):dim(head_array_no_side_curve)[1]


#### 2.5. GPA ####
# Use this to find out who is the closest to the shape sphere centroid
str(head_array)
GPA_head_o <- geomorph::gpagen(A = head_array, curves = as.matrix(curveslide_all), 
                             surfaces = head_surface.lm)

outlier <- plotOutliers_percentile(A = GPA_head_o$coords, percentile = 0.99, save.plot = TRUE)



GPA_head_o2 <- geomorph::gpagen(A = head_array_no_side_curve, curves = as.matrix(curveslide_c), 
                               surfaces = head_surface.lm2)

outlier <- plotOutliers_percentile(A = GPA_head_o2$coords, percentile = 0.99, save.plot = TRUE)
# Looks very much like an outlier. Maybe some errors in the landmarking? Ask Nicholas

# Get rid of it until Nicholas fixes the landmarks
clean_head_array <- head_array[,,-which(dimnames(head_array)[[3]] == row.names(outlier$Proc_d_percentile))]

GPA_head <- geomorph::gpagen(A = clean_head_array, curves = as.matrix(curveslide_all), 
                               surfaces = head_surface.lm)


outliers_c <- plotOutliers_percentile(A = GPA_head$coords, percentile = 0.90, save.plot = FALSE)
# They don't look like outliers to me, just a more severe phenotype
# excellent, keep going

# For the atlas we are going to: (1) perform a GPA, (2) find the specimens the closest to the center of the morphospace
# Pdist <- ShapeDist(GPA_head$coords, GPA_head$consensus)
# find who it is
min(outliers_c$All_Proc_d$`Proc. d. from mean`) # actually I got lazy & did it this easier way, same concept as above though
row.names(outliers_c$All_Proc_d[which(outliers_c$All_Proc_d$`Proc. d. from mean` == min(outliers_c$All_Proc_d$`Proc. d. from mean`)),])
# Ah, very nice! Need to clean this mesh on meshlab now: "chick_ctr_23"

classifiers <- classifiers[-which(row.names(classifiers) == "chick_exp_3"),]
#### 3. ATLAS HEAD & PLOTS ####
# (hypervolume) - which means the smalles Procrustes distance, (3) use the set of LMs & mesh of that specimen
head_mesh <- geomorph::read.ply("./data/atlas/Calgary_Adult_Cranium_Atlas_DS_ascii.ply") # Not sure what is wrong with this mesh
head_mesh <- geomorph::read.ply("./data/atlas/Global_Adult_ONLY_head_Atlas_lowres.ply") # I cleaned and smoothed this one manually, looks good
head_lowres <- vcgQEdecim(head_mesh, percent = 0.15)
atlas_head_lm <- 
  
# Divide the data into type of landmark (we did that before)
head_fixed.lm
head_curves.lm
head_surface.lm


# Plot the mesh with the landmarks, curve semi-landmarks, and surface semi-landmarks
open3d(zoom = 0.75, windowRect = c(0, 0, 700, 700)) # bigger rgl window
shade3d(head_lowres, color = "gray", alpha = 0.8) # alpha for transparency

lateral <- par3d()$userMatrix
frontal <- par3d()$userMatrix
rgl.close()

save(lateral, frontal, file = "./data/atlas/RGL_head_pos.rdata")


load("./data/atlas/RGL_head_pos.rdata")

# Lateral view
open3d(zoom = 0.75, userMatrix = lateral, windowRect = c(0, 0, 1000, 700)) 
rgl::shade3d(head_mesh, color = "gray", alpha =0.9)
rgl::plot3d(atlas_head_lm[head_fixed.lm,], aspect = "iso", type = "s", size=1.2, col = "darkblue", add = T)
# rgl::text3d(x = atlas_head_lm[head_fixed.lm, 1], 
#             y = atlas_head_lm[head_fixed.lm, 2], 
#             z=  atlas_head_lm[head_fixed.lm, 3], 
#             texts = row.names(atlas_head_lm[head_fixed.lm, ]), 
#             cex = 1.5, offset = 0.5, pos = 3)
rgl::plot3d(atlas_head_lm[head_curves.lm,], aspect = "iso", type = "s", size=0.75, col = "orange", add = T)
rgl::plot3d(atlas_head_lm[head_surface.lm,], aspect = "iso", type = "s", size=0.6, col = "turquoise2", add = T)
rgl::rgl.snapshot("./figs/head_LM_lateral.png", top = TRUE)
rgl::rgl.close()

# Frontal view
open3d(zoom = 0.75, userMatrix = frontal, windowRect = c(0, 0, 1000, 700)) 
rgl::shade3d(head_mesh, color = "gray", alpha =0.9)
rgl::plot3d(atlas_head_lm[head_fixed.lm,], aspect = "iso", type = "s", size=1.2, col = "darkblue", add = T)
rgl::plot3d(atlas_head_lm[head_curves.lm,], aspect = "iso", type = "s", size=0.75, col = "orange", add = T)
rgl::plot3d(atlas_head_lm[head_surface.lm,], aspect = "iso", type = "s", size=0.6, col = "turquoise2", add = T)
rgl::rgl.snapshot("./figs/head_LM_frontal.png", top = TRUE)
rgl::rgl.close()

# Combine the two views
side_lm <- image_read("./figs/head_LM_lateral.png")
side_lm <- image_annotate(side_lm, "Lateral", font = "times", location = "+80+120", size = 50)
side_lm <- image_crop(side_lm, "1000x500+0+90")
side_lm

frontal_lm <- image_read("./figs/head_LM_frontal.png")
frontal_lm <- image_annotate(frontal_lm, "Frontal", font = "times", location = "+80+120", size = 50)
frontal_lm <- image_crop(frontal_lm, "1000x500+0+90")
frontal_lm

stack_views <- c(side_lm, frontal_lm)
stacked_images <- image_append(image_scale(stack_views), stack = TRUE)

image_browse(stacked_images)
image_write(stacked_images, path = "./figs/LM_scheme_all_views.png", format = "png")


#### 4. PCA plots HEAD ####
PCA_head <- gm.prcomp(GPA_head$coords)
summary(PCA_head)
str(PCA_head)

# Delete file if it exists
if (file.exists("./output/PCA_head_shape_coords.txt")) {
  file.remove("./output/PCA_head_shape_coords.txt")
}
cat("PCA shape variables raw", capture.output(summary(PCA_head)), 
    file="./output/PCA_head_shape_coords.txt", sep="\n", append=TRUE)

# # Make a different table to get HTML table (or latex)
# eigenvalues <- round(PCA_head$d, digits = 5)
# prop_var <- round(PCA_head$d*100/sum(PCA_head$d), digits = 3)
# cumul_var <- round(cumsum(PCA_head$d)*100/sum(PCA_head$d), digits = 3)
# 
# 
# PCA_head_summary <- as.data.frame(rbind(eigenvalues, prop_var, cumul_var))
# colnames(PCA_head_summary) <- paste0(rep("PC", length(colnames(PCA_head_summary))),
#                                       c(1:length(colnames(PCA_head_summary))))
# PCA_head_summary <- as.data.frame(cbind(c("Eigenvalues", "Proportion of Variance (%)",
#                                            "Cumulative Variance (%)"), PCA_head_summary))
# colnames(PCA_head_summary)[1] <- c(" ")
# 
# PCA_summary_table <- PCA_head_summary %>%
#   gt() %>%
#   tab_header("Principal Component Analysis - head")
# 
# ?gtsave
# gtsave(PCA_summary_table, "./figs/PCA_head_summary.html")


levels(classifiers$treatment)
palette(c("navy", "darkorange"))
classifiers$treatment <- as.character(classifiers$treatment)
classifiers$treatment <- as.factor(classifiers$treatment)

# png("./figs/PCA_head_shape_treatment_raw.png", width = 750, height = 600)
pdf("./figs/PCA_head_shape_treatment_raw.pdf", width = 8.25, height = 6)
plot(PCA_head, pch = 19, col = classifiers$treatment, cex = 1.25)
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
head_mesh <- geomorph::read.ply("./data/atlas/chick_ctr_23.ply") # Make it nice by segmenting in 3DSlicer
# The mesh should only be surface
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

# Just fixed landmarks, no semis
MUT_mesh <- tps3d(head_mesh, as.matrix(atlas_head_lm), MUT_mean_shape, threads = 1)
CTRL_mesh <- tps3d(head_mesh, as.matrix(atlas_head_lm), CTRL_mean_shape, threads = 1)

open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = lateral)
# meshDist(CTRL_mesh, MUT_mesh, rampcolors = diverge_hsv(n = 3), sign = TRUE)
meshDist(CTRL_mesh, MUT_mesh, rampcolors = c("darkblue", "blue", "blue", "white", "red", "red", "darkred"), sign = TRUE)
rgl.snapshot("./figs/Heatmap_CTRL_MUT_head_lateral.png", top = TRUE)

open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = ventral)
meshDist(CTRL_mesh, MUT_mesh, rampcolors = c("darkblue", "blue", "blue", "white", "red", "red", "darkred"), sign = TRUE)
rgl.snapshot("./figs/Heatmap_CTRL_MUT_head_ventral.png", top = TRUE)

open3d(zoom=0.75, windowRect = c(0,0, 1000, 700), userMatrix = dorsal)
meshDist(CTRL_mesh, MUT_mesh, rampcolors = c("darkblue", "blue", "blue", "white", "red", "red", "darkred"), sign = TRUE)
rgl.snapshot("./figs/Heatmap_CTRL_MUT_head_dorsal.png", top = TRUE)

writeWebGL(dir = "webGL", filename = file.path("./figs/Heatmap_CTRL_MUT_head.html"),
           template = system.file(file.path("WebGL", "template.html"), package = "rgl"),
           snapshot = TRUE, width = 1000, height = 700)

# set positions again
lateral <- par3d()$userMatrix
dorsal <- par3d()$userMatrix
ventral <- par3d()$userMatrix
rgl.close()

save(lateral, dorsal, ventral, file = "./figs/RGL_head_heatmaps_pos.rdata")

load("./figs/RGL_head_heatmaps_pos.rdata")


open3d(zoom=0.75, userMatrix = lateral, windowRect = c(0,0,1000,700)) 
shade3d(MUT_mesh, color="gray", alpha=0.9)
rgl.snapshot("./figs/Morph_MUT_lateral_head.png", top = TRUE)
open3d(zoom=0.75, userMatrix = dorsal, windowRect = c(0,0,1000,700)) 
shade3d(MUT_mesh, color="gray", alpha=0.9)
rgl.snapshot("./figs/Morph_MUT_dorsal_head.png", top = TRUE)
open3d(zoom=0.75, userMatrix = ventral, windowRect = c(0,0,1000,700)) 
shade3d(MUT_mesh, color="gray", alpha=0.9)
rgl.snapshot("./figs/Morph_MUT_ventral_head.png", top = TRUE)


open3d(zoom=1, userMatrix = lateral, windowRect = c(0,0,1000,700)) 
shade3d(CTRL_mesh, color = "gray", alpha=0.9)
writePLY("./output/CTRL_mesh.ply")
rgl.close()


open3d(zoom=0.75, userMatrix = lateral, windowRect = c(0,0,1000,700)) 
shade3d(CTRL_mesh, color="gray", alpha=0.9)
rgl.snapshot("./figs/Morph_CTRL_lateral_head.png", top = TRUE)
open3d(zoom=0.75, userMatrix = dorsal, windowRect = c(0,0,1000,700)) 
shade3d(CTRL_mesh, color="gray", alpha=0.9)
rgl.snapshot("./figs/Morph_CTRL_dorsal_head.png", top = TRUE)
open3d(zoom=0.75, userMatrix = ventral, windowRect = c(0,0,1000,700)) 
shade3d(CTRL_mesh, color="gray", alpha=0.9)
rgl.snapshot("./figs/Morph_CTRL_ventral_head.png", top = TRUE)



#Create page with MUT vs CTRL morphs & heatmaps
# Lateral
Morph_CTRL_lateral <- image_read(paste0("./figs/Morph_CTRL_lateral_head.png"))
Morph_MUT_lateral <- image_read(paste0("./figs/Morph_MUT_lateral_head.png"))
Heatmap_lateral <- image_read(paste0("./figs/Heatmap_CTRL_MUT_head_lateral.png"))

Morph_CTRL_lateral <- image_annotate(Morph_CTRL_lateral, "CTRL", font = "times", location = "+400+20", size = 100)
Morph_MUT_lateral <- image_annotate(Morph_MUT_lateral, "MUT", font = "times", location = "+425+20", size = 100)
Heatmap_lateral <- image_annotate(Heatmap_lateral, "Heatmap", font = "times", location = "+350+20", size = 100)

# dorsal
Morph_CTRL_dorsal <- image_read(paste0("./figs/Morph_CTRL_dorsal_head.png"))
Morph_MUT_dorsal <- image_read(paste0("./figs/Morph_MUT_dorsal_head.png"))
Heatmap_dorsal <- image_read(paste0("./figs/Heatmap_CTRL_MUT_head_dorsal.png"))

# ventral
Morph_CTRL_ventral <- image_read(paste0("./figs/Morph_CTRL_ventral_head.png"))
Morph_MUT_ventral <- image_read(paste0("./figs/Morph_MUT_ventral_head.png"))
Heatmap_ventral <- image_read(paste0("./figs/Heatmap_CTRL_MUT_head_ventral.png"))


stack_img <- image_append(c(image_append(image_scale(c(Morph_CTRL_lateral, Morph_MUT_lateral, Heatmap_lateral), "x200"), stack = FALSE), 
                            image_append(image_scale(c(Morph_CTRL_dorsal, Morph_MUT_dorsal, Heatmap_dorsal), "x200"), stack = FALSE),
                            image_append(image_scale(c(Morph_CTRL_ventral, Morph_MUT_ventral, Heatmap_ventral), "x200"), stack = FALSE)), stack = TRUE)
image_browse(stack_img)

image_write(stack_img, path = paste0("./figs/CTRL_MUT_morphs_heatmaps_head.png"), format = "png")


####7.2. HEATMAPS - PC1-10 ####

PC_min <- tps3d(head_mesh, as.matrix(atlas_head_lm), PCA_head$shapes[[1]]$min, threads = 1)
open3d(zoom=0.75, userMatrix = dorsal, windowRect= c(0,0,1000,700))
shade3d(PC_min, color="grey")

PC_max <- tps3d(head_mesh, as.matrix(atlas_head_lm), PCA_head$shapes[[1]]$max, threads = 1)
open3d(zoom=0.75, userMatrix = ventral, windowRect= c(0,0,1000,700))
shade3d(PC_max, color="grey")

# lateral <- par3d()$userMatrix
# dorsal <- par3d()$userMatrix
# ventral <- par3d()$userMatrix
# rgl.close()
# 
# save(lateral, dorsal, ventral, file = "./figs/RGL_head_PCs_pos.rdata")

# Have a look at the heatmap, and change colour if needed
open3d(zoom=0.75, userMatrix = lateral, windowRect= c(0,0,1000,700))
x11(width=1.7, height=8)
meshDist(PC_max, PC_min, rampcolors = c("darkblue", "blue", "blue", "white", "red", "red", "darkred"), sign = TRUE)

load("./figs/RGL_head_PCs_pos.rdata")

n_dimensions <- 10 # number of PCs to include in the figure, decide depending on % of variance explained in PCs
summary(PCA_head)

# This loop automatically positions faces in dorsal, ventral, and lateral views.
for (i in 1:n_dimensions){
  PC_min <- tps3d(head_mesh, as.matrix(atlas_head_lm), PCA_head$shapes[[i]]$min, threads = 1)
  PC_max <- tps3d(head_mesh, as.matrix(atlas_head_lm), PCA_head$shapes[[i]]$max, threads = 1)
  PC <- paste0("PC", i)
  #dorsal views
  open3d(zoom=0.75, userMatrix = dorsal, windowRect= c(0,0,1000,700))
  shade3d(PC_max, color="grey")
  rgl.snapshot(paste0("./figs/pc_morphs/head/",PC,"_","max_dorsal.png"), top = TRUE )
  clear3d()
  rgl.close()
  
  #ventral views
  open3d(zoom=0.75, userMatrix = ventral, windowRect= c(0,0,1000,700))
  shade3d(PC_max, color="grey")
  rgl.snapshot(paste0("./figs/pc_morphs/head/",PC,"_","max_ventral.png"), top = TRUE )
  clear3d()
  rgl.close()
  
  #lateral views
  open3d(zoom=0.75, userMatrix = lateral, windowRect= c(0,0,1000,700)) 
  shade3d(PC_max, color="grey")
  rgl.snapshot(paste0("./figs/pc_morphs/head/",PC,"_","max_lateral.png"), top = TRUE )
  clear3d()
  rgl.close()
  
  #dorsal views
  open3d(zoom=0.75, userMatrix = dorsal, windowRect= c(0,0,1000,700))
  shade3d(PC_min, color="grey")
  rgl.snapshot(paste0("./figs/pc_morphs/head/",PC,"_","min_dorsal.png"), top = TRUE )
  clear3d()
  rgl.close()
  
  #ventral views
  open3d(zoom=0.75, userMatrix = ventral, windowRect= c(0,0,1000,700))
  shade3d(PC_min, color="grey")
  rgl.snapshot(paste0("./figs/pc_morphs/head/",PC,"_","min_ventral.png"), top = TRUE )
  clear3d()
  rgl.close()
  
  #lateral views
  open3d(zoom=0.75, userMatrix = lateral, windowRect= c(0,0,1000,700)) 
  shade3d(PC_min, color="grey")
  rgl.snapshot(paste0("./figs/pc_morphs/head/",PC,"_","min_lateral.png"), top = TRUE )
  clear3d()
  rgl.close()
  
  #Heatmaps
  open3d(zoom=0.75, userMatrix = dorsal, windowRect= c(0,0,1000,700))
  meshDist(PC_max, PC_min, rampcolors = c("darkblue", "blue", "blue", "white", "red", "red", "darkred"), sign = TRUE)
  rgl.snapshot(paste0("./figs/pc_morphs/head/heat_",PC,"_","dorsal.png"), top = TRUE )
  clear3d()
  rgl.close()
  
  open3d(zoom=0.75, userMatrix = ventral, windowRect= c(0,0,1000,700))
  meshDist(PC_max, PC_min, rampcolors = c("darkblue", "blue", "blue", "white", "red", "red", "darkred"), sign = TRUE)
  rgl.snapshot(paste0("./figs/pc_morphs/head/heat_",PC,"_","ventral.png"), top = TRUE )
  clear3d()
  rgl.close()
  
  open3d(zoom=0.75, userMatrix = lateral, windowRect= c(0,0,1000,700))
  x11(width=1.7, height=8)
  meshDist(PC_max, PC_min, rampcolors = c("darkblue", "blue", "blue", "white", "red", "red", "darkred"), sign = TRUE)
  rgl.snapshot(paste0("./figs/pc_morphs/head/heat_",PC,"_","lateral.png"), top = TRUE)
  # dev.copy(pdf, paste0("./figs/pc_morphs/head/heat_",PC,"_scale.pdf"), width=2, height=9.5)
  dev.print(pdf, paste0("./figs/pc_morphs/head/heat_",PC,"_scale.pdf"), width=2, height=9.5)
  dev.off(dev.prev())
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
PC_min_img_sd <- image_read(paste0("./figs/pc_morphs/head/PC",i,"_min_lateral.png"))
PC_max_img_sd <- image_read(paste0("./figs/pc_morphs/head/PC",i,"_max_lateral.png"))
PC_min_img_inf <- image_read(paste0("./figs/pc_morphs/head/PC",i,"_min_ventral.png"))
PC_max_img_inf <- image_read(paste0("./figs/pc_morphs/head/PC",i,"_max_ventral.png"))
PC_HT_img_sup <- image_read(paste0("./figs/pc_morphs/head/heat_PC",i,"_dorsal.png"))
PC_HT_img_sd <- image_read(paste0("./figs/pc_morphs/head/heat_PC",i,"_lateral.png"))
PC_HT_img_inf <- image_read(paste0("./figs/pc_morphs/head/heat_PC",i,"_ventral.png"))

PC_col <- c(PC_min_img_sup,PC_max_img_sup,PC_min_img_inf,PC_max_img_inf,PC_min_img_sd,PC_max_img_sd,PC_HT_img_sup,PC_HT_img_inf,PC_HT_img_sd)

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
imgs <- image_annotate(imgs, paste0("PCmin"), font = "Times",
                       location = "+0+495", size = 32)
imgs <- image_annotate(imgs, paste0("PCmax"), font = "Times",
                       location = "+0+695", size = 32)
imgs <- image_annotate(imgs, paste0("PCmin"), font = "Times",
                       location = "+0+895", size = 32)
imgs <- image_annotate(imgs, paste0("PCmax"), font = "Times",
                       location = "+0+1095", size = 32)
imgs <- image_annotate(imgs, paste0("Heatmap"), font = "Times",
                       location = "+0+1295", size = 32)
imgs <- image_annotate(imgs, paste0("Heatmap"), font = "Times",
                       location = "+0+1495", size = 32)
imgs <- image_annotate(imgs, paste0("Heatmap"), font = "Times",
                       location = "+0+1705", size = 32)

image_browse(imgs)

for (i in 2:n_dimensions){
  PC_min_img_sup <- image_read(paste0("./figs/pc_morphs/head/PC",i,"_min_dorsal.png"))
  PC_max_img_sup <- image_read(paste0("./figs/pc_morphs/head/PC",i,"_max_dorsal.png"))
  PC_min_img_sd <- image_read(paste0("./figs/pc_morphs/head/PC",i,"_min_lateral.png"))
  PC_max_img_sd <- image_read(paste0("./figs/pc_morphs/head/PC",i,"_max_lateral.png"))
  PC_min_img_inf <- image_read(paste0("./figs/pc_morphs/head/PC",i,"_min_ventral.png"))
  PC_max_img_inf <- image_read(paste0("./figs/pc_morphs/head/PC",i,"_max_ventral.png"))
  PC_HT_img_sup <- image_read(paste0("./figs/pc_morphs/head/heat_PC",i,"_dorsal.png"))
  PC_HT_img_sd <- image_read(paste0("./figs/pc_morphs/head/heat_PC",i,"_lateral.png"))
  PC_HT_img_inf <- image_read(paste0("./figs/pc_morphs/head/heat_PC",i,"_ventral.png"))
  
  
  PC_row <- c(PC_min_img_sup,PC_max_img_sup,PC_min_img_inf,PC_max_img_inf,PC_min_img_sd,PC_max_img_sd,PC_HT_img_sup,PC_HT_img_inf,PC_HT_img_sd)
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

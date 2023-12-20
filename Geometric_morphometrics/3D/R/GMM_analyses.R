#### Load R packages ####
# install.packages('devtools')
# library(devtools)
# install_github("marta-vidalgarcia/morpho.tools.GM", force = TRUE)
# install.packages(c('rgl', 'geomorph' 'devtools', 'Morpho', 'Rvcg', 'magick', 'Evomorph', 'ggplot2', 'vegan', 'factoextra', 'gt'))

library(rgl)
library(geomorph)
library(morpho.tools.GM)
library(Morpho)
library(Rvcg)
library(magick)
library(Evomorph)
library(ggplot2)
library(vegan)
library(ggbiplot)
library(factoextra)
library(gt)
library(abind)

library(stringr)
library(dplyr)
library(DescTools)

library(cowplot)


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


get_palette <- function() {
  return(c('#bbbbbb', '#0177bb', '#13783d','#989936', '#882256'))
}


plot_pc_table <- function(PCA_head, classifiers) {
  par(mfrow=c(3,2))
  plot(PCA_head, pch = 16, axis1 = 1, axis2 = 2, col = pal[as.numeric(classifiers$treatment)], cex = 1.25)
  ordiellipse(PCA_head, classifiers$treatment, kind="sd",conf=0.95, border = pal,
              draw = "polygon", alpha = 0, lty = 1)
  legend("topright", pch = 16, col = pal, legend = levels(classifiers$treatment))
  title("PC1 ~ PC2")
  
  plot(PCA_head, pch = 16, axis1 = 1, axis2 = 3, col = pal[as.numeric(classifiers$treatment)], cex = 1.25)
  ordiellipse(PCA_head$x[,c(1,3)], classifiers$treatment, kind="sd",conf=0.95, border = pal,
              draw = "polygon", alpha = 0, lty = 1)
  legend("bottomright", pch = 16, col = pal, legend = levels(classifiers$treatment))
  title("PC1 ~ PC3")
  
  plot(PCA_head, pch = 16, axis1 = 1, axis2 = 4, col = pal[as.numeric(classifiers$treatment)], cex = 1.25)
  ordiellipse(PCA_head$x[,c(1, 4)], classifiers$treatment, kind="sd",conf=0.95, border = pal,
              draw = "polygon", alpha = 0, lty = 1)
  legend("topright", pch = 16, col = pal, legend = levels(classifiers$treatment))
  title("PC1 ~ PC4")
  
  plot(PCA_head, pch = 16, axis1 = 2, axis2 = 3, col = pal[as.numeric(classifiers$treatment)], cex = 1.25)
  ordiellipse(PCA_head$x[,c(2,3)], classifiers$treatment, kind="sd",conf=0.95, border = pal,
              draw = "polygon", alpha = 0, lty = 1)
  legend("bottomleft", pch = 16, col = pal, legend = levels(classifiers$treatment))
  title("PC2 ~ PC3")
  
  plot(PCA_head, pch = 16, axis1 = 2, axis2 = 4, col = pal[as.numeric(classifiers$treatment)], cex = 1.25)
  ordiellipse(PCA_head$x[,c(2,4)], classifiers$treatment, kind="sd",conf=0.95, border = pal,
              draw = "polygon", alpha = 0, lty = 1)
  legend("topright", pch = 16, col = pal, legend = levels(classifiers$treatment))
  title("PC2 ~ PC4")
  
  plot(PCA_head, pch = 16, axis1 = 3, axis2 = 4, col = pal[as.numeric(classifiers$treatment)], cex = 1.25)
  ordiellipse(PCA_head$x[,3:4], classifiers$treatment, kind="sd",conf=0.95, border = pal,
              draw = "polygon", alpha = 0, lty = 1)
  legend("topright", pch = 16, col = pal, legend = levels(classifiers$treatment))
  title("PC3 ~ PC4")
}


plot_pc_heatmap <- function(face_mesh, atlas_head_lm, PCA_head, dimensions) {
  for (i in 1:n_dimensions){
    PC_min <- tps3d(face_mesh, as.matrix(atlas_head_lm), PCA_head$shapes[[i]]$min, threads = 1)
    PC_max <- tps3d(face_mesh, as.matrix(atlas_head_lm), PCA_head$shapes[[i]]$max, threads = 1)
    PC <- paste0("PC", i)
    
    #PCmax
    open3d(zoom=0.75, userMatrix = heatmap_frontal, windowRect= c(0,0,1000,700))
    shade3d(PC_max, color="grey", specular='black')
    rgl.snapshot(paste0("./figs/pc_morphs/",PC,"_","max.png"), top = TRUE )
    clear3d()
    rgl::close3d()
    
    #PCmin
    open3d(zoom=0.75, userMatrix = heatmap_frontal, windowRect= c(0,0,1000,700))
    shade3d(PC_min, color="grey", specular='black')
    rgl.snapshot(paste0("./figs/pc_morphs/",PC,"_","min.png"), top = TRUE )
    clear3d()
    rgl::close3d()
    
    #Heatmaps
    open3d(zoom=0.75, userMatrix = heatmap_frontal, windowRect= c(0,0,1000,700))
    x11(width=1.7, height=8)
    meshDist(PC_min, PC_max, rampcolors = c("darkblue", "blue", "white", "red", "darkred"), sign = TRUE)
    rgl.snapshot(paste0("./figs/pc_morphs/heat_",PC,".png"), top = TRUE)
    # dev.copy(pdf, paste0("./figs/pc_morphs/head/heat_",PC,"_scale.pdf"), width=2, height=9.5)
    dev.print(pdf, paste0("./figs/pc_morphs/heat_",PC,"_scale.pdf"), width=2, height=9.5)
    dev.off()
    clear3d()
    rgl::close3d()
    
    rm(PC_max, PC_min, PC) 
  }
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


#### 1.1 Import landmark data ####
# 1. Import all fcsv files into separate arrays (all specimens for east LM set)
# Traditional landmarks
setwd("./lm_data/Landmarks/") # directory with all the fcsv files in it
LMs <- fcsv2array(string_del = "_Fiducials") # cleanup filenames, from morphotools library

#### 1.2 Load and match classifiers ####
sample_names <- list.files(path = getwd(), pattern="fcsv$") # get all fiducial files
sample_names <- str_remove(sample_names, "_Fiducials.fcsv$") # remove file extension
classifiers_unord <- data.frame(sample_names) # make dataframe
colnames(classifiers_unord) <- "id" # rename first column
# pull treatment from filename
classifiers_unord$treatment <- str_match(classifiers_unord$id, "_(.+)_")[,2] 
# rename treatments
classifiers_unord <- classifiers_unord %>% mutate(treatment = case_when(treatment =='ctr' ~ 'DMSO',
                                                                        treatment =='exp' ~ 'Triple',
                                                                        treatment =='expredo' ~ 'Triple',
                                                                        treatment =='LY' ~ 'LY294002',
                                                                        treatment =='U0' ~ 'U0126',
                                                                        treatment =='U73' ~ 'U73122'))
# note, there shouldn't be any 'exp' anymore, they were all redone and moved away

# make treatment a factor and set level order for graphing
classifiers_unord$treatment <- factor(classifiers_unord$treatment, levels = c('DMSO', 'U0126', 'LY294002', 'U73122', 'Triple'))
row.names(classifiers_unord) <- classifiers_unord$id


# match the landmark sample names to the classifiers csv
classifiers <- classifiers_unord[match(dimnames(LMs)[[3]], row.names(classifiers_unord)),]


# Curve semilandmarks
setwd("../Semi_Curves/")
# finds all files with the pattern in the name
curve_semis_center <- fcsv2array(pattern = "*center*", string_del = "_center_semi-curve")
curve_semis_L <- fcsv2array(pattern = "*_L_semi-curve*", string_del = "_L_semi-curve")
curve_semis_R <- fcsv2array(pattern = "*_R_semi-curve*", string_del = "_R_semi-curve")

# for some reason a bunch of the L and R semi-curve points are in reverse order
# need to reverse the LY, U0, and U73 samples
reverse_ids <- which(classifiers$treatment %in% c("LY294002", 'U0126', 'U73122', 'Triple'))
# reverses the order of the 9 landmarks, but not their xyz order
curve_semis_L_reversed <- apply(curve_semis_L, c(2,3), function (x) Rev(x, 1))
curve_semis_R_reversed <- apply(curve_semis_R, c(2,3), function (x) Rev(x, 1))
# rewrite the original landmarks with reversed order landmarks
curve_semis_L[,,reverse_ids] <- curve_semis_L_reversed[,,reverse_ids]
curve_semis_R[,,reverse_ids] <- curve_semis_R_reversed[,,reverse_ids]

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
# combine LMs, should have an array of 51 landmarks and ?? samples
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
curve_c_left <- c(9, semis_center[c(1,2)])
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
# this csv will be used in the asymmetry analysis!
write.csv(curveslide_all, "./lm_data/curveslide.csv", row.names = FALSE)


#### 1.5 Surface semi-landmarks ####
head_surface.lm <- (dim(LMs)[1]+dim(curve_semis_center)[1]+dim(curve_semis_L)[1]+dim(curve_semis_R)[1]+1):dim(head_array)[1]


#### 2 GPA, PCA, mean shape analysis ####
#### 2.1 GPA of mean shape ####
# Use this to find out who is the closest to the shape sphere centroid
GPA_head <- geomorph::gpagen(A = head_array, curves = as.matrix(curveslide_all), 
                             surfaces = head_surface.lm)
# check for outliers
# looks fine, LY22 I think is a true outlier. Not sure why LY4 is way out there, it looks fine to me
outlier <- plotOutliers_percentile(A = GPA_head$coords, percentile = 0.99, save.plot = FALSE)
min(outlier$All_Proc_d$`Proc. d. from mean`) # actually I got lazy & did it this easier way, same concept as above though
row.names(outlier$All_Proc_d[which(outlier$All_Proc_d$`Proc. d. from mean` == min(outlier$All_Proc_d$`Proc. d. from mean`)),])


test3 <- head_array[,, which(dimnames(head_array)[[3]] %in% c("chick_LY_4", 'chick_ctr_23')), drop=FALSE]
plot_test <- cbind(test2, test3)

open3d(zoom = 0.75, windowRect = c(0, 0, 700, 700))
plotAllSpecimens(A = test3, label=T)
rgl::close3d()

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
rgl::close3d()

save(frontal, file = "./lm_data/RGL_head_pos.rdata")
load("./lm_data/RGL_head_pos.rdata")

# Frontal view
# this figure is in the manuscript
open3d(zoom = 0.75, userMatrix = frontal, windowRect = c(0, 0, 1000, 700)) 
rgl::shade3d(face_mesh, color = "gray", alpha =1, specular = "black" )
rgl::plot3d(atlas_head_lm[head_fixed.lm,], aspect = "iso", type = "s", size=.8, col = "black", add = T)
rgl::plot3d(atlas_head_lm[head_curves.lm,], aspect = "iso", type = "s", size=0.8, col = "orange", add = T)
rgl::plot3d(atlas_head_lm[head_surface.lm,], aspect = "iso", type = "s", size=0.8, col = "turquoise2", add = T)
rgl.snapshot("./figs/head_LM_frontal.png", top = TRUE)
# writeASY(title = "./figs/head_LM_frontal",outtype = 'pdf', prc = FALSE)
rgl::close3d()


#### 2.3. Mean shape PCA plots ####
# perform PCA
PCA_head <- gm.prcomp(GPA_head$coords)


# Save output, delete file if it exists
if (file.exists("./output/PCA_head_shape_coords.txt")) {
  file.remove("./output/PCA_head_shape_coords.txt")
}
if (!dir.exists('./output/')) dir.create('./output/')
cat("PCA shape variables raw", capture.output(summary(PCA_head)), 
    file="./output/PCA_head_shape_coords.txt", sep="\n", append=TRUE)

pal <- get_palette() # get color map used for figures
classifiers$treatment <- as.character(classifiers$treatment)
classifiers$treatment <- as.factor(classifiers$treatment)
classifiers$treatment <- factor(classifiers$treatment, levels = c('DMSO', 'U0126', 'LY294002', 'U73122', 'Triple'))


# Plot PC1 vs PC2
# these are used in the manuscript
# this one generates the histograms on the margins
plot_df <- as.data.frame(PCA_head$x)

pdf("./figs/PCA_mean_shape_treatment_histogram.pdf", width = 8.25, height = 6)

p <- ggplot(plot_df, aes(x = plot_df[,1], y =  plot_df[,2], color = classifiers$treatment)) + geom_point() + scale_color_manual(values = pal) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
px <- ggplot(plot_df, aes(x=plot_df[,1], color=classifiers$treatment)) + geom_density() + scale_color_manual(values = palette()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
py <- ggplot(plot_df, aes(x=plot_df[,2], color=classifiers$treatment)) + geom_density() + scale_color_manual(values = palette()) + coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p %>%
  insert_xaxis_grob(px, grid::unit(1, "in"), position = "top") %>%
  insert_yaxis_grob(py, grid::unit(1, "in"), position = "right") %>%
  ggdraw()

legend("topright", pch = 16, col = palette(), legend = levels(classifiers_filter))
title("PCA of shape coordinates")
dev.off()
dev.off()

# this figure is in the manuscript
pdf("./figs/PCA_mean_shape_treatment.pdf", width = 8.25, height = 6)
plot(PCA_head, pch = 16, col = pal[as.numeric(classifiers$treatment)], cex = 1.25) #, xlim=c(-.1,.15), ylim=c(-.09,.09))
ordiellipse(PCA_head, classifiers$treatment, kind="sd",conf=0.95, border = pal,
            draw = "polygon", alpha = 0, lty = 1)
legend("topleft", pch = 16, col = pal, legend = levels(classifiers$treatment))
# title("PCA of shape coordinates - ctrl vs treatment")
dev.off()

# create a table of PC1-4 comparisons
# PC1 - PC4
pdf("./figs/PCA_head_treatment_PC1-4.pdf", width = 8, height = 12)
pc_table <- plot_pc_table(PCA_head, classifiers)
dev.off()

# plot a scree plot
PCA_comp <- PCA_head
class(PCA_comp) <- "princomp"

pdf("./figs/PCA_head_shape_scree_plot.pdf", height = 5, width = 5)
pca_scree <- fviz_eig(PCA_comp, addlabels=TRUE, hjust = -0.3,
                      barfill="darkgrey", barcolor ="black",
                      linecolor ="blue") + ylim(0, 85) + 
  theme_classic()
print(pca_scree)
dev.off()


#### 2.4 Centroid size, Procruste's distance, allometry ####
Pdist <- ShapeDist(GPA_head$coords, GPA_head$consensus)

# make the dataframe for better ggplot plotting
gdf_head <- geomorph.data.frame(GPA_head, treatment = classifiers$treatment, Pdist = Pdist)

ggplot_df <- as.data.frame(cbind(as.character(gdf_head$Csize), 
                                 as.character(gdf_head$treatment), 
                                 as.character(gdf_head$Pdist)))
colnames(ggplot_df) <- c("Csize", "treatment", "Pdist")
row.names(ggplot_df) <- dimnames(gdf_head$coords)[[3]]

ggplot_df$treatment <- as.factor(ggplot_df$treatment)
ggplot_df$treatment <- factor(ggplot_df$treatment, levels = c('DMSO', 'U0126', 'LY294002', 'U73122', 'Triple'))

ggplot_df$Csize <- as.numeric(as.character(ggplot_df$Csize))
ggplot_df$Pdist <- as.numeric(as.character(ggplot_df$Pdist))

# plot centroid size
pdf("./figs/Csize_boxplot_HEAD.pdf", width = 6.5, height = 6.5)
ggplot(ggplot_df, aes(x=treatment, y=log(Csize), fill=treatment)) + 
  geom_boxplot(alpha = 0.7, outlier.shape=NA) +
  scale_fill_manual(values = pal) +
  geom_jitter(width = 0.1, size = 1.25) +
  theme(plot.title = element_text(size = 20, face = "bold"), axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"), legend.text=element_text(size=12), 
        legend.title=element_text(size=15, face = "bold"))
dev.off()

# plot Procrustes' distance
pdf("./figs/Pdist_treatment.pdf", width = 6.5, height = 6.5)
ggplot(ggplot_df, aes(Pdist, fill = treatment)) +
  scale_fill_manual(values=pal) + geom_density(alpha = 0.25) + 
  ggtitle("Procrustes distances head shape") + xlab("Procrustes distance") + ylab("Relative density") +
  theme(plot.title = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"), legend.text=element_text(size=10), 
        legend.title=element_text(size=12, face = "bold"))
dev.off()

# allometry
allometry_all <- procD.lm(coords ~ Csize, data = gdf_head, iter = 999, RRPP = TRUE)
summary(allometry_all)
# Save output, delete file if it exists
if (file.exists("./output/HEAD_allometry_all.txt")) {
  file.remove("./output/HEAD_allometry_all.txt")
}
cat("Allometry (coords ~ Csize)", capture.output(summary(allometry_all)), 
    file="./output/HEAD_allometry_all.txt", sep="\n", append=TRUE)

shape_residuals <- allometry_all$residuals # let's save this for later. We will run a PCA again on the residuals
write.csv(shape_residuals, "./output/HEAD_shape_residuals.csv")

# allometry by treatment
allometry_treatment <- procD.lm(coords ~ Csize * treatment, data = gdf_head, iter = 999, RRPP = TRUE)
summary(allometry_treatment) # both significant, weak interaction
allometry_ph <- pairwise(allometry_treatment, groups = gdf_head$treatment)
summary(allometry_ph) # no pairwise differences!

# Delete file if it exists
if (file.exists("./output/HEAD_allometry_treatment.txt")) {
  file.remove("./output/HEAD_allometry_treatment.txt")
}

cat("Allometry * treatment (coords ~ Csize * treatment)", capture.output(summary(allometry_treatment)), 
    file="./output/HEAD_allometry_treatment.txt", sep="\n", append=TRUE)


#### 2.5 Mean shape ANOVA and plots ####
# ANOVA of mean shape - this pvalue is in manuscript
treatment <- procD.lm(coords ~ treatment, data = gdf_head, RRPP = TRUE)
summary(treatment)
summary(treatment, test.type = "var")

treatment_ph <- pairwise(treatment, groups = gdf_head$treatment)
summary(treatment_ph)
summary(treatment_ph, test.type = "var", confidence = 0.95, stat.table = TRUE)

# Save and delete file if it exists
if (file.exists("./output/ANOVA_HEAD_shape_treatment.txt")) {
  file.remove("./output/ANOVA_HEAD_shape_treatment.txt")
}
cat("ANOVA shape - treatment (coords ~ treatment)", capture.output(summary(treatment)), 
    file="./output/ANOVA_HEAD_shape_treatment.txt", sep="\n", append=TRUE)

# ANOVA of centroid size by treatment
one_way <- aov(log(gdf_head$Csize) ~ gdf_head$treatment)
summary(one_way) # is this the same as the allometry analysis above? double check w/ Marta

# Save and delete file if it exists
if (file.exists("./output/ANOVA_HEAD_Csize_treatment.txt")) {
  file.remove("./output/ANOVA_HEAD_Csize_treatment.txt")
}
cat("ANOVA Csize - treatment (Csize ~ treatment)", capture.output(summary(one_way)), 
    file="./output/ANOVA_HEAD_Csize_treatment.txt", sep="\n", append=TRUE)


#### 2.5.1 Mean shape heatmaps ####
# load mesh of the specimen closest to the mean shape
setwd('./lm_data/Meshes/')
face_mesh <- get_dec_mesh('face')
setwd("../../")

# get landmarks for each treatment group
triple_coords <- gdf_head$coords[,,which(gdf_head$treat == "Triple")]
ctrl_coords <- gdf_head$coords[,,which(gdf_head$treat == "DMSO")]
U73_coords <- gdf_head$coords[,,which(gdf_head$treat == "U73122")]
U0_coords <- gdf_head$coords[,,which(gdf_head$treat == "U0126")]
LY_coords <- gdf_head$coords[,,which(gdf_head$treat == "LY294002")]

# we are going to calculate the mean shape on the shape coordinates. Everything will be in Procrustes dist.
triple_mean_shape <- mshape(triple_coords)
ctrl_mean_shape <- mshape(ctrl_coords)
U73_mean_shape <- mshape(U73_coords)
U0_mean_shape <- mshape(U0_coords)
LY_mean_shape <- mshape(LY_coords)

# Create morphed meshes
triple_mesh <- tps3d(face_mesh, as.matrix(atlas_head_lm), triple_mean_shape, threads = 1)
ctrl_mesh <- tps3d(face_mesh, as.matrix(atlas_head_lm), ctrl_mean_shape, threads = 1)
U73_mesh <- tps3d(face_mesh, as.matrix(atlas_head_lm), U73_mean_shape, threads = 1)
U0_mesh <- tps3d(face_mesh, as.matrix(atlas_head_lm), U0_mean_shape, threads = 1)
LY_mesh <- tps3d(face_mesh, as.matrix(atlas_head_lm), LY_mean_shape, threads = 1)


# plot heatmap and setup a new reference view
open3d(zoom = 0.75, userMatrix = frontal, windowRect = c(0, 0, 1000, 700)) 
meshDist(ctrl_mesh, triple_mesh, from=-.01, to=0.015, rampcolors = c("blue", "white", "red"), sign = TRUE)
meshDist(ctrl_mesh, U0_mesh, from=-.01, to=0.015, rampcolors = c("blue", "white", "red"), sign = TRUE)
meshDist(ctrl_mesh, LY_mesh, from=-.01, to=0.015, rampcolors = c("blue", "white", "red"), sign = TRUE)
meshDist(ctrl_mesh, U73_mesh, from=-.01, to=0.015, rampcolors = c("blue", "white", "red"), sign = TRUE)
meshDist(triple_mesh, U73_mesh, from=-.01, to=0.015, rampcolors = c("blue", "white", "red"), sign = TRUE)


# get warp for certain PC values rather than mean shapes
PC = PCA_head$x[,1] # get pc1
PC2 = PCA_head$x[,2] # get pc2

PC1_ctrl_extreme <- tps3d(face_mesh, as.matrix(atlas_head_lm), shape.predictor(GPA_head$coords, x= PC, Intercept = FALSE, pred1 = .1)[[1]], threads=1)
PC1_ctrl <- tps3d(face_mesh, as.matrix(atlas_head_lm), shape.predictor(GPA_head$coords, x= PC, Intercept = FALSE, pred1 = .05)[[1]], threads=1)
PC1_trt <- tps3d(face_mesh, as.matrix(atlas_head_lm), shape.predictor(GPA_head$coords, x= PC, Intercept = FALSE, pred1 = -.1)[[1]], threads=1)
PC1_trt_mild <- tps3d(face_mesh, as.matrix(atlas_head_lm), shape.predictor(GPA_head$coords, x= PC, Intercept = FALSE, pred1 = -.05)[[1]], threads=1)
PC1_trt_milder <- tps3d(face_mesh, as.matrix(atlas_head_lm), shape.predictor(GPA_head$coords, x= PC, Intercept = FALSE, pred1 = 0)[[1]], threads=1)

PC2_ctrl <- tps3d(face_mesh, as.matrix(atlas_head_lm), shape.predictor(GPA_head$coords, x= PC2, Intercept = FALSE, pred1 = -.025)[[1]], threads=1)
PC2_trt <- tps3d(face_mesh, as.matrix(atlas_head_lm), shape.predictor(GPA_head$coords, x= PC2, Intercept = FALSE, pred1 = .05)[[1]], threads=1)


open3d(zoom = 0.75, userMatrix = frontal, windowRect = c(0, 0, 1000, 700)) 

meshDist(PC1_ctrl_extreme, PC1_trt, from= -.015, to= 0.03, rampcolors = c("blue", "white", "red"), sign = TRUE)
meshDist(PC1_ctrl, PC1_trt, from= -.015, to= 0.03, rampcolors = c("blue", "white", "red"), sign = TRUE)
meshDist(PC1_ctrl, PC1_trt_mild, from= -.015, to= 0.03, rampcolors = c("blue", "white", "red"), sign = TRUE)
meshDist(PC1_ctrl, PC1_trt_milder, from= -.015, to= 0.03, rampcolors = c("blue", "white", "red"), sign = TRUE)
meshDist(PC1_ctrl, ctrl_mesh, from= -.015, to= 0.03, rampcolors = c("blue", "white", "red"), sign = TRUE)

meshDist(PC2_ctrl, PC2_trt, from=-0.015, to=0.03, rampcolors = c("blue", "white", "red"), sign = TRUE)




# you need to click and move mouse around in the figure until it is en-face with
# the window, then run the line below to save the orientation
heatmap_frontal <- par3d()$userMatrix
rgl::close3d()
save(heatmap_frontal, file = "./lm_data/RGL_heat_head_pos.rdata")

load("./lm_data/RGL_heat_head_pos.rdata")

# plot heatmap - this is used in manuscript
open3d(zoom = 0.75, userMatrix = heatmap_frontal, windowRect = c(0, 0, 1000, 700)) 
pdf("./figs/heatmap_treatment_legend.pdf", width = 2.5, height = 6.5)
meshDist(ctrl_mesh, U0_mesh, from=-.01, to=0.015, rampcolors = c("blue", "white", "red"), sign = TRUE)
rgl.snapshot("./figs/Heatmap_U0.png", top = TRUE) # this one captures 3d output
rgl::close3d() # this one captures the heatmap legend as pdf
dev.off()

# plot the two first PC heatmaps for mean shape
open3d(zoom = 0.75, userMatrix = heatmap_frontal, windowRect = c(0, 0, 1000, 700)) 
pdf("./figs/heatmap_PC1_pt05_to_-0pt1_legend.pdf", width = 2.5, height = 6.5)
meshDist(PC1_ctrl, PC1_trt, from= -.015, to= 0.03, rampcolors = c("blue", "white", "red"), sign = TRUE)
rgl.snapshot("./figs/heatmap_PC1_pt05_to_-0pt1.png", top = TRUE) # this one captures 3d output
rgl::close3d() # this one captures the heatmap legend as pdf
dev.off()

open3d(zoom = 0.75, userMatrix = heatmap_frontal, windowRect = c(0, 0, 1000, 700)) 
pdf("./figs/heatmap_PC2_-0pt025_to_0pt05_legend.pdf", width = 2.5, height = 6.5)
meshDist(PC2_ctrl, PC2_trt, from= -.015, to= 0.03, rampcolors = c("blue", "white", "red"), sign = TRUE)
rgl.snapshot("./figs/heatmap_PC2_-0pt025_to_0pt05.png", top = TRUE) # this one captures 3d output
rgl::close3d() # this one captures the heatmap legend as pdf
dev.off()

# plot the mean triple treated mesh - this is used in manuscript
open3d(zoom=0.75, userMatrix = heatmap_frontal, windowRect = c(0,0,1000,700)) 
shade3d(triple_mesh, color="gray", alpha=1, specular='black')
rgl.snapshot("./figs/mean_shape_triple.png", top = TRUE)
rgl::close3d()

# plot the mean dmso treated mesh - this is used in manuscript
open3d(zoom=0.75, userMatrix = heatmap_frontal, windowRect = c(0,0,1000,700)) 
shade3d(ctrl_mesh, color="gray", alpha=1, specular='black')
rgl.snapshot("./figs/mean_shape_ctrl.png", top = TRUE)
rgl::close3d()

#Create page with triple vs ctrl morphs & heatmaps
Morph_ctrl <- image_read(paste0("./figs/mean_shape_ctrl.png"))
Morph_triple <- image_read(paste0("./figs/mean_shape_triple.png"))
Heatmap <- image_read(paste0("./figs/Heatmap_treatment.png"))

Morph_ctrl <- image_annotate(Morph_ctrl, "control", font = "times", location = "+400+20", size = 100)
Morph_triple <- image_annotate(Morph_triple, "triple", font = "times", location = "+425+20", size = 100)
Heatmap <- image_annotate(Heatmap, "heatmap", font = "times", location = "+350+20", size = 100)

stack_img <- image_append(image_scale(c(Morph_ctrl, Morph_triple, Heatmap), "x200"), stack = TRUE)
image_browse(stack_img)
image_write(stack_img, path = paste0("./figs/ctrl_triple_morphs_heatmaps.png"), format = "png")


#### 2.5.2 Heatmaps - PC1-10 ####
# load mesh of the specimen closest to the mean shape
setwd('./lm_data/Meshes/')
face_mesh <- get_dec_mesh('face')
setwd("../../")

PC_min <- tps3d(face_mesh, as.matrix(atlas_head_lm), PCA_head$shapes[[1]]$min, threads = 1)
open3d(zoom=0.75, userMatrix = heatmap_frontal, windowRect= c(0,0,1000,700))
shade3d(PC_min, color="grey", specular='black')
rgl::close3d()

PC_max <- tps3d(face_mesh, as.matrix(atlas_head_lm), PCA_head$shapes[[1]]$max, threads = 1)
open3d(zoom=0.75, userMatrix = heatmap_frontal, windowRect= c(0,0,1000,700))
shade3d(PC_max, color="grey", specular='black')
rgl::close3d()

# Have a look at the heatmap, and change colour if needed
open3d(zoom=0.75, userMatrix = heatmap_frontal, windowRect= c(0,0,1000,700))
x11(width=1.7, height=8)
meshDist(PC_max, PC_min, rampcolors = c("blue", "white", "red"), sign = FALSE)
rgl::close3d()

# create all the min and max heatmaps for pc1 to n_dimensions
n_dimensions <- 10 # number of PCs to include in the figure, decide depending on % of variance explained in PCs
plot_pc_heatmap(face_mesh, atlas_head_lm, PCA_head, n_dimensions)

#Create page with PCs
i <- 1
PC_min_img_sup <- image_read(paste0("./figs/pc_morphs/PC",i,"_min.png"))
PC_max_img_sup <- image_read(paste0("./figs/pc_morphs/PC",i,"_max.png"))
PC_HT_img_sup <- image_read(paste0("./figs/pc_morphs/heat_PC",i,".png"))

PC_col <- c(PC_min_img_sup,PC_max_img_sup,PC_HT_img_sup)

imgs <- image_append(image_scale(PC_col, "x200"), stack = TRUE)

imgs <- image_border(imgs, "white", "85x15")
imgs <- image_annotate(imgs, paste0("PC",i), font = "Times",
                       location = "+195+0", size = 60)
imgs <- image_annotate(imgs, paste0("PCmin"), font = "Times",
                       location = "+0+95", size = 32)
imgs <- image_annotate(imgs, paste0("PCmax"), font = "Times",
                       location = "+0+295", size = 32)
imgs <- image_annotate(imgs, paste0("Heatmap"), font = "Times",
                       location = "+0+495", size = 32)

for (i in 2:n_dimensions){
  PC_min_img_sup <- image_read(paste0("./figs/pc_morphs/PC",i,"_min.png"))
  PC_max_img_sup <- image_read(paste0("./figs/pc_morphs/PC",i,"_max.png"))
  PC_HT_img_sup <- image_read(paste0("./figs/pc_morphs/heat_PC",i,".png"))
  
  
  PC_row <- c(PC_min_img_sup,PC_max_img_sup,PC_HT_img_sup)
  img2 <-image_append(image_scale(PC_row, "x200"), stack = TRUE)
  img2 <-image_border(img2, "white", "0x15")
  img2 <-image_annotate(img2, paste0("PC",i), font = "Times", 
                       location = "+100+0", size = 60)
  
  imgs <-c(imgs,img2)
  
  imgs <-image_append(image_scale(imgs))
}
image_write(imgs, path = paste0("./figs/mean_treatment_heatmap_morphs_PC1-", i, ".png"), format = "png")

# this analysis continues in the 'Asymmetry.R' script
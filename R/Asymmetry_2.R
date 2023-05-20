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
detect.symmetry(GPA_geomorph$coords, sym.plane = "yz", plot = TRUE)
detect.symmetry(head_array[1:33,,], sym.plane = "yz", plot = TRUE)
detect.symmetry(head_array[,,], sym.plane = "yz", plot = TRUE)


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


# 3. MIRRORING LANDMARKS ####
# Mirror landmarks on both sides and generate two new arrays
# the non-symmetrical landmarks also remain the same

array_side1_mirrored <- GPA_geomorph$coords
array_side2_mirrored <- GPA_geomorph$coords

for (i in 1:dim(head_array)[3]){
  array_side1_mirrored[side.1,,i] <- mirror(GPA_geomorph$coords[side.1,,i], 
                                            v1 = GPA_geomorph$coords[9,,i], 
                                            v2 = GPA_geomorph$coords[10,,i],
                                            v3 = GPA_geomorph$coords[13,,i])
  array_side2_mirrored[side.2,,i] <- mirror(GPA_geomorph$coords[side.2,,i], 
                                            v1 = GPA_geomorph$coords[9,,i], 
                                            v2 = GPA_geomorph$coords[10,,i],
                                            v3 = GPA_geomorph$coords[13,,i])
}

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


# 4. GPA & PCA - both sides mirrored ####
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

# Left side

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





# Right side

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



# 5. Face integration CTRL vs Treatment ####

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


### Need to revise Morpho functions for Mirroring or write my own


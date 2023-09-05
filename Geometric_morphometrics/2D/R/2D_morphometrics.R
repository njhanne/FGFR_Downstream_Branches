library(geomorph)
library(tidyr)
library(dplyr)
library(abind)

# library(devtools)
# install_github("marta-vidalgarcia/morpho.tools.GM", force = TRUE)
library(morpho.tools.GM)
library(Evomorph)
# install_github("marta-vidalgarcia/symmetry")
# library(symmetry)

library(vegan)
library(wesanderson)
library(ggplot2)
library(cowplot)

### LANDMARKS ###
# 1   13  max-LNP junction inside
# 2   12  top nasal
# 3   11  nasal mid (maybe delete)
# 4   10  nasal base
# 5   9   glob base
# 6   8   glob-FNP mid
# 7       FNP midline
# 14  20  maxillary top
# 15  19  max-man junction
# 16  18  mandible
# 17      mandible midline
# 21  25  max-LNP junction outside
# 22  24  eye head junction top
# 23      top head midline
# 26  27  midway between 6-7 and 7-8
# 28  33  midway between 5-6 and 8-9
# 29  32  midway between 6-26 and 8-27
# 30  31  midway between 26-7 and 7-27

setwd("C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/Morphology/2D") # home
setwd("C:/Users/nhanne/Box/FGF_inhibitor_paper_5-26-2020/data/Morphology/2D") # work
setwd("./data/Morphology/2D") # laptop



#### Read and filter landmarks ####
# can just start here if you've run the analysis before
landmarks <- readland.tps('morpho_files/morpho_scaled.tps', specID = "ID")
classifiers <- readRDS('classifiers.Rds')

# This indented area only needs to be run the first time you perform this analysis or if new samples are added
  # read in landmarks. I formatted them as a 'TPS' style
  # but these don't have the 'scale' in them so we will have to get them from the ID names and then re-save with scale
  # landmarks <- readland.tps('morpho_new_tilt.txt', specID = "ID") # old
  landmarks <- readland.tps('morpho_files/morpho.txt', specID = "ID")
  landmarks_build <- readland.tps('morpho_files/U0_morpho.txt', specID = "ID")
  landmarks_build <- abind(landmarks_build, readland.tps('morpho_files/U73_morpho.txt', specID = "ID"))
  landmarks_build <- abind(landmarks_build, readland.tps('morpho_files/LY_morpho.txt', specID = "ID"))
  landmarks_build <- abind(landmarks_build, readland.tps('morpho_files/mix_morpho.txt', specID = "ID"))
  landmarks_build <- abind(landmarks_build, readland.tps('morpho_files/DMSO_morpho.txt', specID = "ID"))

  landmarks_old2 <- landmarks
  classifiers_old2 <- classifiers
  landmarks <- landmarks_build



  # create a list of classifiers. We will start with the 'id' from the TPS file and work from there...
  names <- dimnames(landmarks)[[3]]
  names <- substr(names, 1, nchar(names)-1) # deletes the /t (tab)
  classifiers <- lapply(names, function(x) {gsub("_24hr", "", x)})
  classifiers <- data.frame(x = matrix(unlist(classifiers), nrow=length(classifiers), byrow=TRUE))
  classifiers <- classifiers %>% separate(x, c('treatment', 'stage', 'magnification', 'embryo_num', 'orientation'), sep='_')
  classifiers <- classifiers %>%  mutate(orientation = replace_na(orientation, 'flat'))
  classifiers$treatment <- as.factor(classifiers$treatment)
  classifiers$stage <- as.factor(classifiers$stage)
  classifiers <- classifiers %>% unite(embryo, c("treatment", "stage", "embryo_num"), remove=FALSE)

  # scale from the dissection microscope (um / pixel)
  scale_2 <- 5000 / 1272
  scale_25 <- 4000 / 1569
  scale_32 <- 3000 / 1503

  classifiers <- classifiers %>%  mutate(magnification = case_when(magnification == 2 ~ scale_2,
                                                                   magnification == 2.5 ~ scale_25,
                                                                   magnification == 3.2 ~ scale_32))

  # write the tps file with adjusted scale in there and then reload data
  dimnames(landmarks)[[3]] <- names
  writeland.tps(landmarks, 'morpho_files/morpho_scaled.tps', scale = classifiers$magnification, specID=TRUE)
  saveRDS(classifiers, 'classifiers.Rds')
  landmarks_old <- landmarks
  classifiers_old <- classifiers

  landmarks <- readland.tps('morpho_files/morpho_scaled.tps', specID = "ID")
  classifiers <- readRDS('classifiers.Rds')

# delete x and y data of questionable landmarks
# landmarks_no_nasal <- landmarks[-c(3,11),,]
landmarks_FNP_only <- landmarks[-c(15:19,22:24),,] # this one is used for analysis

# this one isn't used, can rename it without '_alt' to see how it affects results
# landmarks_FNP_only_alt <- landmarks[-c(3,11,15:19,22:24),,]


# delete outlier samples (see below)
# in this case 'LY_st20_2_23' I had noted looked bad - remove
# 'U73_st20_2_1', 'U73_st20_2_4', 'LY_st20_2_25' are the wrong stage - remove
# 'U0_st20_3.2_27' and 'U0_st20_3.2_28' are very blury images - 
# but they aren't coming as outliers so I guess fine
names <- dimnames(landmarks)[[3]]
outlier_remove <- which(names == 'LY_st20_2_23')
outlier_remove <- append(outlier_remove, which(names == 'U73_st20_2_1_flat'))
outlier_remove <- append(outlier_remove, which(names == 'U73_st20_2_4_flat'))
outlier_remove <- append(outlier_remove, which(names == 'LY_st20_2_25'))

landmarks_FNP_only <- landmarks_FNP_only[,,-outlier_remove]
classifiers <- classifiers[-outlier_remove,]

# get the tilted samples so we can remove them
# this isn't needed anymore - I've removed them all from the morpho samples
# tilted <- which(classifiers$orientation == 'tilt')
# landmarks_FNP_only <- landmarks_FNP_only[,,-tilted]
# classifiers <- classifiers[-tilted,]

# get rid of doublemix
# this isn't needed anymore - I've removed them all from the morpho samples
# double <- which(classifiers$treatment == 'DoubleMix')
# landmarks_FNP_only <- landmarks_FNP_only[,,-double]
# classifiers <- classifiers[-double,]

# only look at st20 samples
# based on talking to Diane, I don't think these are actually the wrong stage, the images just are not named well
# stage_20 <- which(classifiers$stage == 'st18')
# landmarks_FNP_st20_only <- landmarks_FNP_only[,,-stage_20]
# classifiers_st20 <- classifiers[-stage_20,]


### Prepare data for GPA

# setup the semilandmarks along curve (26-33)
# for the '_alt' no nasalpit config
# [1,]      5    16     6
# [2,]      6    17     7
# [3,]      4    18     5
# [4,]      5    19    16
# [5,]     16    20     6
# [6,]      6    21    17
# [7,]     17    22     7
# [8,]      7    23     8

# this generates the curveslide csv, doesn't need to be re-run
  # curveslide <- define.sliders(landmarks_FNP_only2[,,1], nsliders=8)
  # curveslide <- as.data.frame(rbind(c(2,3,4), c(9,11,12), c(6,18,7), c(7,19,8), c(5,20,6), c(6,21,18), c(18,22,7), c(7,23,19), c(19,24,8), c(8,25,9)))
  # colnames(curveslide) <- c('before', 'slide', 'after')
  # write.csv(curveslide, 'curveslide.csv', row.names=FALSE, col.names = FALSE)

# curveslide_alt <- read.csv('curveslide_alt.csv') # without nasalpits, for the '_alt' landmark set
curveslide <- read.csv('curveslide.csv') # this one should be used

# setup the paired landmarks
gpa_test <- gpagen(landmarks_FNP_only)
plotAllSpecimens(gpa_test$coords, label=T, plot.param = list(pt.bg = "green", mean.cex=2,txt.pos=3, txt.cex=2))

# side.1 <- c(1:5,13,14,16,18:20)
# side.2 <- c(11:7,12,15,17,23:21)
# pairedLM_alt <- cbind(side.1, side.2)
pairedLM <- cbind(c(1:6,15,16,18,20:22), c(13:8,14,17,19,25:23))


#### Analysis w/ sliding semis ####

## perform GPA - one with all points, one with the points treated as semilandmarks
# will use Procrustes distance to determine semi placement
gpa_FNP_semi <- gpagen(landmarks_FNP_only, curves = curveslide)
# gpa_FNP_semi_20 <- gpagen(landmarks_FNP_st20_only, curves = curveslide)
# gpa_FNP_semi_alt<- gpagen(landmarks_FNP_only, curves = curveslide_alt)
# gpa_FNP_semi_20_alt <- gpagen(landmarks_FNP_st20_only, curves = curveslide_alt)



# check for outliers, they have already been removed above...
# outlier <- plotOutliers_percentile(A = gpa_FNP_semi$coords, percentile = 0.99, save.plot = FALSE)
# outlier <- plotOutliers_percentile(A = gpa_FNP_semi_20$coords, percentile = 0.99, save.plot = FALSE)

## Examine centroid size and Pdistance
PDist <- ShapeDist(gpa_FNP_semi$coords, gpa_FNP_semi$consensus)
# PDist_20 <- ShapeDist(gpa_FNP_semi_20$coords, gpa_FNP_semi_20$consensus)

# need to setup a dataframe for ggplot...
gdf_semi <- geomorph.data.frame(shape = gpa_FNP_semi$coords, treatment = classifiers$treatment, stage = classifiers$stage, cs = gpa_FNP_semi$Csize, pdist = PDist)
# gdf_semi_20 <- geomorph.data.frame(shape = gpa_FNP_semi_20$coords, treatment = classifiers_st20$treatment, stage = classifiers_st20$stage, cs = gpa_FNP_semi_20$Csize, pdist = PDist_20)

gdf_semi_df <- as.data.frame(cbind(gdf_semi$cs, as.character(gdf_semi$treatment), gdf_semi$pdist, gdf_semi$stage))
colnames(gdf_semi_df) <- c("csize", "treatment", "pdist", "stage")
gdf_semi_df$csize <- as.numeric(as.character(gdf_semi_df$csize))
gdf_semi_df$treatment <- as.factor(gdf_semi_df$treatment)
gdf_semi_df$pdist <- as.numeric(as.character(gdf_semi_df$pdist))
gdf_semi_df$stage <- as.factor(gdf_semi_df$stage)

# gdf_semi_20_df <- as.data.frame(cbind(gdf_semi_20$cs, as.character(gdf_semi_20$treatment), gdf_semi_20$pdist, gdf_semi_20$stage))
# colnames(gdf_semi_20_df) <- c("csize", "treatment", "pdist", "stage")
# gdf_semi_20_df$csize <- as.numeric(as.character(gdf_semi_20_df$csize))
# gdf_semi_20_df$treatment <- as.factor(gdf_semi_20_df$treatment)
# gdf_semi_20_df$pdist <- as.numeric(as.character(gdf_semi_20_df$pdist))
# gdf_semi_20_df$stage <- as.factor(gdf_semi_20_df$stage)

# now can plot csize and pdist
# do all samples
ggplot(gdf_semi_df, aes(x=stage, y=log(csize), fill=stage)) +
  geom_boxplot(alpha = 0.7) +
  # scale_fill_manual(values = palette()) +
  geom_jitter(width = 0.1, size = 2, aes(shape = treatment))

ggplot(gdf_semi_df, aes(x=treatment, y=log(csize), fill=treatment)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = palette()) +
  geom_jitter(width = 0.1, size = 2, aes(shape = stage))

ggplot(gdf_semi_df, aes(pdist, fill=treatment)) +
  scale_fill_manual(values = palette()) + geom_density(alpha = 0.65)

# no stage 18
ggplot(gdf_semi_20_df, aes(x=stage, y=log(csize), fill=stage)) +
  geom_boxplot(alpha = 0.7) +
  # scale_fill_manual(values = palette()) +
  geom_jitter(width = 0.1, size = 2, aes(shape = treatment))

ggplot(gdf_semi_20_df, aes(x=treatment, y=log(csize), fill=treatment)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = palette()) +
  geom_jitter(width = 0.1, size = 2, aes(shape = stage))

ggplot(gdf_semi_20_df, aes(pdist, fill=treatment)) +
  scale_fill_manual(values = palette()) + geom_density(alpha = 0.65)

# perform some allometry linear models...
allo_all <- procD.lm(shape ~ cs, data = gdf_semi, iter = 999, RRPP = TRUE)
summary(allo_all)

treatment_allo <- procD.lm(shape ~ cs*treatment, data = gdf_semi, iter = 999, RRPP = TRUE)
summary(treatment_allo)

stage_allo <- procD.lm(shape ~ cs*stage, data = gdf_semi, iter = 999, RRPP = TRUE)
summary(stage_allo) # no effect of centroid size on stage!

plotAllometry(allo_all, size = gdf_semi$cs, method="RegScore")
pcplot <- plotAllometry(allo_all, size = gdf_semi$cs, method="size.shape", col = gdf_semi$treatment, pch = as.integer(gdf_semi$stage))
ordiellipse(pcplot$size.shape.PCA, gdf_semi$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0, lty = 1, border = palette())
legend("topright", pch = 19, col = palette(), legend = levels(gdf_semi$treatment))
legend("bottomright", pch = as.integer(unique(gdf_semi$stage)), legend = as.character(unique(gdf_semi$stage)))


## LM for treatments
# all
### This generates the p-values used in the manuscript ###
treatment_lm <- procD.lm(shape ~ treatment, data = gdf_semi, RRPP = TRUE)
summary(treatment_lm)
summary(treatment_lm, test.type = "var")

# posthoc <- advanced.procD.lm(shape ~ treatment, f2=~1, groups=treatment, data = gdf_semi, RRPP = TRUE) DEPRECATED
treatment_lm_ph <- pairwise(treatment_lm, groups = gdf_semi$treatment)
summary(treatment_lm_ph)
summary(treatment_lm_ph, test.type = "var", confidence = 0.95, stat.table = TRUE)
### End manuscript stats


## let's see if those st18 ones are really a different stage or if they are just not well named images
treatment_stage <- procD.lm(shape ~ treatment*stage, data = gdf_semi, RRPP = TRUE)
summary(treatment_stage)
ph_groups <- interaction(gdf_semi$treatment, gdf_semi$stage)
treatment_stage_ph <- pairwise(treatment_stage, groups = ph_groups)
posthoc3 <- pairwise(treatment_stage, groups = gdf_semi$stage)

summary(posthoc3)
summary(posthoc3, test.type = "var", confidence = 0.95, stat.table = TRUE)
# so at a high level it looks like there's no differences w/ stage, but we should look by treatment
summary(treatment_stage_ph)
# no interaction in LY w/ stage
# no interaction in U0 w/ stage
# U73 has a strong effect, but it's because there are so few st 20
# no interaction in mix w/ stage
# I think the st 18,19,20 are all same stage and just named poorly. I don't think we should filter/remove 'st18' images


# st 20 only
# after analyzing stage I don't think we need to disregard st 18 and st19 named images
# gdf_semi_20 <- geomorph.data.frame(shape = gpa_FNP_semi_20$coords, treatment = classifiers_st20$treatment, stage = classifiers_st20$stage)
# treatment_20 <- procD.lm(shape ~ treatment, data = gdf_semi_20, RRPP = TRUE)
# summary(treatment_20)
# 
# posthoc20 <- pairwise(treatment_20, groups = gdf_semi_20$treatment)
# summary(posthoc20)
# summary(posthoc20, test.type = "var", confidence = 0.95, stat.table = TRUE)
# the results are similar as above, but almost all U73 samples were named st18


# Perform PCA
PCA_FNP_initial <- gm.prcomp(gpa_FNP_semi$coords)
plot(PCA_FNP_initial,col=classifiers$treatment)

# PCA_FNP_20_initial <- gm.prcomp(gpa_FNP_semi_20$coords)
# plot(PCA_FNP_20_initial,col=classifiers$treatment)


#### Plotting ####
# setup the color scheme for plots
pal <- wes_palette('Darjeeling2', n=5, type='discrete')
palette(rev(pal))

# doublemix still shows up as a 'level' of the factor even though we already removed these samples...
classifiers <- classifiers %>% mutate(treatment = recode(treatment, 'DoubleMix' = 'Mix'))
classifiers_filter <- classifiers$treatment

classifiers_st20 <- classifiers_st20 %>% mutate(treatment = recode(treatment, 'DoubleMix' = 'Mix'))
classifiers_filter_20 <- classifiers_st20$treatment

plot_df <- as.data.frame(PCA_FNP_initial$x)

# this one generates the histograms on the margins
png("./figs/PCA_head_shape_treatment_FNP_full_semi_margins_4.png", units = "in", width = 8.25, height = 6, res=300)
pdf("./figs/PCA_head_shape_treatment_FNP_full_semi_margins_4.pdf", width = 8.25, height = 6)
p <- ggplot(plot_df, aes(x = plot_df[,1], y = plot_df[,2], color = classifiers_filter)) + geom_point() + scale_color_manual(values = palette()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
px <- ggplot(plot_df, aes(x=plot_df[,1], color=classifiers_filter)) + geom_density() + scale_color_manual(values = palette()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
py <- ggplot(plot_df, aes(x=plot_df[,2], color=classifiers_filter)) + geom_density() + scale_color_manual(values = palette()) + coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p %>%
  insert_xaxis_grob(px, grid::unit(1, "in"), position = "top") %>%
  insert_yaxis_grob(py, grid::unit(1, "in"), position = "right") %>%
  ggdraw()

legend("topright", pch = 19, col = palette(), legend = levels(classifiers_filter))
title("PCA of shape coordinates")
dev.off()
dev.off()

# this one generates the 95% CI ellipses
png("./figs/PCA_head_shape_treatment_FNP_full_semi_4.png", units = "in", width = 8.25, height = 6, res=300)
pdf("./figs/PCA_head_shape_treatment_FNP_full_semi_4.pdf", width = 8.25, height = 6)
plot(PCA_FNP_initial, pch = 19, col = classifiers_filter, cex = .8, stroke=0)
ordiellipse(PCA_FNP_initial, classifiers_filter, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0, lty = 1, border = palette())
legend("topright", pch = 19, col = palette(), legend = levels(classifiers_filter))
title("PCA of shape coordinates")
dev.off()
dev.off()


#### wireframe plots ####
gdf_semi_20 <- geomorph.data.frame(shape = gpa_FNP_semi_20$coords, treatment = classifiers_st20$treatment, stage = classifiers_st20$stage)
control_20 <- which(gdf_semi_20$treatment == "control")
CTRL_coords_20 <- gdf_semi_20$shape[,,control_20]
CTRL_mean_shape_20 <- mshape(CTRL_coords_20)

control <- which(gdf_semi$treatment == "control")
CTRL_coords <- gdf_semi$shape[,,control]
CTRL_mean_shape <- mshape(CTRL_coords)

# there's something screwy with the way imageJ measures coordinates. To fix the outline you need to invert the y and add the max y
# so in this case make the y negative and add 2048
# now that we've added in scaling, the mean specimen coords are wrong, need to load in the 'old' unscaled coordinates
### Note: these plots are just flipped, not rotated. The left is left side of image, which means right side of face
landmarks_old <- readland.tps('morpho_files/morpho.txt', specID = "ID")

all_ctrl_shapes <- abind(CTRL_coords,CTRL_mean_shape,along=3)
all_ctrl_dist <- as.matrix(dist(two.d.array(all_ctrl_shapes)))
# control_st20_2.5_2 is closest to the 'mean shape' for control samples! Will use this one for the mesh outline

# read in original coordinate data for the specimen chosen above
mean_specimen <- which(dimnames(landmarks_old)[[3]] == 'control_st20_2.5_2_flat\t')
mean_specimen_coords <- landmarks_old[-c(15:19,22:24),,mean_specimen]

# loads in my hand-drawing .txt and applies a shape warp to it with the associated landmarks
ref <- warpRefOutline("st20_DMSO_25x_2_outline.txt", mean_specimen_coords, CTRL_mean_shape)

# outline <- as.matrix(read.table("st20_DMSO_25x_2_outline.txt", header = F))[,1:2]
# plot(ref2$outline, pch = 19, cex = 0.3, main = "Imported outline", asp = T, xlab = "x", ylab = "y")
# points(mean_specimen_coords[,1], mean_specimen_coords[,2], pch = 19, col = "red")


attributes(ref)
write.table(ref$outline, 'actual_mean_specimen_outline.txt', row.names = FALSE, col.names = FALSE)

# warp the mean shape per treatment across ALL PC
mix <- which(gdf_semi$treatment == "Mix")
mix_coords <- gdf_semi$shape[,,mix]
mix_mean_shape <- mshape(mix_coords)

U0 <- which(gdf_semi$treatment == "U0")
U0_coords <- gdf_semi$shape[,,U0]
U0_mean_shape <- mshape(U0_coords)

LY <- which(gdf_semi$treatment == "LY")
LY_coords <- gdf_semi$shape[,,LY]
LY_mean_shape <- mshape(LY_coords)

U73 <- which(gdf_semi$treatment == "U73")
U73_coords <- gdf_semi$shape[,,U73]
U73_mean_shape <- mshape(U73_coords)

png("./figs/mean_shape_mix.png", units = "in", width = 8, height = 8, res=300)
plotRefToTarget(CTRL_mean_shape, mix_mean_shape, outline=ref$outline, method="points", gridPars=gridPar(pt.size = .7, tar.pt.size = .7))
dev.off()

png("./figs/mean_shape_U0.png", units = "in", width = 8, height = 8, res=300)
plotRefToTarget(CTRL_mean_shape, U0_mean_shape, outline=ref$outline, method="points", gridPars=gridPar(pt.size = .7, tar.pt.size = .7))
dev.off()

png("./figs/mean_shape_LY.png", units = "in", width = 8, height = 8, res=300)
plotRefToTarget(CTRL_mean_shape, LY_mean_shape, outline=ref$outline, method="points", gridPars=gridPar(pt.size = .7, tar.pt.size = .7))
dev.off()

png("./figs/mean_shape_U73.png", units = "in", width = 8, height = 8, res=300)
plotRefToTarget(CTRL_mean_shape, U73_mean_shape, outline=ref$outline, method="points", gridPars=gridPar(pt.size = .7, tar.pt.size = .7))
dev.off()

# get the warp at certain values of PC1
PC = PCA_FNP_initial$x[,1]
# preds <- shape.predictor(gpa_FNP_semi$coords, x= PC, Intercept = FALSE,
#                          pred1 = min(PC)*.75, pred2 = max(PC)*.75)
preds <- shape.predictor(gpa_FNP_semi$coords, x= PC, Intercept = FALSE,
                         pred1 = -.05, pred2 = .05)

### Note: these plots are just flipped, not rotated. The left is left side of image, which means right side of face
# png("./figs/PC1_min2.png", units = "in", width = 8, height = 8)
pdf("./figs/PC1_min.pdf", width = 8, height = 8)
plotRefToTarget(CTRL_mean_shape, preds$pred1, outline=ref$outline, method="points", gridPars=gridPar(pt.size = 1, out.col='red', tar.pt.size = .7))
dev.off()

pdf("./figs/PC1_max.pdf", width = 8, height = 8)
plotRefToTarget(CTRL_mean_shape, preds$pred2, outline=ref$outline, method="points", gridPars=gridPar(pt.size = 0, out.col='white', tar.pt.size = .7))
dev.off()



PC = PCA_FNP_initial$x[,2]
pred2 <- shape.predictor(gpa_FNP_semi$coords, x= PC, Intercept = FALSE,
                         pred1 = -.05, pred2 = .05)
plotRefToTarget(CTRL_mean_shape, preds$pred1, outline=ref$outline, method="points")
plotRefToTarget(CTRL_mean_shape, pred2$pred2, outline=ref$outline, method="points")

# we are going to calculate the mean shape on the shape coordinates. Everything will be in Procrustes dist.

MUT_mean_shape <- mshape(MUT_coords)

MUT_coords <- gdf_head$coords[,,who_is_MUT]
dim(MUT_coords)


dim(CTRL_coords)

dim(gdf_head$coords)
plot(gpa_FNP_semi)



#### Asymmetry ####
# try GPA with bilateral object symmetry - see above for landmark pairs
gpa_symmetry <- bilat.symmetry(landmarks_FNP_only, side = NULL, replicate = NULL, object.sym = TRUE, curves=curveslide,
                          ind = dimnames(landmarks_FNP_only)[[3]], land.pairs = pairedLM, iter = 999, seed = NULL, RRPP = TRUE)


# gpa_symmetry_20 <- bilat.symmetry(landmarks_FNP_st20_only2, side = NULL, replicate = NULL, object.sym = TRUE, curves=curveslide2,
                          # ind = dimnames(landmarks_FNP_st20_only2)[[3]], land.pairs = pairedLM2, iter = 999, seed = NULL, RRPP = TRUE)

summary(gpa_symmetry)
plot(gpa_symmetry)

str(gpa_symmetry$FA.component)
str(gpa_symmetry$DA.component)

# check for outliers, they are removed above...
# outlier <- plotOutliers_percentile(A = gpa_symmetry$coords, percentile = 0.99, save.plot = FALSE)


## Analyze symmetric component
ANOVA_ALL_sym <- procD.lm(gpa_symmetry$symm.shape ~ classifiers$treatment,
                          iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_ALL_sym)
ANOVA_ALL_sym_ph <- pairwise(ANOVA_ALL_sym, groups = classifiers$treatment)
# this gives p values used in manuscript
summary(ANOVA_ALL_sym_ph)


# Perform PCA
PCA_SYM <- gm.prcomp(gpa_symmetry$symm.shape)


# PLOTTING
# setup the color scheme for plots
pal <- wes_palette('Darjeeling2', n=5, type='discrete')
palette(pal)

pdf("./figs/PCA_2D_symmetric_component_FGF.pdf", width = 7.5, height = 6)
plot(PCA_SYM, pch = 19, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_SYM, classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0, lty = 1, border = palette())
legend("bottomleft", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA SYMMETRIC COMPONENT")
dev.off()

PC = PCA_SYM$x[,1]
preds <- shape.predictor(gpa_symmetry$symm.shape, x= PC, Intercept = FALSE,
                         pred1 = -.05, pred2 = .05)
plotRefToTarget(CTRL_mean_shape, preds$pred1, outline=ref$outline, method="points", gridPars=gridPar(pt.size = 0, out.col='white', tar.pt.size = .7))
plotRefToTarget(CTRL_mean_shape, preds$pred2, outline=ref$outline, method="points", gridPars=gridPar(pt.size = 0, out.col='white', tar.pt.size = .7))



# temp_filter <- which(classifiers$treatment != 'LY')
# classifiers_filter <- classifiers$stage[-temp_filter]
# 
# pdf("./figs/PCA_head_shape_stage_raw_FNP_only.pdf", width = 8.25, height = 6)
# plot(PCA_FNP_LY_initial, pch = 19, col = classifiers_filter, cex = 1.25)
# # ordiellipse(PCA_FNP_LY_initial, classifiers_filter[1], kind="sd",conf=0.95, col = palette(),
# #             draw = "polygon", alpha = 0, lty = 1, border = palette())
# legend("topright", pch = 19, col = palette(), legend = levels(classifiers_filter))
# title("PCA of shape coordinates - CTRL vs treatment")
# dev.off()


## Analyze asymmetric component
ANOVA_ALL_asym <- procD.lm(gpa_symmetry$asymm.shape ~ classifiers$treatment,
                          iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_ALL_asym)
ANOVA_ALL_asym_ph <- pairwise(ANOVA_ALL_asym, groups = classifiers$treatment)
# this gives pvalues used in manuscript
summary(ANOVA_ALL_asym_ph)


# Perform PCA
PCA_ASYM <- gm.prcomp(gpa_symmetry$asymm.shape)


# PLOTTING
# setup the color scheme for plots
pal <- wes_palette('Darjeeling2', n=5, type='discrete')
palette(pal)

pdf("./figs/PCA_2D_asymmetric_component_FGF.pdf", width = 7.5, height = 6)
plot(PCA_ASYM, pch = 19, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_ASYM, classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0, lty = 1, border = palette())
legend("bottomleft", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA ASYMMETRIC COMPONENT")
dev.off()

PC = PCA_ASYM$x[,2]
preds <- shape.predictor(gpa_symmetry$asymm.shape, x= PC, Intercept = FALSE,
                         pred1 = 0, pred2 = .05)
plotRefToTarget(CTRL_mean_shape, preds$pred1, outline=ref$outline, method="points", gridPars=gridPar(pt.size = 0, out.col='white', tar.pt.size = .7))
plotRefToTarget(CTRL_mean_shape, preds$pred2, outline=ref$outline, method="points", gridPars=gridPar(pt.size = 0, out.col='white', tar.pt.size = .7))

PC = PCA_ASYM$x[,2]
preds <- shape.predictor(gpa_symmetry$asymm.shape, x= PC, Intercept = FALSE,
                         pred1 = min(PC), pred2 = max(PC))
plotRefToTarget(CTRL_mean_shape, preds$pred1, outline=ref2$outline, method="points")
plotRefToTarget(CTRL_mean_shape, preds$pred2, outline=ref2$outline, method="points")


## Analyze fluctuating asymmetric component
ANOVA_ALL_asym_fluct <- procD.lm(gpa_symmetry$FA.component ~ classifiers_st20$treatment,
                          iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_ALL_asym_fluct)
ANOVA_ALL_asym_fluct_ph <- pairwise(ANOVA_ALL_asym_fluct, groups = classifiers_st20$treatment)
summary(ANOVA_ALL_asym_fluct_ph)


# Perform PCA
PCA_ASYM_fluct <- gm.prcomp(gpa_symmetry$FA.component)


# PLOTTING
# setup the color scheme for plots
pal <- wes_palette('Darjeeling2', n=5, type='discrete')
palette(pal)

pdf("./figs/PCA_asymmetric_component_FGF.pdf", width = 7.5, height = 6)
plot(PCA_ASYM_fluct, pch = 19, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_ASYM_fluct, classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0, lty = 1, border = palette())
legend("bottomleft", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PCA ASYMMETRIC COMPONENT")
dev.off()




png("./figs/PCA_head_treatment_PC1-4.png", width = 750, height = 1200)
pdf("./figs/PCA_head_treatment_PC1-4.pdf", width = 8, height = 12)
par(mfrow=c(3,2))
plot(PCA_FNP_initial, pch = 19, axis1 = 1, axis2 = 2, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_FNP_initial, classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0, lty = 1, border = palette())
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PC1 ~ PC2")

plot(PCA_FNP_initial, pch = 19, axis1 = 1, axis2 = 3, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_FNP_initial$x[,c(1,3)], classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0, lty = 1, border = palette())
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PC1 ~ PC3")

plot(PCA_FNP_initial, pch = 19, axis1 = 1, axis2 = 4, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_FNP_initial$x[,c(1, 4)], classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0, lty = 1, border = palette())
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PC1 ~ PC4")

plot(PCA_FNP_initial, pch = 19, axis1 = 2, axis2 = 3, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_FNP_initial$x[,c(2,3)], classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0, lty = 1, border = palette())
legend("bottomleft", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PC2 ~ PC3")

plot(PCA_FNP_initial, pch = 19, axis1 = 2, axis2 = 4, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_FNP_initial$x[,c(2,4)], classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0, lty = 1, border = palette())
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PC2 ~ PC4")

plot(PCA_FNP_initial, pch = 19, axis1 = 3, axis2 = 4, col = classifiers$treatment, cex = 1.25)
ordiellipse(PCA_FNP_initial$x[,3:4], classifiers$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0, lty = 1, border = palette())
legend("topright", pch = 19, col = palette(), legend = levels(classifiers$treatment))
title("PC3 ~ PC4")

dev.off()
dev.off()

par(mfrow=c(1,1))

PCA_comp <- PCA_FNP_initial
class(PCA_comp) <- "princomp"

png("./figs/PCA_head_shape_scree_plot.png", width = 300, height = 300)
pdf("./figs/PCA_head_shape_scree_plot.pdf", height = 5, width = 5)
pca_scree <- fviz_eig(PCA_comp, addlabels=TRUE, hjust = -0.3,
                      barfill="darkgrey", barcolor ="black",
                      linecolor ="blue") + ylim(0, 85) +
  theme_classic()

print(pca_scree)
dev.off()

# I think unneeded but don't want to delete yet
wireframe_links <- define.links(CTRL_mean_shape, ptsize = 1, links = NULL)
# [,1] [,2]
# sel    1   14
# sel    1    2
# sel    2    3
# sel    3    4
# sel    4   18
# sel    5   18
# sel    5   19
# sel   16   19
# sel   16   20
# sel    6   20
# sel    6   21
# sel   17   21
# sel   17   22
# sel    7   22
# sel    7   23
# sel    8   23
# sel    8    9
# sel    9   10
# sel   10   11
# sel   11   15
# sel    1   13
# sel   11   12
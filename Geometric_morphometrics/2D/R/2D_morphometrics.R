# install.packages('devtools')
# library(devtools)
# install_github("marta-vidalgarcia/morpho.tools.GM", force = TRUE)
# install.packages(c('devtools', 'geomorph', 'tidyr', 'dplyr', 'abind', 'Evomorph', 'vegan', 'ggplot2', 'cowplot', 'factoextra'))

library(tidyr)
library(dplyr)
library(abind)

library(geomorph)
library(morpho.tools.GM)
library(Evomorph)

library(vegan)
library(ggplot2)
library(cowplot)
library(factoextra)


#### 0 Helper functions ####
#### 0.1 Landmark loading helpers ####
compile_landmark_files <- function() {
  # read in landmarks. I formatted them as a 'TPS' style
  # but these don't have the 'scale' in them so we will have to get them from the ID names and then re-save with scale
  landmarks_build <- readland.tps('morpho_files/U0_morpho.txt', specID = "ID")
  landmarks_build <- abind(landmarks_build, readland.tps('morpho_files/U73_morpho.txt', specID = "ID"))
  landmarks_build <- abind(landmarks_build, readland.tps('morpho_files/LY_morpho.txt', specID = "ID"))
  landmarks_build <- abind(landmarks_build, readland.tps('morpho_files/mix_morpho.txt', specID = "ID"))
  landmarks_build <- abind(landmarks_build, readland.tps('morpho_files/DMSO_morpho.txt', specID = "ID"))
  landmarks <- landmarks_build
  return(landmarks)
}


initialize_classifiers <- function(landmarks, names) {
  # create a list of classifiers. We will start with the 'id' from the TPS file and work from there...
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
  return(classifiers)
}


#### 0.2 GPA helpers ####
define_curveslide<- function(landmarks) {
  curveslide <- define.sliders(landmarks[,,1], nsliders=8)
  curveslide <- as.data.frame(rbind(c(2,3,4), c(9,11,12), c(6,18,7), c(7,19,8), c(5,20,6), c(6,21,18), c(18,22,7), c(7,23,19), c(19,24,8), c(8,25,9)))
  colnames(curveslide) <- c('before', 'slide', 'after')
  return(curveslide)
}


gpa_ggplot_df <- function(gpa_in, classifiers, pdist) {
  # need to setup a dataframe for ggplot...
  gm_df <- geomorph.data.frame(shape = gpa_in$coords, treatment = classifiers$treatment, stage = classifiers$stage, cs = gpa_in$Csize, pdist = pdist)

  gg_df <- as.data.frame(cbind(gm_df$cs, as.character(gm_df$treatment), gm_df$pdist, gm_df$stage))
  colnames(gg_df) <- c("csize", "treatment", "pdist", "stage")
  gg_df$csize <- as.numeric(as.character(gg_df$csize))
  gg_df$treatment <- as.factor(gg_df$treatment)
  gg_df$pdist <- as.numeric(as.character(gg_df$pdist))
  gg_df$stage <- as.factor(gg_df$stage)
  return(gg_df)
}


#### 0.3 Plot helpers ####
define_pal <- function() {
  return(c('lightgrey', '#d55c00', '#0072b2','#009e74', '#cc79a7'))
}


plot_symmetry_PC <- function(PC_data, classifiers, pal, tit) {
  plot(PC_data, pch = 19, col = classifiers$treatment, cex = 1.25)
  ordiellipse(PCA_SYM, classifiers$treatment, kind="sd",conf=0.95, col = pal,
              draw = "polygon", alpha = 0, lty = 1, border = pal)
  legend("bottomleft", pch = 19, col = pal, legend = levels(classifiers$treatment))
  title(title())
}

#### Main ####

#### 1 Load landmarks and setup dataframe ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd() #check our working directory
setwd("../../data/morphology/2D")
getwd() #check our working directory


# this checks if you've run the analysis before and skips over the pre-processing
if (file.exists('morpho_files/morpho_scaled.tps')){
  landmarks <- readland.tps('morpho_files/morpho_scaled.tps', specID = "ID")
  classifiers <- readRDS('classifiers.Rds')
} else {
  landmarks <- compile_landmark_files()
  names <- dimnames(landmarks)[[3]]
  names <- substr(names, 1, nchar(names)-1) # deletes the /t (tab)
  classifiers <- initialize_classifiers(landmarks, names)
  
  # write the tps file with adjusted scale in there and then reload data
  dimnames(landmarks)[[3]] <- names
  writeland.tps(landmarks, 'morpho_files/morpho_scaled.tps', scale = classifiers$magnification, specID=TRUE)
  saveRDS(classifiers, 'classifiers.Rds')
  landmarks_old <- landmarks
  classifiers_old <- classifiers
  
  landmarks <- readland.tps('morpho_files/morpho_scaled.tps', specID = "ID")
  classifiers <- readRDS('classifiers.Rds')
}


#### 2 Clean landmark data ####
#### 2.1 Landmarks definition ####
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


#### 2.2 Landmark cleanup ####
# delete x and y data of questionable landmarks
landmarks_clean <- landmarks[-c(15:19,22:24),,] # this one is used for analysis


#### 2.3 Outlier cleanup ####
# this is an iterative process, but there are some outliers I noted before GPA 
# was performed that can be removed from the outset

# in this case 'LY_st20_2_23' I had noted looked bad - remove
# 'U73_st20_2_1', 'U73_st20_2_4', 'LY_st20_2_25' are the wrong stage - remove
# 'U0_st20_3.2_27' and 'U0_st20_3.2_28' are very blury images - 
# but they aren't coming as outliers so I guess fine
names <- dimnames(landmarks)[[3]]
# create list of names to remove
outlier_remove <- which(names == 'LY_st20_2_23')
outlier_remove <- append(outlier_remove, which(names == 'U73_st20_2_1_flat'))
outlier_remove <- append(outlier_remove, which(names == 'U73_st20_2_4_flat'))
outlier_remove <- append(outlier_remove, which(names == 'LY_st20_2_25'))

# and remove them from the dataframe and classifiers list
landmarks_clean <- landmarks_clean[,,-outlier_remove]
classifiers <- classifiers[-outlier_remove,]


#### 2.3 Remove wrong stage ####
# remove stage HH18 samples
# based on talking to Diane, I don't think these are actually st18,
# the images just are not named well/correctly
# stage_20 <- which(classifiers$stage == 'st18')
# landmarks_FNP_st20_only <- landmarks_FNP_only[,,-stage_20]
# classifiers_st20 <- classifiers[-stage_20,]


#### 3 Generalized Procrustes analysis - overall shape ####
#### 3.1 Define and load curveslide for semilandmarks #### 
# this defines which landmarks are semilandmarks that 'slide' between other landmarks
# it's defined in a function in the helpers, but it's easy to just load them from csv
# after they've been defined once before
if (file.exists('curveslide.csv')) {
  curveslide <- read.csv('curveslide.csv')
} else {
  # this generates the curveslide csv, doesn't need to be re-run
  curveslide <-define_curveslide(landmarks_clean)
  write.csv(curveslide, 'curveslide.csv', row.names=FALSE, col.names = FALSE)
}

# setup the paired landmarks
# we'll just do a quick GPA so we can plot them all and get an idea of where the landmarks are
# also it allows at a glance to notice outliers
gpa_test <- gpagen(landmarks_clean)
plotAllSpecimens(gpa_test$coords, label=T, plot.param = list(pt.bg = "green", mean.cex=2,txt.pos=3, txt.cex=2))

# create a list of landmarks that are paired, these will be used for 
# analyzing symmetry later
pairedLM <- cbind(c(1:6,15,16,18,20:22), c(13:8,14,17,19,25:23))


#### 3.2 GPA w/ sliding semilandmarks ####

## perform GPA - one with all points, one with the points treated as semilandmarks
# will use Procrustes distance to determine semi placement
gpa_FNP_semi <- gpagen(landmarks_clean, curves = curveslide)

# check for outliers, many have already been removed above...
# if there are any more this is an easy way to find them, but I think I got them all
# LY_st20_2_27 and Mix_st20_2_45 are very malformed from treatment and should stay
outlier <- plotOutliers_percentile(A = gpa_FNP_semi$coords, percentile = 0.99, save.plot = FALSE)

## Generate a dataframe that can be used by ggplot2
PDist <- ShapeDist(gpa_FNP_semi$coords, gpa_FNP_semi$consensus)
gpa_FNP_semi_df <- gpa_ggplot_df(gpa_FNP_semi, classifiers, PDist)
gpa_FNP_semi_df$treatment <- factor(gpa_FNP_semi_df$treatment, levels = c('control', 'U0', 'LY', 'U73', 'Mix'))
gpa_FNP_semi_df <- gpa_FNP_semi_df %>% mutate(treatment = recode(treatment, 'Mix' = 'Triple', 'control' = 'DMSO'))
gpa_FNP_semi_df <- gpa_FNP_semi_df %>% mutate(stage = recode(stage, '1' = '18', '2' = '19', '3' = '20'))


#### 3.2 Centroid sizec Procruste's distance, allometry ####
# make a nice color palette 
# https://jfly.uni-koeln.de/color/
pal <- define_pal()

# now can plot centroid size and procrustes distance
ggplot(gpa_FNP_semi_df, aes(x=stage, y=log(csize), fill=stage)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2, aes(shape = treatment))

ggplot(gpa_FNP_semi_df, aes(x=treatment, y=log(csize), fill=treatment)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = pal) +
  geom_jitter(width = 0.1, size = 2, aes(shape = stage))

ggplot(gpa_FNP_semi_df, aes(pdist, fill=treatment)) +
  scale_fill_manual(values = pal) + geom_density(alpha = 0.65)

# perform some allometry linear models...
gdf_semi <- geomorph.data.frame(shape = gpa_FNP_semi$coords, treatment = classifiers$treatment, stage = classifiers$stage, cs = gpa_FNP_semi$Csize, pdist = PDist)
allo_all <- procD.lm(coords ~ Csize, data = gdf_semi, iter = 999, RRPP = TRUE)
summary(allo_all)

treatment_allo <- procD.lm(shape ~ cs*treatment, data = gdf_semi, iter = 999, RRPP = TRUE)
summary(treatment_allo)

stage_allo <- procD.lm(shape ~ cs*stage, data = gdf_semi, iter = 999, RRPP = TRUE)
summary(stage_allo)

plotAllometry(allo_all, size = gdf_semi$cs, method="RegScore")
pcplot <- plotAllometry(allo_all, size = gdf_semi$cs, method="size.shape", col = gdf_semi$treatment, pch = as.integer(gdf_semi$stage))
ordiellipse(pcplot$size.shape.PCA, gdf_semi$treatment, kind="sd",conf=0.95, col = palette(),
            draw = "polygon", alpha = 0, lty = 1, border = palette())
legend("topright", pch = 19, col = palette(), legend = levels(gdf_semi$treatment))
legend("bottomright", pch = as.integer(unique(gdf_semi$stage)), legend = as.character(unique(gdf_semi$stage)))


#### 3.3 ANOVA for mean shape GPA ####
# all
### This generates the p-values used in the manuscript ###
treatment_lm <- procD.lm(shape ~ treatment, data = gdf_semi, RRPP = TRUE)
summary(treatment_lm)
summary(treatment_lm, test.type = "var")

treatment_lm_ph <- pairwise(treatment_lm, groups = gdf_semi$treatment)
summary(treatment_lm_ph)
summary(treatment_lm_ph, test.type = "var", confidence = 0.95, stat.table = TRUE)
### End manuscript stats


#### 3.4 Treatment by stage ANOVA ####
## let's see if those st18 ones are really a different stage or if they are just not well named images
treatment_stage <- procD.lm(shape ~ treatment*stage, data = gdf_semi, RRPP = TRUE)
summary(treatment_stage)
ph_groups <- interaction(gdf_semi$treatment, gdf_semi$stage)
# full post-hoc treatment*stage. lots of interactions we don't care about...
# none of the 'main interactions' look significant
treatment_stage_ph <- pairwise(treatment_stage, groups = ph_groups)
summary(treatment_stage_ph)
# just compare stages to one another, no treatment
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


#### 4.0 Mean shape PCA ####
PCA_FNP_initial <- gm.prcomp(gpa_FNP_semi$coords)
plot(PCA_FNP_initial, col=classifiers$treatment)
# another way to check for some outliers
# mix_st20_3.2_7 is also just a strongly affected sample
text(PCA_FNP_initial$x, PCA_FNP_initial$y, dimnames(gpa_FNP_semi$coords)[[3]])

## scree plot
PCA_comp <- PCA_FNP_initial
class(PCA_comp) <- "princomp"

pdf("./figs/PCA_head_shape_scree_plot.pdf", height = 5, width = 5)
pca_scree <- fviz_eig(PCA_comp, addlabels=TRUE, hjust = -0.3,
                      barfill="darkgrey", barcolor ="black",
                      linecolor ="blue") + ylim(0, 85) +
  theme_classic()

print(pca_scree)
dev.off()


#### 4.1 PC Plotting ####
## these appear in the manuscript
# setup the color scheme for plots
pal <- define_pal()

# setup the df for better ordering
classifiers$treatment <- factor(classifiers$treatment, levels = c('control', 'U0', 'LY', 'U73', 'Mix'))
classifiers <- classifiers %>% mutate(treatment = recode(treatment, 'Mix' = 'Triple', 'control' = 'DMSO'))
classifiers_filter <- classifiers$treatment

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


#### 4.2 wireframe plots ####
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

# Double check that weird outlier one more time
weird_specimen <- which(dimnames(landmarks_old)[[3]] == 'Mix_st20_3.2_7\t')
weird_specimen_coords <- landmarks_old[-c(15:19,22:24),,weird_specimen]

# yea, it's not wrong, it's just got a bad truncation
ref <- warpRefOutline("st20_DMSO_25x_2_outline.txt", weird_specimen_coords, CTRL_mean_shape)


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



#### 5.0 Asymmetry ####
#### 5.1 Bilateral symmetry GPA ####
# try GPA with bilateral object symmetry - see above for landmark pairs
gpa_symmetry <- bilat.symmetry(landmarks_clean, side = NULL, replicate = NULL, object.sym = TRUE, curves=curveslide,
                          ind = dimnames(landmarks_clean)[[3]], land.pairs = pairedLM, iter = 999, seed = NULL, RRPP = TRUE)
summary(gpa_symmetry)
plot(gpa_symmetry)

str(gpa_symmetry$FA.component)
str(gpa_symmetry$DA.component)


#### 5.2 Bilateral symmetry symmetric component ####
ANOVA_ALL_sym <- procD.lm(gpa_symmetry$symm.shape ~ classifiers$treatment,
                          iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_ALL_sym)
ANOVA_ALL_sym_ph <- pairwise(ANOVA_ALL_sym, groups = classifiers$treatment)
# this gives p values used in manuscript
summary(ANOVA_ALL_sym_ph)


#### 5.3 Bilateral symmetry symmetric PCA and plots ####
PCA_SYM <- gm.prcomp(gpa_symmetry$symm.shape)

# PLOTTING
# setup the color scheme for plots
pal <- define_pal()
plot_symmetry_PC(PCA_SYM, classifiers, pal, "PCA SYMMETRIC COMPONENT")
pdf("./figs/PCA_2D_symmetric_component_FGF.pdf", width = 7.5, height = 6)

dev.off()

PC = PCA_SYM$x[,1]
preds <- shape.predictor(gpa_symmetry$symm.shape, x= PC, Intercept = FALSE,
                         pred1 = -.05, pred2 = .05)
plotRefToTarget(CTRL_mean_shape, preds$pred1, outline=ref$outline, method="points", gridPars=gridPar(pt.size = 0, out.col='white', tar.pt.size = .7))
plotRefToTarget(CTRL_mean_shape, preds$pred2, outline=ref$outline, method="points", gridPars=gridPar(pt.size = 0, out.col='white', tar.pt.size = .7))


####  5.4 Bilateral symmetry asymmetric component ####
ANOVA_ALL_asym <- procD.lm(gpa_symmetry$asymm.shape ~ classifiers$treatment,
                          iter=999, RRPP=TRUE, print.progress = FALSE)
summary(ANOVA_ALL_asym)
ANOVA_ALL_asym_ph <- pairwise(ANOVA_ALL_asym, groups = classifiers$treatment)
# this gives pvalues used in manuscript
summary(ANOVA_ALL_asym_ph)

####  5.5 Bilateral symmetry asymmetric PCA and plots ####
PCA_ASYM <- gm.prcomp(gpa_symmetry$asymm.shape)

# PLOTTING
pal <- define_pal()

pdf("./figs/PCA_2D_asymmetric_component_FGF.pdf", width = 7.5, height = 6)
plot_symmetry_PC(PCA_ASYM, classifiers, pal, "PCA ASYMMETRIC COMPONENT")
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
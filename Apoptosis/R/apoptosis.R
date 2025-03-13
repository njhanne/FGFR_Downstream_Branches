library(dplyr)
library(tidyr)
library(ggplot2)


# Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd() #check our working directory
setwd("./../../data/apoptosis")
getwd() #check our working directory


# Load data
df_cellpose <- read.csv('results.csv')

## Clean and combine data
# add percent to the cellpose df
df_cellpose <- df_cellpose %>% mutate(apoptosis = tunel_count / nuclei_count * 100)

df_cellpose$side <- factor(df_cellpose$side, labels = c("treated", "contralateral"))
df_cellpose$treatment <- factor(df_cellpose$treatment, levels = c("DMSO", "U0126", "LY294002", "U73122", 'mix'))

## Compare treated-control side and groups (anova)
# average the section results
df_cellpose <- df_cellpose %>% group_by(treatment, side, sample) %>% summarise_at(c('nuclei_count', 'tunel_count', 'apoptosis'), mean)

cellpose_summarise <- summarise(group_by(df_cellpose, treatment,side), mean=mean(apoptosis),sd=sd(apoptosis))
cellpose_summarise <- cellpose_summarise %>% unite('treatment_side', c('treatment', 'side'), remove=FALSE)


# Plot
pdf("./figs/cellpose_tunel.pdf", width = 7.5, height = 6)
p <- ggplot(data=df_cellpose, aes(fill = side, y=apoptosis, x = treatment)) +
     geom_bar(stat = "summary", position = "dodge", fun = mean) +
     geom_point(position=position_jitterdodge(dodge.width=0.9)) +
     geom_errorbar(data = cellpose_summarise, aes(y=mean, x=treatment, ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))
# file_name <- paste("cellpose_pHH3.png")
# ggsave(filename=file_name, p, width = 15, heigh = 25, units='cm')
print(p)
dev.off()


# Analysis
## Paired T-test
df_cellpose_pair <- df_cellpose %>% filter(! treatment %in% c('DMSO', 'mix'))
df_cellpose_pair <- df_cellpose_pair %>% filter(! sample %in% c('U0_1', 'U73_9'))
by(df_cellpose_pair, df_cellpose_pair$treatment, function(x) t.test(x$apoptosis ~ x$side, paired=TRUE, data=x))
# these are nearly significant so I don't think we should combine the sides in the plots

## ANOVA
### main effects only
test <- aov(apoptosis ~ treatment, data= df_cellpose)
summary(test)
TukeyHSD(test)

### treatment by side
test <- aov(apoptosis ~ treatment*side, data= df_cellpose)
summary(test)
TukeyHSD(test)
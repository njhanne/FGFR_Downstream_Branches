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

### Clean and combine data
# add percent to the cellpose df
df_cellpose <- df_cellpose %>% mutate(apoptosis = tunel_count / nuclei_count * 100)

df_cellpose$side <- factor(df_cellpose$side, labels = c("treated", "contralateral"))
df_cellpose$treatment <- factor(df_cellpose$treatment, levels = c("DMSO", "U0126", "LY294002", "U73122", 'mix'))

### Compare treated-control side
# Should use Student's paired t-test since we want difference within each sample
# First we need to average the section results or it won't work...
df_cellpose <- df_cellpose %>% group_by(treatment, side, sample) %>% summarise_at(c('nuclei_count', 'tunel_count', 'apoptosis'), mean)

cellpose_summarise <- summarise(group_by(df_cellpose, treatment,side), mean=mean(apoptosis),sd=sd(apoptosis))
cellpose_summarise <- cellpose_summarise %>% unite('treatment_side', c('treatment', 'side'), remove=FALSE)

pdf("./figs/cellpose_pHH3.pdf", width = 7.5, height = 6)
p <- ggplot() + geom_bar(data = cellpose_summarise, aes(y=mean, x = side), stat="identity") +
  geom_jitter(data = df_cellpose, aes(x = side, y = apoptosis), shape=16) +
  geom_errorbar(data = cellpose_summarise, aes(y=mean,x=side,ymin=mean-sd,ymax=mean+sd)) +
  facet_wrap(~ treatment, nrow=1)
# file_name <- paste("cellpose_pHH3.png")
# ggsave(filename=file_name, p, width = 15, heigh = 25, units='cm')
print(p)
dev.off()


# ANOVA
### main effects only
test <- aov(apoptosis ~ treatment, data= df_cellpose)
summary(test)
TukeyHSD(test)

### side
test <- aov(apoptosis ~ treatment*side, data= df_cellpose)
summary(test)
TukeyHSD(test)
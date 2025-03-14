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

df_apoptosis <- df_cellpose
df_apoptosis$side <- factor(df_apoptosis$side, labels = c("contralateral", 'treated'))
df_apoptosis$treatment <- factor(df_apoptosis$treatment, levels = c("DMSO", "LY294002", 'mix', "U0126", "U73122"))

### Compare treated-control side
# Should use Student's paired t-test since we want difference within each sample
# First we need to average the section results or it won't work...
df_apoptosis <- df_apoptosis %>% group_by(treatment, side, sample) %>% summarise_at(c('nuclei_count', 'tunel_count', 'apoptosis'), mean)

apoptosis_summarise <- summarise(group_by(df_apoptosis, treatment,side), mean=mean(apoptosis),sd=sd(apoptosis))
apoptosis_summarise <- apoptosis_summarise %>% unite('treatment_side', c('treatment', 'side'), remove=FALSE)

pdf("./figs/cellpose_tunel.pdf", width = 7.5, height = 6)
p <- ggplot(data=df_apoptosis, aes(fill = side, y=apoptosis, x = treatment)) +
     geom_bar(stat = "summary", position = "dodge", fun = mean) +
     geom_point(position=position_jitterdodge(dodge.width=0.9)) +
     geom_errorbar(data = apoptosis_summarise, aes(y=mean, x=treatment, ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))
print(p)
dev.off()


# ANOVA
### main effects only
test <- aov(apoptosis ~ treatment, data= df_apoptosis)
summary(test)
TukeyHSD(test)

### side
test <- aov(apoptosis ~ treatment*side, data= df_apoptosis)
summary(test)
TukeyHSD(test)


## Paired t-tests
# Then we gotta get rid of U73_21 since it only has one side and can't be used for the t-test
df_apoptosis_t <- df_apoptosis %>% filter(! sample %in% c('U0_1', 'U73_9'))

# run the t-tests - this gives pvalues used in manuscript
df_apoptosis_t <- df_apoptosis_t %>% filter(! treatment %in% c('DMSO', 'mix'))
by(df_apoptosis_t, df_apoptosis_t$treatment, function(x) t.test(x$apoptosis ~ x$side, paired=TRUE, data=x))

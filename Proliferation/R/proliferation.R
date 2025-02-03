library(dplyr)
library(tidyr)
library(ggplot2)

# Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd() #check our working directory
setwd("./../../data/proliferation")
getwd() #check our working directory

# Load data
df_cellpose <- read.csv('results.csv')

### Clean and combine data
# add percent to the cellpose df
df_cellpose <- df_cellpose %>% mutate(proliferation = phh3_count / nuclei_count * 100)

df_cellpose$side <- factor(df_cellpose$side, levels = c("treated", "control"))
df_cellpose$side <- factor(df_cellpose$side, labels = c("treated", "contralateral"))
df_cellpose <- rename(df_cellpose, treatment = group)
df_cellpose$treatment <- factor(df_cellpose$treatment, levels = c("DMSO", "U0126", "LY294002", "U73122"))

### Compare treated-control side
# Should use Student's paired t-test since we want difference within each sample
# First we need to average the section results or it won't work...
df_cellpose <- df_cellpose %>% group_by(treatment, side, sample) %>% summarise_at(c('nuclei_count', 'phh3_count', 'proliferation'), mean)

# Then we gotta get rid of U73_21 since it only has one side and can't be used for the t-test
df_cellpose <- df_cellpose %>% filter(sample != 'U73_21')

# run the t-tests - this gives pvalues used in manuscript
df_cellpose_no_DMSO <- df_cellpose %>% filter(treatment != 'DMSO')
by(df_cellpose_no_DMSO, df_cellpose_no_DMSO$treatment, function(x) t.test(x$proliferation ~ x$side, paired=TRUE, data=x))

cellpose_summarise <- summarise(group_by(df_cellpose, treatment,side), mean=mean(proliferation),sd=sd(proliferation))
cellpose_summarise <- cellpose_summarise %>% unite('treatment_side', c('treatment', 'side'), remove=FALSE)

pdf("./figs/cellpose_pHH3.pdf", width = 7.5, height = 6)
p <- ggplot() + geom_bar(data = cellpose_summarise, aes(y=mean, x = side), stat="identity") +
  geom_jitter(data = df_cellpose, aes(x = side, y = proliferation), shape=16) +
  geom_errorbar(data = cellpose_summarise, aes(y=mean,x=side,ymin=mean-sd,ymax=mean+sd)) +
  facet_wrap(~ treatment, nrow=1)
# file_name <- paste("cellpose_pHH3.png")
# ggsave(filename=file_name, p, width = 15, heigh = 25, units='cm')
print(p)
dev.off()


# ANOVA
### main effects only
test <- aov(proliferation ~ treatment, data= df_cellpose)
summary(test)
TukeyHSD(test)

### side
test <- aov(proliferation ~ treatment*side, data= df_cellpose)
summary(test)
TukeyHSD(test)


### Now that we have DMSO I don't think we need these comparisons anymore
### Looks like the U0126 don't have significant difference but the control side seems 'affected'
# Will use Welch's t test as the sample sizes are different
df_cellpose_temp <- df_cellpose %>% mutate(welch_sort = case_when(treatment == 'U0126' ~ 'U0126',
                                              treatment != 'U0126' ~ 'notU0126'))
df_cellpose_temp <- df_cellpose_temp %>% filter(side == 'contralateral')
# df_cellpose_temp <- df_cellpose_temp %>% filter(treatment != 'U73122')

# test normality
# dplyr note: 'pull' gets the values while 'select' gets it as a dataframe
shapiro.test(df_cellpose_temp %>% filter(welch_sort == 'U0126') %>% pull(proliferation))
shapiro.test(df_cellpose_temp %>% filter(welch_sort == 'notU0126') %>% pull(proliferation))

# Welch test
t.test(proliferation ~ welch_sort, data=df_cellpose_temp)

ggplot(df_cellpose_temp, aes(x = welch_sort, y = proliferation)) +
  geom_boxplot(aes(fill = welch_sort), alpha = .2) +
  geom_point(size = 2)

# Lets just do absolute mean of U73 vs contralateral side for LY and contralateral U0
shapiro.test(df_cellpose %>% filter(treatment == 'U73122') %>% pull(proliferation))
shapiro.test(df_cellpose %>% filter(treatment == 'U0126', side == 'contralateral') %>% pull(proliferation))
shapiro.test(df_cellpose %>% filter(treatment == 'LY294002', side == 'contralateral') %>% pull(proliferation))

temp_sort <- df_cellpose %>% filter(treatment == 'U73122')
temp_sort <- bind_rows(temp_sort, df_cellpose %>% filter(treatment == 'U0126', side == 'contralateral'))
t.test(proliferation ~ treatment, data=temp_sort)

temp_sort <- df_cellpose %>% filter(treatment == 'U73122')
temp_sort <- bind_rows(temp_sort, df_cellpose %>% filter(treatment == 'LY294002', side == 'contralateral'))
t.test(proliferation ~ treatment, data=temp_sort)
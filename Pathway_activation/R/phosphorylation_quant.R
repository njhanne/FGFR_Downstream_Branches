library(dplyr)
library(tidyr)
library(ggplot2)

# Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd() #check our working directory
setwd("./../../data/Branches_activation")
getwd() #check our working directory

# Load data
df_activation <- read.csv('samples.csv')


### Clean and combine data
# remove unwanted cols
df_clean <- df_activation %>% group_by(AB, treatment, time, new_id, side) %>%
  select(area, mean, thresh_area, thresh_mean) %>% drop_na()

# make variables into factors
df_clean$side <- factor(df_clean$side, levels = c("treat", "contra"))
df_clean$AB <- factor(df_clean$AB, levels = c("AKT", "plcy", 'ERK'))
df_clean$treatment <- factor(df_clean$treatment, levels = c("DMSO", "U0126", "LY294002", "U73122", 'mix'))
df_clean$time <- factor(df_clean$time, levels = c("24hr", "6hr"))

# combine samples with more than one pair of images
df_clean <- df_clean %>% group_by(AB, treatment, time, new_id, side) %>%summarise_at(c('area', 'mean', 'thresh_area', 'thresh_mean'), mean)

# calculate total area divided by threshold area
df_clean <- df_clean %>% mutate(rel_area = thresh_area / area)

plot_df <- df_clean %>% group_by(AB, treatment, time, new_id) %>% summarise(rel_area = rel_area[side == "treat"] / rel_area[side == "contra"])
plot_df_summary <- summarise(group_by(plot_df, treatment, AB, time), mean=mean(rel_area),sd=sd(rel_area))

plot_df_summary <- summarise(group_by(df_clean, treatment, AB, time, side), mean=mean(rel_area),sd=sd(rel_area))
ggplot(data = plot_df_summary %>% filter(AB == 'AKT', time=='24hr'), aes(fill = side, y=mean, x = treatment), stat='identity') +
  geom_bar(position="dodge", stat="identity") + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9))



p <- ggplot() + geom_bar(data = plot_df_summary %>% filter(AB == 'AKT', time=='24hr'), aes(y=mean, x = treatment), stat='identity') +
  geom_jitter(data = plot_df %>% filter(AB == 'AKT', time=='24hr'), aes(x = treatment, y = rel_area), shape=16) +
  geom_errorbar(data = plot_df_summary %>% filter(AB == 'AKT', time=='24hr'), aes(y=mean,x=treatment,ymin=mean-sd,ymax=mean+sd))


test <- aov(rel_area ~ treatment, data= plot_df %>% filter(AB == 'AKT', time=='24hr'))
summary(test)
TukeyHSD(test)
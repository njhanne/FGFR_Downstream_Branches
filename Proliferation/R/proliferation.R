library(dplyr)
library(tidyr)
library(ggplot2)

# Directory
getwd() #check our working directory
setwd("./data/proliferation")

# Load data
df_handcount <- read.csv("pHH3_handcount_combined.csv")
df_cellpose <- read.csv('results.csv')

### Clean and combine data
# add percent to the cellpose df
df_cellpose <- df_cellpose %>% mutate(proliferation = stain_count / nuclei_count * 100)

# combine handcount and cellpose data into one df
df_cellpose <- df_cellpose %>% mutate(measurer = 'cellpose')
df_cellpose <- rename(df_cellpose, treatment = group)
df <- df_cellpose %>% select(sample, side, treatment, proliferation, measurer, nuclei_count, stain_count)

df_temp <- df_handcount %>% select(!side)
df_temp <- df_temp %>% pivot_longer(cols= (starts_with('control') | starts_with('treated')), names_to = c('side', '.value'), names_sep = '_')
df_temp <- df_temp %>% mutate(measurer = 'handcount')
df_temp <- rename(df_temp, proliferation = ratio, nuclei_count = total, stain_count = phh3)
df <- bind_rows(df, df_temp)


### Compare handcount vs cellpose
# Should use Student's paired t-test since we want difference within each sample
df$measurer <- as.factor(df$measurer)
df <- df %>% unite(sample_side, c('sample', 'side'), remove=FALSE)
t_test_temp = group_by(df, sample_side) %>% filter(n() != 1) %>% arrange(measurer, sample_side) %>% ungroup()

with(t_test_temp, t.test(proliferation ~ measurer, paired=TRUE))

difference_nuc <- t_test_temp %>% group_by(sample_side) %>%
  mutate(diff = -diff(nuclei_count))
difference_phh3 <- t_test_temp %>% group_by(sample_side) %>%
  mutate(diff = -diff(stain_count))

# they are different... we should graph them to compare

p <- ggplot(t_test_temp, aes(x = measurer, y = proliferation)) +
    geom_boxplot(aes(fill = measurer), alpha = .2) +
    geom_line(aes(group = sample_side)) +
    geom_point(size = 2) +
    facet_wrap(~ treatment)
file_name <- paste("hand-vs-cellpose.png")
ggsave(filename=file_name, p, width = 15, heigh = 25, units='cm')
print(p)

# to me the cellpose seems consistently a similar amount higher, but there doesn't appear much interaction
# For now let's just use both and see if they give similar results
df <- group_by(df, sample_side) %>% arrange(sample_side) %>% ungroup()
df$treatment <- as.factor(df$treatment)
df_cellpose <- t_test_temp %>% filter(measurer == 'cellpose')
df_cellpose <- df %>% filter(measurer == 'cellpose') %>% group_by(sample) %>% filter(n() != 1)
df_handcount <- df %>% filter(measurer == 'handcount')

### Compare treated-control side
# Should use Student's paired t-test since we want difference within each sample
by(df_handcount, df_handcount$treatment, function(x) t.test(x$proliferation ~ x$side, paired=TRUE, data=x))

df_handcount$side <- factor(df_handcount$side, levels = c("treated", "control"))
df_handcount$treatment <- factor(df_handcount$treatment, levels = c("U0126", "U73122", "LY294002"))

df_cellpose$side <- factor(df_cellpose$side, levels = c("treated", "control"))
df_cellpose$treatment <- factor(df_cellpose$treatment, levels = c("U0126", "U73122", "LY294002"))



handcount_summarise <- summarise(group_by(df_handcount, treatment,side), mean=mean(proliferation),sd=sd(proliferation))
p <- ggplot() + geom_bar(data = handcount_summarise, aes(y=mean, x = side), stat="identity") +
  geom_jitter(data = df_handcount, aes(x = side, y = proliferation)) +
  geom_errorbar(data = handcount_summarise, aes(y=mean,x=side,ymin=mean-sd,ymax=mean+sd)) +
  facet_wrap(~ treatment)
file_name <- paste("handcount_pHH3.png")
ggsave(filename=file_name, p, width = 15, heigh = 25, units='cm')
print(p)


by(df_cellpose, df_cellpose$treatment, function(x) t.test(x$proliferation ~ x$side, paired=TRUE, data=x))

df_cellpose$side <- factor(df_cellpose$side, levels = c("treated", "control"))
df_cellpose$treatment <- factor(df_cellpose$treatment, levels = c("U0126", "U73122", "LY294002"))

cellpose_summarise <- summarise(group_by(df_cellpose, treatment,side), mean=mean(proliferation),sd=sd(proliferation))
p <- ggplot() + geom_bar(data = cellpose_summarise, aes(y=mean, x = side), stat="identity") +
  geom_jitter(data = df_cellpose, aes(x = side, y = proliferation)) +
  geom_errorbar(data = cellpose_summarise, aes(y=mean,x=side,ymin=mean-sd,ymax=mean+sd)) +
  facet_wrap(~ treatment)
file_name <- paste("cellpose_pHH3.png")
ggsave(filename=file_name, p, width = 15, heigh = 25, units='cm')
print(p)

### Looks like the U0126 don't have significant difference but the control side seems 'affected'
# Will use Welch's t test as the sample sizes are different
df_cellpose_temp <- df_cellpose %>% mutate(welch_sort = case_when(treatment == 'U0126' ~ 'U0126',
                                              treatment != 'U0126' ~ 'notU0126'))
df_cellpose_temp <- df_cellpose_temp %>% filter(side == 'control')
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

# Lets just do absolute mean of U73 vs control side for LY and control U0
shapiro.test(df_cellpose %>% filter(treatment == 'U73122') %>% pull(proliferation))
shapiro.test(df_cellpose %>% filter(treatment == 'U0126', side == 'control') %>% pull(proliferation))
shapiro.test(df_cellpose %>% filter(treatment == 'LY294002', side == 'control') %>% pull(proliferation))

temp_sort <- df_cellpose %>% filter(treatment == 'U73122')
temp_sort <- bind_rows(temp_sort, df_cellpose %>% filter(treatment == 'U0126', side == 'control'))
t.test(proliferation ~ treatment, data=temp_sort)

temp_sort <- df_cellpose %>% filter(treatment == 'U73122')
temp_sort <- bind_rows(temp_sort, df_cellpose %>% filter(treatment == 'LY294002', side == 'control'))
t.test(proliferation ~ treatment, data=temp_sort)
library(dplyr)
library(tidyr)
library(ggplot2)

# Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd() #check our working directory
setwd("../../../data/Branches_activation")
getwd() #check our working directory

# Load data
df_activation <- read.csv('fluor_results.csv')


### Clean and combine data
# remove unwanted cols
df_clean <- df_activation %>% group_by(pathway, treatment, sample_id, section, sample_name) %>%
  select(total_cells, pos_cells, total_masked_cells, pos_masked_cells) %>% drop_na()

# make variables into factors
df_clean$pathway <- factor(df_clean$pathway, levels = c("P-AKT", 'P-ERK', 'P-PLCr'))
df_clean$treatment <- factor(df_clean$treatment, levels = c("DMSO", 'mix'))

# combine samples with more than one pair of images
df_clean <- df_clean %>% group_by(pathway, treatment, sample_id) %>% summarise_at(c('total_cells', 'pos_cells', 'total_masked_cells', 'pos_masked_cells'), mean)

# calculate total area divided by threshold area
df_clean <- df_clean %>% mutate(rel_positive = pos_cells / total_cells)
df_clean <- df_clean %>% mutate(rel_positive_masked = pos_masked_cells / total_masked_cells)


plot_df_summary <- summarise(group_by(df_clean, treatment, pathway), mean=mean(rel_positive_masked),sd=sd(rel_positive_masked), n=n())

pdf("./analysis/figs/masked_activated.pdf", width = 7.5, height = 6)
p <- ggplot(data=df_clean, aes(fill = treatment, y=rel_positive_masked, x = pathway)) +
     geom_bar(stat = "summary", position = "dodge", fun = mean) +
     geom_point(position=position_jitterdodge(dodge.width=0.9)) +
     geom_errorbar(data = plot_df_summary, aes(y=mean, x=pathway, ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))
print(p)
dev.off()


for (pathway_i in 1:length(levels(df_clean$pathway))) {
  temp_df <- df_clean %>% filter(pathway == levels(df_clean$pathway)[pathway_i])
  print(levels(df_clean$pathway)[pathway_i])
  print(summary(aov(rel_positive_masked ~ treatment, data=temp_df)))
}

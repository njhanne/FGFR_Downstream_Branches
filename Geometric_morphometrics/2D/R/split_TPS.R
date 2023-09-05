library(geomorph)
library(dplyr)
setwd("~/Box/Desktop/UCSF/Branches/Morphology")

df<- readland.tps('morpho_new_tilt.txt', specID="ID") #load tps file
df_simple <- two.d.array(df) # convert from geomorph 3D to normal dataframe


# make code to break up new groups
df_treated <- df_simple[,c(1:14,33:46,51:52,55:60)]
df_untreated <- df_simple[,c(25:26,23:24,21:22,19:20,17:18,15:16,13:14,
                             33:34,31:32,29:30,27:28,
                             49:50,47:48,45:46,
                             53:54,
                             65:66,63:64,61:62)]

row_names <- gsub('.{1}$', '', row.names(df_simple)) # regex deletes /t at end of string
row.names(df_treated) <- paste(row_names, 'treated', sep='_')
row.names(df_untreated) <- paste(row_names, 'untreat', sep='_')


df_sides <- rbind(df_treated, df_untreated)
df_sides <- arrayspecs(df_sides, 18, 2) # convert back to geomorph 3D
writeland.tps(df_sides, 'morpho_new_tilt_sides', specID = TRUE)

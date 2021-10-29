##### Determines proportion of met and abx genes per total genes and
##### and displays in a jitter plot with box plot overlay

rm(list = ls())

# example installation
#install.packages("ggplot2")

library(readxl)
library(ggplot2) # ggplot() for plotting
library(dplyr) # data reformatting
library(tidyr) # data reformatting
library(stringr) # string manipulation
library(cowplot)
library(BoutrosLab.plotting.general)
library(tidyverse) 

# set working directory
wd_pth <- "~/Desktop/scripts_and_datasets/"
setwd(wd_pth)
data_pth <- "~/Desktop/scripts_and_datasets/" 
data_pth <- paste(data_pth,"TableS11.xlsx",sep="")
info <- read_excel(data_pth)

# save figure as PDF
figure_pth = paste(wd_pth,"fig2D.pdf",sep="")
pdf(file = figure_pth,
    width = 5, # The width of the plot in inches
    height = 4)

# plotting
gscat <- ggplot(info, aes(x=ST_PREV, y=genes_per_total, fill = GENE_TYPE)) + 
  # scatter plot with jitter
  geom_point(mapping = aes(color = GENE_TYPE), 
             position = position_jitterdodge(jitter.height=0.1, 
                                             jitter.width=0.1), dodge.width = 0.1) +
  # set axis properties
  theme_bw(base_line_size = 2, base_rect_size = 2) +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  # box plot
  geom_boxplot(outlier.size=0) +
  scale_y_continuous(trans='log10') +
  # change the fill (for boxes) and color (for scatter points)
  scale_fill_manual(values = c("#01B2FF","#FF193F")) +
  scale_color_manual(values = c("#95DAFF","#FFADB3")) +
  # x and y labels
  xlab("ST") + ylab("# genes per total genes") +
  # change legend label
  labs(fill = "Gene type") +
  guides(color = F)


# plot figure
plot(gscat)

# generate figure as file
dev.off()
plot(gscat)

# run ANOVA
res.aov <- aov(genes_per_total ~ ST_PREV + GENE_TYPE + ST_PREV:GENE_TYPE, data = info)
summary(res.aov)


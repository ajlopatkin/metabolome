##### Determines percentage of plasmids per mobility type (top graph)
##### and number of plasmids (bottom graph). Both per prevalent ST.

# example installation
#install.packages("ggtree")


rm(list = ls())

install.packages("ggtree")
library(ggtree)
library(readxl)
library(ape)
library(ggplot2) # ggplot() for plotting
library(dplyr) # data reformatting
library(tidyr) # data reformatting
library(stringr) # string manipulation
library(cowplot)
library(gcookbook) # Load gcookbook for the cabbage_exp data set

set.seed(2017-02-16)


# set working directory

wd_pth <- "~/Desktop/scripts_and_datasets/"
setwd(wd_pth)
data_pth <- "~/Desktop/scripts_and_datasets/" 
data_pth <-paste(data_pth,"TableS3.xlsx",sep="") 
info <- read_excel(data_pth)



# save figure file
figure_pth = paste(wd_pth,"fig1B.pdf",sep="")
pdf(file = figure_pth,
    width = 4, # The width of the plot in inches
    height = 6)

# select columns from the info dataset to make a submatrix
mobility <- select(info,PredictedMobility)
ST <- select(info,ST_PREV)
strain <- select(info,BioSample)
data <- data.frame(mobility,ST,strain)

# data for the bar graph
data_bar <- data %>% count(PredictedMobility,ST_PREV)

# start defining the g_bar
gbar <- ggplot(data_bar, aes(x = ST_PREV, y = n, fill = PredictedMobility)) +
  geom_bar(position="fill", stat="identity", width=c(0.7,0.7,0.7,0.7,0.7)) + 
  theme_bw(base_line_size = .5, base_rect_size = 1) +
  scale_fill_discrete(name = element_blank()) +
  scale_fill_grey() + 
  theme(legend.position="bottom", 
        legend.box = "horizontal",
        legend.box.background = element_rect(fill="NA", colour="black")) +
  labs(fill = "")
# save gbar as a figure
gbar_legend <- get_legend(gbar)

# now remove the legend to save the bar graph alone as a figure
gbar <- gbar +   theme(axis.ticks.x = element_blank(), 
                       axis.text.x = element_blank(),
                       legend.position = "None",
                       axis.text = element_text(size=18),
                       axis.title = element_text(size=14),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_cartesian(ylim=c(0, 1)) +
  xlab(element_blank()) +
  ylab("Percent plasmids") +
  scale_y_continuous(labels=c("0.00" = "0", "0.25" = "",
                            "0.50" = "50", "0.75" = "", "1.00" = "100"),expand = c(0,0))

# data for the scatter
data_scatter <- data %>% count(ST_PREV,BioSample)
my_colors <- RColorBrewer::brewer.pal(6, "Pastel2")
my_colors <- c(my_colors[1:4],my_colors[6])

# plot data for scatter
gscat <- ggplot(data_scatter, aes(x=ST_PREV, y=n, fill = ST_PREV)) + 
  # set axis properties
  theme_bw(base_line_size = .5, base_rect_size = 1) +
  theme(legend.position = "None",
        axis.text = element_text(size=18),
        axis.title = element_text(size=18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # scatter plot with jitter
  geom_point(mapping = aes(color = ST_PREV), 
             position = position_jitterdodge(jitter.height=0.1, 
                                             jitter.width=0.1), 
             alpha = 0.5) +
  # boxplot
  geom_boxplot(outlier.size=0, width=0.5/length(unique(info$ST_PREV))) +
  coord_cartesian(ylim=c(0, 12)) +
  scale_fill_manual(values= my_colors) +
  scale_color_manual(values= my_colors) +
  scale_y_continuous(breaks=seq(0, 12, 3),expand = c(0,0)) +
  xlab("ST") + ylab("Number of plasmids")

# combine together
pg1 <- plot_grid(gbar, NULL, gscat,labels = NA, ncol = 1,align = "v", rel_heights = c(1,0.1, 2))
pg2 <- plot_grid(gbar_legend,pg1,ncol = 1, rel_heights = c(1, 10))
plot(pg2)


# generate figure as file
dev.off()


plot(pg2)

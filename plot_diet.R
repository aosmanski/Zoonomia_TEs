library(magrittr)
library(dplyr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(ggdist)
library(tidybayes)
library(ggplot2)
library(cowplot)
library(rstan)
library(brms)
library(ggrepel)
library(RColorBrewer)
library(gganimate)
library(posterior)

##remove prior data
rm(list=ls()) 

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##get data
load("TE_diet.RData")

##much better plotting of bayes results
##https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html
##set theme
theme_set(theme_tidybayes() + panel_border())

##plot differences for big data set
mamm<-m.die %>%
  spread_draws(r_Diet[condition,]) %>%
  compare_levels(r_Diet, by = condition, comparison = list(c("carnivore", "herbivore"), c("carnivore", "omnivore"), c("herbivore", "omnivore"))) %>%
  ungroup() %>%
#  mutate(condition = reorder(condition, r_Diet)) %>%
  ggplot(aes(y = r_Diet, x = condition, fill=condition)) +
  stat_halfeye() +
  geom_hline(yintercept = 0, linetype = "dashed")  + labs(y="Difference in Genomic Proportion\nof Recent DNA Transposon\nAccumulation", x="Mammals") + scale_fill_brewer(palette = "Set1") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position="none", legend.title=element_blank()) + scale_x_discrete(position = "top")

##save 
  ggsave("mamm_diet.pdf", h=3, w=2.5)

arti<-m.dia %>%
  spread_draws(r_Diet[condition,]) %>%
  compare_levels(r_Diet, by = condition, comparison = list(c("carnivore", "herbivore"), c("carnivore", "omnivore"), c("herbivore", "omnivore"))) %>%
  ungroup() %>%
#  mutate(condition = reorder(condition, r_Diet)) %>%
  ggplot(aes(y = r_Diet, x = condition, fill=condition)) +
  stat_halfeye() +
  geom_hline(yintercept = 0, linetype = "dashed")  + labs(x="Artiodactyla") + scale_fill_brewer(palette = "Set1") + theme(axis.title.y=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position="bottom", legend.title=element_blank()) + scale_x_discrete(position = "top")

##save
  ggsave("arti_diet.pdf", h=4, w=2)

chir<-m.dic %>%
  spread_draws(r_Diet[condition,]) %>%
  compare_levels(r_Diet, by = condition, comparison = list(c("carnivore", "herbivore"))) %>%
  ungroup() %>%
#  mutate(condition = reorder(condition, r_Diet)) %>%
  ggplot(aes(y = r_Diet, x = condition, fill=condition)) +
  stat_halfeye() +
  geom_hline(yintercept = 0, linetype = "dashed")  + labs(x="Chiroptera") + scale_fill_brewer(palette = "Set1") + theme(axis.title.y=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position="none") + scale_x_discrete(position = "top")

##save
  ggsave("chir_diet.pdf", h=3, w=1.5)

##make a palette
pal<-brewer.pal(3, "Set1")

prim<-m.dip %>%
  spread_draws(r_Diet[condition,]) %>%
  compare_levels(r_Diet, by = condition, comparison = list(c("herbivore", "omnivore"))) %>%
  ungroup() %>%
#  mutate(condition = reorder(condition, r_Diet)) %>%
  ggplot(aes(y = r_Diet, x = condition, fill=condition)) +
  stat_halfeye() +
  geom_hline(yintercept = 0, linetype = "dashed")  + labs(x="Primates") + scale_fill_manual(values=pal[3]) + theme(axis.title.y=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position="none") + scale_x_discrete(position = "top")

##save
  ggsave("prim_diet.pdf", h=3, w=1.5)

rode<-m.dir %>%
  spread_draws(r_Diet[condition,]) %>%
  compare_levels(r_Diet, by = condition, comparison = list(c("herbivore", "omnivore"))) %>%
  ungroup() %>%
#  mutate(condition = reorder(condition, r_Diet)) %>%
  ggplot(aes(y = r_Diet, x = condition, fill=condition)) +
  stat_halfeye() +
  geom_hline(yintercept = 0, linetype = "dashed")  + labs(x="Rodentia") + scale_fill_manual(values=pal[3]) + theme(axis.title.y=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position="none") + scale_x_discrete(position = "top")

##save
    ggsave("rode_diet.pdf", h=3, w=1.5)

  
##print
pdf("diet_plot.pdf", h=3, w=8)
multiplot(mamm, arti, chir, prim, rode, cols=5)
dev.off()

  

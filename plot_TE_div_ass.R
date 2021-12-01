##libraries
library(ggplot2)
library(ggthemes)
library(brms)

##get data
load("TE_diversity.RData")

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

##original conversion
##asse$size<-scale(log10(asse$Assembly_size))

##back convert values
a<-sd(log10(asse$Assembly_size))
b<-mean(log10(asse$Assembly_size))

##values in exponents are given by this
c<-seq(-2,4, by=2)*a+b

##extract model data
dat <- m.ass$data

##work with what you have
pma <- conditional_effects(m.ass) 
p1 <- plot(pma, plot = FALSE)[[1]]

assp<-p1 + 
  geom_point(
    aes(x = size, y = yDNAp), 
    # this is the key!
    data = dat, 
    shape = 21,
    # This tells it to ignore the ymin and ymax settings used elsewhere
    inherit.aes = FALSE
  ) + theme_few() + labs(x = "Genome size", y="Proportion of the genome\nattributed to TEs") + scale_x_continuous(breaks = seq(-2,4, by=2), label = c(expression(paste("2 x", 10^9)), expression(paste("3 x", 10^9)), expression(paste("4 x", 10^9)), expression(paste("5 x", 10^9)))) + scale_y_continuous(breaks = seq(0.1,0.7, by=.15))

##save plot
ggsave("assembly_plot.pdf", h=3, w=6)

##extract model data
dat <- m.sha$data

##work with what you have
pms <- conditional_effects(m.sha) 
p2 <- plot(pms, plot = FALSE)[[1]]

shap<-p2 + 
  geom_point(
    aes(x = SHANNON, y = yDNAp), 
    # this is the key!
    data = dat, 
    shape = 21,
    # This tells it to ignore the ymin and ymax settings used elsewhere
    inherit.aes = FALSE
  ) + theme_few() + labs(x = "Shannon Diversity Index", y="Proportion of genome attributed\nto recently accumulated TEs")

##save plot
ggsave("shannon_plot.pdf", h=3, w=3.15)

##extract model data
dat <- m.pie$data

##work with what you have
pmp <- conditional_effects(m.pie) 
p3 <- plot(pmp, plot = FALSE)[[1]]

piep<-p3 + 
  geom_point(
    aes(x = pielou, y = yDNAp), 
    # this is the key!
    data = dat, 
    shape = 21,
    # This tells it to ignore the ymin and ymax settings used elsewhere
    inherit.aes = FALSE
  ) + theme_few() + labs(x = "Pielou's Evenness Index") + theme(axis.title.y=element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

##save plot
ggsave("pielou_plot.pdf", h=3, w=2.5)

##print
pdf("diversity_plot.pdf", h=3, w=6)
multiplot(shap, piep, cols=2)
dev.off()

sink("assembly_diversity.txt")
print(summary(m.ass))
print(summary(m.sha))
print(summary(m.pie))
sink()

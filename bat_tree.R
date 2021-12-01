library(ggtree)
library(treeio)
library(ggplot2)
library(viridis)
library(ggtreeExtra)

setwd("C:/Users/austi/Dropbox/250_mammals_TE_validations/R_trees")

my_tree <- read.tree("C:/Users/austi/Dropbox/250_mammals_TE_validations/R_trees/bat_tree.nwk")
young_bat_data = read.table("C:/Users/austi/Dropbox/250_mammals_TE_validations/July2021_analyses/bat_data.txt",header=T,stringsAsFactors = T)

# Just a stacked bar chart
young.g.bar = ggplot(young_bat_data,aes(x=SPECIES,y=PROP,fill=TE)) + geom_bar(stat="identity",alpha=0.5)
young.g.bar

# Just the tree, with tip labels offset to make room for the bar plot
p.tree = ggtree(my_tree, branch.length=0.005) + geom_tiplab(size=2.5, align=TRUE, linesize=.5, offset=1.5) + xlim_expand(c(0,15))

p.tree.bar = p.tree + geom_fruit(data=young_bat_data, geom=geom_bar, mapping=aes(x=PROP,fill=TE,y=SPECIES), orientation="y", stat="identity", offset=.005) + scale_fill_manual("TE Class", values = c("DNA" = 'darkred', "LINE" = 'dodgerblue', "LTR" = 'maroon1', "RC" = 'green', "SINE" = 'gold'))
p.tree.bar

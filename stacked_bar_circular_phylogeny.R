library(ggtree)
library(treeio)
library(ggplot2)
library(viridis)
library(ggtreeExtra)

# Read in tree and data
my_tree <- read.tree("C:/Users/austi/Dropbox/250_mammals_TE_validations/R_trees/species_tree_v2.nwk")
total_data = read.table("C:/Users/austi/Dropbox/250_mammals_TE_validations/figures_7oct2021/total_te_data.txt",header=T,stringsAsFactors = T)
young_data = read.table("C:/Users/austi/Dropbox/250_mammals_TE_validations/figures_7oct2021/young_te_data.txt",header=T,stringsAsFactors = T)

# Just a stacked bar chart
total.g.bar = ggplot(total_data,aes(x=SPECIES,y=PROP,fill=TE)) + geom_bar(stat="identity",alpha=0.5)
total.g.bar

young.g.bar = ggplot(young_data,aes(x=SPECIES,y=PROP,fill=TE)) + geom_bar(stat="identity",alpha=0.5)
young.g.bar

#Set up dominant TE colors
my_info <- read.delim("C:/Users/austi/Dropbox/250_mammals_TE_validations/July2021_analyses/shannon_data.txt")
cols <- c(LINE='dodgerblue', SINE='gold', LTR='maroon1', DNA='darkred', RC='green')

# Just the tree, with tip labels offset to make room for the bar plot
p.tree = ggtree(my_tree, layout='circular') %<+% my_info + geom_tippoint(aes(color=DOM_TE), size=1.5) + scale_color_manual(values=cols) + geom_tiplab(size=2.5, align=TRUE, linesize=.5, offset=9.4)

p.tree.bar = p.tree + geom_fruit(data=total_data, geom=geom_bar, offset=0.5, mapping=aes(x=PROP,fill=TE,y=SPECIES), orientation="y", stat="identity") + geom_fruit(data=young_data, geom=geom_bar, offset=0.5, mapping=aes(x=PROP,fill=TE,y=SPECIES), orientation="y", stat="identity") + scale_fill_manual("TE Class", values = c("DNA" = 'darkred', "LINE" = 'dodgerblue', "LTR" = 'maroon1', "RC" = 'green', "SINE" = 'gold'))
p.tree.bar

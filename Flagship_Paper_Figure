#Circle Tree for flagship paper

library(ggplot2)
library(ggtree)
library(ggnewscale)

setwd("~/")

my_tree <- read.tree("~/species_tree_v2.nwk")
my_info <- read.delim("~/total_data.txt")

p <- ggtree(my_tree, layout='circular') %<+% my_info + geom_tiplab2(size=2.5, offset=4)

total_data <- read.delim("~/total_data.txt", row.names=1)

pgp1 <- gheatmap(p, total_data[, "TE_proportion", drop=F], offset=0, width=.05, colnames_angle = 0, colnames_offset_x = -0.5, font.size=3, low='green', high='red', color='white')) #+ scale_fill_viridis_c(option="B", name = "Total TE")
pgp2 <- pgp1 + new_scale_fill()
pgp3 <- gheatmap(pgp2, total_data[, "LINE_proportion", drop=F], offset = 1, width = 0.05, colnames_angle = 0, colnames_offset_x = 1.1, font.size = 3, low='green', high='red', color='white')) #+ scale_fill_viridis_c(option = "B", name="LINE") 
pgp4 <- pgp3 + new_scale_fill()
pgp5 <- gheatmap(pgp4, total_data[, "SINE_proportion", drop=F], offset = 1, width = 0.05, colnames_angle = 0, colnames_offset_x = 1.1, font.size = 3, low='green', high='red', color='white')) #+ scale_fill_viridis_c(option = "B", name="SINE") 
pgp6 <- pgp5 + new_scale_fill()
pgp7 <- gheatmap(pgp6, total_data[, "LTR_proportion", drop=F], offset = 1, width = 0.05, colnames_angle = 0, colnames_offset_x = 1.1, font.size = 3, low='green', high='red', color='white')) #+ scale_fill_viridis_c(option = "B", name="LTR") 
pgp8 <- pgp7 + new_scale_fill()
pgp9 <- gheatmap(pgp8, total_data[, "DNA_proportion", drop=F], offset = 1, width = 0.05, colnames_angle = 0, colnames_offset_x = 1.1, font.size = 3, low='green', high='red', color='white')) #+ scale_fill_viridis_c(option = "B", name="DNA") 
pgp10 <- pgp9 + new_scale_fill()
pgp_all <- gheatmap(pgp10, total_data[, "DNA_proportion", drop=F], offset = 2, width = 0.05, colnames_angle = 0, colnames_offset_x = 1.1, font.size = 3, low='green', high='red', color='white')) #+ scale_fill_viridis_c(option = "B", name="RC") 
pgp_all

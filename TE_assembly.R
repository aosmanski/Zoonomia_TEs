##libraries
library(phytools)
library(geiger)
library(brms)

##clear out cache
rm(list=ls())

##get data
asse<-read.table("fig1_data_liliana_v2.txt", header=T)

##for treedata
rownames(asse)<-asse$Species

##tree
tree<-read.tree("ChrX_topology_241_species_PhyloP_20210208/Zoonomia_ChrX_lessGC40_241species_30Consensus.tree")

##add two tips
##find river dolphins
n1<-findMRCA(tree, tips=c("Inia_geoffrensis", "Lipotes_vexillifer"))

##now bind using branch lengths
tre1<-bind.tip(tree,"Pontoporia_blainvillei",where=n1, position= 0.5*tree$edge.length[which(tree$edge[,2]==n1)], edge.length= 0.5*tree$edge.length[which(tree$edge[,2]==n1)] + tree$edge.length[which(tree$edge[,2]==getDescendants(tree,n1)[1])])

##find sloths
n2<-findMRCA(tre1, tips=c("Choloepus_hoffmanni", "Choloepus_didactylus"))

##now bind using branch lengths
tre2<-bind.tip(tre1,"Bradypus_variegatus",where=n2, position= 0.5*tre1$edge.length[which(tre1$edge[,2]==n2)], edge.length= 0.5*tre1$edge.length[which(tre1$edge[,2]==n2)] + tre1$edge.length[which(tre1$edge[,2]==getDescendants(tre1,n2)[1])])

##find and bind Rhinolophus_sinicus
n3<- which(tre2$tip.label=="Rhinolophus_sinicus")

tre3<-bind.tip(tre2,"Rhinolophus_ferrumequinum",where=n3, position= 0.5*tre2$edge.length[which(tre2$edge[,2]==n3)], edge.length= 0.5*tre2$edge.length[which(tre2$edge[,2]==n3)] + tre2$edge.length[which(tre2$edge[,2]==getDescendants(tre2,n3)[1])])

##find and bind Molossus_molossus
n4<- which(tre3$tip.label=="Tadarida_brasiliensis")

tre4<-bind.tip(tre3,"Molossus_molossus",where=n4, position= 0.5*tre3$edge.length[which(tre3$edge[,2]==n4)], edge.length= 0.5*tre3$edge.length[which(tre3$edge[,2]==n4)] + tre3$edge.length[which(tre3$edge[,2]==getDescendants(tre3,n4)[1])])

##find and bind Pipistrellus_kuhlii
n5<- which(tre4$tip.label=="Pipistrellus_pipistrellus")

tre5<-bind.tip(tre4,"Pipistrellus_kuhlii",where=n5, position= 0.5*tre4$edge.length[which(tre4$edge[,2]==n5)], edge.length= 0.5*tre4$edge.length[which(tre4$edge[,2]==n5)] + tre4$edge.length[which(tre4$edge[,2]==getDescendants(tre4,n5)[1])])

##find and bind Phyllostomus_discolor
n6<- which(tre5$tip.label=="Tonatia_saurophila")

tre6<-bind.tip(tre5,"Phyllostomus_discolor",where=n6, position= 0.5*tre5$edge.length[which(tre5$edge[,2]==n6)], edge.length= 0.5*tre5$edge.length[which(tre5$edge[,2]==n6)] + tre5$edge.length[which(tre5$edge[,2]==getDescendants(tre5,n6)[1])])

##find and bind Rhizomys_pruinosus
n7<- which(tre6$tip.label=="Nannospalax_galili")

tre7<-bind.tip(tre6,"Rhizomys_pruinosus",where=n7, position= 0.5*tre6$edge.length[which(tre6$edge[,2]==n7)], edge.length= 0.5*tre6$edge.length[which(tre6$edge[,2]==n7)] + tre6$edge.length[which(tre6$edge[,2]==getDescendants(tre6,n7)[1])])

##make sure names match and prune tree
tre7<-treedata(tre7,asse)$phy

##make ultrametric
time<-chronos(tre7)

##clean up
attr(time,"order") <- NULL
attr(time,"call") <- NULL
attr(time,"ploglik") <- NULL
attr(time,"rates") <- NULL
attr(time,"message") <- NULL
attr(time,"PHIIC") <- NULL
attr(time,"class") <- "phylo"

##get phylogenetic matrix
A <- ape::vcv.phylo(time)

##names brms doesn't like
colnames(asse)[5]<-"yDNAp"

asse$size<-scale(log10(asse$Assembly_size))

##young DNA as a function of diet
##geiger also has a bf function
pro.s <- brms:::bf(yDNAp ~ size + (1|gr(Species, cov = A)))

##run uni models together at a time
m.ass <- brm(pro.s, 
            family = "beta", 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 2000, thin = 10, 
                  control = list(adapt_delta = 0.9999, max_treedepth=20),
                  data = asse,
                  cores = 4, chains = 4)

##save
save.image("TE_assembly.RData")

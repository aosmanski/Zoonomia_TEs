##libraries
library(phytools)
library(geiger)
library(brms)

##clear out cache
rm(list=ls())

##get data
diet<-read.table("diet_data_liliana_v2.txt", header=T)

##for treedata
rownames(diet)<-diet$Species

##tree
tree<-read.tree("ChrX_topology_241_species_PhyloP_20210208/Zoonomia_ChrX_lessGC40_241species_30Consensus.tree")

##fix so all species in data are included
##replace sinicus with ferrumequinum on tree
tree$tip.label[tree$tip.label=="Rhinolophus_sinicus"]<-"Rhinolophus_ferrumequinum"

##add two tips
##find river dolphins
n1<-findMRCA(tree, tips=c("Inia_geoffrensis", "Lipotes_vexillifer"))

##now bind using branch lengths
tre1<-bind.tip(tree,"Pontoporia_blainvillei",where=n1, position= 0.5*tree$edge.length[which(tree$edge[,2]==n1)], edge.length= 0.5*tree$edge.length[which(tree$edge[,2]==n1)] + tree$edge.length[which(tree$edge[,2]==getDescendants(tree,n1)[1])])

##find sloths
n2<-findMRCA(tree, tips=c("Choloepus_hoffmanni", "Choloepus_didactylus"))

##now bind using branch lengths
tre2<-bind.tip(tre1,"Bradypus_variegatus",where=n2, position= 0.5*tre1$edge.length[which(tre1$edge[,2]==n2)], edge.length= 0.5*tre1$edge.length[which(tre1$edge[,2]==n2)] + tre1$edge.length[which(tre1$edge[,2]==getDescendants(tre1,n2)[1])])

##make ultrametric
time<-chronos(tre2)

##clean up tree
attr(time,"order") <- NULL
attr(time,"call") <- NULL
attr(time,"ploglik") <- NULL
attr(time,"rates") <- NULL
attr(time,"message") <- NULL
attr(time,"PHIIC") <- NULL
attr(time,"class") <- "phylo"

##make sure names match and prune tree
tre3<-treedata(time,diet)$phy

##get phylogenetic matrix
A <- ape::vcv.phylo(tre3)

##names brms doesn't like
colnames(diet)[5]<-"yDNAp"

##young DNA as a function of diet
##geiger also has a bf function
pro.d <- brms:::bf(yDNAp ~ -1 + (1|gr(Diet)) + (1|gr(Species, cov = A)))

##run uni model
m.die <- brm(pro.d, 
            family = "beta", 
                  data2 = list(A = A),
                  iter  =5000, warmup = 1000, thin = 10, 
                  control = list(adapt_delta = 0.9999, max_treedepth=20),
                  data = diet,
                  cores = 4, chains = 4)

##save
save.image("TE_diet.RData")

##other groups B) Cetartiodactyla (KW, p<0.001); C) Chiroptera (W, p=0.0002); D) Primates (W, p=0.609); E) Rodentia

##artiodactyla
dieA<-subset(diet, Order=="Artiodactyla")

##make sure names match and prune tree
treA<-treedata(tre3,dieA)$phy

##get phylogenetic matrix
Aa <- ape::vcv.phylo(treA)

##young DNA as a function of diet
##geiger also has a bf function
pro.a <- brms:::bf(yDNAp ~ -1 + (1|gr(Diet)) + (1|gr(Species, cov = Aa)))

##run uni model
m.dia <- brm(pro.a, 
            family = "beta", 
                  data2 = list(Aa = Aa),
                  iter  =5000, warmup = 1000, thin = 10, 
                  control = list(adapt_delta = 0.9999, max_treedepth=20),
                  data = dieA,
                  cores = 4, chains = 4)

##save
save.image("TE_diet.RData")

##chiroptera
dieC<-subset(diet, Order=="Chiroptera")

##make sure names match and prune tree
treC<-treedata(tre3,dieC)$phy

##get phylogenetic matrix
Ac <- ape::vcv.phylo(treC)

##young DNA as a function of diet
##geiger also has a bf function
pro.c <- brms:::bf(yDNAp ~ -1 + (1|gr(Diet)) + (1|gr(Species, cov = Ac)))

##run uni model
m.dic <- brm(pro.c, 
            family = "beta", 
                  data2 = list(Ac = Ac),
                  iter  =5000, warmup = 1000, thin = 10, 
                  control = list(adapt_delta = 0.9999, max_treedepth=20),
                  data = dieC,
                  cores = 4, chains = 4)

##save
save.image("TE_diet.RData")

##primates
dieP<-subset(diet, Order=="Primates")

##make sure names match and prune tree
treP<-treedata(tre3,dieP)$phy

##get phylogenetic matrix
Ap <- ape::vcv.phylo(treP)

##young DNA as a function of diet
##geiger also has a bf function
pro.p <- brms:::bf(yDNAp ~ -1 + (1|gr(Diet)) + (1|gr(Species, cov = Ap)))

##run uni model
m.dip <- brm(pro.p, 
            family = "beta", 
                  data2 = list(Ap = Ap),
                  iter  =5000, warmup = 1000, thin = 10, 
                  control = list(adapt_delta = 0.9999, max_treedepth=20),
                  data = dieP,
                  cores = 4, chains = 4)

##save
save.image("TE_diet.RData")

##rodents
dieR<-subset(diet, Order=="Rodentia")

##make sure names match and prune tree
treR<-treedata(tre3,dieR)$phy

##get phylogenetic matrix
Ar <- ape::vcv.phylo(treR)

##young DNA as a function of diet
##geiger also has a bf function
pro.r <- brms:::bf(yDNAp ~ -1 + (1|gr(Diet)) + (1|gr(Species, cov = Ar)))

##run uni model
m.dir <- brm(pro.r, 
            family = "beta", 
                  data2 = list(Ar = Ar),
                  iter  =5000, warmup = 1000, thin = 10, 
                  control = list(adapt_delta = 0.9999, max_treedepth=20),
                  data = dieR,
                  cores = 4, chains = 4)

##save
save.image("TE_diet.RData")

##clear out cache
rm(list=ls())

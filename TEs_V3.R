##libraries
library(MCMCglmm)
library(geiger)
library(MCMCvis)

##clear out cache
rm(list=ls())

##data
dat<-read.csv("Mammal_TEcountdata_Liliana_v2.csv")

##tree
tre<-read.tree("ChrX_topology_241_species_PhyloP_20210208/Zoonomia_ChrX_lessGC40_241species_30Consensus.tree")

##make ultrametric
tim<-chronos(tre)

##clean up
attr(tim,"order") <- NULL
attr(tim,"call") <- NULL
attr(tim,"ploglik") <- NULL
attr(tim,"rates") <- NULL
attr(tim,"message") <- NULL
attr(tim,"PHIIC") <- NULL
attr(tim,"class") <- "phylo"

##make matrix
inv.phylo <- inverseA(tim, nodes = "TIPS", scale = T)


##priors
##one random effect
prior <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002)))

##two random effects
prior2 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002)))

##three random effects
prior3 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002)))

##two random effects us structure
prior4 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = diag(5) * 0.02, nu = 4)))

##just family in random
m1 <- MCMCglmm(Number_of_Hits ~ TE_Type + log(Genome_Size), random = ~Species + Family, pr = TRUE,
                   family = "poisson", ginverse = list(Species = inv.phylo$Ainv),
                   prior=prior2, data = dat, nitt = 50000, burnin = 10000, thin = 100, verbose=F)

##save
save.image("TE_v3.RData")

##does the identity of the family explain different amounts of variation depending on te type?

##model te type within family as in (1|sex:dam)
m2 <- MCMCglmm(Number_of_Hits ~ TE_Type + log(Genome_Size), random = ~Species + TE_Type:Family, pr = TRUE,
                   family = "poisson", ginverse = list(Species = inv.phylo$Ainv),
                   prior=prior2, data = dat, nitt = 50000, burnin = 10000, thin = 100, verbose=F)

##save
save.image("TE_v3.RData")

##model te type within family as in (1|dam)+(1|sex:dam) 
m3 <- MCMCglmm(Number_of_Hits ~ TE_Type + log(Genome_Size), random = ~Species + Family + TE_Type:Family, pr = TRUE,
                   family = "poisson", ginverse = list(Species = inv.phylo$Ainv),
                   prior=prior3, data = dat, nitt = 50000, burnin = 10000, thin = 100, verbose=F)

##save
save.image("TE_v3.RData")                   

##model te type within family as in (sex-1|dam) or us structure
m4 <- MCMCglmm(Number_of_Hits ~ TE_Type + log(Genome_Size), random = ~Species + us(TE_Type):Family, pr = TRUE,
                   family = "poisson", ginverse = list(Species = inv.phylo$Ainv),
                   prior=prior4, data = dat, nitt = 500000, burnin = 100000, thin = 500, verbose=F)

##save
save.image("TE_v3.RData")                   

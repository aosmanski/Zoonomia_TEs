##libraries
library(phytools)
library(geiger)
library(brms)

##clear out cache
rm(list=ls())

##get data
dive<-read.table("fig5_data_liliana.txt", header=T)

##same spp missing as in the fig 1 dataset
load("TE_assembly.RData")

##names brms doesn't like
colnames(dive)[5]<-"pielou"
colnames(dive)[6]<-"yDNAp"

##young DNA as a function of diversity
##geiger also has a bf function
pro.h <- brms:::bf(yDNAp ~ SHANNON + (1|gr(Species, cov = A)))

##run uni models together at a time
m.sha <- brm(pro.h, 
            family = "beta", 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 2000, thin = 10, 
                  control = list(adapt_delta = 0.9999, max_treedepth=20),
                  data = dive,
                  cores = 4, chains = 4)

##save
save.image("TE_diversity.RData")

pro.p <- brms:::bf(yDNAp ~ pielou + (1|gr(Species, cov = A)))

##run uni models together at a time
m.pie <- brm(pro.p, 
            family = "beta", 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 2000, thin = 10, 
                  control = list(adapt_delta = 0.9999, max_treedepth=20),
                  data = dive,
                  cores = 4, chains = 4)

##save
save.image("TE_diversity.RData")

##clear out cache
rm(list=ls())

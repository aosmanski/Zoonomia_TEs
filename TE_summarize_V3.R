##libraries
library(MCMCglmm)
library(MCMCvis)
library(bayesplot)
library(tidyverse)

##load results
load("TE_v3.RData")

##print comparisons across models
sink("model_DIC.txt")
print(" TE type and genome size in fixed, TE and family no interaction")
print(summary(m1)$DIC)

print("TE type and genome size in fixed, TE type within family as in (1|type:family)")
print(summary(m2)$DIC)

print("TE type and genome size in fixed, TE type within family as in (1|family)+(1|type:family)")
print(summary(m3)$DIC)

print("TE type and genome size in fixed, TE type within family as in (type-1|family) or us structure")
print(summary(m4)$DIC)

sink()

###model 2 is the best

##print family and type:family coefficients 
sink("model_coefs.txt")
#print(MCMCsummary(m2$Sol, params="Family", exact=F, excl="TE_Type:Family."))
print("sample wide effects")
print(summary(m2))
print("cluster specific effects")
print(MCMCsummary(m2$Sol, params="TE_Type:Family", exact=F))
sink()

##plot family and type:family coefficients 
#pdf("family_coefs.pdf", h=12, w=4)
#MCMCplot(m2$Sol, params="Family", exact=F, excl="TE_Type:Family.", ref_ovl=T, sz_labels = 0.5)
#dev.off()

pdf("te_family_coefs.pdf", h=48, w=4)
MCMCplot(m2$Sol, params="TE_Type:Family", exact=F, ref_ovl=T, sz_labels = 0.5)
dev.off()

pdf("te_family_coefs_break.pdf", h=12, w=4)
MCMCplot(m2$Sol, params="TE_Type:Family.DNA", exact=F, ref_ovl=T, sz_labels = 0.5)
MCMCplot(m2$Sol, params="TE_Type:Family.LINE", exact=F, ref_ovl=T, sz_labels = 0.5)
MCMCplot(m2$Sol, params="TE_Type:Family.SINE", exact=F, ref_ovl=T, sz_labels = 0.5)
MCMCplot(m2$Sol, params="TE_Type:Family.LTR", exact=F, ref_ovl=T, sz_labels = 0.5)
MCMCplot(m2$Sol, params="TE_Type:Family.RC", exact=F, ref_ovl=T, sz_labels = 0.5)
dev.off()

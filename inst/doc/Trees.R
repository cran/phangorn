###################################################
### chunk number 1: 
###################################################
options(width=70)
foo <- packageDescription("phangorn")


###################################################
### chunk number 2: 
###################################################
library(phangorn)
require("multicore")
primates = read.phyDat("primates.dna", format = "phylip", type = "DNA")


###################################################
### chunk number 3: 
###################################################
dm = dist.dna(as.DNAbin(primates))
treeUPGMA = upgma(dm)
treeNJ = NJ(dm)


###################################################
### chunk number 4: plotNJ
###################################################
par(mfrow =c(1,2), mar = c(1,1,4,1))
plot(treeUPGMA, main="UPGMA")
plot(treeNJ, "unrooted", main="NJ")


###################################################
### chunk number 5: figNJ
###################################################
par(mfrow =c(1,2), mar = c(1,1,4,1))
plot(treeUPGMA, main="UPGMA")
plot(treeNJ, "unrooted", main="NJ")


###################################################
### chunk number 6: 
###################################################
parsimony(treeUPGMA, primates)
parsimony(treeNJ, primates)


###################################################
### chunk number 7: 
###################################################
treePars = optim.parsimony(treeUPGMA, primates)
parsimony(treePars, primates)


###################################################
### chunk number 8: 
###################################################
fit = pml(treeNJ, data=primates)
fit


###################################################
### chunk number 9: 
###################################################
methods(class="pml")


###################################################
### chunk number 10: 
###################################################
fitJC = optim.pml(fit, TRUE)
logLik(fitJC)


###################################################
### chunk number 11: 
###################################################
fitGTR = update(fit, k=4, inv=0.2) 
fitGTR = optim.pml(fitGTR, TRUE,TRUE, TRUE, TRUE, TRUE)


###################################################
### chunk number 12: 
###################################################
fitGTR 


###################################################
### chunk number 13: 
###################################################
anova(fitJC, fitGTR) 


###################################################
### chunk number 14: 
###################################################
AIC(fitGTR) 
AIC(fitJC)


###################################################
### chunk number 15: 
###################################################
SH.test(fitGTR, fitJC) 


###################################################
### chunk number 16: 
###################################################
bs = bootstrap.pml(fitJC, bs=100, optNni=TRUE)


###################################################
### chunk number 17: plotBS
###################################################
par(mar=c(.1,.1,.1,.1))
plotBS(fitJC$tree, bs)


###################################################
### chunk number 18: figBS
###################################################
par(mar=c(.1,.1,.1,.1))
plotBS(fitJC$tree, bs)


###################################################
### chunk number 19: 
###################################################
options(prompt=" ")
options(continue="  ")


###################################################
### chunk number 20:  eval=FALSE
###################################################
## library(multicore) # supports parallel computing
## library(phangorn)
## file="myfile"
## dat = read.phyDat(file)
## dm = dist.ml(dat)
## tree = NJ(dm)
## fitNJ = pml(tree, dat, k=4, inv=.2)
## fit = optim.pml(fitNJ,TRUE,TRUE,TRUE,TRUE,TRUE)
## fit
## bs = bootstrap.pml(fit, bs=100, optNni=TRUE)


###################################################
### chunk number 21:  eval=FALSE
###################################################
## library(multicore) # supports parallel computing
## library(phangorn)
## file="myfile"
## dat = read.phyDat(file, type = "AA")
## dm = dist.ml(dat, model="JTT")
## tree = NJ(dm)
## fitNJ = pml(tree, dat, model="JTT", k=4, inv=.2)
## fit = optim.pml(fitNJ, optNni=TRUE, optInv=TRUE, optGamma=TRUE)
## fit
## bs = bootstrap.pml(fit, bs=100, optNni=TRUE)


###################################################
### chunk number 22: 
###################################################
toLatex(sessionInfo())



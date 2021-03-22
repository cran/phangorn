## ----setup, echo=FALSE--------------------------------------------------------
# set global chunk options: images will be bigger
knitr::opts_chunk$set(fig.width=6, fig.height=4)
#, global.par=TRUE
options(digits = 4)

## -----------------------------------------------------------------------------
library(phangorn)
data = matrix(c("r","a","y","g","g","a","c","-","c","t","c","g",
    "a","a","t","g","g","a","t","-","c","t","c","a",
    "a","a","t","-","g","a","c","c","c","t","?","g"),
    dimnames = list(c("t1", "t2", "t3"),NULL), nrow=3, byrow=TRUE)
data

## -----------------------------------------------------------------------------
gapsdata1 = phyDat(data)
gapsdata1

## -----------------------------------------------------------------------------
gapsdata2 = phyDat(data, type="USER", levels=c("a","c","g","t","-"),
    ambiguity = c("?", "n"))
gapsdata2

## -----------------------------------------------------------------------------
contrast = matrix(data = c(1,0,0,0,0,
    0,1,0,0,0,
    0,0,1,0,0,
    0,0,0,1,0,
    1,0,1,0,0,
    0,1,0,1,0,
    0,0,0,0,1,
    1,1,1,1,0,
    1,1,1,1,1),
    ncol = 5, byrow = TRUE)
dimnames(contrast) = list(c("a","c","g","t","r","y","-","n","?"),
    c("a", "c", "g", "t", "-"))
contrast
gapsdata3 = phyDat(data, type="USER", contrast=contrast)
gapsdata3

## -----------------------------------------------------------------------------
library(ape)
tree = unroot(rtree(3))
fit = pml(tree, gapsdata3)
fit = optim.pml(fit, optQ=TRUE, subs=c(1,0,1,2,1,0,2,1,2,2),
    control=pml.control(trace=0))
fit

## -----------------------------------------------------------------------------
library(phangorn)
fdir <- system.file("extdata/trees", package = "phangorn")
primates <- read.phyDat(file.path(fdir, "primates.dna"),
                        format = "interleaved")
tree <- NJ(dist.ml(primates))
dat <- dna2codon(primates)
fit <- pml(tree, dat, bf="F3x4")
fit0 <- optim.pml(fit, model="codon0", control=pml.control(trace=0))
fit1 <- optim.pml(fit, model="codon1", control=pml.control(trace=0))
fit2 <- optim.pml(fit, model="codon2", control=pml.control(trace=0))
fit3 <- optim.pml(fit, model="codon3", control=pml.control(trace=0))
anova(fit0, fit2, fit3, fit1)

## -----------------------------------------------------------------------------
trees <- allTrees(5)
par(mfrow=c(3,5), mar=rep(0,4))
for(i in 1:15)plot(trees[[i]], cex=1, type="u")

## -----------------------------------------------------------------------------
nni(trees[[1]])

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()


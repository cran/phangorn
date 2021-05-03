## -----------------------------------------------------------------------------
library(ape)
library(phangorn)
fdir <- system.file("extdata/trees", package = "phangorn")
primates <- read.phyDat(file.path(fdir, "primates.dna"),
                        format = "interleaved")

## -----------------------------------------------------------------------------
dm  <- dist.ml(primates)
treeUPGMA  <- upgma(dm)
treeNJ  <- NJ(dm)

## ----plot1, fig.cap="Rooted UPGMA tree.", echo=TRUE---------------------------
plot(treeUPGMA, main="UPGMA")

## ----plot2, fig.cap="Unrooted NJ tree.", echo=TRUE----------------------------
plot(treeNJ, "unrooted", main="NJ")

## -----------------------------------------------------------------------------
parsimony(treeUPGMA, primates)
parsimony(treeNJ, primates)

## -----------------------------------------------------------------------------
treePars  <- optim.parsimony(treeUPGMA, primates)
treeRatchet  <- pratchet(primates, trace = 0)
parsimony(c(treePars, treeRatchet), primates)

## -----------------------------------------------------------------------------
(trees <- bab(subset(primates,1:10)))

## -----------------------------------------------------------------------------
fit = pml(treeNJ, data=primates)
fit

## -----------------------------------------------------------------------------
methods(class="pml")

## -----------------------------------------------------------------------------
fitJC  <- optim.pml(fit, TRUE)
logLik(fitJC)

## -----------------------------------------------------------------------------
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
    rearrangement = "NNI", control = pml.control(trace = 0))
fitGTR

## -----------------------------------------------------------------------------
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
    rearrangement = "stochastic", control = pml.control(trace = 0))
fitGTR

## -----------------------------------------------------------------------------
anova(fitJC, fitGTR)

## -----------------------------------------------------------------------------
SH.test(fitGTR, fitJC)

## -----------------------------------------------------------------------------
AIC(fitJC)
AIC(fitGTR)
AICc(fitGTR)
BIC(fitGTR)

## ---- echo=FALSE--------------------------------------------------------------
load("Trees.RData")

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  mt = modelTest(primates)

## ---- echo=FALSE--------------------------------------------------------------
library(knitr)
kable(mt, digits=2)

## ---- echo=TRUE---------------------------------------------------------------
env <- attr(mt, "env")
ls(envir=env)
(fit <- eval(get("HKY+G+I", env), env))

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  bs = bootstrap.pml(fitJC, bs=100, optNni=TRUE,
#      control = pml.control(trace = 0))

## ----plotBS, fig.cap="Tree with bootstrap support. Unrooted tree (midpoint rooted) with bootstrap support values.", echo=TRUE----
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")

## ----ConsensusNet, fig.cap="ConsensusNet from the bootstrap sample.", echo=TRUE----
cnet <- consensusNet(bs, p=0.2)
plot(cnet, show.edge.label=TRUE)

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  library(phangorn)
#  file <- "myfile"
#  dat <- read.phyDat(file)
#  dm <- dist.ml(dat, "F81")
#  tree <- NJ(dm)
#  # as alternative for a starting tree:
#  tree <- pratchet(dat)          # parsimony tree
#  tree <- nnls.phylo(tree, dm)   # need edge weights
#  
#  
#  # 1. alternative: quick and dirty: GTR + G
#  fitStart <- pml(tree, dat, k=4)
#  fit <- optim.pml(fitStart, model="GTR", optGamma=TRUE, rearrangement="stochastic")
#  
#  # 2. alternative: preper with modelTest
#  mt <- modelTest(dat, tree=tree, multicore=TRUE)
#  mt[order(mt$AICc),]
#  # choose best model from the table according to AICc
#  bestmodel <- mt$Model[which.min(mt$AICc)]
#  
#  env <- attr(mt, "env")
#  fitStart <- eval(get("GTR+G+I", env), env)
#  
#  # or let R search the table
#  fitStart <- eval(get(bestmodel, env), env)
#  # equivalent to:   fitStart = eval(get("GTR+G+I", env), env)
#  fit <- optim.pml(fitStart, rearrangement = "stochastic",
#      optGamma=TRUE, optInv=TRUE, model="GTR")
#  bs <- bootstrap.pml(fit, bs=100, optNni=TRUE, multicore=TRUE)

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  library(phangorn)
#  file <- "myfile"
#  dat <- read.phyDat(file, type = "AA")
#  dm <- dist.ml(dat, model="JTT")
#  tree <- NJ(dm)
#  
#  # parallel will only work safely from command line
#  # and not at all windows
#  (mt <- modelTest(dat, model=c("JTT", "LG", "WAG"),
#      multicore=TRUE))
#  # run all available amino acid models
#  (mt <- modelTest(dat, model="all", multicore=TRUE))
#  
#  env <- attr(mt, "env")
#  fitStart <- eval(get(mt$Model[which.min(mt$BIC)], env), env)
#  
#  fitNJ <- pml(tree, dat, model="JTT", k=4, inv=.2)
#  fit <- optim.pml(fitNJ, rearrangement = "stochastic",
#      optInv=TRUE, optGamma=TRUE)
#  fit
#  
#  bs <- bootstrap.pml(fit, bs=100, optNni=TRUE, multicore=TRUE)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()


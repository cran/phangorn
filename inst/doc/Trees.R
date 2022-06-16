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

## ----bootstrap_dist, echo=TRUE------------------------------------------------
fun <- function(x) upgma(dist.ml(x))
bs_upgma <- bootstrap.phyDat(primates,  fun)

## ----bootstrap_dist_new, echo=TRUE, eval=FALSE--------------------------------
#  bs_upgma <- bootstrap.phyDat(primates,  \(x){dist.ml(x) |> upgma})

## ----plot_bs, fig.cap="Rooted UPGMA tree.", echo=TRUE-------------------------
plotBS(treeUPGMA, bs_upgma, main="UPGMA")

## -----------------------------------------------------------------------------
parsimony(treeUPGMA, primates)
parsimony(treeNJ, primates)

## -----------------------------------------------------------------------------
treePars  <- optim.parsimony(treeUPGMA, primates)

## ----pratchet-----------------------------------------------------------------
treeRatchet  <- pratchet(primates, trace = 0)
parsimony(c(treePars, treeRatchet), primates)

## -----------------------------------------------------------------------------
treeRatchet  <- acctran(treeRatchet, primates)

## ----midpoint-----------------------------------------------------------------
plotBS(midpoint(treeRatchet), type="phylogram")

## ----bab----------------------------------------------------------------------
(trees <- bab(subset(primates,1:10)))

## ----pml----------------------------------------------------------------------
fit <- pml(treeNJ, data=primates)
fit

## -----------------------------------------------------------------------------
methods(class="pml")

## -----------------------------------------------------------------------------
fitJC  <- optim.pml(fit, rearrangement="NNI")
logLik(fitJC)

## ----GTR+G+I------------------------------------------------------------------
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
    rearrangement = "NNI", control = pml.control(trace = 0))
fitGTR

## ----stochastic---------------------------------------------------------------
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
    rearrangement = "stochastic", control = pml.control(trace = 0))
fitGTR

## ----anova--------------------------------------------------------------------
anova(fitJC, fitGTR)

## ----SH_test------------------------------------------------------------------
SH.test(fitGTR, fitJC)

## ----AIC----------------------------------------------------------------------
AIC(fitJC)
AIC(fitGTR)
AICc(fitGTR)
BIC(fitGTR)

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  mt <- modelTest(primates)

## ---- echo=FALSE--------------------------------------------------------------
load("Trees.RData")

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  mt <- modelTest(primates, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR"))

## ---- echo=FALSE--------------------------------------------------------------
library(knitr)
kable(mt, digits=2)

## ----as.pml, echo=TRUE--------------------------------------------------------
(fit <- as.pml(mt, "HKY+G+I"))
(fit <- as.pml(mt, "BIC"))

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE,
#      control = pml.control(trace = 0))

## ----plotBS, fig.cap="Tree with bootstrap support. Unrooted tree (midpoint rooted) with bootstrap support values.", echo=TRUE----
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")

## ----ConsensusNet, fig.cap="ConsensusNet from the bootstrap sample.", echo=TRUE----
cnet <- consensusNet(bs, p=0.2)
plot(cnet, show.edge.label=TRUE)

## ---- echo=TRUE, cache=TRUE---------------------------------------------------
fit_strict <- pml(treeUPGMA, data=primates, k=4, bf=baseFreq(primates))
fit_strict <- optim.pml(fit_strict, model="GTR", optRooted = TRUE, 
                        rearrangement = "NNI", optGamma = TRUE, optInv = TRUE, 
                        control = pml.control(trace = 0))

## -----------------------------------------------------------------------------
plot(fit_strict)

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
#  # 2. alternative: prepare with modelTest
#  mt <- modelTest(dat, tree=tree, multicore=TRUE)
#  mt[order(mt$AICc),]
#  # choose best model from the table according to AICc
#  fitStart <- as.pml(mt, "AICc")
#  
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
#  fitStart <- as.pml(mt)
#  
#  fitNJ <- pml(tree, dat, model="JTT", k=4, inv=.2)
#  fit <- optim.pml(fitNJ, rearrangement = "stochastic",
#      optInv=TRUE, optGamma=TRUE)
#  fit
#  
#  bs <- bootstrap.pml(fit, bs=100, optNni=TRUE, multicore=TRUE)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()


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

## ----bootstrap_dist_new, echo=TRUE, eval=FALSE, cache=TRUE--------------------
#  bs_upgma <- bootstrap.phyDat(primates,  \(x){dist.ml(x) |> upgma})

## ----plot_bs, fig.cap="Rooted UPGMA tree.", echo=TRUE-------------------------
plotBS(treeUPGMA, bs_upgma, main="UPGMA")

## -----------------------------------------------------------------------------
parsimony(treeUPGMA, primates)
parsimony(treeNJ, primates)

## ----pratchet-----------------------------------------------------------------
treeRatchet  <- pratchet(primates, trace = 0, minit=100)
parsimony(treeRatchet, primates)

## ----acctran------------------------------------------------------------------
treeRatchet  <- acctran(treeRatchet, primates)

## ----di2multi-----------------------------------------------------------------
treeRatchet  <- di2multi(treeRatchet)

## ----unique trees-------------------------------------------------------------
if(inherits(treeRatchet, "multiPhylo")){
  treeRatchet <- unique(treeRatchet)
}

## ----midpoint-----------------------------------------------------------------
plotBS(midpoint(treeRatchet), type="phylogram")
add.scale.bar()

## -----------------------------------------------------------------------------
treeRA <- random.addition(primates)
treeSPR  <- optim.parsimony(treeRA, primates)
parsimony(c(treeRA, treeSPR), primates)

## ----bab----------------------------------------------------------------------
(trees <- bab(primates[1:10,]))

## ----pml----------------------------------------------------------------------
fit <- pml(treeNJ, data=primates)
fit

## -----------------------------------------------------------------------------
methods(class="pml")

## -----------------------------------------------------------------------------
fitJC  <- optim.pml(fit, rearrangement="NNI")
logLik(fitJC)

## ----F81+G+I, cache=TRUE------------------------------------------------------
fitF81 <- update(fitJC, k=4, inv=0.2, bf=baseFreq(primates))
fitF81

## ----GTR+G+I, cache=TRUE------------------------------------------------------
fitGTR <- optim.pml(fitF81, model="GTR", optInv=TRUE, optGamma=TRUE,
    rearrangement = "NNI", control = pml.control(trace = 0))
fitGTR

## ----stochastic, cache=TRUE---------------------------------------------------
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
    rearrangement = "stochastic", control = pml.control(trace = 0))
fitGTR

## ----pml_bb, cache=TRUE-------------------------------------------------------
fitGTR <- pml_bb(primates, model="GTR+G(4)+I", control = pml.control(trace = 0))
fitGTR

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  mt <- modelTest(primates)

## ---- echo=FALSE--------------------------------------------------------------
load("Trees.RData")

## ---- echo=TRUE, eval=TRUE, cache=TRUE----------------------------------------
mt <- modelTest(primates, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR"), 
                control = pml.control(trace = 0))

## ---- echo=FALSE--------------------------------------------------------------
library(knitr)
kable(mt, digits=2)

## ----as.pml, echo=TRUE--------------------------------------------------------
(fit <- as.pml(mt, "HKY+G(4)+I"))
(fit <- as.pml(mt, "BIC"))

## ----pml_bb_modelTest, cache=TRUE---------------------------------------------
fit_mt <- pml_bb(mt, control = pml.control(trace = 0))
fit_mt

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
#  bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE,
#      control = pml.control(trace = 0))

## ----plotBS, fig.cap="Tree with bootstrap support. Unrooted tree (midpoint rooted) with bootstrap support values.", echo=TRUE----
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")

## ----ConsensusNet, fig.cap="ConsensusNet from the bootstrap sample.", echo=TRUE----
cnet <- consensusNet(bs, p=0.2)
plot(cnet, show.edge.label=TRUE)

## ---- echo=TRUE, cache=TRUE---------------------------------------------------
fit_strict <- pml_bb(primates, model="HKY+G(4)", method="ultrametric",
                     control = pml.control(trace = 0))

## -----------------------------------------------------------------------------
plot(fit_strict)

## ---- tipdated data-----------------------------------------------------------
fdir <- system.file("extdata/trees", package = "phangorn")
tmp <- read.csv(file.path(fdir,"H3N2_NA_20.csv"))
H3N2 <- read.phyDat(file.path(fdir,"H3N2_NA_20.fasta"), format="fasta")

## ---- tipdated processing-----------------------------------------------------
dates <- setNames(tmp$numdate_given, tmp$name)
head(dates)

## ---- tipdated fit------------------------------------------------------------
fit_td <- pml_bb(H3N2, model="GTR+G(4)", method="tipdated", tip.dates=dates, 
               rearrangement="NNI", control = pml.control(trace = 0))

## ---- tipdated plot-----------------------------------------------------------
tree_td <- fit_td$tree
root_time <- max(dates) - max(node.depth.edgelength(tree_td))
plot(tree_td, show.tip.label = FALSE)
axisPhylo(root.time = root_time, backward = FALSE)

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  library(phangorn)
#  file <- "myfile"
#  dat <- read.phyDat(file, format="fasta")
#  dm <- dist.ml(dat, "F81")
#  tree <- NJ(dm)
#  # as alternative for a starting tree:
#  tree <- pratchet(dat)          # parsimony tree
#  tree <- nnls.phylo(tree, dm)   # need edge weights
#  
#  # 1. alternative: quick and dirty: GTR + G
#  (fit <- pml_bb(dat, model="GTR+G"))
#  
#  # 2. alternative: choose with modelTest
#  mt <- modelTest(dat, multicore=TRUE)
#  mt[order(mt$BIC),]
#  # chooses best model from the table according to BIC (default)
#  fit <- pml_bb(mt)
#  bs <- bootstrap.pml(fit, bs=100, optNni=TRUE, multicore=TRUE)
#  
#  tree_ml <- plotBS(fit$tree, bs)
#  # export the tree
#  write.tree(tree_ml, file="tree_ml.nwk")

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  library(phangorn)
#  file <- "myfile"
#  dat <- read.phyDat(file, type = "AA", format="fasta")
#  
#  # compute a neighbor joining tree
#  dm <- dist.ml(dat, model="JTT")
#  tree <- NJ(dm)
#  
#  # parallel will only work safely from command line and not at all windows
#  # Select specific models.
#  (mt <- modelTest(dat, model=c("JTT", "LG", "WAG"), multicore=TRUE))
#  # run all available amino acid models
#  (mt <- modelTest(dat, model="all", multicore=TRUE))
#  
#  (fit <- pml_bb(mt))
#  # non-parametric bootstrap
#  bs <- bootstrap.pml(fit, bs=100, optNni=TRUE, multicore=TRUE)
#  
#  tree_ml <- plotBS(fit$tree, bs)
#  # export the tree
#  write.tree(tree_ml, file="tree_ml.nwk")

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()


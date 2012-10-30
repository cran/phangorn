fnodesNew <- function (tree, weight, nr, external = FALSE) 
{
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]

    pars <- integer(nr)
    root <- as.integer(node[!match(node, edge, 0)][1])
    nTips = length(tree$tip)   

    n = length(node)
    m= as.integer(max(tree$edge)+1L)

    tmp <- .C("fnhelp", as.integer(node), as.integer(edge), as.integer(n), as.integer(m), root, 
        integer(2L*n), integer(2L*n), integer(2L*n))
    pc <- tmp[[8]]
    edge2 <- tmp[[6]]
    node2 <- tmp[[7]]

    if(!external){
        ind <- which(node2>nTips)
        pc <- pc[ind]
        edge2 <- edge2[ind]   
        node2 <- node2[ind]
    }

    m2 = length(node2)
    .Call("FNALL3", as.integer(nr), as.integer(node), as.integer(edge), as.integer(node2), as.integer(edge2),
        as.integer(length(edge)), as.double(weight), as.integer(m), as.integer(m2), as.integer(nTips),  
        as.integer(pc) )
}   


random.addition <- function(data, method="fitch") 
{
    label <- names(data)
    nTips <- as.integer(length(label))
    remaining <- as.integer(sample(nTips))  
    tree <- structure(list(edge = structure(c(rep(nTips+1L, 3), remaining[1:3]), .Dim = c(3L, 2L)), 
    tip.label = label, Nnode = 1L), .Names = c("edge", "tip.label", "Nnode"), class = "phylo", order = "postorder")
    remaining <- remaining[-c(1:3)]
    
    if(nTips==3L) return(tree)
 
    nr <- attr(data, "nr")
    storage.mode(nr) <- "integer"
    n <- length(data) #- 1L
     
    data <- subset(data,,order(attr(data, "weight"), decreasing=TRUE))   
    data <- prepareDataFitch(data) 
    weight <- attr(data, "weight")

    m = nr*(2L*nTips - 2L)

    on.exit(.C("fitch_free"))
    .C("fitch_init", as.integer(data), as.integer(nTips*nr), as.integer(m) )
 
    storage.mode(weight) <- "double"

    for (i in remaining) {   

        dat = fnodesNew(tree, weight, nr, TRUE) 
        score = (dat[[3]] + dat[[4]]) 
        edge = tree$edge[,2]

        res <- .Call("FITCHTRIP", data[,i], as.integer(nr), as.integer(edge), as.double(weight))    

        score = score[edge] + res              

        res = min(score) 
        nt = which.min(score)  
        tree = addOne(tree, i, nt) 
        }
    attr(tree, "pscore") = res
    tree 
}


# weight, nr auch in C ??
fast.fitch <- function (tree, weight, nr, ps = TRUE) 
{
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    m = max(tree$edge) 
    .Call("FITCH345", as.integer(nr), as.integer(node), as.integer(edge), as.integer(length(edge)), 
        as.double(weight), as.integer(m), as.integer(ps), PACKAGE = "phangorn")
}


fitch.spr <- function(tree, data){
    nTips = as.integer(length(tree$tip))
    nr = attr(data, "nr")
    weight <- attr(data, "weight")
    minp = fast.fitch(tree, weight, nr, TRUE)

    for(i in 1:nTips){
        treetmp = dropTip(tree, i)          
        dat = fnodesNew(treetmp, weight, nr, TRUE)
        score = (dat[[3]] + dat[[4]])
 
        edge = treetmp$edge[,2] 

        res <- .Call("FITCHTRIP", data[,i], as.integer(nr), as.integer(edge), as.double(weight)) 
        score = score[edge] + res   
        if(min(score)<minp){
            nt = which.min(score)
            tree = addOne(treetmp, i, nt) 
            minp=min(score)
            print(minp)
        }
    }
    tree
}


           
indexNNI2 <- function(tree){
    parent = tree$edge[, 1]
    child = tree$edge[, 2]
 
    ind = which(child %in% parent)
    Nnode = tree$Nnode
    edgeMatrix = matrix(0L, 6, length(ind))

    pvector <- numeric(max(parent))
    pvector[child] <- parent
    cvector <- allChildren(tree)  

    k=0
    for(i in ind){        
            p1 = parent[i]          
            p2 = child[i]
            e34 = cvector[[p2]]
            ind1 = cvector[[p1]]
            e12 = ind1[ind1 != p2]
            if(pvector[p1]) edgeMatrix[, k+1] = c(p1,e12, e34, p2, 1L) 
            else edgeMatrix[, k+1] = c(e12, e34, p2, 0L)
            k=k+1
    } 
    cbind(edgeMatrix[c(1,3,2,4,5,6),], edgeMatrix[c(1,4,2,3,5,6),])
}
       


fitch.nni <- function (tree, data, ...) 
{
    nTips = as.integer(length(tree$tip))
    INDEX <- indexNNI2(tree)    
    nr = attr(data, "nr")
    weight <- attr(data, "weight")

    tmp = fnodesNew(tree, weight, nr, TRUE)
    p0 = tmp[[1]]
    m <- dim(INDEX)[2]
    pscore <- .C("fitchQuartet", as.integer(INDEX), as.integer(m), as.integer(nr), as.double(tmp[[3]]), as.double(tmp[[4]]),
          as.double(weight), double(m))[[7]]    

    swap <- 0
    candidates <- pscore < p0

    while (any(candidates)) {

        ind = which.min(pscore)
        pscore[ind] = Inf
        tree2 <- changeEdge(tree, INDEX[c(2,3), ind])

        test <- fast.fitch(tree2, weight, nr)

        if (test >= p0) 
            candidates[ind] = FALSE

        if (test < p0) {
            p0 <- test
            swap = swap + 1
            tree <- tree2
            indi <- which(INDEX[5,] %in% INDEX[1:5, ind])

            candidates[indi] <- FALSE
            pscore[indi] <- Inf
        }
    }
    list(tree = tree, pscore = p0, swap = swap)
}


optim.fitch <- function(tree, data, trace=1, rearrangements = "SPR", ...) {
    if(class(tree)!="phylo") stop("tree must be of class phylo") 
    if(!is.binary.tree(tree)) tree <- multi2di(tree)
    if(is.rooted(tree))tree <- unroot(tree)
    if(is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") tree <- reorder(tree, "postorder")  
    if (class(data)[1] != "phyDat") stop("data must be of class phyDat")

    rt = FALSE
    
    nTips = as.integer(length(tree$tip))
    nr = attr(data, "nr")
   
    data <- subset(data,tree$tip,order(attr(data, "weight"), decreasing=TRUE))   
    dat <- prepareDataFitch(data) 
    weight <- attr(data, "weight")

    m = nr*(2L*nTips - 2L)
    on.exit(.C("fitch_free"))
    .C("fitch_init", as.integer(dat), as.integer(nTips*nr), as.integer(m) )

    tree$edge.length=NULL
    swap = 0
    iter = TRUE
    pscore <- fast.fitch(tree, weight, nr)
    
    while (iter) {
        res <- fitch.nni(tree, dat, ...)
        tree <- res$tree
        if(trace>1)cat("optimize topology: ", pscore , "-->", res$pscore, 
            "\n")
        pscore = res$pscore
        swap = swap + res$swap
        if (res$swap == 0){
            if(rearrangements=="SPR"){
                tree <- fitch.spr(tree, dat)
                psc <- fast.fitch(tree, weight, nr)
                if(trace>1)cat("optimize topology (SPR): ", pscore , "-->", psc, "\n")
                if(pscore < psc+1e-6) iter=FALSE
                pscore <- psc
            } 
            else iter = FALSE
        }
    }
    if(trace>0)cat("Final p-score",pscore,"after ",swap, "nni operations \n") 
    if(rt)tree <- ptree(tree, data)  
    attr(tree, "pscore") = pscore
    tree
}



getOrder <- function (x) 
{
    label = names(x)
    dm = as.matrix(dist.hamming(x, FALSE))
    ind = as.vector(which(dm == max(dm), arr.ind = TRUE)[1, ])
    nTips = length(label)
    added = ind
    remaining <- c(1:nTips)[-ind]


    tree <- structure(list(edge = structure(c(rep(nTips+1L, 3), c(ind, 0L)), .Dim = c(3L, 2L)), 
        tip.label = label, Nnode = 1L), .Names = c("edge", "tip.label", "Nnode"), class = "phylo", order = "postorder")      

    l = length(remaining)
    res = numeric(l)

    nr <- attr(x, "nr")
    storage.mode(nr) <- "integer"
    n <- length(x) #- 1L
      
    data <- prepareDataFitch(x) 
    weight <- attr(data, "weight")
    storage.mode(weight) <- "double"

    m = nr*(2L*nTips - 2L)

    on.exit(.C("fitch_free"))
    .C("fitch_init", as.integer(data), as.integer(nTips*nr), as.integer(m) )

    for(i in 1:length(remaining)){
        tree$edge[3,2]= remaining[i]
        res[i] = fast.fitch(tree, weight, nr) 
    }
    tmp = which.max(res)
    added = c(added, remaining[tmp])
    remaining <- remaining[-tmp]
    tree$edge[,2]= added

    for (i in 4:(nTips - 1L)) {   

        dat = fnodesNew(tree, weight, nr, TRUE) 
        score0 = (dat[[3]] + dat[[4]]) 
        edge = tree$edge[,2]
        l = length(remaining)
        res = numeric(l)
        nt = numeric(l)
        k = length(added)+1L
        for(j in 1:l){
            psc <- .Call("FITCHTRIP", data[,remaining[j]], as.integer(nr), as.integer(edge), as.double(weight)) 
            score = score0[edge] + psc
            res[j] = min(score) 
            nt[j] = which.min(score)
        }
        tmp = which.max(res)
        added = c(added, remaining[tmp])        
        tree = addOne(tree, remaining[tmp], nt[tmp])
        remaining <- remaining[-tmp]  
    }
    added = c(added, remaining) 
    added 
}


addOne2 <- function (edge, tip, i){
    parent = edge[,1]
    l = dim(edge)[1]
    m = max(edge)+1L 
    p = edge[i,1]
    k = edge[i,2] 
    edge[i, 2] = m
    ind = match(p, parent)
    if(ind==1) edge = rbind(matrix(c(m,m,k,tip), 2, 2), edge)
    else edge = rbind(edge[1:(ind-1), ], matrix(c(m,m,k,tip), 2, 2), edge[ind:l, ])  
    edge
}   


bab <- function (data, tree = NULL, trace = 1, ...) 
{
    o = order(attr(data, "weight"), decreasing = TRUE)
    data = subset(data, , o)
    nr <- attr(data, "nr")
    pis <- parsinfo(data)
    p0 <- sum(attr(data, "weight")[pis[, 1]] * pis[, 2])
    if (length(pis) > 0) 
        data <- getRows(data, c(1:nr)[-pis[, 1]], TRUE)
    tree <- pratchet(data, start = tree, trace = trace - 1, ...)
    data <- subset(data, tree$tip.label) 
    nr <- as.integer(attr(data, "nr"))
    inord <- getOrder(data)
    lb <- lowerBound(data)
    nTips <- m <- length(data)
    
    nr <- as.integer(attr(data, "nr"))
    TMP <- matrix(0, m, nr)
    for (i in 4:m) {
        TMP[i, ] = lowerBound(subset(data, inord[1:i]))
    }

    weight <- as.double(attr(data, "weight"))
    data <- prepareDataFitch(data)
    m = nr*(2L*nTips - 2L)
    on.exit(.C("fitch_free"))
    .C("fitch_init", as.integer(data), as.integer(nTips*nr), as.integer(m) )

    mmsAmb = 0
    mmsAmb = TMP %*% weight  
    mmsAmb = mmsAmb[nTips] - mmsAmb
    mms0 = 0 
    mms0 = mms0 + mmsAmb

    minPars = mms0[1]
    kPars = 0

    if (trace) 
        print(paste("lower bound:", p0 + mms0[1]))
    bound <- fast.fitch(tree, weight, nr)
    if (trace) 
        print(paste("upper bound:", bound + p0))

    startTree <- structure(list(edge = structure(c(rep(nTips+1L, 3), as.integer(inord)[1:3]), .Dim = c(3L, 2L)), 
        tip.label = tree$tip.label, Nnode = 1L), .Names = c("edge", "tip.label", "Nnode"), class = "phylo", order = "postorder")

    trees <- vector("list", nTips)
    trees[[3]] <- list(startTree)
    pscores <- vector("list", nTips)
    pscores[[3]] <- fast.fitch(startTree, weight, nr)
    result <- list()

    k = 4L
    Nnode = 1L

    nTrees <- sapply(trees, length)
    result <- list() 
    while (any(nTrees) > 0) {

        a = max(which(nTrees>0))
        b = which.min(pscores[[a]])
 
        tmpTree <- trees[[a]][[b]]
        trees[[a]][[b]] <- NULL
        pscores[[a]] <- pscores[[a]][-b]

        dat = fnodesNew(tmpTree, weight, nr, TRUE) 
        score = (dat[[3]] + dat[[4]]) 
        edge = tmpTree$edge[,2]
        res <- .Call("FITCHTRIP", data[,inord[a+1L]], as.integer(nr), as.integer(edge), as.double(weight))    
        score = score[edge] + res + mms0[a+1L]              
        ms = min(score)
        if(ms<=bound){
            if((a+1L)<nTips){
                ind = which(score<=bound)
                trees[[a+1]] <- vector("list", length(ind))
                for(i in 1:length(ind))trees[[a+1]][[i]] <- addOne(tmpTree, inord[a+1L], ind[i])
                pscores[[a+1]] <- score[ind]  
            }
            else{
                ind = which(score==ms) 
                tmp <- vector("list", length(ind))
                for(i in 1:length(ind))tmp[[i]] <- addOne(tmpTree, inord[a+1L], ind[i])  
                if(ms < bound){
                     bound = ms
                     if(trace)print(bound) 
                     result = tmp    
                     for(i in 4:(nTips-1)){
                         dind = which(pscores[[i]] < bound)    
                         if(length(dind)>0){
                              pscores[[i]] <- pscores[[i]][-dind] 
                              trees[[i]][dind] <- NULL
                         } 
                     } 
                }
                else result = c(result, tmp)  
            }
        }    
        nTrees <- sapply(trees, length)
    }
    class(result) <- "multiPhylo"
    return(result)
}



.packageName <- "phangorn"


.onLoad  <- function(libname, pkgname) {
    library.dynam("phangorn", pkgname, libname)
}


#
# some general functions
#

# return vector of positions
dec2Bin = function (x) 
{
    res = NULL
    i = 1L
    while (x > 0) {
        if (x%%2L) 
            res = c(res, i)
            x = x%/%2L
            i = i + 1L
        }
    res
}

# returns binary (0, 1) vector of length k
dec2bin <- function (x, k=ceiling(log2(x))) 
{
    i = 1L
    res = integer(k)
    while (x > 0) {
        if (x%%2L) 
            res[i] = 1L
        x = x%/%2L
        i = i + 1L
    }
    res
}

# double factorial: log version
"ldfactorial" <- function(x){
    x = (x+1)/2
    res = lgamma(2*x)-(lgamma(x)+(x-1)*log(2))
    res
}

# double factorial
"dfactorial" <- function(x){exp(ldfactorial(x))}


#
# Hadamard Conjugation
#

hadamard <- function(x){
    res=1
    while(x>0){
        res=rbind(cbind(res,res),cbind(res,-res))
        x=x-1
    }
    res
}


fhm <- function(v){
    n = length(v)
    n = log2(n)
    res = .C("fhm", v = as.double(v), n = as.integer(n), PACKAGE = "phangorn")$v
    res
}


seq2split = function(s){
    n=length(s)
    res= fhm(log(fhm(s)))/n
    res
}


split2seq = function(q){
    n=length(q)
    res= fhm(exp(fhm(q)))/n
    res
}


distanceHadamard <- function (dm, eps = 0.001) 
{
    if (class(dm) == "dist") {
        n <- attr(dm, "Size")
        Labels = attr(dm, "Labels")
    }
    if (class(dm) == "matrix") {
        n <- dim(dm)[1]
        Labels <- colnames(dm)
        dm <- dm[lower.tri(dm)]
    }
    ns <- 2^(n - 1)
    if (n > 23) 
        stop("Hadamard conjugation works only efficient for n < 24")
    result <- .Call("dist2spectra", dm, as.integer(n), as.integer(ns), 
        PACKAGE = "phangorn")
    weights = -fhm(result)/2^(n - 2)    
    
    if(eps>0){
        weights = weights[-1]
        ind2 = which(weights>eps)
        n2 = length(ind2)
        splits = vector("list", n2)
        for(i in 1:n2)splits[[i]] = dec2Bin(ind2[i])
        attr(splits, "weights") = weights[ind2]
        attr(splits, "labels") = Labels
        attr(splits, 'dm') = dm
        class(splits)='splits'
        return(splits)      
    }  
    res <- data.frame(distance = result, edges = weights, index = 0:(ns - 1))
    attr(res, "Labels") <- Labels
    res
}


# phyDat, DNAbin support
h4st = function(obj, levels=c('a','c','g','t')){
    if (is.matrix(obj)) 
        obj = as.data.frame(t(obj))
    if (class(obj) == "phyDat") 
        obj = as.data.frame(t(as.character(obj)))	
#    if(is.matrix(obj)) obj = as.data.frame(t(obj))
#    DNA = as.data.frame(obj)
#    DNA = t(as.character(obj))

    n = dim(obj)[1]
    p = dim(obj)[2]

    if(p>11) stop("4-state Hadamard conjugation works only efficient for n < 12")

    DNAX = matrix(0,n,p)
    DNAY = matrix(0,n,p)

    DNAX[obj==levels[1]]=0
    DNAX[obj==levels[2]]=1
    DNAX[obj==levels[3]]=1
    DNAX[obj==levels[4]]=0

    DNAY[obj==levels[1]]=0
    DNAY[obj==levels[2]]=1
    DNAY[obj==levels[3]]=0
    DNAY[obj==levels[4]]=1

    DNAY = DNAY - DNAY[,p]
    DNAX = DNAX - DNAX[,p]

    DNAY = abs(DNAY[,-p])
    DNAX = abs(DNAX[,-p])
    dy = DNAY %*% (2^(0:(p-2))) 
    dx = DNAX %*% (2^(0:(p-2))) 

    INDEX =  dx + 2^(p-1) * dy
    blub = table(INDEX)
    index = as.numeric(rownames(blub)) + 1
    sv = numeric(4^(p-1))
    sv[index] = blub
    qv = matrix(seq2split(sv),2^(p-1),2^(p-1))
    sv = matrix(sv,2^(p-1),2^(p-1))

#    q = cbind(transversion = qv[-1,1], transition.1 = diag(qv)[-1], transition.2 = qv[1,-1])

    transversion <- transition.1 <- transition.2 <- allSplits(p, colnames(obj)) 
    attr(transversion,"weights") = qv[-1,1]
    attr(transition.1,"weights") = diag(qv)[-1]
    attr(transition.2,"weights") = qv[1,-1]

#    result = list(q = q, qv = qv, sv=sv, n=sum(sv), names=names(obj))
    result = list(transversion = transversion, transition.1=transition.1, transition.2 = transition.2, 
        qv = qv, sv=sv, n=sum(sv), names=names(obj))
    result
}


h2st <- function (obj, eps=0.001) 
{
    if (class(obj) != "phyDat") stop("Error") 
    if (attr(obj,"nc") != 2)stop("Error")
    nr = attr(obj, "nr") #n
    p = length(obj) #p
    weight = attr(obj, "weight")
    if (p > 23) 
        stop("Hadamard conjugation works only efficient for n < 24")
    DNAX = matrix(0, nr, p-1)
    for(i in 1:(p-1)) DNAX[,i] = obj[[i]]-1
    DNAX[obj[[p]]==2,] = 1 - DNAX[obj[[p]]==2,]

    index = DNAX %*% (2^(0:(p - 2))) + 1
    sv = numeric(2^(p - 1))
    for(i in 1:nr)sv[index[i]] = sv[index[i]]+ weight[i]
    qv = seq2split(sv)
    
    if(eps>0){
	    qv = qv[-1]
        ind2 = which(qv>eps)
        indT= c(2L^(0:(p-2)), 2L^(p-1)-1) 
        ind2 = union(ind2, indT)
        n2 = length(ind2)
        splits = vector("list", n2)
        for(i in 1:n2)splits[[i]] = dec2Bin(ind2[i])
        attr(splits, "weights") = qv[ind2]
        attr(splits, "labels") = names(obj)
        class(splits)='splits'
        return(splits)    
        }
    result = data.frame(edges = qv, splits = sv, index = 0:(2^(p - 
        1) - 1))
    attr(result, "Labels") = names(obj)
    result
}


#
# tree distance functions
#
treedist <- function (tree1, tree2) 
{
    tree1 = unroot(tree1)
    tree2 = unroot(tree2)
    tree1 = reorderPruning(tree1)
    tree2 = reorderPruning(tree2)
    symmetric.difference = NULL
    branch.score.difference = NULL
    path.difference = NULL
    quadratic.path.difference = NULL
    if(!is.binary.tree(tree1) | !is.binary.tree(tree2))warning("Trees are not binary!")

    o1 = order(tree1$tip.label)
    o2 = order(tree2$tip.label)
    ll = length(o1)
    p1 = bipartition(tree1)
    p2 = bipartition(tree2)
    p = dim(p1)[1]
    pa = dim(p1)[1]
    pb = dim(p2)[1]
    M1 = p1[, o1]
    M2 = p2[, o2]
    if (!is.null(tree1$edge.length) & !is.null(tree2$edge.length)) {
        v1 = tree1$edge.length
        v2 = tree2$edge.length

        dv1 = crossprod(M1 * v1, 1-M1) 
        dv1 = dv1 + t(dv1)
        dv2 = crossprod(M2 * v2, 1-M2) 
        dv2 = dv2 + t(dv2)
        ind= lower.tri(dv1)
        quadratic.path.difference = sqrt(sum((dv1[ind] - dv2[ind])^2))
    }
    R = M1 %*% t(M2) + (1 - M1) %*% t(1 - M2)
    R = (R%%ll == 0)
    r1 = rowSums(R) > 0
    r2 = colSums(R) > 0
    symmetric.difference = 2 * (p - sum(r1))
    if (!is.null(tree1$edge.length) & !is.null(tree2$edge.length)) {
        v1 = tree1$edge.length
        v2 = tree2$edge.length
        ind1 <- (1:p)[r1]
        ind2 <- unlist(apply(R, 1, which, TRUE))
        s1 = sum((v2[ind2] - v1[ind1])^2)
        zaehler = abs(v2[ind2] - v1[ind1])
        nenner = (v2[ind2] + v1[ind1])/2
        difference = matrix(0, sum(r1), 4)
        difference[, 1] = zaehler
        difference[, 2] = nenner
        difference[, 3] = ind1
        difference[, 4] = ind2
        s2 = sum((v1[(1:p)[!r1]])^2)
        s3 = sum((v2[(1:p)[!r2]])^2)
        branch.score.difference = sqrt(s1 + s2 + s3)
    }
    M1[M1 == 0] = -1
    M2[M2 == 0] = -1
    dt1 = (pa - t(M1) %*% M1)/2
    dt2 = (pb - t(M2) %*% M2)/2
    ind = lower.tri(dt1)
    path.difference = sqrt(sum((dt1[ind] - dt2[ind])^2))
    result = c(symmetric.difference = symmetric.difference, 
        branch.score.difference = branch.score.difference, path.difference = path.difference, 
        quadratic.path.difference = quadratic.path.difference)
    result              
}


RF.dist <- function (tree1, tree2, check.labels = TRUE)
{
   r1 = is.rooted(tree1)
   r2 = is.rooted(tree2)
   if(r1 != r2){
       warning("one tree is unrooted, unrooted both")
   }
   if (check.labels) {
       ind <- match(tree1$tip.label, tree2$tip.label)
       if (any(is.na(ind)) | length(tree1$tip.label) !=
       length(tree1$tip.label))
           stop("trees have different labels")
       tree2$tip.label <- tree2$tip.label[ind]
#       tree2$edge[match(ind, tree2$edge[, 2]), 2] <- 1:length(ind)
       ind2 <- match(1:length(ind), tree2$edge[, 2])
       tree2$edge[ind2, 2] <- order(ind)
   }

   if(!r1 | !r2){
       if(r1) tree1 = unroot(tree1)
       if(r2) tree2 = unroot(tree2)
       ref1 <- Ancestors(tree1, 1, "parent")
       tree1 <- reroot(tree1, ref1)
       ref2 <- Ancestors(tree2, 1, "parent")
       tree2 <- reroot(tree2, ref2)
   }
   if(!is.binary.tree(tree1) | !is.binary.tree(tree2))warning("Trees are not binary!")
   p1 = bipart(tree1)
   p2 = bipart(tree2)
   ind <- sum(p1 %in% p2)
   l = length(tree1$tip)
   l = l - 2 + is.rooted(tree1)
   RF = 2 * (l-ind)
   RF
}

#
# UPGMA, NJ and UNJ
#
"upgma" <- function(D,method="average",...){
    DD=as.dist(D)
    hc = hclust(DD,method=method,...)
    result = as.phylo(hc)
    result = reorderPruning(result)
    result
}


"wpgma" <- function(D,method="mcquitty",...){
    DD=as.dist(D)
    hc = hclust(DD,method=method,...)
    result = as.phylo(hc)
    result = reorderPruning(result)
    result
}


NJ <- function (x) 
{
    x = as.matrix(x)
    labels <- attr(x, "Labels")[[1]]
    edge.length = NULL
    edge = NULL
    d = as.matrix(x)
    if (is.null(labels)) 
        labels = colnames(d)
    l = dim(d)[1]
    m = l - 2
    nam = 1L:l
    k = 2L * l - 2L
    while (l > 2) {
        r = rowSums(d)/(l - 2)
        i = 0
        j = 0
        tmp <- .C("out", as.double(d), as.double(r), as.integer(l), 
            as.integer(i), as.integer(j), PACKAGE = "phangorn")
        e2 = tmp[[5]]
        e1 = tmp[[4]]
        l1 = d[e1, e2]/2 + (r[e1] - r[e2])/(2)
        l2 = d[e1, e2] - l1
        edge.length = c(l1, l2, edge.length)
        edge = rbind(c(k, nam[e2]), edge)
        edge = rbind(c(k, nam[e1]), edge)
        nam = c(nam[c(-e1, -e2)], k)
        dnew = (d[e1, ] + d[e2, ] - d[e1, e2])/2
        d = cbind(d, dnew)
        d = rbind(d, c(dnew, 0))
        d = d[-c(e1, e2), -c(e1, e2)]
        k = k - 1L
        l = l - 1L
    }
    edge.length = c(d[2, 1], edge.length)
    attr(edge.length,"names") = NULL
    result = list(edge = rbind(c(nam[2], nam[1]), edge), edge.length = edge.length,
     tip.label = labels, Nnode = m)
    class(result) <- "phylo" 
    reorderPruning(result)
}


UNJ = function (x) 
{
    x = as.matrix(x)
    labels <- attr(x, "Labels")[[1]]
    edge.length = NULL
    edge = NULL
    d = as.matrix(x)
    if (is.null(labels)) 
        labels = colnames(d)
    l = dim(d)[1]
    n = l
    nam = as.character(1:l)
    m=l-2
    nam = 1:l
    k = 2*l-2       
    w = rep(1,l)
    while (l > 2) {
        r = rowSums(d)/(l - 2)
        i = 0
        j = 0
        tmp <- .C("out", as.double(d), as.double(r), as.integer(l), as.integer(i), as.integer(j), PACKAGE = "phangorn")
        e2 = tmp[[5]]
        e1 = tmp[[4]]
        l1 = d[e1, e2]/2 + sum((d[e1,-c(e1,e2)] - d[e2,-c(e1,e2)])*w[-c(e1,e2)])/(2*(n-w[e1]-w[e2]))
        l2 = d[e1, e2]/2 + sum((d[e2,-c(e1,e2)] - d[e1,-c(e1,e2)])*w[-c(e1,e2)])/(2*(n-w[e1]-w[e2]))
        edge.length = c(l1, l2, edge.length)
        edge = rbind(c(k, nam[e2]), edge)
        edge = rbind(c(k, nam[e1]), edge)
        nam = c(nam[c(-e1, -e2)], k)
      
        dnew = (w[e1]*d[e1, ] + w[e2]*d[e2, ] - w[e1]*l1 - w[e2]*l2)/(w[e1] + w[e2])
        d = cbind(d, dnew)
        d = rbind(d, c(dnew, 0))
        d = d[-c(e1, e2), -c(e1, e2)]
        w = c(w, w[e1] + w[e2])
        w = w[-c(e1, e2)]
        k = k - 1
        l = l - 1
    }
    edge.length=c(d[2,1],edge.length)
    result = list(edge = rbind(c(nam[2], nam[1]), edge), 
    edge.length=edge.length, tip.label = labels, Nnode=m)
    class(result) <- "phylo"
    reorderPruning(result)  
}


PNJ <- function (data) 
{
    q <- l <- r <- length(data)
    weight <- attr(data,"weight")
        
    height = NULL    
    parentNodes <- NULL
    childNodes <- NULL    
    nam <- names(data)
    tip.label <- nam
    edge = 1:q
    
    z = 0
    D = matrix(0, q, q)
    
    for (i in 1:(l - 1)) {
        for (j in (i + 1):l) {
            w = (data[[i]] * data[[j]]) %*% c(1, 1, 1, 1)
            D[i, j] = sum(weight[w==0])
        }
    }

    while (l > 1) {
        l = l - 1
        z = z + 1
        d = D + t(D)
        if(l>1) r = rowSums(d)/(l-1)
        if(l==1) r = rowSums(d)
        M = d - outer(r,r,"+")
        diag(M) = Inf

        e=which.min(M)
        e0=e%%length(r)
        e1 = ifelse(e0==0, length(r), e0)
        e2= ifelse(e0==0, e%/%length(r), e%/%length(r) + 1)
        
        ind = c(e1,e2)       
        len = d[e]/2
        nam = c(nam[-ind], as.character(-l))
           
        parentNodes = c(parentNodes,-l,-l)            
        childNodes = c(childNodes,edge[e1],edge[e2])        
        
        height = c(height, len, len)
        edge = c(edge[-ind], -l)
        w = (data[[e1]] * data[[e2]]) %*% c(1, 1, 1, 1)
        w = which(w == 0)
        newDat = data[[e1]] * data[[e2]]
        newDat[w, ] = data[[e1]][w, ] + data[[e2]][w, ]   
        data = data[-c(e1,e2)]
        data[[l]] = newDat 
        if (l > 1) {
            D = as.matrix(D[, -ind])
            D = D[-ind, ]
            dv = numeric(l - 1)
            for (i in 1:(l - 1)) {
                w = (data[[i]] * data[[l]]) %*% c(1, 1, 1, 1)
                dv[i] = sum(weight[w==0])
            }
            D = cbind(D, dv)
            D = rbind(D, 0)
        }
    }
    tree <- list(edge = cbind(as.character(parentNodes),as.character(childNodes)),tip.label=tip.label) 
    class(tree) <- "phylo"
    tree <- old2new.phylo(tree)   
    reorderPruning(tree)    
}


dist.hamming <- function (x, ratio = TRUE) 
{
    if (class(x) != "phyDat") 
        stop("x has to be element of class phyDat")
    l = length(x)
    weight <- attr(x, "weight")  
    nc <- attr(x, "nc")
    d = numeric((l * (l - 1))/2)
    if(nc > 31){
        contrast <- attr(x, "contrast")
        k = 1
        for (i in 1:(l - 1)) {
            X = contrast[x[[i]], , drop = FALSE]
            for (j in (i + 1):l) {
                d[k] = sum(weight * (rowSums(X * contrast[x[[j]], , drop = FALSE]) == 0))
                k = k + 1
            }
        }
    } # end if  
    else{
        nr <- attr(x, "nr")
        x <- prepareDataFitch(x) 
        res <- .C("distHamming", as.integer(x), as.double(weight), as.integer(nr), as.integer(l), as.double(d))
        d <- res[[5]]
    }     

    if (ratio) 
        d = d/sum(weight)
    attr(d, "Size") <- l
    if (is.list(x)) 
        attr(d, "Labels") <- names(x)
    else attr(d, "Labels") <- colnames(x)
    attr(d, "Diag") <- FALSE
    attr(d, "Upper") <- FALSE
    attr(d, "call") <- match.call()
    attr(d, "method") <- "hamming"
    class(d) <- "dist"
    return(d)
}


dist.tip <- function (x1, x2, el, eig = edQt(), bf) 
{
    d = crossprod(x1, x2)
    el = 1 - (sum(diag(d))/sum(d))
    m = dim(d)[1]
    if(el> (m-1) / m){ 
        el <- 10
        ll1 <- sum((d * log(getP(el, eig)[[1]] * bf))[d>0])
        return(c(10, 10, ll1, ll1))   
        } 
    else el = min(10, -((m-1)/m)*log(1-(m/(m-1))*el))   
    eps = 1  
    if (el == 0) {
        ll1 <- sum((d * log(getP(el, eig)[[1]] * bf))[d>0])
        return(c(0, 0, ll1, ll1))
    }
    ll0 <- sum((d * log(getP(el, eig)[[1]] * bf))[d>0])
    while (eps > 1e-06) {
        f = getP(el, eig)[[1]] * bf
        df = getdP(el, eig)[[1]] * bf
        d2f = getd2P(el, eig)[[1]] * bf
        dl = sum(d * (df/f))
        d2l = sum(d * ((f * d2f) - df^2)/(f^2))
        el1 = log(el) - dl/d2l
        el1 = exp(el1)
        if(el1 > 100) el1 = 100
        eps = abs(el1 - el)
        el = el1
    }
    df = getdP2(el, eig)[[1]] * bf
    d2f = getd2P2(el, eig)[[1]] * bf
    dl = sum(d * (df/f))
    d2l = sum(d * ((f * d2f) - df^2)/(f^2))
    ll1 <- sum((d * log(getP(el, eig)[[1]] * bf))[d>0])
    c(el, 1/-d2l, ll0, ll1)
}

# ohne 2. Ableitung
dist.tip2 = function (x1, x2, el, eig = edQt(), bf) 
{
    d = crossprod(x1, x2)
    el = 1 - (sum(diag(d))/sum(d))  
    ll0 <- sum((d * log(getP(el, eig)[[1]] * bf))[d>0])
    eps = 1
    if (el == 0) {
        return(c(0, 0, ll0, ll0))
    }
    while (eps > 1e-06) {
        f = getP(el, eig)[[1]] * bf
        df = (getdP(el, eig)[[1]] * bf) / f
        dl = sum(d * df)
        d2l = sum(d * df^2)
        el1 = log(el) + dl/d2l
        el1 = exp(el1)
        if(el1 > 100) el1 = 100
        eps = abs(el1 - el)
        el = el1
    }
    df = getdP2(el, eig)[[1]] * bf
    d2f = getd2P2(el, eig)[[1]] * bf
    dl = sum(d * (df/f))
    d2l = sum(d * ((f * d2f) - df^2)/(f^2))
    ll1 <- sum((d * log(getP(el, eig)[[1]] * bf))[d>0])
    c(el, 1/-d2l, ll0, ll1)
}



dist.ml <- function(x, model="JC69", pairwise.deletion = FALSE, bf=NULL, Q=NULL, ...){
    if (class(x) != "phyDat") 
        stop("x has to be element of class phyDat")
    l = length(x)
    contrast <- attr(x, 'contrast')
    nc <- attr(x, "nc")
    w=1
    g=1
    d = numeric((l * (l - 1))/2)
    v = numeric((l * (l - 1))/2) 
    type <- attr(x,"type")
    con = rowSums(contrast>0)<2

    if(!pairwise.deletion){
        index = con[x[[1]]]
        for(i in 2:l) index = index & con[x[[i]]]
        index = which(index)
        x = subset(x,,index)
    }
    
    weight <- attr(x, "weight")
    ll.0 <- numeric(length(weight))

    model <- match.arg(model, c("JC69", "WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24"))
    if(!is.na(match(model, c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24"))))
            getModelAA(model, bf=is.null(bf), Q=is.null(Q))   

    if(type == "DNA" & model == "JC69"){
        if(is.null(bf))bf = c(.25,.25,.25,.25)
        if(is.null(Q))Q = rep(1,6)
    }
    if(model == "JC69"){
        if(is.null(bf))bf = rep(1, nc) / nc
        if(is.null(Q))Q = rep(1,(nc-1)*nc/2L)
    }
    eig <- edQt(Q=Q, bf=bf)
    old.el <- 0.1  # start with dist.hamming ??
    k = 1
    if(!pairwise.deletion){
        for (i in 1:(l - 1)) {
            Xi = weight * contrast[x[[i]], , drop=FALSE] 
            for (j in (i + 1):l) {
                tmp <- dist.tip(Xi, contrast[x[[j]], , drop=FALSE], old.el, eig=eig, bf)
                d[k] <- tmp[1] 
                v[k] <- tmp[2]  
                k = k + 1
            }
        }
    }
    if(pairwise.deletion){
        for (i in 1:(l - 1)) {
            Xi = weight * contrast[x[[i]], , drop=FALSE] 
            ind = con[x[[i]]]
            for (j in (i + 1):l) {
                index = which(ind & con[x[[j]]]) 
                tmp <- dist.tip(Xi[index, , drop=FALSE], contrast[x[[j]], , drop=FALSE][index, , drop=FALSE], old.el, eig=eig, bf)
                d[k] <- tmp[1] 
                v[k] <- tmp[2]  
                k = k + 1
            }
        }
    }
    attr(d, "Size") <- l
    if (is.list(x)) 
        attr(d, "Labels") <- names(x)
    else attr(d, "Labels") <- colnames(x)
    attr(d, "Diag") <- FALSE
    attr(d, "Upper") <- FALSE
    attr(d, "call") <- match.call()
    attr(d, "method") <- model
    attr(d, "var") <- v
    class(d) <- "dist"
    return(d)   
}


dist.logDet = function (x) 
{
    if (class(x) != "phyDat") 
        stop("x has to be element of class phyDat")
    weight <- attr(x, "weight")
    contrast <- attr(x, 'contrast')
    r <- attr(x, "nc")
    l = length(x)
    d = numeric((l * (l - 1))/2)
    k = 1
    for (i in 1:(l - 1)) {
        Xi = weight * contrast[x[[i]], , drop=FALSE]
        for (j in (i + 1):l) {
            tmp = crossprod(Xi, contrast[x[[j]], , drop=FALSE])
            class(tmp) = "matrix"
            z = determinant.matrix(tmp, logarithm=TRUE)  
            res = z$sign*z$modulus
            if (is.nan(res)) {
                d[k] = 10
            }
            else d[k] = (-res + sum(log(rowSums(tmp) * colSums(tmp)))/2)/r
            k = k + 1
        }
    }
    attr(d, "Size") <- l
    if (is.list(x)) 
        attr(d, "Labels") <- names(x)
    else attr(d, "Labels") <- colnames(x)
    attr(d, "Diag") <- FALSE
    attr(d, "Upper") <- FALSE
    attr(d, "call") <- match.call()
    attr(d, "method") <- "logDet"
    class(d) <- "dist"
    return(d)
}


#
# splits
#
splitsNetwork <- function(dm, gamma=.1, lambda=1e-6, weight=NULL){
    dm = as.matrix(dm)
    k = dim(dm)[1]
    X2 = designAll(k)
    X=X2[[1]]
    
    y = dm[lower.tri(dm)]
    ind = c(2^(0:(k-2)),2^(k-1)-1)
    
    y2 = lm(y~X[,ind]-1)$res
    n = dim(X)[2]

    ridge <- lambda * diag(n) 
    ridge[ind,ind] <- 0
    if(!is.null(weight)) Dmat <- crossprod(X * sqrt(weight)) + ridge
    else Dmat <- crossprod(X) + ridge
    if(!is.null(weight)) dvec <- crossprod(X * sqrt(weight),y * sqrt(weight))
    else dvec <- crossprod(X, y)

    
    ind1       <- rep(1,n)
    ind1[ind]  <- 0 

    Amat       <- cbind(ind1,diag(n)) 
    bvec       <- c(gamma, rep(0,n))

    solution <- quadprog::solve.QP(Dmat,dvec,Amat,bvec=bvec, meq=1)$sol   
    
    ind2 <- which(solution>1e-8)
    n2 <- length(ind2)

    ind3 = which(duplicated(c(ind2, ind), fromLast = TRUE)[1:n2])
    ridge2 <- lambda * diag(n2) 
    ridge2[ind3,ind3] <- 0
    
    if(!is.null(weight)) Dmat <- crossprod(X[, ind2] * sqrt(weight)) + ridge2
    else Dmat <- crossprod(X[, ind2]) + ridge2
    if(!is.null(weight)) dvec <- crossprod(X[, ind2] * sqrt(weight),y * sqrt(weight))
    else dvec <- crossprod(X[, ind2], y)
    
    Amat2 <- diag(n2)
    bvec2 <- rep(0, n2)
    solution2  <- quadprog::solve.QP(Dmat, dvec, Amat2)$sol
    
    RSS1 = sum((y-X[,ind2]%*%solution[ind2])^2)
    RSS2 = sum((y-X[,ind2]%*%solution2)^2)
    
    
    splits = vector("list", length(ind2))
    for(i in 1:length(ind2))splits[[i]] = which(X2[[2]][ind2[i],]==1)
    attr(splits, "weights") = solution[ind2]
    attr(splits, "unrestricted") = solution2
    attr(splits, "stats") = c(df=n2, RSS_p = RSS1, RSS_u=RSS2)
    attr(splits,"labels") =dimnames(dm)[[1]]
    class(splits)='splits'
    return(splits)           
}


allSplits = function(k, labels=NULL){
    result <- lapply(1:(2^(k-1)-1),dec2Bin)
    if(is.null(labels)) labels=(as.character(1:k))
    attr(result, 'labels') =labels
    class(result)='splits'
    result
    }   

    
splits2design <- function(obj, weight=NULL){
    labels= attr(obj,'labels')
    m = length(labels)
    n=length(obj)
    l = 1:m 
    sl = sapply(obj, length)
    p0 = sl * (m-sl)
    p = c(0,cumsum(p0))
    i = numeric()
    for(k in 1:n){
        sp = obj[[k]]
        if(p0[k]!=0) i[(p[k]+1):p[k+1]] = getIndex(sp, l[-sp], m) 
    }
    dims=c(m*(m-1)/2,n)
    sparseMatrix(i=i, p=p, dims=dims) 
    }


#
# Data structures for ML and MP
# 
fast.table <- function (data)                                                            
{                                                                                 
    if(!is.data.frame(data)) 
        data = as.data.frame(data, stringsAsFactors = FALSE)                    
    da = do.call("paste", c(data, sep = "\r"))                                             
    ind = !duplicated(da)                                                                  
    levels = da[ind]                                                                       
    cat <- factor(da,levels = levels)                                                      
    nl <- length(levels(cat))                                                        
    bin <- (as.integer(cat) - 1)                                                           
    pd <- nl                                                                               
    bin <- bin[!is.na(bin)]                                                                
    if (length(bin)) bin <- bin + 1                                                        
    y <- tabulate(bin, pd)                                                                 
    result=list(index = bin, weights = y, data = data[ind,])                                                                                  
    result                                                                                 
}                                                                                        


phyDat.default <- function (data, levels = NULL, return.index = TRUE, contrast = NULL, 
    ambiguity = "?", compress=TRUE, ...) 
{
    if (is.matrix(data)) 
        nam = row.names(data)
    else nam = names(data)
    if (class(data) == "DNAbin") 
        data = as.character(data)
    if (is.matrix(data)) 
        data = as.data.frame(t(data), stringsAsFactors = FALSE)
    else data = as.data.frame(data, stringsAsFactors = FALSE)
    if(compress){
        ddd = fast.table(data)
        data = ddd$data
        weight = ddd$weight
        index = ddd$index
    }
    else{
        p = length(data[[1]])
        weight = rep(1, p)
        index = 1:p
    }
    q = length(data)
    p = length(data[[1]])
    tmp <- vector("list", q)
    if (!is.null(contrast)) {
        levels = colnames(contrast)
        all.levels = rownames(contrast)
        rownames(contrast) = NULL
    }
    else {
        if (is.null(levels)) 
            stop("Either argument levels or contrast has to be supplied")
        l = length(levels)
        contrast = diag(l)
        all.levels = levels
        if (!is.null(ambiguity)) {
            all.levels = c(all.levels, ambiguity)
            k = length(ambiguity)
            if (k > 0) 
                contrast = rbind(contrast, matrix(1, k, l))
        }
    }
    d = dim(data)
    att = attributes(data) 
    data = match(unlist(data), all.levels)
    attr(data, "dim") = d
    data = as.data.frame(data, stringsAsFactors=FALSE)  
    attributes(data) = att

    row.names(data) = as.character(1:p)
    data = na.omit(data)
   
    aaa = match(index, attr(data, "na.action"))
    index = index[is.na(aaa)] 
    index = match(index, unique(index))
    rn = as.numeric(rownames(data))
    attr(data, "na.action") = NULL  
        
    weight = weight[rn] 
    p = dim(data)[1]
    names(data) = nam
    attr(data, "row.names") = NULL
    attr(data, "weight") = weight
    attr(data, "nr") = p
    attr(data, "nc") = length(levels)
    if (return.index) 
        attr(data, "index") = index
    attr(data, "levels") = levels
    attr(data, "allLevels") = all.levels
    attr(data, "type") = "USER"
    attr(data, "contrast") = contrast
    class(data) = "phyDat"
    data
}

 
phyDat.DNA = function (data, return.index = TRUE) 
{
    if (is.matrix(data)) 
        nam = row.names(data)
    else nam = names(data)
    if (class(data) == "DNAbin") 
        data = as.character(data)
    if (is.matrix(data)) 
        data = as.data.frame(t(data), stringsAsFactors = FALSE)
    else data = as.data.frame(data, stringsAsFactors = FALSE)

    data = data.frame(tolower(as.matrix(data)), stringsAsFactors = FALSE)
 
    ac = c("a", "c", "g", "t", "u", "m", "r", "w", "s", "y", 
        "k", "v", "h", "d", "b", "n", "?", "-")
    AC = matrix(c(c(1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 
        0, 1, 1, 1), c(0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 
        0, 1, 1, 1, 1), c(0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 
        0, 1, 1, 1, 1, 1), c(0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 
        0, 1, 1, 1, 1, 1, 1)), 18, 4, dimnames = list(NULL, c("a", 
        "c", "g", "t")))
    
    ddd = fast.table(data)
    data = ddd$data
    index = ddd$index
    q = length(data)
    p = length(data[[1]])
    d = dim(data)
    att = attributes(data) 
    data = match(unlist(data), ac)
    attr(data, "dim") = d
    data = as.data.frame(data, stringsAsFactors=FALSE)
    attributes(data) = att

    row.names(data) = as.character(1:p)
    data = na.omit(data)
    rn = as.numeric(rownames(data))

    aaa = match(index, attr(data, "na.action"))
    index = index[is.na(aaa)] 
    index = match(index, unique(index))
    rn = as.numeric(rownames(data))
    attr(data, "na.action") = NULL
    
    weight = ddd$weight[rn]
    p = dim(data)[1]
    names(data) = nam
    attr(data, "row.names") = NULL 
    attr(data, "weight") = weight
    attr(data, "nr") = p
    attr(data, "nc") = 4
    if (return.index) 
        attr(data, "index") = index
    attr(data, "levels") = c("a", "c", "g", "t")
    attr(data, "allLevels") = ac
    attr(data, "type") = "DNA"
    attr(data, "contrast") = AC
    class(data) = "phyDat"
    data
}


phyDat.AA <- function (data, return.index = TRUE) 
{
    if(is.matrix(data)) nam = row.names(data)
    else nam = names(data)  
    if (class(data) == "DNAbin") 
        data = as.character(data)
    if (is.matrix(data)) 
        data = as.data.frame(t(data), stringsAsFactors = FALSE)
    else data = as.data.frame(data, stringsAsFactors = FALSE)
  
    data = data.frame(tolower(as.matrix(data)), stringsAsFactors = FALSE)

    aa <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", 
        "l", "k", "m", "f", "p", "s", "t", "w", "y", "v")
    aa2 <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", 
        "l", "k", "m", "f", "p", "s", "t", "w", "y", "v", "b", 
        "z", "x", "-", "?")
    AA <- diag(20)
    AA <- rbind(AA, matrix(0, 5, 20))
    AA[21, 3] <- AA[21, 4] <- 1 # Aspartate or Asparagine
    AA[22, 6] <- AA[22, 7] <- 1 #
    AA[23:25, ] = 1
    dimnames(AA) <- list(aa2, aa)
    
    ddd = fast.table(data)
    data = ddd$data
    index = ddd$index
    q = length(data)
    p = length(data[[1]])
    tmp <- vector("list", q)

    d = dim(data)
    att = attributes(data) 
    data = match(unlist(data), aa2)
    attr(data, "dim") = d
    data = as.data.frame(data, stringsAsFactors=FALSE)
    attributes(data) = att

    row.names(data) = as.character(1:p)
    data = na.omit(data)
    rn = as.numeric(rownames(data))

    aaa = match(index, attr(data, "na.action"))
    index = index[is.na(aaa)] 
    index = match(index, unique(index))
    rn = as.numeric(rownames(data))
    attr(data, "na.action") = NULL
        
    weight = ddd$weight[rn]
    p = dim(data)[1]
    names(data) = nam
    attr(data, "row.names") = NULL
    attr(data, "weight") = weight
    attr(data, "nr") = p
    attr(data, "nc") = 20
    if (return.index) 
        attr(data, "index") = index
    attr(data, "levels") = aa
    attr(data, "allLevels") = aa2
    attr(data, "type") = "AA"
    attr(data, "contrast") = AA    
    class(data) = "phyDat"
    data
}



phyDat.codon <- function (data, return.index = TRUE) 
{
    if(is.matrix(data)) nam = row.names(data)
    else nam = names(data)  
    if (class(data) == "DNAbin") 
        data = as.character(data)

    if (is.matrix(data)) 
        data = as.data.frame(t(data), stringsAsFactors = FALSE)
    else data = as.data.frame(data, stringsAsFactors = FALSE)
    
    data = data.frame(tolower(as.matrix(data)), stringsAsFactors = FALSE)

    data[data=="u"] = "t" 

    splseq = function (seq, frame = 0) 
    {
        starts <- seq(from = frame + 1, to = length(seq), by = 3L)
        sapply(starts, function(x) paste(seq[x:(x + 2L)], collapse=""))
    } 
 
    data = sapply(data, splseq)
    
    ddd = fast.table(data)
    codon = c("aaa", "aac", "aag", "aat", "aca", "acc", "acg", "act", 
      "aga", "agc", "agg", "agt", "ata", "atc", "atg", "att", 
      "caa", "cac", "cag", "cat", "cca", "ccc", "ccg", "cct", "cga", 
      "cgc", "cgg", "cgt", "cta", "ctc", "ctg", "ctt", "gaa", "gac", 
      "gag", "gat", "gca", "gcc", "gcg", "gct", "gga", "ggc", "ggg", 
      "ggt", "gta", "gtc", "gtg", "gtt", "tac", "tat", 
      "tca", "tcc", "tcg", "tct", "tgc", "tgg", "tgt", "tta", 
      "ttc", "ttg", "ttt")
# ohne Stopcodons "taa", "tag", "tga",     

    CODON <- diag(61)
    dimnames(CODON) <- list(codon, codon)

    data = ddd$data
    index = ddd$index
    q = length(data)
    p = length(data[[1]])
    tmp <- vector("list", q)

    d = dim(data)
    att = attributes(data) 
    data = match(unlist(data), codon)
    attr(data, "dim") = d
    data = as.data.frame(data, stringsAsFactors=FALSE)
    attributes(data) = att

    row.names(data) = as.character(1:p)
    data = na.omit(data)
    rn = as.numeric(rownames(data))

    aaa = match(index, attr(data, "na.action"))
    index = index[is.na(aaa)] 
    index = match(index, unique(index))
    rn = as.numeric(rownames(data))
    attr(data, "na.action") = NULL
        
    weight = ddd$weight[rn]
    p = dim(data)[1]
    names(data) = nam
    attr(data, "row.names") = NULL
    attr(data, "weight") = weight
    attr(data, "nr") = p
    attr(data, "nc") = 61
    if (return.index) 
        attr(data, "index") = index
    attr(data, "levels") = codon
    attr(data, "allLevels") = codon
    attr(data, "type") = "CODON"
    attr(data, "contrast") = CODON    
    class(data) = "phyDat"
    data
}


as.phyDat <- function (x, ...){
    if (class(x) == "phyDat") return(x)
    UseMethod("as.phyDat")
}


as.phyDat.DNAbin <- function(x,...) phyDat.DNA(x,...)


as.phyDat.alignment <- function (x, type="DNA",...) 
{
    x$seq <- tolower(x$seq)
    data <- sapply(x$seq, strsplit, "")
    names(data) <- x$nam
    if(type=="DNA") dat <- phyDat.DNA(data,...)
    if(type=="AA") dat <- phyDat.AA(data, ...)
    if(type=="CODON") dat <- phyDat.codon(data, ...)
    if(type=="USER") dat <- phyDat.default(data, ...)
    dat
}


as.alignment.phyDat <- function(x, ...) ape:::as.alignment(as.character(x))


as.phyDat.matrix <- function (x, ...) phyDat(data=x, ...)


as.phyDat.data.frame <- function (x, ...) phyDat(data=x, ...)
 

acgt2ry <- function(obj){
   ac = c("a", "c", "g", "t", "u", "m", "r", "w", "s", "y", 
        "k", "v", "h", "d", "b", "n", "?", "-")
   AC = matrix(c(c(1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 
        0, 1, 1, 1), c(0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 
        0, 1, 1, 1, 1), c(0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 
        0, 1, 1, 1, 1, 1), c(0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 
        0, 1, 1, 1, 1, 1, 1)), 18, 4, dimnames = list(NULL, c("a", 
        "c", "g", "t")))
   ry = AC[c(7,10),]
   RY = AC %*% t(ry)
   RY[RY==2] = 1
   dimnames(RY) = list(NULL, c("r", "y"))
   attr(obj, "levels") = c("r", "y")
   attr(obj, "nc") = 2
   attr(obj, "type") = "USER"
   attr(obj, "contrast") = RY
   obj=phyDat.default(as.character(obj, allLevels=FALSE), levels = c("r", "y"), ambiguity = NULL)
   obj  
}


as.character.phyDat <- function (x, allLevels=TRUE, ...) 
{
    nr <- attr(x, "nr")
    nc <- attr(x, "nc")
    type <- attr(x, "type")
    if (type == "DNA") {
        labels <- c("a", "c", "g", "t", "u", "m", "r", "w", "s", 
            "y", "k", "v", "h", "d", "b", "n", "?", "-")
    }
    if (type == "AA") {
        labels <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", 
            "i", "l", "k", "m", "f", "p", "s", "t", "w", "y", 
            "v", "b", "z", "x", "-", "?")
    }
    if (type == "USER") {
        #levels
        if(allLevels)labels = attr(x, "allLevels")
        else{
            tmp = attr(x, "levels")
            contrast = attr(x, "contrast") # contrast=AC
            contrast[contrast>0] = 1
            ind = which(rowSums(contrast)==1)
            contrast[rowSums(contrast)>1,] = 0 
            labels = rep(NA, length(attr(x, "allLevels")))
            labels[ind] = tmp[contrast%*%c(1:length(tmp))]
            }
    }
    result = matrix(NA, nrow = length(x), ncol = nr)
    for (i in 1:length(x)) result[i, ] <- labels[x[[i]]]
    if (is.null(attr(x, "index"))) 
        index = rep(1:nr, attr(x, "weight"))
    else {
        index = attr(x, "index")
        if (is.data.frame(index)) 
            index <- index[, 1]
    }
    result = result[, index, drop = FALSE]
    rownames(result) = names(x)
    result
}


as.data.frame.phyDat <- function(x, ...){
   data.frame(t(as.character(x, ...)), stringsAsFactors=FALSE)
}


as.DNAbin.phyDat <- function(x,...) {
   if(attr(x, "type")=="DNA") return(as.DNAbin(as.character(x, ...)))
   else stop("x must be a nucleotide sequence")
}

 
phyDat <- function (data, type="DNA", levels=NULL, return.index = TRUE,...) 
{
    if (class(data) == "DNAbin") type <- "DNA"
    pt <- match.arg(type, c("DNA", "AA", "CODON", "USER"))  
    if(pt=="DNA") dat <- phyDat.DNA(data, return.index=return.index,...)
    if(pt=="AA") dat <- phyDat.AA(data, return.index=return.index, ...)
    if(pt=="CODON") dat <- phyDat.codon(data, return.index=return.index, ...)
    if(pt=="USER") dat <- phyDat.default(data, levels = levels, return.index=return.index, ...)
    dat
}


print.phyDat = function (x, ...) 
{
    cat(length(x), "sequences with",sum(attr(x,"weight")), "character and",attr(x,"nr"),"different site patterns.\n")
    cat("The states are",attr(x,"levels"), "\n")
}


c.phyDat <- function(...){
    object <- as.list(substitute(list(...)))[-1]    
    x <- list(...)
    n <- length(x) 
    match.names <- function(a,b){
        if(any(!(a %in% b)))stop("names do not match previous names") 
        }
    if (n == 1) 
        return(x[[1]])
    type <- attr(x[[1]], "type")
    nr = numeric(n)
    nr[1] <- sum(attr(x[[1]], "weight"))
    levels <- attr(x[[1]], "levels")
    snames <- names(x[[1]])
    objNames<-as.character(object)
    if(any(duplicated(objNames))) objNames <- paste(objNames,1:n,sep="")
    tmp <- as.character(x[[1]])
    for(i in 2:n){
        match.names(snames,names(x[[i]]))
        x[[i]] <- getCols(x[[i]],snames)
        nr[i] <- sum(attr(x[[i]], "weight"))
        tmp <- cbind(tmp, as.character(x[[i]]))
    }
    if (type == "DNA") 
        dat <- phyDat.DNA(tmp, return.index = TRUE)
    if (type == "AA") 
        dat <- phyDat.AA(tmp, return.index = TRUE)
    if (type == "USER") 
        dat <- phyDat.default(tmp, levels = levels, return.index = TRUE)
     if (type == "CODON") 
        dat <- phyDat.codon(tmp, return.index = TRUE)       
    attr(dat,"index") <- data.frame(index=attr(dat,"index"), genes=rep(objNames, nr))   
    dat
}


cbind.phyDat <- function(..., gaps="-"){
    object <- as.list(substitute(list(...)))[-1]    
    x <- list(...)
    n <- length(x) 
    if (n == 1) 
        return(x[[1]])
    type <- attr(x[[1]], "type")
    nr = numeric(n)
    nr[1] <- sum(attr(x[[1]], "weight"))
    levels <- attr(x[[1]], "levels")
    snames <- vector("list", n)  # names(x[[1]])
    vec = numeric(n+1)
    objNames<-as.character(object)
    if(any(duplicated(objNames))) objNames <- paste(objNames,1:n,sep="")
    tmp <- as.character(x[[1]])
    for(i in 1:n){
        snames[[i]] = names(x[[i]]) #match.names(snames,names(x[[i]]))
        nr[i] <- sum(attr(x[[i]], "weight")) 
        vec[i+1] = sum(attr(x[[i]], "weight"))
    }
    vec = cumsum(vec)
    snames = unique(unlist(snames))

    tmp = matrix(gaps, length(snames), vec[n+1], dimnames = list(snames, NULL))

    for(i in 1:n){
        nam = names(x[[i]])
        tmp[nam,(vec[i]+1):vec[i+1] ] <- as.character(x[[i]])
    }
    if (type == "DNA") 
        dat <- phyDat.DNA(tmp, return.index = TRUE)
    if (type == "AA") 
        dat <- phyDat.AA(tmp, return.index = TRUE)
    if (type == "USER") 
        dat <- phyDat.default(tmp, levels = levels, 
            return.index = TRUE)
    if (type == "CODON") 
        dat <- phyDat.codon(tmp, return.index = TRUE)            
    attr(dat,"index") <- data.frame(index=attr(dat,"index"), genes=rep(objNames, nr))   
    dat
}


write.phyDat <- function(x, file, format="phylip",...){
    if(format=="fasta") write.dna(as.character(x), file, format="fasta", ...)
    if(format=="phylip") write.dna(as.character(x), file, format="sequential", ...)    
    if(format=="nexus"){   
         type = attr(x, "type")
         if(type=="DNA") write.nexus.data(as.character(x), file, format = "dna",...)
         else write.nexus.data(as.list(as.data.frame(x)), file, format = "protein", ...)
         }
    }


read.phyDat <- function(file, format="phylip", type="DNA", ...){
    if(format=="nexus") data=read.nexus.data(file, ...)
    else {
        if(format=="phylip")format="interleaved"  #"sequential"
        if (type == "DNA" || type == "CODON"){ 
            data = read.dna(file, format, as.character = TRUE, ...)
        }
        if (type == "AA") data = read.aa(file, format=format, ...)
    }
    phyDat(data, type, return.index = TRUE)
}

 
baseFreq <- function(dat){
    if (class(dat) != "phyDat") 
        stop("data must be of class phyDat")
    levels <- attr(dat,"levels")
    weight <- attr(dat,"weight")
    n <- length(dat)    
    res <- numeric(length(levels))  
    dat <- new2old.phyDat(dat)   
    for(i in 1:n)res <- res+colSums(dat[[i]]*weight)    
    res <- res/sum(res)
    names(res) <- levels    
    res    
}


phylo <- function(edge, tip, edge.length=NULL){
    res <- list(edge=edge, tip.label=tip, edge.length=edge.length)
    class(res)="phylo"
    res
    }


getCols <- function (data, cols) 
{
    attrib = attributes(data)
    attr(data, "class") <- "list"
    data = data[cols]
    if (is.character(cols)) 
        attrib$names = cols
    else attrib$names = attrib$names[cols]
    attributes(data) = attrib
    attr(data, "class") <- "phyDat" 
    data
}


getRows <- function (data, rows, site.pattern = TRUE) 
{              
	for (i in 1:length(data)){ 
        if(is.matrix(data[[i]]))data[[i]] = data[[i]][rows,]
        else data[[i]] = data[[i]][rows]
        }  
    if(site.pattern) attr(data, "weight") = attr(data, "weight")[rows]
    else attr(data, "weight") = rep(1, length(rows))
    attr(data, "nr") = length(rows)
    attr(data, "index") = NULL
    data
}


subset.phyDat <- function (x, subset, select, site.pattern = TRUE,...) 
{  
     
    if (!missing(subset)) x <- getCols(x, subset)
    if (!missing(select)){
         if(!site.pattern)select <- attr(x, "index")[select] 
         if(any(is.na(select))) return(NULL) 
         x <- getRows(x, select, site.pattern=site.pattern)
    }    
    x 
}


unique.phyDat <- function(x, incomparables=FALSE, ...) getCols(x, !duplicated(x))


allSitePattern <- function(n,levels=c("a","c","g","t"), names=NULL){
    l=length(levels)
    X=matrix(0, l^n,n)
    for(i in 1:n)
    X[, i] = rep(rep(c(1:l), each=l^(i-1)),l^(n-i))
    for(i in 1:l)X[X==i] = levels[i]
    if(is.null(names))colnames(X) = paste("t",1:n, sep="")
    else colnames(X)=names
    phyDat.default(t(X), levels)
} 


write.phylip <- function(data, weight, file=""){
        n = sum(weight)
        m = dim(data)[2]
        cat(m,n,"\n",file = file)
        for(i in 1:m)
        cat(colnames(data)[i],"   ",toupper(rep(data[,i],weight)),"\n", sep="", file=file, append=TRUE)
}


read.aa <- function (file, format = "interleaved", skip = 0, nlines = 0, 
    comment.char = "#", seq.names = NULL) 
{
    getTaxaNames <- function(x) {
        x <- sub("^ +", "", x)
        x <- sub(" +$", "", x)
        x <- sub("^['\"]", "", x)
        x <- sub("['\"]$", "", x)
        x
    }
    format <- match.arg(format, c("interleaved", "sequential", "fasta"))
    phylip <- if (format %in% c("interleaved", "sequential")) 
        TRUE
    else FALSE
    X <- scan(file = file, what = character(), sep = "\n", quiet = TRUE, 
        skip = skip, nlines = nlines, comment.char = comment.char)      
           
    if (phylip) {
        fl <- X[1]
        oop <- options(warn = -1)
        fl.num <- as.numeric(unlist(strsplit(gsub("^ +", "", fl), " +")))
        options(oop)
        if (all(is.na(fl.num))) 
            stop("the first line of the file must contain the dimensions of the data")
        if (length(fl.num) != 2) 
            stop("the first line of the file must contain TWO numbers")
        else {
            n <- fl.num[1]
            s <- fl.num[2]
        }
        X <- X[-1]
        obj <- vector("character", n * s)
        dim(obj) <- c(n, s)
    }
    if (format == "interleaved") {
        fl <- X[1]
        fl <- unlist(strsplit(fl, NULL))
        bases <- grep("[-AaRrNnDdCcQqEeGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzXx?]", fl)        
        z <- diff(bases)
        for (i in 1:length(z)) if (all(z[i:(i + 8)] == 1)) 
            break
        start.seq <- bases[i]
        if (is.null(seq.names)) 
            seq.names <- getTaxaNames(substr(X[1:n], 1, start.seq - 1))
        X[1:n] <- substr(X[1:n], start.seq, nchar(X[1:n]))
        X <- gsub(" ", "", X)
        nl <- length(X)
        for (i in 1:n) obj[i, ] <- unlist(strsplit(X[seq(i, nl, n)], NULL))
    }
    if (format == "sequential") {
        fl <- X[1]
        taxa <- character(n)
        j <- 1
        for (i in 1:n) {
            bases <- grep("[-AaRrNnDdCcQqEeGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzXx?]", 
                unlist(strsplit(X[j], NULL)))
            z <- diff(bases)
            for (k in 1:length(z)) if (all(z[k:(k + 8)] == 1)) 
                break
            start.seq <- bases[k]
            taxa[i] <- substr(X[j], 1, start.seq - 1)
            sequ <- substr(X[j], start.seq, nchar(X[j]))
            sequ <- gsub(" ", "", sequ)
            j <- j + 1
            while (nchar(sequ) < s) {
                sequ <- paste(sequ, gsub(" ", "", X[j]), sep = "")
                j <- j + 1
            }
            obj[i, ] <- unlist(strsplit(sequ, NULL))
        }
        if (is.null(seq.names)) 
            seq.names <- getTaxaNames(taxa)
    }
    if (format == "fasta") {
        start <- grep("^ {0,}>", X)
        taxa <- X[start]
        n <- length(taxa)
        obj <- vector("list", n)
        if (is.null(seq.names)) {
            taxa <- sub("^ {0,}>", "", taxa)
            seq.names <- getTaxaNames(taxa)
        }
        start <- c(start, length(X) + 1)
        for (i in 1:n) obj[[i]] <- unlist(strsplit(gsub(" ", 
            "", X[(start[i] + 1):(start[i + 1] - 1)]), NULL))
    }
    if (phylip) {
        rownames(obj) <- seq.names
        obj <- tolower(obj)
    }
    else {
        names(obj) <- seq.names
        obj <- lapply(obj, tolower)
    }
    obj   
}


#
# tree manipulation
# 
getRoot = function(tree){
    res = unique(tree$edge[,1][!match(tree$edge[,1], tree$edge[,2], 0)])
    if(length(res)==1) return(res)
    else stop("There are apparently two root edges in your tree")
}


# renames root
reroot <-  function (tree, node) 
{
    anc = Ancestors(tree, node, "all")
    l = length(anc)
    if(is.na(match(node,tree$edge[,1])))stop("node not in tree")
    if(l==0)return(tree)
    ind = match(c(node, anc[-l]), tree$edge[, 2])
    tree$edge[ind, c(1, 2)] = tree$edge[ind, c(2, 1)]
    root = anc[l]
    tree$edge[tree$edge == root] = 0L
    tree$edge[tree$edge == node] = root
    tree$edge[tree$edge == 0L] = node
# needed for unrooted trees    
    tree <- collapse.singles(tree)
    reorderPruning(tree)
}


reroot2 <- function(tree, node) {
    if(node==getRoot(tree)) return(tree)
    anc = Ancestors(tree, node, "all")
    l = length(anc)
    ind = match(c(node, anc[-l]), tree$edge[, 2])
    tree$edge[ind, c(1, 2)] = tree$edge[ind, c(2, 1)]
    reorderPruning(tree) 
}    


midpoint <- function(tree){
    dm = cophenetic(tree)
    tree = unroot(tree)
    rn = max(tree$edge)+1
    maxdm = max(dm)
    ind =  which(dm==maxdm,arr.ind=TRUE)[1,]
    tmproot = Ancestors(tree, ind[1], "parent")
    tree = reroot(tree, tmproot)
    edge = tree$edge
    el = tree$edge.length
    children = tree$edge[,2]
    left = match(ind[1], children)
    tmp = Ancestors(tree, ind[2], "all")
    tmp= c(ind[2], tmp[-length(tmp)]) 
    right = match(tmp, children)
    if(el[left]>= (maxdm/2)){
         edge = rbind(edge, c(rn, ind[1]))       
         edge[left,2] = rn 
         el[left] = el[left] - (maxdm/2)
         el = c(el, maxdm/2) 
    }
    else{
        sel = cumsum(el[right]) 
        i = which(sel>(maxdm/2))[1]
        edge = rbind(edge, c(rn, tmp[i]))       
        edge[right[i],2] = rn  
        eltmp =  sel[i] - (maxdm/2)
        el = c(el, el[right[i]] - eltmp)
        el[right[i]] = eltmp
    }
    tree$edge.length = el
    tree$edge=edge
    tree$Nnode  = tree$Nnode+1
    reorderPruning(reroot(tree, rn))
}


pruneTree = function(tree, ..., FUN = ">="){
     if(is.null(tree$node)) stop("no node labels")
     if(is.rooted(tree)) tree = unroot(tree)
     m = max(tree$edge)
     nTips = length(tree$tip)
     bs = rep(TRUE, m)
     bs[ (nTips+1) : m] = sapply(as.numeric(as.character(tree$node)), FUN,...)    
     tree$edge.length[!bs[tree$edge[,2]]] = 0
     reorderPruning(di2multi(tree))
}


reorderPruning <- function (x, ...) 
{
    parents <- as.integer(x$edge[, 1])
    child <- as.integer(x$edge[, 2])
    root <- as.integer(parents[!match(parents, child, 0)][1])  # unique out
    if (length(root) > 2) 
        stop("more than 1 root found")
    n = length(parents)    
    m = max(x$edge)  # edge  parents 
    neworder = .C("reorder", parents, child, as.integer(n), as.integer(m), integer(n), as.integer(root-1L), DUP=FALSE)[[5]]    
    x$edge = x$edge[neworder,]
    x$edge.length = x$edge.length[neworder]
    attr(x, "order") <- "pruningwise"
    x
}


add.tip <- function(phy, n, edgeLength=NULL, tip=""){ 
     ind <- which(phy$edge[,2] == n)
     phy <- new2old.phylo(phy) 
     edge <- matrix(as.numeric(phy$edge),ncol=2)
     k <- min(edge)
     l <- max(edge)
     phy$edge <- rbind(phy$edge, c(k-1,phy$edge[ind,2]))
     phy$edge <- rbind(phy$edge, c(k-1,l+1))
     phy$edge[ind,2] = k-1 
     phy$edge.length[ind] = edgeLength[1]
     phy$edge.length <- c(phy$edge.length, edgeLength[-1])
     phy$tip.label <- c(phy$tip.label, tip) 
     phy <- old2new.phylo(phy)
     phy <- reorderPruning(phy) 
     phy
}


nnin <- function (tree, n) 
{
    tree1 = tree
    tree2 = tree
    edge = matrix(tree$edge, ncol = 2)
    parent = edge[, 1]
    child = tree$edge[, 2]
    k = min(parent) - 1
    ind = which(child > k)[n]
    p1 = parent[ind]
    p2 = child[ind]
    ind1 = which(parent == p1)
    ind1 = ind1[ind1 != ind][1]
    ind2 = which(parent == p2)
    e1 = child[ind1]
    e2 = child[ind2[1]]
    e3 = child[ind2[2]]
    tree1$edge[ind1, 2] = e2
    tree1$edge[ind2[1], 2] = e1
    tree2$edge[ind1, 2] = e3
    tree2$edge[ind2[2], 2] = e1
    if(!is.null(tree$edge.length)){
        tree1$edge.length[c(ind1, ind2[1])] = tree$edge.length[c(ind2[1] ,ind1)]
        tree2$edge.length[c(ind1, ind2[2])] = tree$edge.length[c(ind2[2] ,ind1)]
        }
    tree1 <- reorderPruning(tree1) 
    tree2 <- reorderPruning(tree2) 
    result = list(tree1, tree2)
    result
} 


nni <- function (tree) 
{
    tip.label <- tree$tip.label
    tree$tip.label <- NULL
    k = min(tree$edge[, 1]) - 1
    n = sum(tree$edge[, 2] > k)
    result = vector("list", 2*n)
    l=1
    for (i in 1:n) {
          result[c(l, l+1)] = nnin(tree, i)
          l = l + 2
          }
    attr(result, "TipLabel") <- tip.label
    class(result) <- "multiPhylo"
    result
}


allTrees <- function (n, rooted = FALSE, tip.label = NULL) 
{
	n <- as.integer(n)  
    nt <- as.integer(round(dfactorial(2 * (n + rooted) - 5))) 
    if ((n + rooted) > 10) {
        nt <- dfactorial(2 * (n + rooted) - 5)
        stop("That would generate ", round(nt), " trees, and take up more than ", 
            round(nt/1000), " MB of memory!")
    }
    if (n < 2) {
        stop("A tree must have at least two taxa.")
    }
    if (!rooted && n == 2) {
        stop("An unrooted tree must have at least three taxa.")
    }

    if (rooted) {
        edge <- matrix(NA, 2*n-2, 2)
        edge[1:2,] <- c(n+1L, n+1L, 1L, 2L)
    }
    else {
        edge <- matrix(NA,  2*n-3, 2)
        edge[1:3,] <- c(n+1L, n+1L, n+1L, 1L, 2L, 3L)
    }
    edges <- list()
    edges[[1]] <- edge

    m=1     
    nedge = 1
    trees <- vector("list", nt)
    if ((n + rooted) > 3) {
        i = 3L  + (!rooted)    
        pa = n + 2L
        nr = 2L + (!rooted)
        while(i < (n+1L)){
            nedge = nedge+2
            m2 = m * nedge 
            
            newedges <- vector("list", m2)
            for (j in 1:m) {
                edge <- edges[[j]]
                l <- nr  # nrow(edge)

                    edgeA <- edge
                    edgeB <- edge

                    for (k in 1L:l) {
                       edge = edgeA
                       node <- edge[k, 1]
                       edge[k, 1] <- pa             
                       edge[l + 1, ] <- c(pa, i)
                       edge[l + 2, ] <- c(node, pa)

                       newedges[[(j - 1) * (l + rooted) + k]] <- edge
                       }

                if(rooted) { 
                  edgeB[] <- as.integer(sub(n+1L, pa, edgeB))
                  edge = edgeB
                  edge[l + 1, ] <- c(n+1L, i)
                  edge[l + 2, ] <- c(n+1L, pa) 
                  newedges[[j * (l + 1)]] <- edge
                }
            } # end for 
            edges <- newedges
            m=m2
            i = i + 1L
            pa = pa + 1L  
            nr = nr + 2L 
        } # end for m
    } # end if
    for (x in 1:m) {
        tree <- list(edge = edges[[x]])
        tree$Nnode <- n - 2L + rooted
        class(tree) <- "phylo"       
        trees[[x]] <- reorderPruning(tree)
    }
    attr(trees, "TipLabel") <- if (is.null(tip.label)) 
        paste("t", 1:n, sep = "")
    else tip.label
    class(trees) <- "multiPhylo"
    trees
}

 
# dn sehr viel schneller
# dist.nodes ape replacement   
bipi <- function (x)
{
   if (is.null(attr(x, "order")) || attr(x, "order") == "cladewise")
       x = reorderPruning(x)
   nNode = x$Nnode
   nTips = length(x$tip)
   parent <- as.integer(x$edge[, 1])
   child <- as.integer(x$edge[, 2])
   res = vector("list", max(x$edge))
   p = parent[1]
   tmp = p
   for (i in 1:nTips) res[[i]] = i
   for (i in 1:length(parent)) {
       pi = parent[i]
       ci = child[i]
       if (pi == p) {
           if (ci < (nTips + 1))
               tmp = cisort(tmp, ci)
           else tmp = cisort(tmp, res[[ci]])
       }
       else {
           res[[p]] = (tmp)
           if (ci < (nTips + 1))
               tmp = c(ci, pi)
           else tmp = cisort(pi, res[[ci]])
           p = pi
       }
   }
   res[[p]] = (tmp)
   res
}


all.dist <- function (tree, edge.length = FALSE)
{
   lab = tree$tip.label
   l = dim(tree$edge)[1]
   bp = bipi(tree)
   m = length(bp)
   b1 = matrix(0L, m, m)
   for (i in 1:m) b1[i, bp[[i]]] = 1L
   b2 = (1 - b1)
   res1 = crossprod(b1, b2)
   res2 = crossprod(b2, b1)
   res1 + res2
}


dn <- function (x){
    if (!is.binary.tree(x) ) 
        x <- multi2di(x, random = FALSE)  
    x = reroot2(x, 1)       
    n <- length(x$tip.label)
    n.node <- x$Nnode
    N <- n + n.node
    x <- reorderPruning(x)
    res <- matrix(NA, N, N)
    res[cbind(1:N, 1:N)] <- 0
    res[x$edge] <- res[x$edge[, 2:1]] <- 1
    for (i in seq(from = 1, by = 2, length.out = n.node)) {
        j <- i + 1
        anc <- x$edge[i, 1]
        des1 <- x$edge[i, 2]
        des2 <- x$edge[j, 2]
        if (des1 > n) 
            des1 <- which(!is.na(res[des1, ]))
        if (des2 > n) 
            des2 <- which(!is.na(res[des2, ]))
        for (y in des1) res[y, des2] <- res[des2, y] <- res[anc, 
            y] + res[anc, des2]
        if (anc != 1) {
            ind <- which(x$edge[, 2] == anc)
            nod <- x$edge[ind, 1]
            l <- length(ind)
            res[des2, nod] <- res[nod, des2] <- res[anc, des2] + 
                l
            res[des1, nod] <- res[nod, des1] <- res[anc, des1] + 
                l
        }
    }
    dimnames(res)[1:2] <- list(1:N)
    res
}


rSPR = function (tree, moves = 1, n = 1, k=NULL) 
{
    if (n == 1) {
        trees = tree
        for (i in 1:moves) trees = kSPR(trees, k=k)
    }
    else {
        trees = vector("list", n)
        for (j in 1:n) {
            tmp = tree
            for (i in 1:moves) tmp = kSPR(tmp, k=k)
            tmp$tip.label = NULL
            trees[[j]] = tmp
        }
        attr(trees, "TipLabel") <- tree$tip.label
        class(trees) <- "multiPhylo"
    }
    trees
}


kSPR = function(tree, k=NULL){  
    l <- length(tree$tip.label)
    root= getRoot(tree)
    distN = dn(tree)[-c(1:l), -c(1:l)]
    distN[upper.tri(distN)]=Inf
    dN = distN[lower.tri(distN)]
    tab = table(dN) 
    tab[1] = tab[1] * 2 
    tab[-1] = tab[-1] * 8   
    if(is.null(k)) k = 1:length(tab)
    k = na.omit((1:length(tab))[k])
    if(length(k)>1)k = sample((1:length(tab))[k], 1, prob=tab[k] / sum(tab[k]) )
    if(k==1) return(rNNI(tree, 1, 1))
    index = which(distN==k, arr.ind=TRUE) + l
    m = dim(index)[1]
    if(m==0)stop("k is chosen too big")
    ind = index[sample(m, 1),]
    s1 = sample(c(1,2),1) 
    if(s1==1)res = (oneOf4(tree, ind[1], ind[2], sample(c(1,2),1), sample(c(1,2),1)))
    if(s1==2)res = (oneOf4(tree, ind[2], ind[1], sample(c(1,2),1), sample(c(1,2),1))) 
    res=reroot2(res, root)
    reorderPruning(res)    
}


oneOf4 = function(tree, ind1, ind2, from=1, to=1){
    if (!is.binary.tree(tree)) 
        stop("Sorry, trees must be binary!")        
    tree=reroot2(tree, ind2)
    trees = vector('list', 8)
    kids1 = Children(tree, ind1)
    anc = Ancestors(tree, ind1, "all")
    l = length(anc)
    kids2 = Children(tree, ind2)
    kids2 = kids2[kids2!=anc[l-1]]

    child = tree$edge[,2]
    tmp = numeric(max(tree$edge))
    tmp[child] = 1:length(child)

    edge = tree$edge
    edge[tmp[kids1[-from]],1] = Ancestors(tree, ind1, "parent")
    edge[tmp[kids2[to]],1] = ind1
    edge[tmp[ind1]] = ind2
    tree$edge=edge
    reorderPruning(tree)   
}


# faster than kSPR
rSPR_Old <- function(tree, moves=1, n=1){
    k=length(tree$edge[,1])
    if(n==1){
        trees = tree
        for(i in 1:moves) trees = sprMove(trees,sample(k,1))  
    }  
    else{
        trees = vector("list", n)
        for(j in 1:n){
            tmp = tree 
            for(i in 1:moves) tmp = sprMove(tmp,sample(k,1))
            tmp$tip.label=NULL
            trees[[j]] = tmp
        }
        attr(trees, "TipLabel") <- tree$tip.label
        class(trees) <- "multiPhylo"   
    }
    trees
}


sprMove <- function(tree, m){
    if (is.rooted(tree)) tree <- unroot(tree)
    #stop("Sorry trees must be unrooted")
    if(!is.binary.tree(tree))stop("Sorry trees must be binary!")

    reroot2 <- function(tree, node){
        anc = Ancestors(tree, node, "all")
        l = length(anc)
        ind = match(c(node, anc[-l]), tree$edge[,2])
        tree$edge[ind,c(1,2)] = tree$edge[ind,c(2,1)]
        tree    
    }    
    changeEdge <- function(tree, new, old){
        tree$edge[tree$edge==old] = 0L
        tree$edge[tree$edge==new] = old
        tree$edge[tree$edge==0L] = new
    # needed for unrooted trees
        tree <- collapse.singles(tree)
        tree          
    }

    edge = tree$edge    
    k = max(edge)
    nTips = length(tree$tip)
    nEdges = 2*nTips-3 
    if(m > nEdges) stop("m to big")

    parent = edge[,1]
    child = edge[,2]
    pv = integer(k)      
    pv[child] = parent
    cv = list()
    for(i in unique(parent)) cv[[i]] = child[parent==i]
    bp = bip(tree)
    root <- parent[!match(parent, child, 0)][1]    
       
    ch = child[m]
    pa = parent[m] 

    candidates = !logical(k)
    candidates[root] = FALSE     
    candidates[cv[[ch]]] = FALSE
    candidates[cv[[pa]]] = FALSE
    candidates[pv[pa]] = FALSE
    candidates[pa] = FALSE

    ind = which(candidates)
    l = sample(ind,1)

    cr=FALSE 

    if(!any(is.na(match(bp[[l]], bp[[ch]]))) ){
        
        newroot = cv[[ch]] #[ 1]
        newroot = newroot[newroot>nTips][1]
        tree <- reroot2(tree, newroot)
        edge = tree$edge
        parent = tree$edge[,1]
        child = tree$edge[,2]
        pv = integer(k)      
        pv[child] = parent
        cv = list()
        for(i in unique(parent)) cv[[i]] = child[parent==i]
        
        tmp = pa
        pa=ch
        ch=tmp
        cr = TRUE
    }

    if(pa==root){
        cp = cv[[pa]]
        newroot = cp[cp!=ch]
        
        newroot = newroot[newroot>nTips][1]
        if(length(newroot)==0)browser()
        #!newroot = cp[cp>nTips][1]
        tree = reroot2(tree, newroot)
        edge = tree$edge
        parent = tree$edge[,1]
        child = tree$edge[,2]
        pv = integer(k)      
        pv[child] = parent
        cv = list()
        for(i in unique(parent)) cv[[i]] = child[parent==i]
        
        cr = TRUE 
    }

    el = tree$edge.length
    cp = cv[[pa]]
    sib = cp[cp!=ch]

    edge[child==l,1] = pa
    edge[child==pa,1] = pv[l]  
    edge[child==sib,1] = pv[pa]

    el[child==sib] = el[child==sib] + el[child==pa]
    el[child==l] = el[child==l] / 2
    el[child==pa] = el[child==l]   

    tree$edge=edge
    tree$edge.length = el
    if(cr) tree <- changeEdge(tree,root,newroot)    
    tree <- reorderPruning(tree) 
    tree    
}
 

rNNI <- function(tree, moves=1, n=1){   
    k = length(na.omit(match(tree$edge[,2], tree$edge[,1])))   
    if(n==1){
        trees = tree
        for(i in 1:moves) trees = nnin(trees,sample(k,1))[[sample(2,1)]] 
    }  
    else{
        trees = vector("list", n)
        for(j in 1:n){
            tmp = tree 
            for(i in 1:moves) tmp = nnin(tmp, sample(k,1))[[sample(2,1)]]
            tmp$tip.label=NULL
            trees[[j]] = tmp
        }
        attr(trees, "TipLabel") <- tree$tip.label
        class(trees) <- "multiPhylo"   
    }
    trees
}


#
# Maximum Parsimony 
#
rowMin = function(X){
    d=dim(X)
    .Call("rowMin", X, as.integer(d[1]), as.integer(d[2]), PACKAGE = "phangorn") 
}


sankoff.quartet <- function (dat, cost, p, l, weight) 
{
    erg <- .Call("sankoffQuartet", sdat = dat, sn = p, scost = cost, 
        sk = l, PACKAGE = "phangorn")
    sum(weight * erg)
}


parsimony <- function(tree, data, method='fitch', ...){
    if (class(data)[1] != "phyDat") stop("data must be of class phyDat")
    if(method=='sankoff') result <- sankoff(tree, data, ...)
    if(method=='fitch') result <- fitch(tree, data, ...)
    result 
}


ancestral.pars <- function (tree, data, type = c("MPR", "ACCTRAN")) 
{
    call <- match.call()
    type <- match.arg(type)
    if (type == "ACCTRAN") 
        res = ptree(tree, data, retData = TRUE)[[2]]
    if (type == "MPR") 
        res = mpr(tree, data)
    l = length(tree$tip)
    
    x = attributes(data)
    m = dim(res)[2]
    label = as.character(1:m)
    nam = tree$tip.label
    label[1:length(nam)] = nam
    x[["names"]] = label

    nc = attr(data, "nc")
    result = vector("list", m) 
    Z = unique(as.vector(res))
    tmp = t(sapply(Z, function(x, k=4)dec2bin(x, nc)))
    tmp = tmp / rowSums(tmp)
    rownames(tmp) = Z
   
    for(i in 1:m){ 
#        tmp = t(sapply(res[,i], function(x, k=4)dec2bin(x, nc)))
#        result[[i]] = tmp / rowSums(tmp) no indices
         test = match(res[,i], Z)
         result[[i]] = tmp[as.character(res[,i]),]
        }

    attributes(result) = x
    attr(result, "call") <- call
    result
}


pace <- ancestral.pars


MinMaxPars <- function(tree, data){
    edge = tree$edge[,2]
    node = tree$edge[,1]
    nr = as.integer(attr(data, "nr"))
    weight = attr(data, "weight")
    dat = mpr(tree, data)
    d = length(edge)
    ext = integer(d)
    ext[is.na(match(edge, node))] = 1L
    res = matrix(0L, d, 2) 
    for(i in 1:d){
        res[i, ] = .C("countMPR", numeric(2), dat[,edge[i]], dat[,node[i]], nr, as.double(weight), ext[i])[[1]]
    }
    res 
}


mpr <- function (tree, data, returnData = FALSE) 
{
    data = prepareDataFitch(data)
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorderPruning(tree)
#    if (!is.binary.tree(tree)) 
#        warning("tree is not binary, parsimony score may be wrong!!!")
    nr = attr(data, "nr")
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    weight = as.double(attr(data, "weight"))
    m = length(edge) + 1
    dat = integer(m * nr)
    attr(dat, "dim") <- c(nr, m)
    q = length(tree$tip)
    dat[, 1:q] = data[, tree$tip.label]
    pars <- integer(nr)
    root <- as.integer(node[!match(node, edge, 0)][1])
    nTips = length(tree$tip)
    node0 = node[node != root]
    node0 = unique(node0[length(node0):1L])
    res = integer(2 * length(node0))
    sibs = Siblings(tree, node0)
    anc = Ancestors(tree, node0, type = "parent")
    k = 1
    for (i in 1:length(node0)) {
        tmp = anc[i]
        res[k] = sibs[[i]][1]
        if (tmp == root) 
            res[k + 1] = sibs[[i]][2]
        else res[k + 1] = tmp
        k = k + 2
    }
    edge2 = res
    node2 = rep(node0, each = 2)
    m2 = length(res)
    pc = rep(c(0L, 1L), length = m2)
    pc[node2 %in% Descendants(tree, root, "children")] = 0L
    na = which(is.na(res))
    if (any(is.na(res))) {
        edge2 = edge2[-na]
        node2 = node2[-na]
        pc = pc[-na]
        m2 = length(node2)
    }
    dat[(nr * (root - 1L) + 1L):(nr * root)] = 0L
    tmp <- .Call("FNALL", dat, as.integer(nr), as.integer(node), 
        as.integer(edge), as.integer(node2), as.integer(edge2), 
        as.integer(length(edge)), as.double(weight), as.integer(length(edge) + 
            1L), as.integer(m2), as.integer(length(tree$tip)), 
        as.integer(pc))

    if (!is.rooted2(tree)) {
        root = getRoot(tree)
        ind = edge[node == root]
        rSeq = .C("fitchTriplet", integer(nr), tmp[[3]][, ind[1]], 
            tmp[[3]][, ind[2]], tmp[[3]][, ind[3]], as.integer(nr))
        tmp[[3]][, root] = rSeq[[1]]
        tmp[[4]][, root] = rSeq[[1]]
    }
    result = tmp[[3]]
    for (i in node0) {
        ind = Children(tree, i)
        result[, i] = .C("fitchTriplet", integer(nr), tmp[[3]][, 
            ind[1]], tmp[[3]][, ind[2]], tmp[[4]][, i], as.integer(nr))[[1]]
    }
    attr(result, "dim") = c(nr, m)
    row.names = node0
    result
}


plotAnc <- function (tree, data, i = 1, col=NULL, ...)
{
   y = subset(data, , i)
   args <- list(...)
   CEX <- if ("cex" %in% names(args))
       args$cex
   else par("cex")
   xrad <- CEX * diff(par("usr")[1:2])/50
   levels = attr(data, "levels")
   nc = attr(data, "nc")
   y = matrix(unlist(y[]), ncol = nc, byrow = TRUE)
   l = dim(y)[1]
   dat = matrix(0, l, nc)
   for (i in 1:l) dat[i, ] = y[[i]]
   plot(tree, label.offset = 1.1 * xrad, plot = FALSE, ...)
   lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
   XX <- lastPP$xx
   YY <- lastPP$yy
   xrad <- CEX * diff(lastPP$x.lim * 1.1)/50
   par(new = TRUE)
   plot(tree, label.offset = 1.1 * xrad, plot = TRUE, ...)
   if(is.null(col)) col = rainbow(nc)
   if(length(col)!=nc) warning("Length of color vector differs from number of levels!")
   BOTHlabels(pie = y, XX = XX, YY = YY, adj = c(0.5, 0.5),
       frame = "rect", pch = NULL, sel = 1:length(XX), thermo = NULL,
       piecol = col, col = "black", bg = "lightblue", horiz = FALSE,
       width = NULL, height = NULL)
   legend("bottomright", levels, text.col = col)
}


prepareDataFitch <- function (data) 
{
    lev <- attr(data, "levels")
    l <- length(lev)
    nr <- attr(data, "nr")
    nc <- length(data)
    contrast <- attr(data, "contrast")
    tmp = contrast %*% 2L^c(0L:(l - 1L))
    tmp = as.integer(tmp)
    attrData <- attributes(data)
    nam <- attrData$names
    attrData$names <- NULL
    data = unlist(data, FALSE, FALSE)
    X = tmp[data]  
    attributes(X) <- attrData
    attr(X, "dim") <- c(nr, nc)
    dimnames(X) <- list(NULL, nam)
    X
}


compressSites <- function(data){
    attrData <- attributes(data)
    lev <- attr(data, "levels")  
    LEV <- attr(data,"allLevels")
    l <- length(lev)
    nr <- attr(data, "nr")
    nc <- length(data)

    data = unlist(data, FALSE, FALSE)

    attr(data, "dim") <- c(nr, nc)
    uni <- match(lev, LEV)    
    fun = function(x, uni) {
        u = unique.default(x)
        res=  if(any(is.na(match(u, uni)))) return(x)
        match(x, u)
    }
    data = t(apply(data, 1, fun, uni))         
    ddd = fast.table(data)
    data = ddd$data
    class(data) = "list"
    attrData$weight = tapply(attrData$weight,ddd$index, sum)
    attrData$index=NULL
    attrData$nr <- length(attrData$weight)
    attrData$compressed <- TRUE
    attributes(data) <- attrData
    data
}


  
fitch <- function (tree, data, site="pscore") 
{ 
    if (class(data) != "phyDat") 
        stop("data must be of class phyDat")
    levels <- attr(data, "levels")
    if (!is.null(tree$TipLabel)) data = subset(data, tree$TipLabel[[1]])
    data <- prepareDataFitch(data)
    d = attributes(data)
    data <- as.integer(data)
    attributes(data) <- d
    if(class(tree)=="phylo") return(fit.fitch(tree, data, site))
    if(class(tree)=="multiPhylo"){
        if(is.null(tree$TipLabel)){
            tree = unclass(tree)
            return(sapply(tree, fit.fitch, data, site))
        }    
        else{
#            data <- data[, tree$TipLabel]
            tree = unclass(tree)
            site = ifelse(site == "pscore", 1L, 0L) #ifelse(site=="pscore", 1L, 0L)
            return(sapply(tree, fast.fitch, data, site)) # 50% schneller??? 
        }       
    }
}


fit.fitch = function (tree, data, returnData = c("pscore", "site", "data")) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorderPruning(tree)
    returnData <- match.arg(returnData)
    nr = attr(data, "nr")
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    weight = attr(data, "weight")
    m = length(edge) + 1
    q = length(tree$tip)
    result <- .Call("FITCH", data[, tree$tip.label], as.integer(nr), as.integer(node), as.integer(edge), as.integer(length(edge)), as.double(weight), as.integer(m), as.integer(q), PACKAGE = "phangorn")
    if (returnData == "site") return(result[[2]])
    pscore <- result[[1]]
    res = pscore
    if (returnData == "data") 
        res <- list(pscore = pscore, dat = result[[3]], site = result[[2]])
    res
}   



fast.fitch <- function (tree, data, ps = TRUE) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorderPruning(tree)
    nr = attr(data, "nr")
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    weight = attr(data, "weight")
    m = length(edge) + 1L
    q = dim(data)[2]
    .Call("FITCH123", data, as.integer(nr), as.integer(node), as.integer(edge), as.integer(length(edge)), 
        as.double(weight), as.integer(m), as.integer(q), as.integer(ps), PACKAGE = "phangorn")
}

     
fnodes <- function (tree, data, external = FALSE) 
{
    if(is.rooted(tree)){
         tree = unroot(tree)
         warning("tree is rooted, I unrooted the tree!")
    }
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorderPruning(tree)
    if (!is.binary.tree(tree)) 
        warning("tree is not binary, parsimony score may be wrong!!!")
    nr = attr(data, "nr")
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    weight = attr(data, "weight")
    m = length(edge) + 1
    dat = integer(m * nr)
    attr(dat, "dim") <- c(nr, m)
    q = length(tree$tip)
    dat[, 1:q] = data[, tree$tip.label]
    pars <- integer(nr)
    root <- as.integer(node[!match(node, edge, 0)][1])
    nTips = length(tree$tip)   
    if(external) node0 = edge
    else node0 = node[node != root]
    node0 = unique(node0[length(node0):1L])  
    res = integer(2*length(node0))
    sibs = Siblings(tree, node0)
    anc = Ancestors(tree, node0, type = "parent")
# fuer nicht binaere Baeume programmieren    
    k = 1
    for (i in 1:length(node0)) {
        tmp = anc[i]
        res[k] = sibs[[i]][1]
        if (tmp == root) 
            res[k + 1] = sibs[[i]][2]
        else res[k + 1] = tmp
        k = k + 2
    }
    node2 = rep(node0, each = 2)
    edge2 = res
    m2 = length(res)

    dat[(nr * (root - 1L) + 1L):(nr * root)] = 0L

    pc = rep(c(0L, 1L), length = m2)   
    pc[node2 %in% Descendants(tree, root, "children")] = 0L #res[, 1]
    
    .Call("FNALL", dat, as.integer(nr), as.integer(node), as.integer(edge), as.integer(node2), as.integer(edge2),
        as.integer(length(edge)), as.double(weight), as.integer(length(edge) + 1L), as.integer(m2), as.integer(length(tree$tip)),  
        as.integer(pc) )
}     


fnodesFast = function (tree, data, external = FALSE)  {
    if (is.rooted(tree)) {
        tree = unroot(tree)
        warning("tree is rooted, I unrooted the tree!")
    }
    if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") 
        tree <- reorderPruning(tree)
    nr = attr(data, "nr")
    m = max(tree$edge)
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    weight = attr(data, "weight")
    m = length(edge) + 1
    dat = integer(m * nr)
    attr(dat, "dim") <- c(nr, m)
    q = length(tree$tip)
    dat[, 1:q] = data[, tree$tip.label]
    pars <- integer(nr)
    root <- as.integer(node[!match(node, edge, 0)][1])
    nTips = length(tree$tip)
    if (external) 
        node0 = edge
    else node0 = node[node != root]
    node0 = unique(node0[length(node0):1L])
    res = integer(2 * length(node0))

    pvector <- integer(m)
    pvector[edge] <- node
    ch= vector("list", m)
    for(i in 1:length(node)) ch[[node[i]]] = c(ch[[node[i]]], edge[i])
    k=1 
    sibs=list(length(node0))
    for (i in node0) {
        if (i != root) {
        tmp <- ch[[pvector[i]]]
        sibs[[k]] = tmp[tmp != i]
        }
        k = k + 1
    }
    anc = pvector[node0] 
    k = 1
    for (i in 1:length(node0)) {
        tmp = anc[i]
        res[k] = sibs[[i]][1]
        if (tmp == root) 
            res[k + 1] = sibs[[i]][2]
        else res[k + 1] = tmp
        k = k + 2
    }
    node2 = rep(node0, each = 2)
    edge2 = res
    m2 = length(res)
    dat[(nr * (root - 1L) + 1L):(nr * root)] = 0L
    pc = rep(c(0L, 1L), length = m2)
    pc[pvector[node2] == root] = 0L # 1:6
    .Call("FNALL", dat, as.integer(nr), as.integer(node), as.integer(edge), 
        as.integer(node2), as.integer(edge2), as.integer(length(edge)), 
        as.double(weight), as.integer(length(edge) + 1L), as.integer(m2), 
        as.integer(length(tree$tip)), as.integer(pc))
}


is.rooted2 = function(tree){
    length(tree$edge[, 1][!match(tree$edge[, 1], tree$edge[, 2], 0)]) < 3
}


optim.fitch <- function(tree, data, trace=1, ...) {
    if(class(tree)!="phylo") stop("tree must be of class phylo") 
    if(!is.binary.tree(tree)) tree <- multi2di(tree)
    if(is.rooted(tree))tree <- unroot(tree)
    if(is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") tree <- reorderPruning(tree)  
    if (class(data)[1] != "phyDat") stop("data must be of class phyDat")
    if(is.null(attr(data, "compressed")) || attr(data, "compressed") == FALSE) data <- compressSites(data)

    rt = FALSE
    dat <- prepareDataFitch(data)
    l <- attr(data, "nc")

    tree$edge.length=NULL
    swap = 0
    iter = TRUE
    pscore <- fit.fitch(tree, dat)
    while (iter) {
        res <- fitch.nni(tree, dat, ...)
        tree <- res$tree
        if(trace>1)cat("optimize topology: ", pscore , "-->", res$pscore, 
            "\n")
        pscore = res$pscore
        swap = swap + res$swap
        if (res$swap == 0) iter = FALSE
        }
    if(trace>0)cat("Final p-score",pscore,"after ",swap, "nni operations \n") 
    if(rt)tree <- ptree(tree, data)  
    attr(tree, "pscore") = pscore
    tree
}


fitch.nni <- function (tree, data, ...) 
{
    if (is.rooted(tree)) 
        tree <- reorderPruning(unroot(tree))
    INDEX <- indexNNI(tree)
    rootEdges <- attr(INDEX, "root")
    if (class(data) != "phyDat") 
        stop("data must be of class phyDat")
    levels <- attr(data, "levels")
    l = length(levels)
    weight = attr(data, "weight")
    p = attr(data, "nr")  
    kl = TRUE
    i = 1
    tmp = fnodes(tree, data)
    p0 = tmp[[1]]
    datf = tmp[[3]] # no matrix needed
    datp = tmp[[4]] #matrix(tmp[[2]], nrow = p)
    parent = tree$edge[, 1]
    child = tree$edge[, 2]
    m <- dim(INDEX)[1]
    k = min(parent)
    pscore = numeric(m * 2)
    datn = integer(p*6L)
    psc = numeric(4) 
    for (i in 1:m) {
        ei = INDEX[i, ]
        psc <- tmp[[5]][ei[1:4]]
        datn[1:(4L*p)] <- datf[, ei[1:4]]        
        if (!(ei[5] %in% rootEdges)){
            datn[1:p] <- datp[, ei[1]]
            psc[1] = tmp[[6]][ei[1]]
            }     
        tmp2 <- p0 - sum(psc)
        pscore[(2 * i) - 1] <- fitch.quartet(datn, p, weight, edge = c(1L, 3L, 2L, 4L, 6L)) - tmp2
        pscore[(2 * i)] <- fitch.quartet(datn, p, weight, edge = c(1L, 4L, 3L, 2L, 6L)) - tmp2
    }
    swap <- 0
    candidates <- pscore < 0

    while (any(candidates)) {
        ind = which.min(pscore)
        pscore[ind] = Inf
        if (ind%%2) 
            swap.edge = c(2, 3)
        else swap.edge = c(2, 4)
        tree2 <- changeEdge(tree, INDEX[(ind + 1)%/%2, swap.edge])
#         tree2 <- changeEdge(tree, INDEX2[ind, swap.edge])
        test <- fit.fitch(tree2, data) # fast.fitch
        if (test >= p0) 
            candidates[ind] = FALSE
        if (test < p0) {
            p0 <- test
            swap = swap + 1
            tree <- tree2
            candidates[ind] = FALSE          
            indi <- unique(which(INDEX %in% (INDEX[(ind + 1)%/%2, ]))%%m)
            indi[indi == 0] <- m
            indi = c((2 * indi) - 1, 2 * indi)
            candidates[indi] <- FALSE
            pscore[indi] <- Inf
        }
    }
    list(tree = tree, pscore = p0, swap = swap)
}


fitch.quartet = function(data, nr, weight, node = c(6L, 6L, 5L, 5L, 5L), edge = c(1L, 2L, 3L, 4L, 6L)) 
{
    .Call("FITCH", data, as.integer(nr), node, edge, 5L, as.double(weight), 6L , 4L)[[1]]
}

####################
# Branch and bound #
####################

getorder <- function (x) 
{
    label = names(x)
    dm = as.matrix(dist.hamming(x, FALSE))
    ind = which(dm == max(dm), arr.ind = TRUE)[1, ]
    added = label[ind]
    remaining <- label[-ind]

    tree <- stree(n = 3, tip.label = label[1:3])

    l = length(remaining)
    res = numeric(l)
    for(i in 1:l){
        tree$tip.label= c(added,remaining[i])
        res[i] = fitch(tree, x) 
    }
    tmp = which.max(res)
    added = c(added, remaining[tmp])
    remaining <- remaining[-tmp]
    tree$tip.label= added

    n = length(x) - 1L
    x = prepareDataFitch(x) 

    for (i in 4:n) {   
        temp <- add.everywhere(tree, remaining[1])
        l = length(remaining)
        res = numeric(l)
        nt = numeric(l)
        k = length(added)+1L
        for(j in 1:l){
            ind = c(added, remaining[j])
            datTmp = x[, ind]
            mostattributes(datTmp) = attributes(x)
            attr(datTmp, "dim") = c(attr(datTmp, "nr"), k)
            score = vapply(temp, fast.fitch, 0, datTmp)
            res[j] = min(score) 
            nt[j] = which.min(score)
        }
        tmp = which.max(res)
        added = c(added, remaining[tmp])
        remaining <- remaining[-tmp]
        tree = temp[[nt[tmp]]]
        tree$tip.label= added
    }
    temp <- add.everywhere(tree, remaining[1])
    added = c(added, remaining) 
    added 
}

#
# still too slow 
#
fitch.triplet = function(dat1, dat2, dat3, nr, weight) 
{
    .Call("FITCH", c(dat1, dat2, dat3), as.integer(nr), node = c(4L, 4L, 4L), edge = c(1L, 2L, 3L), 3L, as.double(weight), 4L , 3L)[[1]]
}


fitch.trip <- function(dat1, dat2, dat3, nr, weight){ 
    .C("fitch3", integer(nr), dat1, dat2, dat3, as.integer(nr), integer(nr), as.double(weight), double(1L), DUP=FALSE)[[8]]
}


random.addition <- function(x, type="fitch") 
{
    label = sample(names(x))
    ind = 1:3
    added = label[ind]
    remaining <- label[-ind]

    x = subset(x, label)
    
    tree <- stree(n = 3, tip.label = label[1:3])
    tree$tip.label= added
     
    nr = attr(x, "nr")
    storage.mode(nr) <- "integer"
    n = length(x) #- 1L
    x = prepareDataFitch(x) 

    weight = attr(x, "weight")
    storage.mode(weight) <- "double"
    m = 3
    for (i in 4:n) {   
        
        dat = fnodes(tree, x, TRUE)
        score = numeric(m+1)
        score = (dat[[5]] + dat[[6]]) 
        score[i] = Inf
        edge = tree$edge[,2]
        for(j in edge){       
# ohne function 
# as.double  
# .C("fitch3", integer(nr), dat[[3]][,j], dat[[4]][,j], x[,i], nr, integer(nr), weight, double(1L), DUP=FALSE)[[8]]
             res = fitch.triplet(dat[[3]][,j], dat[[4]][,j], x[,i], nr, weight)
             score[j] = score[j] + res
        }              

        m=m+2
        l = length(remaining)
        k = length(added)+1L
        added = c(added, remaining[1])

        res = min(score) 
        nt = which.min(score)
        nt = match(nt, edge)
        
        remaining <- remaining[-1]
        tree = add.one(tree, label[i], nt) 
        tree$tip.label= added
    }
    attr(tree, "pscore") = res
    tree 
}


parsinfo <- function (x) 
{
    low = lowerBound(x)
    up = upperBound(x)
    ind = which(low==up)
    cbind(ind, low[ind])
}


lowerBound <- function(x, cost=NULL){
    tip <- names(x)
    att = attributes(x)
    nc = attr(x, "nc")
    nr = attr(x, "nr")
    contrast = attr(x, "contrast")
    rownames(contrast) = attr(x, "allLevels")
    colnames(contrast) = attr(x, "levels")
    attr(x, "weight") = rep(1, nr)
    attr(x, "index") = NULL
 
    y <- as.character(x)
    states <- apply(y, 2, unique.default)
    if(nr==1) nst <- length(states)   
    else nst <- sapply(states, length)

    res = numeric(nr)
    ust = sort(unique(nst))
    if(any(ust>1)){ 
        ust = ust[ust>1]
        m <- max(ust)    
        tips = paste("t", 1:m, sep="") 
        if(is.null(cost))cost <- 1 - diag(nc) 
        for(i in ust){
            dat = matrix(unlist(states[nst==i]), nrow=i, dimnames=list(tips[1:i], NULL))
            dat = phyDat(dat, type="USER", contrast=contrast)      
            tree = stree(i)
            res[nst==i] = sankoff(tree, dat, cost=cost, site="site")[attr(dat, "index")]
        }
    }
    res
}


upperBound <- function(x, cost=NULL){
    tree = stree(length(x), tip.label=names(x))
    if(is.null(cost))cost <- 1 - diag(attr(x, "nc")) 
    sankoff(tree, x, cost=cost, site="site")
}


CI <- function (tree, data){
    pscore = sankoff(tree, data)
    weight = attr(data, "weight")
    data = subset(data, tree$tip.label) 
    m = lowerBound(data)    
    sum(m * weight)/pscore
}


RI <- function (tree, data){
    pscore = sankoff(tree, data)
    data = subset(data, tree$tip.label) 
    weight = attr(data, "weight")
    m = lowerBound(data)
    m = sum(m * weight)
    g = upperBound(data)
    g = sum(g * weight)
    (g - pscore)/(g - m)
}


add.everywhere <- function(tree,tip.name, rooted = FALSE){
       if(class(tree)!="phylo") stop("tree should be an object of class 'phylo.'")
#       if(!rooted)tree = unroot(tree)
       nTips = length(tree$tip)
       tmpedge = tree$edge 
       m = max(tmpedge)
       l = nrow(tmpedge)
       trees<-vector("list", l) 
       tmp = tree
       tmp$tip.label = c(tree$tip.label, tip.name)
       tmpedge[tmpedge>nTips] <- tmpedge[tmpedge>nTips] + 1L

       tmp$Nnode = tmp$Nnode + 1L
       tmp$edge.length <- NULL
       tmpedge = rbind(tmpedge, matrix(c(m+2L, m+2L, 0L, nTips+1L),2,2))
       for(i in 1:l){
            edge = tmpedge
            edge[l+1L,2] <- edge[i,2]
            edge[i, 2] <- m+2L
            neworder = .C("reorder", edge[,1], edge[,2], as.integer(l+2L), as.integer(m+2L), 
                integer(l+2L), as.integer(nTips+1L), DUP = FALSE)[[5]]
            tmp$edge <- edge[neworder,]
            trees[[i]] = tmp
       }
       return(trees)
}


add.one <- function (tree, tip.name, i){
    if (class(tree) != "phylo") 
        stop("tree should be an object of class 'phylo.'")
    nTips = length(tree$tip)
    tmpedge = tree$edge
    m = max(tmpedge)
    l = nrow(tmpedge)
    trees <- vector("list", l)
    tmp = tree
    tmp$tip.label = c(tree$tip.label, tip.name)
    tmpedge[tmpedge > nTips] <- tmpedge[tmpedge > nTips] + 1L
    tmp$Nnode = tmp$Nnode + 1L
    tmp$edge.length <- NULL
    tmpedge = rbind(tmpedge, matrix(c(m + 2L, m + 2L, 0L, nTips + 1L), 2, 2))
    edge = tmpedge
    edge[l + 1L, 2] <- edge[i, 2]
    edge[i, 2] <- m + 2L
    neworder = .C("reorder", edge[, 1], edge[, 2], as.integer(l + 
           2L), as.integer(m + 2L), integer(l + 2L), as.integer(nTips + 
           1L), DUP = FALSE)[[5]]
    tmp$edge <- edge[neworder, ]
    tmp
}


mmsNew0 <- function (x, Y) 
{
    w <- attr(x, "weight")
    names(w) = NULL
    m = length(x)
    data <- matrix(unlist(x[1:m]), ncol = m)
    l = nrow(data)
    v = Y[,1] + 1L
#    v = numeric(l)
#    for (i in 1:l) v[i] = length(.Internal(unique(data[i, ], 
#        FALSE, FALSE)))
    result = matrix(NA, sum(w), 6)
    Res = matrix(NA, sum(w), m)
    Res2 = matrix(NA, sum(w), m)
    j = 1
    res = 0
    bin = as.integer(2L^c(0L:30L))
    for (i in 1:(l - 1L)) {
        if (w[i] > 0) {
            v2 = v[i] + v[(i + 1L):l] - 2L
            v3 = integer(l - i)
            ind = which(w[(i + 1):l] > 0)
            V3 = matrix(NA, m, l - i)
            k = length(ind)
            V3[, ind] <- t(data[ind + i, , drop = FALSE]) + 100L * 
                data[i, ]
            v3[ind] <- apply(V3[, ind, drop = FALSE], 2, function(x) {
                length(.Internal(unique(x, FALSE, FALSE))) - 
                  1L })
            r = v3 - v2
            while (any(r > 0) && w[i] > 0) {
                a = which.max(r)
                w0 = min(w[i], w[i + a])
                if (w0 == 0) {
                  r[a] = 0
                }
                else {
                  res = res + w0 * v3[a]
                  w[i] = w[i] - w0
                  w[i + a] = w[i + a] - w0
                  result[j, ] = c(i, a + i, w0, r[a], v3[a], 
                    v2[a])
                  abc = V3[, a]
                  Res[j, ] = bin[match(abc, unique(abc))]
                  Res2[j, ] = match(abc, unique(abc))
                  r[a] = 0
                  j = j + 1
                }
            }
        }
    }
    result = na.omit(result)
    mm = max(result[, 5])
    Res = na.omit(Res)
    Res2 = na.omit(Res2)
    maxr = max(Res2)
    resm = apply(Res2, 1, function(x) {
            length(.Internal(unique(x, FALSE, FALSE))) - 1L })
    Res2 = t(Res2)
    Res2 = phyDat(Res2, type="USER", levels=1:maxr)
    names(Res2) = as.character(1:m)
    resm = lowerBound(Res2)
    ind = which(w > 0)
#    data = data[ind, ]

    tmp = matrix(0, attr(Res2, "nr"), m)
    for (i in 4:m) {
        tmp[, i] = resm - upperBound(subset(Res2, 1:i))
    }
    tmp = tmp[attr(Res2, "index"), , drop=FALSE]
    tmp2 = Y[result[,1],] + Y[result[,2],]
    tmp3 = pmax(tmp, tmp2) 
#    Res = rbind(Res, data[ind, ])
    tmp = rbind(tmp3, Y[ind, ])
    weight = c(result[, 3], w[ind])
    res = t(tmp) %*% weight
    #res[m] - res
    res
}


bab <- function (data, tree = NULL, trace = 1, ...) 
{
    if (is.null(attr(data, "compressed")) || attr(data, "compressed") == 
        FALSE) 
        data <- phangorn:::compressSites(data)
    o = order(attr(data, "weight"), decreasing = TRUE)
    data = subset(data, , o)
    nr <- attr(data, "nr")
    pis <- phangorn:::parsinfo(data)
    p0 <- sum(attr(data, "weight")[pis[, 1]] * pis[, 2])
    if (length(pis) > 0) 
        data <- phangorn:::getRows(data, c(1:nr)[-pis[, 1]], 
            TRUE)
    tree = pratchet(data, start = tree, trace = trace - 1, ...)
    nr <- as.integer(attr(data, "nr"))
    label = phangorn:::getorder(data)
    data = subset(data, label)
    lb = lowerBound(data)
    m = length(data)
    lb2 = apply(matrix(unlist(data[1:m]), ncol = m), 1, function(x) length(.Internal(unique(x, 
        FALSE, FALSE))) - 1L)
    Y=NULL
    
    ind = which(lb2 > lb, useNames = FALSE)
    nr <- as.integer(attr(data, "nr"))
    TMP = matrix(0, m, nr)
    for (i in 4:m) {
        TMP[i, ] = lowerBound(subset(data, 1:i))
    }
    Y = t(TMP) #[-ind, ]
    Y = Y[, m] - Y

    mmsAmb = 0
    mmsAmb = TMP %*% attr(data, "weight") # TMP[, ind] %*% attr(data, "weight")[ind] 
    mmsAmb = mmsAmb[m] - mmsAmb
    if(any(lb2 > lb)){
        Y = Y[-ind, ]
        data3 = subset(data, , (1:nr)[-ind])
    }
    else data3 = data
    mms0 = 0 
#    mms0 = mmsNew0(data3, Y)
#    if (any(lb2 > lb)) 
    mms0 = mms0 + mmsAmb
    data <- phangorn:::prepareDataFitch(data)
    minPars = mms0[1]
    kPars = 0
    weight = as.double(attr(data, "weight"))
    if (trace) 
        print(paste("lower bound:", p0 + mms0[1]))
    bound <- phangorn:::fit.fitch(tree, data)
    if (trace) 
        print(paste("upper bound:", bound + p0))
    blub <- list(stree(n = 3, tip.label = label[1:3]))
    attr(blub[[1]], "order") = "pruningwise"
    new = list(blub[[1]]$edge)
    added <- label[1:3]
    remaining <- setdiff(label, added)
    k = 4L
    Nnode = 1L
    while (length(remaining) > 0) {
        old <- new
        new <- vector("list", length(old) * (2 * k - 5))
        new.tip <- remaining[1]
        pscores <- numeric(length(old) * (2 * k - 5))
        added <- c(added, new.tip)
        datTmp = data[, added]
        mostattributes(datTmp) = attributes(data)
        attr(datTmp, "dim") = c(attr(datTmp, "nr"), k)
        l = 0
        mmsi = mms0[k]
        for (i in 1:length(old)) {
            temp <- phangorn:::add.every(old[[i]], Nnode)
            score = vapply(temp, fast.fitch3, 0, data=datTmp, weight=weight, nr=nr, USE.NAMES = FALSE)
            score = score + mmsi
            ind = score <= bound
            tmp2 = sum(ind)
            if (tmp2 > 0) {
                new[(l + 1):(l + tmp2)] = temp[ind]
                pscores[(l + 1):(l + tmp2)] <- score[ind]
            }
            l = l + tmp2
        }
        new = new[1:l]
        pscores = pscores[1:l]
        if (trace) {
            print(paste(length(added), "species added;", length(new), 
                "trees retained", collapse = ""))
            print(paste("lower bound:", min(pscores) + p0))
        }
        remaining <- setdiff(label, added)
        k = k + 1L
        Nnode = Nnode + 1L
    }
    edges <- new[pscores == min(pscores)]
    trees = vector("list", length(edges))
    for (i in 1:length(trees)) {
        trees[[i]]$edge = edges[[i]]
        attr(trees[[i]], "pscore") <- min(pscores) + p0
        trees[[i]]$tip.label = label
        trees[[i]]$Nnode = Nnode
        class(trees[[i]]) = "phylo"
    }
    class(trees) <- "multiPhylo"
    return(trees)
}


add.every <- function (tmpedge, Nnode) 
{
    m = max(tmpedge)
    nTips = as.integer(m - Nnode)
    l = nrow(tmpedge)
    result <- vector("list", l)
    tmpedge[tmpedge > nTips] <- tmpedge[tmpedge > nTips] + 1L
    newE = matrix(c(m + 2L, m + 2L, 0L, nTips + 1L), 2, 2)
    edge0 = matrix(0L, l+2L, 2L)
    l2= l-3L
    ind = c(rep( 0L : (l2/2L),  each = 2)*2L, l2)
    for (i in 1:l) {
        edge1 <- newE
        edge2 <- tmpedge
        edge1[1,2] <- tmpedge[i, 2]
        edge2[i,2] <- m + 2L
        k = ind[i]    
        edge = edge0      
        if(k>0) edge[1:k,] = edge2[1:k,]
        edge[(k+1L):(k+2L),] = edge1
        edge[(k+3L):(l+2L),] = edge2[(k+1):l,]
        result[[i]] <- edge
    }
    return(result)
}


allAncestors <- function(x){
    x = reorderPruning(x)
    parents <- x$edge[, 1]
    child <- x$edge[, 2]
    l = length(parents)
    res <- vector("list",max(x$edge))
    for(i in l:1){
          pa = parents[i]  
          res[[child[i]]] = c(pa, res[[pa]])
     } 
     res
}


fast.fitch3 <- function (edge, data, ps = TRUE, weight, nr){
    l = dim(edge)[1] 
    m = l+1L
    q = dim(data)[2]
    .Call("FITCH123", data, as.integer(nr), as.integer(edge[, 1]), 
        as.integer(edge[, 2]), as.integer(l), as.double(weight), 
        as.integer(m), as.integer(q), as.integer(ps), PACKAGE = "phangorn")
}


###########
# Sankoff #
###########

old2new.phyDat <- function(data){}


new2old.phyDat <- function(data){
    contrast = attr(data, "contrast")
    for(i in 1:length(data)) data[[i]] = contrast[data[[i]],,drop=FALSE]
    data
    }


prepareDataSankoff <- function(data){
    contrast = attr(data, "contrast")
    contrast[contrast == 0] = 1e+06
    contrast[contrast == 1] <- 0
    for (i in 1:length(data)) data[[i]] = contrast[data[[i]], , drop = FALSE]
    data
}


sankoff <- function (tree, data, cost = NULL, site = 'pscore') 
{
    if (class(data) != "phyDat") 
        stop("data must be of class phyDat")
    data <- prepareDataSankoff(data)
    levels <- attr(data, "levels")
    l = length(levels)  

    if (is.null(cost)) {
        cost <- matrix(1, l, l)
        cost <- cost - diag(l)
    }   
    for (i in 1:length(data)) storage.mode(data[[i]]) = "double"
    if(class(tree)=="phylo") return(fit.sankoff(tree, data, cost, returnData =site))
    if(class(tree)=="multiPhylo"){
	    if(is.null(tree$TipLabel))tree = unclass(tree)
	    return(sapply(tree, fit.sankoff, data, cost, site))
    }    
}


fit.sankoff <- function (tree, data, cost, returnData = c("pscore", "site", "data")) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorderPruning(tree)
    returnData <- match.arg(returnData) 
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    weight = attr(data, "weight")
    nr = p = attr(data, "nr")
    q = length(tree$tip.label)
    nc = l = attr(data, "nc")
    m = length(edge) + 1
    dat = vector(mode = "list", length = m)
    dat[1:q] = data[tree$tip.label]
    node = as.integer(node - 1)
    edge = as.integer(edge - 1)
    nTips = as.integer(length(tree$tip))
    mNodes = as.integer(max(node) + 1)
    tips = as.integer((1:length(tree$tip))-1)
    res <- .Call("sankoff3", dat, as.numeric(cost), as.integer(nr),as.integer(nc),
         node, edge, mNodes, tips, PACKAGE="phangorn")  
    root <- getRoot(tree) 
    erg <- .Call("rowMin", res[[root]], as.integer(nr), as.integer(nc), PACKAGE = "phangorn")
    if (returnData=='site') return(erg)
    pscore <- sum(weight * erg)
    result = pscore
    if (returnData=="data"){ 
        #if (is.null(attr(data, "index"))) 
        #    index = rep(1:nr, attr(data, "weight"))
        #else index = attr(data, "index")
        #    result <- res[[root]][index,]#
        result <- list(pscore = pscore, dat = res)
        }
    result
}


pnodes <- function (tree, data, cost) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorderPruning(tree)
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    nr = nrow(data[[1]])
    nc = ncol(data[[1]])
    node = as.integer(node - 1)
    edge = as.integer(edge - 1)  
    .Call("pNodes", data, as.numeric(cost), as.integer(nr),as.integer(nc),
         node, edge, PACKAGE="phangorn")
}

           
indexNNI <- function(tree){
    parent = tree$edge[, 1]
    child = tree$edge[, 2]
 
    ind = which(child %in% parent)
    Nnode = tree$Nnode
    edgeMatrix = matrix(0,(Nnode-1),5)

    pvector <- numeric(max(parent))
    pvector[child] <- parent
    tips  <- !logical(max(parent))
    tips[parent] <-  FALSE
    cvector <- vector("list",max(parent))   
    for(i in 1:length(parent))  cvector[[parent[i]]] <- c(cvector[[parent[i]]], child[i]) 
    k=0
    for(i in ind){        
            p1 = parent[i]
            p2 = child[i]
            e34 = cvector[[p2]]
            ind1 = cvector[[p1]]
            e12 = ind1[ind1 != p2]
            if(pvector[p1])e12=c(p1,e12)
            edgeMatrix[k+1, ] = c(e12,e34,p2)
            k=k+1
    } 
    attr(edgeMatrix, 'root') <-cvector[[min(parent)]]  
    edgeMatrix
}
                   
        
changeEdge <- function(tree, swap, edge=NULL, edge.length=NULL){
    child <- tree$edge[,2]
    tmp = numeric(max(child))
    tmp[child] = 1:length(child)
    tree$edge[tmp[swap[1]], 2] = swap[2]
    tree$edge[tmp[swap[2]], 2] = swap[1]
    if(!is.null(edge)){
        tree$edge.length[tmp[edge]] = edge.length
    }
    reorderPruning(tree) 
}


sankoff.nni = function (tree, data, cost, ...) 
{   
    if(is.rooted(tree))tree<- reorderPruning(unroot(tree))     
    INDEX <-  indexNNI(tree)
    rootEdges <- attr(INDEX,"root")
    if (class(data) != "phyDat") 
        stop("data must be of class phyDat")
    levels <- attr(data, "levels")
    l = length(levels)
    weight = attr(data, "weight")
    p = attr(data, "nr")
    kl = TRUE
    i = 1
    tmp = fit.sankoff(tree, data, cost, returnData='data')
    p0 = tmp[[1]]
    datf = tmp[[2]]
    datp = pnodes(tree, datf, cost) 
    
    parent = tree$edge[,1]
    child = tree$edge[,2]
    m <- dim(INDEX)[1]
    k = min(parent)
    pscore = numeric(2*m)

    for(i in 1:m){
        ei = INDEX[i,]
        datn <- datf[ei[1:4]]
        if (!(ei[5] %in% rootEdges)) datn[1] = datp[ei[1]]
        pscore[(2*i)-1] <- sankoff.quartet(datn[ c(1, 3, 2, 4)], 
            cost, p, l, weight)
        pscore[(2*i)] <- sankoff.quartet(datn[ c(1, 4, 3, 2)], 
            cost, p, l, weight)    
    }
    swap <- 0
    candidates <- pscore < p0
    while(any(candidates)){
    
        ind = which.min(pscore)
        pscore[ind]=Inf
        if( ind %% 2 ) swap.edge = c(2,3)
        else swap.edge = c(2,4)

        tree2 <- changeEdge(tree, INDEX[(ind+1)%/%2,swap.edge])
        test <- fit.sankoff(tree2, data, cost, 'pscore')

        if(test >= p0) candidates[ind] = FALSE
        if(test < p0) {
            p0 <- test
            swap=swap+1
            tree <- tree2
            candidates[ind] = FALSE
            indi <- which(rep(colSums(apply(INDEX,1,match,INDEX[(ind+1)%/%2,],nomatch=0))>0,each=2))
            candidates[indi] <- FALSE
            pscore[indi] <- Inf
        }
    }
    list(tree = tree, pscore = p0, swap = swap)
}


optim.parsimony <- function(tree,data, method='fitch', cost=NULL, trace=1, ...){
    if(method=='fitch') result <- optim.fitch(tree=tree, data=data, trace=trace, ...) 
    if(method=='sankoff') result <- optim.sankoff(tree=tree, data=data, cost=cost, trace=trace, ...)
    result 
}


pratchet <- function (data, start = NULL, k = 20, np = 1, trace = 1, all=FALSE, method="fitch", multicore=FALSE,  ...) 
{
    eps = 1e-08
    if(method=="fitch" && (is.null(attr(data, "compressed")) || attr(data, "compressed") == FALSE)) data <- compressSites(data)
    trace = trace - 1
    uniquetree <- function(trees) {
        k = 1
        res = trees[[1]]
        result = list()
        result[[1]]=res
        k=2
        trees = trees[-1]
        while (length(trees) > 0) {
            rf = sapply(trees, RF.dist, res, FALSE) # added FALSE 
            if(any(rf==0))trees = trees[-which(rf == 0)]
            if (length(trees) > 0) {
                res = trees[[1]]
                result[[k]] = res
                k=k+1 
                trees = trees[-1]
            }
        }
        result
    }
    if (is.null(start)) 
        start = optim.parsimony(nj(dist.hamming(data)), data, trace = trace, method=method, ...)
    tree = start
    subset(data, tree$tip.label) 
    attr(tree, "pscore") = parsimony(tree, data, method=method, ...)
    if (trace >= 0) 
        print(paste("Best pscore so far:",attr(tree, "pscore")))
        FUN = function(data, tree, method=method, ...) {
            optim.parsimony(tree, data = data, method=method, ...)
    }
    result = list()
    result[[1]] = tree
    for (i in 1:k) {
        bstrees <- bootstrap.phyDat(data, FUN, tree = tree, bs = np, 
            trace = trace, method=method, ...)
        eval.success <- FALSE
        if (!eval.success & multicore) {
            if (!require(parallel) || .Platform$GUI!="X11") {
                warning("package 'parallel' not found or GUI is used, 
                analysis is performed in serial")
            } else {       
                trees <- mclapply(bstrees, optim.parsimony, data, trace = trace, method=method, ...)
                eval.success <- TRUE
            } 
        }
        if (!eval.success) trees <- lapply(bstrees, optim.parsimony, data, trace = trace, method=method, ...)
        if(class(result)=="phylo")m=1
        else m = length(result)
        if(m>0){
              trees[(np+1) : (np+m)] = result[1:m]
        }
        pscores <- sapply(trees, function(data) attr(data, "pscore"))

        mp = min(pscores)
        if (trace >= 0) 
            print(paste("Best pscore so far:",mp))
        ind = which(pscores < mp + eps)
        if (length(ind) == 1) {
            result = trees[ind]
            tree = result[[1]]
        }
        else {
            result = uniquetree(trees[ind])
            l = length(result)
            tree = result[[sample(l, 1)]]
        }
    }
    if(!all) return(tree)
    if(length(result)==1) return(result[[1]])
    class(result) = "multiPhylo"
    result
}


optim.sankoff <- function(tree, data, cost=NULL, trace=1, ...) {
    if(class(tree)!="phylo") stop("tree must be of class phylo") 
    if(is.rooted(tree))tree <- unroot(tree)
    if(is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") tree <- reorderPruning(tree)
    if (class(data)[1] != "phyDat") stop("data must be of class phyDat")
    
    rt = FALSE
    dat <- prepareDataSankoff(data)
    l <- attr(dat, "nc")
    if (is.null(cost)) {
        cost <- matrix(1, l, l)
        cost <- cost - diag(l)
        #       rt = TRUE
    }
    tree$edge.length=NULL
    swap = 0
    iter = TRUE
    pscore <- fit.sankoff(tree,dat,cost,'pscore')
    while (iter) {
        res <- sankoff.nni(tree,dat,cost,...)
        tree <- res$tree
        if(trace>1)cat("optimize topology: ", pscore , "-->", res$pscore, 
            "\n")
        pscore = res$pscore
        swap = swap + res$swap
        if (res$swap == 0) iter = FALSE
        }
    if(trace>0)cat("Final p-score",pscore,"after ",swap, "nni operations \n") 
    if(rt)tree <- ptree(tree, data)  
    attr(tree, "pscore") = pscore
    tree
}


#
# ACCTRAN
#
ptree <- function (tree, data, type = "ACCTRAN", retData = FALSE) 
{
    if (class(data) != "phyDat") 
        stop("data must be of class phyDat")
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- ape:::reorder.phylo(tree, "pruningwise") #reorderPruning??
 #   if (!is.binary.tree(tree)) 
 #       stop("Tree must be binary!")
    tmp = fitch(tree, data, site = "data")
    nr = attr(data, "nr")
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    weight = attr(data, "weight")
    m = length(edge) + 1
    q = length(tree$tip)
    l = as.integer(length(edge))
    nTips = length(tree$tip)
    dat = tmp[[2]]
    if (!is.rooted2(tree)) {
        root = getRoot(tree)
        ind = edge[node == root]
        rSeq = .C("fitchTriplet", integer(nr), dat[, ind[1]], 
            dat[, ind[2]], dat[, ind[3]], as.integer(nr))
        dat[, root] = rSeq[[1]]
    }
    result <- .C("ACCTRAN2", dat, as.integer(nr), numeric(nr), 
        as.integer(node[l:1L]), as.integer(edge[l:1L]), l, as.double(weight), 
        numeric(l), as.integer(nTips))
    el = result[[8]][l:1L]
    if (!is.rooted2(tree)) {
        ind2 = which(node[] == root)
        dat = matrix(result[[1]], nr, max(node))
        result <- .C("ACCTRAN3", result[[1]], as.integer(nr), 
            numeric(nr), as.integer(node[(l - 3L):1L]), as.integer(edge[(l - 
                3L):1L]), l - 3L, as.double(weight), numeric(l), 
            as.integer(nTips))
        el = result[[8]][(l - 3L):1L]
        pars = .C("fitchTripletACC4", dat[, root], dat[, ind[1]], 
            dat[, ind[2]], dat[, ind[3]], as.integer(nr), numeric(1), 
            numeric(1), numeric(1), as.double(weight), numeric(nr), 
            integer(nr))
        el[ind2[1]] = pars[[6]]
        el[ind2[2]] = pars[[7]]
        el[ind2[3]] = pars[[8]]
    }
    else {
        result <- .C("ACCTRAN3", result[[1]], as.integer(nr), 
            numeric(nr), as.integer(node[l:1L]), as.integer(edge[l:1L]), 
            l, as.double(weight), numeric(l), as.integer(nTips))
        el = result[[8]][l:1L]
    }
    tree$edge.length = el
    if (retData) 
        return(list(tree, matrix(result[[1]], nr, max(node))))
    tree
}


acctran <- function(tree, data) ptree(tree, data, type="ACCTRAN", retData=FALSE)


parsimony.plot <- function(tree, ...){
   x = numeric(max(tree$edge))
   x[tree$edge[,2]] = tree$edge.length 
   plot(tree, ...)
   ind <- get("last_plot.phylo", envir = .PlotPhyloEnv)$edge[, 2]
   edgelabels(prettyNum(x[ind]), frame = "none")
}



#
# Maximum likelihood estimation
#
discrete.gamma <- function (alpha, k) 
{
    if (k == 1) return(1)
    quants <- qgamma((1:(k - 1))/k, shape = alpha, rate = alpha)
    diff( c(0, pgamma(quants * alpha, alpha + 1),1)) * k
}


optimQ <- function (tree, data, Q=rep(1,6), subs=rep(1,length(Q)), trace = 0, ...) 
{
    m = length(Q)
    n = max(subs)
    ab = numeric(n)
    for(i in 1:n) ab[i]=log(Q[which(subs==i)[1]])
    fn = function(ab, tree, data, m, n, subs,...) {
        Q = numeric(m)
        for(i in 1:n)Q[subs==i] = ab[i]
        pml2(tree, data, Q = exp(Q),...)# Q^2, ...)
    }
    res = optim(par = ab, fn = fn, gr = NULL, method = "L-BFGS-B", 
        lower = -Inf, upper = Inf, control = list(fnscale = -1, 
        maxit = 25, trace = trace), tree = tree, data = data, m=m, n=n, subs=subs,...)
    Q = rep(1, m)
    for(i in 1:n) Q[subs==i] = exp(res[[1]][i])
    res[[1]] = Q
    res
}    

  
optimCodon <- function (tree, data, Q=rep(1,1830), subs=rep(1,length(Q)), syn = rep(0, length(Q)), trace = 0L, ab = c(0,0), optK=TRUE, optW=TRUE, ...) 
{
    m = length(Q)
    n = 1L # max(subs)

    fn = function(ab, tree, data, m, n, subs, syn, optK, optW, ...) {
        Q = numeric(m)
        Q[subs==1] = 0 # transversion
        if(optK) Q[subs==2] = ab[1] # transition
        else Q[subs==2] = 0
        if(optW) Q[syn==1] = Q[syn==1] + ab[2] # ab[n+1] dnds
        Q[syn<0] = -Inf
        phangorn:::pml2(tree, data, Q = exp(Q),...)# Q^2, ...)
    }
    res = optim(par = ab, fn = fn, gr = NULL, method = "L-BFGS-B", 
        lower = -Inf, upper = Inf, control = list(fnscale = -1, 
        maxit = 25, trace = trace), tree = tree, data = data, m=m, n=n, 
        subs=subs, syn=syn, optK=optK, optW=optW, ...)
    ab = exp(res[[1]])
    Q[subs==1] = 1 # transversion
    if(optK) Q[subs==2] = ab[1] # transition
    else{ 
        Q[subs==2] = 1
        ab[1] = 1 
        }  
    if(optW) Q[syn==1] = Q[syn==1] * ab[2] # dnds
    else ab[2] = 1
    Q[syn<0] = 0
    res[[5]] = ab
    res[[1]] = Q
    res
} 


subsChoice <- function(type=c("JC", "F81", "K80", "HKY", "TrNe", "TrN", "TPM1", "K81", "TPM1u", "TPM2", "TPM2u", "TPM3", "TPM3u", "TIM1e", "TIM1", "TIM2e", "TIM2", "TIM3e", "TIM3", "TVMe", "TVM", "SYM", "GTR")){
    type = match.arg(type)
    switch(type,
         JC = list(optQ=FALSE, optBf=FALSE,   subs=c(0, 0, 0, 0, 0, 0)),
         F81 = list(optQ=FALSE, optBf=TRUE,   subs=c(0, 0, 0, 0, 0, 0)),
         K80 = list(optQ=TRUE, optBf=FALSE,   subs=c(0, 1, 0, 0, 1, 0)),
         HKY = list(optQ=TRUE, optBf=TRUE,    subs=c(0, 1, 0, 0, 1, 0)),
         TrNe = list(optQ=TRUE, optBf=FALSE,  subs=c(0, 1, 0, 0, 2, 0)),
         TrN = list(optQ=TRUE, optBf=TRUE,    subs=c(0, 1, 0, 0, 2, 0)),
         TPM1 = list(optQ=TRUE, optBf=FALSE,  subs=c(0, 1, 2, 2, 1, 0)),
         K81 = list(optQ=TRUE, optBf=FALSE,   subs=c(0, 1, 2, 2, 1, 0)),
         TPM1u = list(optQ=TRUE, optBf=TRUE,  subs=c(0, 1, 2, 2, 1, 0)),
         TPM2 = list(optQ=TRUE, optBf=FALSE,  subs=c(1, 2, 1, 0, 2, 0)),
         TPM2u = list(optQ=TRUE, optBf=TRUE,  subs=c(1, 2, 1, 0, 2, 0)),
         TPM3 = list(optQ=TRUE, optBf=FALSE,  subs=c(1, 2, 0, 1, 2, 0)),
         TPM3u = list(optQ=TRUE, optBf=TRUE,  subs=c(1, 2, 0, 1, 2, 0)),
         TIM1e = list(optQ=TRUE, optBf=FALSE, subs=c(0, 1, 2, 2, 3, 0)),
         TIM1 = list(optQ=TRUE, optBf=TRUE,   subs=c(0, 1, 2, 2, 3, 0)),
         TIM2e = list(optQ=TRUE, optBf=FALSE, subs=c(1, 2, 1, 0, 3, 0)),
         TIM2 = list(optQ=TRUE, optBf=TRUE,   subs=c(1, 2, 1, 0, 3, 0)),
         TIM3e = list(optQ=TRUE, optBf=FALSE, subs=c(1, 2, 0, 1, 3, 0)),
         TIM3 = list(optQ=TRUE, optBf=TRUE,   subs=c(1, 2, 0, 1, 3, 0)),
         TVMe = list(optQ=TRUE, optBf=FALSE,  subs=c(1, 2, 3, 4, 2, 0)),
         TVM = list(optQ=TRUE, optBf=TRUE,    subs=c(1, 2, 3, 4, 2, 0)),
         SYM = list(optQ=TRUE, optBf=FALSE,   subs=c(1, 2, 3, 4, 5, 0)),
         GTR = list(optQ=TRUE, optBf=TRUE,    subs=c(1, 2, 3, 4, 5, 0))
         )
}


modelTest <- function (object, tree = NULL, model = c("JC", "F81", "K80", 
    "HKY", "SYM", "GTR"), G = TRUE, I = TRUE, k = 4, control = pml.control(epsilon = 1e-08, 
    maxit = 3, trace = 1), multicore = FALSE) 
{	
    if (class(object) == "phyDat") 
        data = object
    if (class(object) == "pml") {
        data = object$data
        if (is.null(tree)) 
            tree = object$tree
    }

    if(attr(data, "type")=="DNA") type = c("JC", "F81", "K80", "HKY", "TrNe", "TrN", "TPM1", 
        "K81", "TPM1u", "TPM2", "TPM2u", "TPM3", "TPM3u", "TIM1e", 
        "TIM1", "TIM2e", "TIM2", "TIM3e", "TIM3", "TVMe", "TVM", 
        "SYM", "GTR")
    if(attr(data, "type")=="AA") type = c("WAG", "JTT", "Dayhoff", "LG", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24")    
    model = match.arg(model, type, TRUE)

    env = new.env()
    assign("data", data, envir=env)

    if (is.null(tree)) 
        tree = NJ(dist.hamming(data))
    trace <- control$trace
    control$trace = trace - 1
    fit = pml(tree, data)
    fit = optim.pml(fit, control = control)
    l = length(model)
    n = 1L + sum(I + G + (G & I))
    nseq = sum(attr(data, "weight"))
    fitPar = function(model, fit, G, I, k) {
        m = 1
        res = matrix(NA, n, 5)
        res = as.data.frame(res)
        colnames(res) = c("Model", "df", "logLik", "AIC", "BIC")
        data.frame(c("Model", "df", "logLik", "AIC", "BIC"))
        calls = vector("list", n)
        trees = vector("list", n)
        fittmp = optim.pml(fit, model = model, control = control)
        res[m, 1] = model
        res[m, 2] = fittmp$df
        res[m, 3] = fittmp$logLik
        res[m, 4] = AIC(fittmp)
        res[m, 5] = AIC(fittmp, k = log(nseq))
        calls[[m]] = fittmp$call

        trees[[m]] = fittmp$tree
        m = m + 1
        if (I) {
            fitI = optim.pml(fittmp, model = model, optInv = TRUE, 
                control = control)
            res[m, 1] = paste(model, "+I", sep = "")
            res[m, 2] = fitI$df
            res[m, 3] = fitI$logLik
            res[m, 4] = AIC(fitI)
            res[m, 5] = AIC(fitI, k = log(nseq))
            calls[[m]] = fitI$call
            trees[[m]] = fitI$tree
            m = m + 1
        }
        if (G) {
            fitG = update(fittmp, k = k)
            fitG = optim.pml(fitG, model = model, optGamma = TRUE, 
                control = control)
            res[m, 1] = paste(model, "+G", sep = "")
            res[m, 2] = fitG$df
            res[m, 3] = fitG$logLik
            res[m, 4] = AIC(fitG)
            res[m, 5] = AIC(fitG, k = log(nseq))
            calls[[m]] = fitG$call
            trees[[m]] = fitG$tree
            m = m + 1
        }
        if (G & I) {
            fitGI = optim.pml(fitG, model = model, optGamma = TRUE, 
                optInv = TRUE, control = control)
            res[m, 1] = paste(model, "+G+I", sep = "")
            res[m, 2] = fitGI$df
            res[m, 3] = fitGI$logLik
            res[m, 4] = AIC(fitGI)
            res[m, 5] = AIC(fitGI, k = log(nseq))
            calls[[m]] = fitGI$call
            trees[[m]] = fitGI$tree
            m = m + 1
        }
        list(res, trees, calls)
    }
    eval.success <- FALSE
    if (!eval.success & multicore) {
        if (!require(parallel) || .Platform$GUI != "X11") {
            warning("package 'parallel' not found or GUI is used, \n      analysis is performed in serial")
        }
        else {
            RES <- mclapply(model, fitPar, fit, G, I, k)
            eval.success <- TRUE
        }
    }
    if (!eval.success) 
        res <- RES <- lapply(model, fitPar, fit, G, I, k)

    RESULT = matrix(NA, n * l, 5)
    RESULT = as.data.frame(RESULT)
    colnames(RESULT) = c("Model", "df", "logLik", "AIC", "BIC")
    for (i in 1:l) RESULT[((i - 1) * n + 1):(n * i), ] = RES[[i]][[1]]
    for(i in 1:l){
        for(j in 1:n){
            mo = RES[[i]][[1]][j,1]
            tname = paste("tree_", mo, sep = "")
            tmpmod = RES[[i]][[3]][[j]]
            tmpmod["tree"] = call(tname)
            if(!is.null(tmpmod[["k"]]))tmpmod["k"] = k
            if(attr(data, "type")=="AA") tmpmod["model"] = RES[[i]][[1]][1,1]          
    	    assign(tname, RES[[i]][[2]][[j]], envir=env)
            assign(mo, tmpmod, envir=env) 
        }
    }
    attr(RESULT, "env") = env 
    RESULT
}

  
optimGamma = function(tree, data, shape=1, k=4,...){
    fn = function(shape, tree, data, k,...)pml2(tree, data, shape=shape, k=k,...)
    res = optimize(f=fn, interval = c(0,100), lower = 0, upper = 100, maximum = TRUE,
        tol = .01, tree=tree, data=data, k=k,...)
    res
    }
    
 
optimInv = function(tree, data, inv=0.01, INV=NULL, ll.0=NULL,...){
    fn = function(inv, tree, data,...)pml2(tree, data, inv=inv, INV=INV, ll.0=NULL,...)
    res = optimize(f=fn, interval = c(0,1), lower = 0, upper = 1, maximum = TRUE,
         tol = .0001, tree=tree, data=data,...)
    res
    }
  

# changed to c(-10,10) from c(-5,5)
optimRate <- function(tree, data, rate=1, ...){
    fn <- function(rate, tree, data, ...) pml2(tree, data, rate=exp(rate), ...)
    res <- optimize(f = fn, interval = c(-10, 10), tree = tree, data = data, ..., maximum = TRUE)
    res[[1]] <- exp(res[[1]])
    res
}
    

optimBf = function(tree, data, bf=c(.25,.25,.25,.25), trace=0,...){
    l=length(bf)
    nenner = 1/bf[l]
    lbf = log(bf * nenner)
    lbf = lbf[-l]
    fn = function(lbf, tree, data,...){
        bf = exp(c(lbf,0))
        bf = bf/sum(bf)
        pml2(tree, data, bf=bf, ...)
        }
    res = optim(par=lbf, fn=fn, gr=NULL, method="Nelder-Mead", control=list(fnscale=-1, maxit=500, trace=trace),tree=tree, data=data,...)
    bf = exp(c(res[[1]],0))
    bf = bf/sum(bf)
    result = list(bf=bf, loglik = res[[2]])
    result
    }


optimW = function(fit,...){
    w = fit$w
    g = fit$g
    siteLik = fit$siteLik
    k = length(w)
    l = dim(siteLik[[1]])[1]
    x=matrix(0,l,k)
    for(i in 1:k)x[,i] = rowSums(siteLik[[i]])
    weight = fit$weight
    nenner = 1/w[k]
    eta = log(w * nenner)
    eta = eta[-k]
    fn = function(eta,x,g,weight){
        eta = c(eta,0)
        p = exp(eta)/sum(exp(eta))
        res = x%*%p
        res = sum(weight*log(res))  * (1 + abs(sum(p*g) - 1))
        res
    }  
    res = optim(eta, fn = fn, method = "Nelder-Mead", control=list(fnscale=-1, reltol = 1e-12),gr=NULL, x=x,g=g, weight=weight)
    p = exp(c(res$par,0))
    p = p/sum(p)
    result = list(par = p, value = res$value)
    result    
}


#predict.pml <- function(object, newdata,...) sum(object$site * newdata)


logLik.pml <- function(object,...){
    res <- object$logLik
    attr(res,"df") <- object$df
    class(res) <- "logLik"
    res
}

anova.pml <- function (object, ...) 
{
    X <- c(list(object), list(...))
    df <- sapply(X, "[[", "df")
    ll <- sapply(X, "[[", "logLik")
    dev <- c(NA, 2 * diff(ll)) 
    ddf <- c(NA, diff(df))
    table <- data.frame(ll, df, ddf, dev, pchisq(dev, ddf, lower.tail = FALSE))
    dimnames(table) <- list(1:length(X), c("Log lik.", "Df", 
        "Df change", "Diff log lik.", "Pr(>|Chi|)"))
    structure(table, heading = "Likelihood Ratio Test Table", 
        class = c("anova", "data.frame"))
}
    
    
#vcov.pml <- function(object, obs=FALSE,...){
#    if(obs) FI = score4(object)[[2]]
#    else FI = score(object,FALSE)[[2]]
#    l = dim(FI)[1]
#    res = try(solve(FI))
#    if(class(res) == "try-error"){
#        cat("Covariance is ill-conditioned !! \n")
#        res = solve(FI + diag(l)* 1e-8)
#        }
#    res
#}
                             
vcov.pml <- function(object, ...){
    FI = score(object,FALSE)[[2]]
    l = dim(FI)[1]
    res = try(solve(FI))
    if(class(res) == "try-error"){
        cat("Covariance is ill-conditioned !! \n")
        res = solve(FI + diag(l)* 1e-8)
        }
    res
}


getd2P <- function(el, eig=edQt(), g=1.0){
    n <- length(eig$values)    
    res <- .Call("getd2PM",eig,as.integer(n),as.double(el),as.double(g), PACKAGE = "phangorn")
    attr(res,"dim") <- c(length(g),length(el))
    res
}


getdP <- function(el, eig=edQt(), g=1.0){
    n <- length(eig$values)    
    res <- .Call("getdPM",eig,as.integer(n),as.double(el),as.double(g), PACKAGE = "phangorn")
    attr(res,"dim") <- c(length(g),length(el))
    res
}


# version without transformation (used for vcov)
getdP2 <- function(el, eig=edQt(), g=1.0){
    n <- length(eig$values)    
    res <- .Call("getdPM2",eig,as.integer(n),as.double(el),as.double(g), PACKAGE = "phangorn")
    attr(res,"dim") <- c(length(g),length(el))
    res
}


# version without transformation 
getd2P2 <- function(el, eig=edQt(), g=1.0){
    n <- length(eig$values)    
    res <- .Call("getd2PM2",eig,as.integer(n),as.double(el),as.double(g), PACKAGE = "phangorn")
    attr(res,"dim") <- c(length(g),length(el))
    res
}


getP <- function(el, eig=edQt(), g=1.0){
    #if(el<0)stop("need positive edge length")
    n <- length(eig$values)
    res <- .Call("getPM", eig, as.integer(n), as.double(el), as.double(g), PACKAGE = "phangorn")
    attr(res, "dim") <- c(length(g), length(el)) 
    res
}


lli = function (data, tree, ...) 
{
    contrast = attr(data, "contrast")
    nr = attr(data, "nr")
    nc = attr(data, "nc")
    nco = as.integer(dim(contrast)[1])
    .Call("invSites", data[tree$tip.label], as.integer(nr), as.integer(nc), contrast, as.integer(nco), PACKAGE = "phangorn")    
}


edQt <- function (Q = c(1, 1, 1, 1, 1, 1), bf = c(0.25, 0.25, 0.25, 0.25)) 
{
    l = length(bf)
    res = matrix(0, l, l)
    res[lower.tri(res)] = Q
    res = res + t(res)
    res = res * bf
    res2 = res * rep(bf, each = l)    
    diag(res) = -colSums(res)
    res = res/sum(res2)
    e = eigen(res, FALSE)
    e$inv = solve(e$vec)
    e
}


edQ <- function(Q=c(1,1,1,1,1,1), bf=c(0.25,.25,.25,.25)){
    l=length(bf)
    res = matrix(0, l, l)
    res[lower.tri(res)] = Q
    res = res+t(res)
    res = res * rep(bf,each=l)
    diag(res) = -rowSums(res)
    res2 = res * rep(bf,l)
    diag(res2)=0 
    res = res/sum(res2)
    e = eigen(res, FALSE)
    e$inv = solve(e$vec)
    e
}


# this should replace pml
pmlScale <- function (tree, data, bf = NULL, Q = NULL, inv = 0, k = 1, shape = 1, 
    rate = 1, model = NULL, ret="pml", ...) 
{
    call <- match.call()
    extras <- match.call(expand.dots = FALSE)$...
    pmla <- c("wMix", "llMix")
    existing <- match(pmla, names(extras))
    ret  <- match.arg(ret, c("pml", "logLik"))
    wMix <- ifelse(is.na(existing[1]), 0, eval(extras[[existing[1]]], 
        parent.frame()))
    llMix <- ifelse(is.na(existing[2]), 0, eval(extras[[existing[2]]], 
        parent.frame()))
    if (class(tree) != "phylo") 
        stop("tree must be of class phylo")
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorderPruning(tree)
    if (class(data)[1] != "phyDat") 
        stop("data must be of class phyDat")
    if (is.null(tree$edge.length)) 
        stop("tree must have edge weights")
    if (any(is.na(match(tree$tip, attr(data, "names"))))) 
        stop("tip labels are not in data")
    levels <- attr(data, "levels")
    weight <- attr(data, "weight")
    nr <- attr(data, "nr")
    type <- attr(data, "type")
    if (type == "AA" & !is.null(model)) {
        model <- match.arg(model, c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24"))
        getModelAA(model, bf=is.null(bf), Q=is.null(Q))
    }
    if (is.null(bf)) 
        bf <- rep(1/length(levels), length(levels))
    if (is.null(Q)) 
        Q <- rep(1, length(levels) * (length(levels) - 1)/2)
    m <- 1
    eig <- edQt(bf = bf, Q = Q)
    w <- rep(1/k, k)
    if (inv > 0) 
        w <- (1 - inv) * w
    if (wMix > 0) 
        w <- wMix * w
    g <- discrete.gamma(shape, k)
    if (inv > 0) 
        g <- g/(1 - inv)
    g <- rate * g
    INV <- Matrix(lli(data, tree), sparse=TRUE)
    ll.0 <- as.matrix(INV %*% (bf * inv))
    if (wMix > 0) 
        ll.0 <- ll.0 + llMix
    resll <- matrix(0, nr, k)
    lll <- matrix(0, nr, k)
    while (m <= k) {
        tmp = ll3(data, tree, bf = bf, g = g[m], Q = Q, eig = eig, assign.dat = FALSE, ...)
        resll[, m] = tmp[[1]]
        lll[, m] = tmp[[2]] 
        m = m + 1
    }
    sca = .Call("rowMax", resll, length(weight), k) # + 1 
    resll = resll - sca 

    lll <- exp(resll + log(lll)) 
    siteLik <- (lll%*%w)
    siteLik <- log(siteLik) + sca

    loglik = sum(weight * siteLik)
    if(ret == "logLik")return(loglik) 
    df = length(tree$edge.length) + (k > 1) + (inv > 0) + length(unique(bf)) - 
        1 + length(unique(Q)) - 1
    result = list(logLik = loglik, inv = inv, k = k, shape = shape, 
        Q = Q, bf = bf, rate = rate, siteLik = siteLik, weight = weight, 
        g = g, w = w, eig = eig, data = data, model = model, 
        INV = INV, ll.0 = ll.0, tree = tree, lv = resll, call = call, 
        df = df, wMix = wMix, llMix = llMix)
    class(result) = "pml"
    result
}


ll <- function (dat1, tree, bf = c(0.25, 0.25, 0.25, 0.25), g = 1, 
    Q = c(1, 1, 1, 1, 1, 1), eig = NULL, assign.dat = FALSE, ...) 
{
    q = length(tree$tip.label) 
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    m = length(edge) + 1
    if (is.null(eig)) eig = edQt(bf = bf, Q = Q)
    el <- tree$edge.length
    P <- getP(el, eig, g)  
    nr <- as.integer(attr(dat1,"nr"))   
    nc <- as.integer(attr(dat1,"nc"))
    node = as.integer(node-min(node))
    edge = as.integer(edge-1L)
    nTips = as.integer(length(tree$tip))
    mNodes = as.integer(max(node) + 1L)
    contrast = attr(dat1, "contrast")
    nco = as.integer(dim(contrast)[1])
    res <- .Call("LogLik2", dat1[tree$tip.label], P, nr, nc, node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")
    result = res[[1]] %*% bf  # root statt 1
    if (assign.dat){
        dat = vector(mode = "list", length = m)
        dat[(q+1):m] <- res
        attr(dat, "names") = c(tree$tip.label, as.character((q + 1):m))
        assign("asdf", dat, envir = parent.frame(n = 1))
        }
    result
}


# scaled version, needs more speeding up  
ll3 <- function (dat1, tree, bf = c(0.25, 0.25, 0.25, 0.25), g = 1, 
    Q = c(1, 1, 1, 1, 1, 1), eig = NULL, assign.dat = FALSE, 
    ...) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorderPruning(tree)
    q = length(tree$tip.label)
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    m = length(edge) + 1
    if (is.null(eig)) 
        eig = edQt(bf = bf, Q = Q)
    el <- tree$edge.length
    P <- getP(el, eig, g)
    nr <- as.integer(attr(dat1, "nr"))
    nc <- as.integer(attr(dat1, "nc"))
    node = as.integer(node - min(node))
    edge = as.integer(edge - 1)
    nTips = as.integer(length(tree$tip))
    mNodes = as.integer(max(node) + 1)
    contrast = attr(dat1, "contrast")
    nco = as.integer(dim(contrast)[1])
    res <- .Call("LogLik4", dat1[tree$tip.label], P, nr, nc, node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")
   
    if (assign.dat) {
        dat = vector(mode = "list", length = m)
        dat[(q + 1):m] <- res
        attr(dat, "names") = c(tree$tip.label, as.character((q + 
            1):m))
        assign("asdf", dat, envir = parent.frame(n = 1))
    }
    list(res[[2]][[1]],res[[1]][[1]] %*% bf)
}


# raus
fn.quartet <- function(old.el, eig, bf, dat,  g=1, w=1, weight, ll.0) {
    l= length(dat[,1]) 
    ll = ll.0
    res = vector("list", 2*l)
    tmp1 = NULL
    tmp2 = NULL
    attr(res,"dim") = c(l,2)
    for(j in 1:l){
            P = getP(old.el, eig, g[j])
            tmp1 = (dat[[j,1]] %*% P[[1]]) *(dat[[j,2]] %*% P[[2]])
            tmp2 = (dat[[j,3]] %*% P[[3]]) * (dat[[j,4]] %*% P[[4]])
            res[[j,1]] = tmp1 * (tmp2 %*% P[[5]])
            res[[j,2]] = tmp2
            ll = ll +  res[[j,1]] %*% (w[j]*bf)
        } 
    l0 = sum(weight * log(ll))
    list(ll=l0,res=res)
}


fn.quartet2 <- function (old.el, eig, bf, dat1, dat2, dat3, dat4, g = 1, w = 1, 
    weight, ll.0, contrast, ext) 
{
    l = length(w)
    ll = ll.0
    res = vector("list", 2 * l)
    tmp1 = NULL
    tmp2 = NULL
    attr(res, "dim") = c(l, 2)
    for (j in 1:l) {
        P = getP(old.el, eig, g[j])
        if (ext[1] == FALSE && ext[2] == FALSE) 
            tmp1 = (dat1[[j]] %*% P[[1]]) * (dat2[[j]] %*% P[[2]])
        if (ext[1] == FALSE && ext[2] == TRUE) 
            tmp1 = (dat1[[j]] %*% P[[1]]) * (contrast %*% P[[2]])[dat2, ]
        if (ext[1] == TRUE && ext[2] == FALSE) 
            tmp1 = (contrast %*% P[[1]])[dat1, ] * (dat2[[j]] %*% P[[2]])
        if (ext[1] == TRUE && ext[2] == TRUE) 
            tmp1 = (contrast %*% P[[1]])[dat1, ] * (contrast %*% P[[2]])[dat2, ]
        if (ext[3] == FALSE && ext[4] == FALSE) 
            tmp2 = (dat3[[j]] %*% P[[3]]) * (dat4[[j]] %*% P[[4]])
        if (ext[3] == FALSE && ext[4] == TRUE) 
            tmp2 = (dat3[[j]] %*% P[[3]]) * (contrast %*% P[[4]])[dat4, ]
        if (ext[3] == TRUE && ext[4] == FALSE) 
            tmp2 = (contrast %*% P[[3]])[dat3, ] * (dat4[[j]] %*% P[[4]])
        if (ext[3] == TRUE && ext[4] == TRUE) 
            tmp2 = (contrast %*% P[[3]])[dat3, ] * (contrast %*% P[[4]])[dat4, ]
        res[[j, 1]] = tmp1 * (tmp2 %*% P[[5]])
        res[[j, 2]] = tmp2
        ll = ll + res[[j, 1]] %*% (w[j] * bf)
    }
    l0 = sum(weight * log(ll))
    list(ll = l0, res = res)
}


optim.quartet2 <- function (old.el, eig, bf, dat1, dat2, dat3, dat4, g = 1, w = 1, 
    weight, ll.0 = weight * 0, control = list(eps = 1e-08, maxit = 5, 
        trace = 0), llcomp = -Inf, evi, contrast, contrast2, 
    ext = c(FALSE, FALSE, FALSE, FALSE)) 
{
    eps = 1
    iter = 0
    while (eps > control$eps && iter < control$maxit) {
        tmp <- fn.quartet2(old.el = old.el, eig = eig, bf = bf, 
            dat1 = dat1, dat2 = dat2, dat3 = dat3, dat4 = dat4, 
            g = g, w = w, weight = weight, ll.0 = ll.0, contrast=contrast, ext = ext)
        old.ll = tmp$ll
   
        el1 <- fs3(old.el[1], eig, tmp$res[, 1], dat1, weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, contrast=contrast, contrast2=contrast2, evi=evi, ext = ext[1], getA=TRUE, getB=FALSE)
        el2 <- fs3(old.el[2], eig, el1[[2]], dat2, weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, contrast=contrast, contrast2=contrast2, evi=evi, ext = ext[2], getA=TRUE, getB=FALSE)
        el5 <- fs3(old.el[5], eig, el2[[2]], tmp$res[, 2], weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, contrast=contrast, contrast2=contrast2, evi=evi, ext = 0L, getA=FALSE, getB=TRUE)
        el3 <- fs3(old.el[3], eig, el5[[3]], dat3, weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, contrast=contrast, contrast2=contrast2, evi=evi, ext = ext[3], getA=TRUE, getB=FALSE)
        el4 <- fs3(old.el[4], eig, el3[[2]], dat4, weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, contrast=contrast, contrast2=contrast2, evi=evi, ext = ext[4], getA=FALSE, getB=FALSE)
        old.el[1] = el1[[1]]
        old.el[2] = el2[[1]]
        old.el[3] = el3[[1]]
        old.el[4] = el4[[1]]
        old.el[5] = el5[[1]]
        iter = iter + 1
        ll = el4[[4]]
        eps = (old.ll - ll)/ll
        if (ll < llcomp) 
            return(list(old.el, ll))
        old.ll = ll
    }
    list(old.el, ll)
}



pml.nni <- function (tree, data, w, g, eig, bf, ll.0, ll, ...) 
{        
    k = length(w)
    INDEX <-  indexNNI(tree)
    rootEdges <- attr(INDEX,"root")
    .dat <- NULL

    data = getCols(data, tree$tip)

    parent = tree$edge[,1]
    child = tree$edge[,2]
    weight = attr(data, "weight")
#    datp = rnodes2(tree, data, w, g, eig, bf) #eval('rnodes', fit) eval(rnodes, object, parent.frame())
    datp = rnodes3(tree, data, w, g, eig, bf) 
# rnodes2 ersetzen    
    contrast <- attr(data, "contrast")
    contrast2 <- contrast %*% eig[[2]] 
    evi = (t(eig[[3]]) * bf)

    nTips = length(tree$tip.label)
    evector <- numeric(max(parent)) 
    evector[child] <- tree$edge.length
    m <- dim(INDEX)[1]
    loglik = numeric(2*m)
    edgeMatrix <- matrix(0, 2*m, 5)
    l = length(datp[, 1])
    for(i in 1:m){
        ei = INDEX[i,]
        el0 = evector[INDEX[i,]]
 #       datn = vector("list", 4 * l)
 #       attr(datn, "dim") = c(l, 4)
 #       datn <- .dat[, ei[1:4], drop = FALSE]

        ext = ei[1:4] < nTips+1L
        if (!(ei[5] %in% rootEdges)) dat1 = datp[, ei[1], drop = FALSE]
        else{ if(ext[1]) dat1 = data[[ ei[1] ]]
             else dat1 = .dat[, ei[1], drop=FALSE]
        } 
        if(ext[2]) dat2 = data[[ ei[2] ]]
             else dat2 = .dat[, ei[2], drop=FALSE] 
        if(ext[3]) dat3 = data[[ ei[3] ]]
             else dat3 = .dat[, ei[3], drop=FALSE]
        if(ext[4]) dat4 = data[[ ei[4] ]]
             else dat4 = .dat[, ei[4], drop=FALSE]

        new1 <- optim.quartet2(el0[c(1, 3, 2, 4, 5)], eig, bf, 
            dat1, dat3, dat2, dat4, g, w, weight, ll.0, llcomp=ll, evi=evi, contrast=contrast, contrast2=contrast2, ext=ext[c(1, 3, 2, 4)])
        new2 <- optim.quartet2(el0[c(1, 4, 3, 2, 5)], eig, bf,  
            dat1, dat4, dat3, dat2, g, w, weight, ll.0, llcomp=ll, evi=evi, contrast=contrast, contrast2=contrast2, ext=ext[c(1, 4, 3, 2)])


        loglik[(2*i)-1]=new1[[2]]
        loglik[(2*i)]=new2[[2]] 
        edgeMatrix[(2*i)-1,]=new1[[1]]
        edgeMatrix[(2*i),]=new2[[1]]           
    }
    swap <- 0
    eps0 <- 1e-6
    candidates <- loglik > ll + eps0
    while(any(candidates)){     
        ind = which.max(loglik)
        loglik[ind]=-Inf
        if( ind %% 2 ) swap.edge = c(2,3)
        else swap.edge = c(2,4)
        tree2 <- changeEdge(tree, INDEX[(ind+1)%/%2,swap.edge], INDEX[(ind+1)%/%2,], edgeMatrix[ind,])
        
        test <- pml2(tree2, data, bf = bf, k=k, g=g, w=w, eig=eig, ll.0=ll.0, ...) # maybe INV
        if(test <= ll + eps0) candidates[ind] = FALSE
        if(test > ll + eps0) {
            ll = test 
            swap=swap+1
            tree <- tree2
            indi <- which(rep(colSums(apply(INDEX,1,match,INDEX[(ind+1)%/%2,],nomatch=0))>0,each=2))
            candidates[indi] <- FALSE
            loglik[indi] <- -Inf
        }
    } 
    list(tree=tree, ll=ll, swap=swap)     
}


rnodes <- function (fit)  
{
    tree = fit$tree 
    data = getCols(fit$data, tree$tip) 
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorderPruning(tree)
    q = length(tree$tip.label)
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    m = length(edge) + 1  # max(edge)
    w = fit$w
    g = fit$g
    l = length(w)        
    dat = vector(mode = "list", length = m*l)
    dim(dat) <- c(l,m)
    
    tmp = length(data)
    for(i in 1:length(w))dat[i,1:tmp]=new2old.phyDat(data) 
    
    eig = fit$eig

    bf = fit$bf
    el <- tree$edge.length
    P <- getP(el, eig, g)
    nr <- as.integer(attr(data, "nr"))
    nc <- as.integer(attr(data, "nc"))
    node = as.integer(node - min(node))
    edge = as.integer(edge - 1)
    nTips = as.integer(length(tree$tip))
    mNodes = as.integer(max(node) + 1)
    contrast = attr(data, "contrast")
    nco = as.integer(dim(contrast)[1])
    for(i in 1:l)dat[i,(q + 1):m] <- .Call("LogLik2", data, P[i,], nr, nc, node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")

    parent <- tree$edge[, 1]
    child <- tree$edge[, 2]
    nTips = min(parent) - 1
    datp = vector("list", m)   
    dat2 = vector("list", m * l)

    dim(dat2) <- c(l,m)

    for(i in 1:l){     
      datp[(nTips + 1)] = dat[i,(nTips + 1)]
      for (j in (m - 1):1) {
          if (child[j] > nTips){
             tmp2 = (datp[[parent[j]]]/(dat[[i,child[j]]] %*% P[[i,j]]))
             datp[[child[j]]] = (tmp2 %*% P[[i,j]]) * dat[[i,child[j]]]  
             dat2[[i, child[j]]] = tmp2
             }
       }
    }
    assign(".dat", dat, envir = parent.frame(n = 1))
    dat2
}


rnodes3 <- function (tree, data, w, g, eig, bf) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorderPruning(tree)
    data = getCols(data, tree$tip) 
    q = length(tree$tip.label)
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    m = length(edge) + 1  # max(edge)
    l = length(w)        
    dat = vector(mode = "list", length = m*l)
    dim(dat) <- c(l,m)
    tmp = length(data)
#    for(i in 1:length(w))dat[i,1:tmp]=new2old.phyDat(data) #
#    dat[1,1:tmp] <- data  vielleicht gebraucht
    el <- tree$edge.length
    P <- getP(el, eig, g)
    nr <- as.integer(attr(data, "nr"))
    nc <- as.integer(attr(data, "nc"))
    node = as.integer(node - min(node))
    edge = as.integer(edge - 1)
    nTips = as.integer(length(tree$tip))
    mNodes = as.integer(max(node) + 1)
    contrast = attr(data, "contrast")
    nco = as.integer(dim(contrast)[1])
    for(i in 1:l)dat[i,(q + 1):m] <- .Call("LogLik2", data, P[i,], nr, nc, node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")
    parent <- tree$edge[, 1]
    child <- tree$edge[, 2]
    nTips = min(parent) - 1
    datp = vector("list", m)   
    dat2 = vector("list", m * l)
    dim(dat2) <- c(l,m)
# prep???
    for(i in 1:l){     
      datp[(nTips + 1)] = dat[i,(nTips + 1)]
      for (j in (m - 1):1) {
          if (child[j] > nTips){
             tmp2 = (datp[[parent[j]]]/(dat[[i,child[j]]] %*% P[[i,j]]))
             datp[[child[j]]] = (tmp2 %*% P[[i,j]]) * dat[[i,child[j]]]  
             dat2[[i, child[j]]] = tmp2
             }
       }
    }
    assign(".dat", dat, envir = parent.frame(n = 1))
    dat2
}


score <- function (fit, transform=TRUE) 
{
    tree = fit$tree
    child <- tree$edge[, 2]
    l = length(child)
    sc = numeric(l)
    weight = as.numeric(fit$weight)
    f <- drop(exp(fit$site))
    dl = dl(fit, transform)
    dl = dl/f
    sc = colSums(weight * dl)
    F = crossprod(dl*weight,dl) 
    names(sc) = child
    dimnames(F) = list(child, child) 
    result = list(sc = sc, F = F)
    result
}


# wird noch in partition models verwendet
optim.quartet <- function (old.el, eig, bf, dat, g = 1, w = 1, weight, ll.0 = weight * 
    0, control = list(eps = 1e-08, maxit = 5, trace = 0), llcomp=-Inf) 
{
    eps = 1
    iter = 0
    evi = (t(eig[[3]]) * bf)
    while (eps > control$eps && iter < control$maxit) {
        tmp <- fn.quartet(old.el = old.el, eig = eig, bf = bf, dat = dat, 
            g = g, w = w, weight = weight, ll.0 = ll.0)
        old.ll = tmp$ll 
        el1 <- fs(old.el[1], eig, tmp$res[, 1], dat[, 1], weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, evi, getA=TRUE, getB=FALSE)
        el2 <- fs(old.el[2], eig, el1[[2]], dat[, 2], weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, evi, getA=TRUE, getB=FALSE)
        el5 <- fs(old.el[5], eig, el2[[2]], tmp$res[, 2], weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, evi, getA=FALSE, getB=TRUE)
        el3 <- fs(old.el[3], eig, el5[[3]], dat[, 3], weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, evi, getA=TRUE, getB=FALSE)
        el4 <- fs(old.el[4], eig, el3[[2]], dat[, 4], weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, evi, getA=FALSE, getB=FALSE)
        old.el[1] = el1[[1]]
        old.el[2] = el2[[1]]
        old.el[3] = el3[[1]]
        old.el[4] = el4[[1]]
        old.el[5] = el5[[1]]
        iter = iter + 1
        ll = el4[[4]]
        eps = (old.ll - ll) / ll
        if(ll<llcomp)return(list(old.el, ll))  
        old.ll = ll
    }
    list(old.el, ll)
}


optim.quartet2 <- function (old.el, eig, bf, dat1, dat2, dat3, dat4, g = 1, w = 1, weight, ll.0 = weight * 
    0, control = list(eps = 1e-08, maxit = 5, trace = 0), llcomp=-Inf, evi, contrast, contrast2, ext=c(FALSE, FALSE, FALSE, FALSE)) 
{
    eps = 1
    iter = 0
    while (eps > control$eps && iter < control$maxit) {
        tmp <- fn.quartet2(old.el = old.el, eig = eig, bf = bf, dat1 = dat1, dat2 = dat2, dat3 = dat3, dat4 = dat4,
            g = g, w = w, weight = weight, ll.0 = ll.0, contrast=contrast, ext=ext)
        old.ll = tmp$ll 
        el1 <- fs3(old.el[1], eig, tmp$res[, 1], dat1, weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, contrast=contrast, contrast2=contrast2, evi=evi, ext = ext[1], getA=TRUE, getB=FALSE)
        el2 <- fs3(old.el[2], eig, el1[[2]], dat2, weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, contrast=contrast, contrast2=contrast2, evi=evi, ext = ext[2], getA=TRUE, getB=FALSE)
        el5 <- fs3(old.el[5], eig, el2[[2]], tmp$res[, 2], weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, contrast=contrast, contrast2=contrast2, evi=evi, ext = 0L, getA=FALSE, getB=TRUE)
        el3 <- fs3(old.el[3], eig, el5[[3]], dat3, weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, contrast=contrast, contrast2=contrast2, evi=evi, ext = ext[3], getA=TRUE, getB=FALSE)
        el4 <- fs3(old.el[4], eig, el3[[2]], dat4, weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, contrast=contrast, contrast2=contrast2, evi=evi, ext = ext[4], getA=FALSE, getB=FALSE)
        old.el[1] = el1[[1]]
        old.el[2] = el2[[1]]
        old.el[3] = el3[[1]]
        old.el[4] = el4[[1]]
        old.el[5] = el5[[1]]
        iter = iter + 1
        ll = el4[[4]]
        eps = (old.ll - ll) / ll
        if(ll<llcomp)return(list(old.el, ll))  
        old.ll = ll
    }
    list(old.el, ll)
}


plot.pml<-function(x,...)plot.phylo(x$tree,...)


phangornParseFormula <- function(model){

    parseSide <- function(model) {
        model.vars <- list()
        while (length(model) == 3 && model[[1]] == as.name("+")) {
            model.vars <- c(model.vars, model[[3]])
            model <- model[[2]]
        }
        unlist(rev(c(model.vars, model)))

    } 

    if (!inherits(model, "formula")) 
        stop("model must be a formula object")
    l <- length(model)
    varsLHS <- NULL       
    if(l==3){        
        modelLHS <- model[[2]]
        modelRHS <- model[[3]]
        varsRHS <- parseSide(modelRHS)
        varsRHS <- unlist(lapply(varsRHS,as.character))
        varsLHS <- parseSide(modelLHS)
        varsLHS <- unlist(lapply(varsLHS,as.character))
    }
    if(l==2){
       modelRHS <- model[[2]]
       varsRHS <- parseSide(modelRHS)
       varsRHS <- unlist(lapply(varsRHS,as.character))
    }
    list(left=varsLHS, right=varsRHS)
}


pml.control <- function (epsilon = 1e-08, maxit = 10, trace = 1) 
{
    if (!is.numeric(epsilon) || epsilon <= 0) 
        stop("value of 'epsilon' must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    list(epsilon = epsilon, maxit = maxit, trace = trace)
}


optim.pml <- function (object, optNni = FALSE, optBf = FALSE, optQ = FALSE, 
    optInv = FALSE, optGamma = FALSE, optEdge = TRUE, optRate = FALSE, optRooted=FALSE, 
    control = pml.control(epsilon = 1e-8, maxit = 10, trace = 1L), 
    model = NULL, subs = NULL, ...) 
{
    extras <- match.call(expand.dots = FALSE)$...
    pmla <- c("wMix", "llMix")
    wMix <- object$wMix
    llMix <- object$llMix
    if(is.null(llMix)) llMix=0
    if (!is.null(extras)) {
        names(extras) <- pmla[pmatch(names(extras), pmla)]
        existing <- match(pmla, names(extras))
        if (!is.na(existing[1])) 
            wMix <- eval(extras[[existing[1]]], parent.frame())
        if (!is.na(existing[2])) 
            llMix <- eval(extras[[existing[2]]], parent.frame())
    }
    tree = object$tree
    call = object$call
    if(optNni) {
        if(!is.binary.tree(tree)) 
            tree = multi2di(tree)
        optEdge = TRUE     
    }
    if(is.rooted(tree)) {
        if(optRooted==FALSE && optEdge==TRUE){
            tree = unroot(tree)
            tree = reorderPruning(tree) # may in wrong order for pml when using reorder
            warning("I unrooted the tree (rooted trees are not yet supported)", 
                call. = FALSE)
        }    
    }
    if(is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorderPruning(tree)
    if(any(tree$edge.length < 1e-08)) {
        tree$edge.length[tree$edge.length < 1e-08] <- 1e-08
        object <- update.pml(object, tree = tree)
    }
    if(optEdge & optRate) {
        warning("You can't optimise edges and rates at the same time, only edges are optimised!", call. = FALSE)
        optRate = FALSE
    }
    if(optRooted){
	    if(optNni)warning("Cannot perform tree rearrangements for rooted trees!", call. = FALSE)
        optNni = FALSE 
        optEdge = FALSE
        if(!is.rooted(tree)) stop("Tree must be rooted!")
        if(!is.ultrametric(tree)) stop("Tree must be ultrametric!")
	}
    trace <- control$trace
    
    data = object$data
    type <- attr(data, "type")
    if (type == "AA" & !is.null(model)){
        object = update(object, model=model)  
#        getModelAA(model)
    }     
    if (type == "CODON") {
        dnds <- object$dnds 
        tstv <- object$tstv
        if(!is.null(model)){
            if(model == "codon0") optQ = FALSE
            else  optQ = TRUE
        }
    }       
    Q = object$Q
    if(is.null(subs)) subs = c(1:(length(Q) - 1), 0)
    bf = object$bf
    eig = object$eig
    inv = object$inv
    k = object$k
    if(k==1 & optGamma){
        optGamma = FALSE
        warning('only one rate class, ignored optGamma')
    }
    shape = object$shape
    w = object$w
    g = object$g
    if (type == "DNA" & !is.null(model)) {
        tmp = subsChoice(model)
        optQ = tmp$optQ
        if (!optQ) 
            Q = rep(1, 6)
        optBf = tmp$optBf
        if (!optBf) 
            bf = c(0.25, 0.25, 0.25, 0.25)
        subs = tmp$subs
    }   
    ll0 <- object$logLik
    INV <- object$INV
    ll.0 <- object$ll.0
    rate <- object$rate
    ll = ll0
    ll1 = ll0
    opti = TRUE
    if (optEdge) {
         res <- optimEdge(tree, data, eig=eig, Q=Q, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0,
              control = pml.control(epsilon = 1e-07, maxit = 5, trace=trace - 1)) 
        if (trace > 0) 
            cat("optimize edge weights: ", ll, "-->", res[[2]], "\n")
        if (res[[2]] > ll){  
           ll <- res[[2]]
           tree <- res[[1]]
        }
    }
    if(optRooted){
	    res <- optimRooted(object, control = pml.control(epsilon = 1e-07, maxit = 25, trace = trace-1))
	    if (res[[2]] > ll){  
           ll <- res[[2]]
           tree <- res[[1]]
        }     
	}
    rounds = 1
    while (opti) {
        if (optBf) {
            res = optimBf(tree, data, bf = bf, inv = inv, Q = Q, 
                w = w, g = g, INV = INV, rate = rate, k = k, 
                llMix = llMix)
            bf = res[[1]]
            eig = edQt(Q = Q, bf = bf)
            if (inv > 0) 
                ll.0 <- as.matrix(INV %*% (bf * inv))
            if (wMix > 0) 
                ll.0 <- ll.0 + llMix
            if (trace > 0) 
                cat("optimize base frequencies: ", ll, "-->", 
                  res[[2]], "\n")
            ll = res[[2]]
        }
        if (optQ) {
            if(type=="CODON"){
                 if(is.null(model)) model <- "codon1"
                 model <- match.arg(model, c("codon0", "codon1", "codon2", "codon3"))
                 ab <- c(tstv, dnds)
                 res <- switch(model, 
                     codon1 = optimCodon(tree,data, Q=rep(1,1830), subs=.sub, syn=.syn, 
                         bf = bf, w = w, g = g, inv = inv, INV = INV, ll.0 = ll.0, rate = rate, k = k, ab=log(ab),
                         optK=TRUE, optW = TRUE),  
                     codon2 = optimCodon(tree,data, Q=rep(1,1830), subs=.sub, syn=.syn, 
                         bf = bf, w = w, g = g, inv = inv, INV = INV, ll.0 = ll.0, rate = rate, k = k, ab=log(ab), 
                         optK=FALSE, optW = TRUE),
                     codon3 = optimCodon(tree,data, Q=rep(1,1830), subs=.sub, syn=.syn, 
                         bf = bf, w = w, g = g, inv = inv, INV = INV, ll.0 = ll.0, rate = rate, k = k, ab=log(ab),
                         optK=TRUE, optW = FALSE))
                 tmp <- res[[5]]
                 m = length(tmp)
                 dnds = tmp[m]
                   
                 if(m>1) tstv <- tmp[1]
            }
            else
            res = optimQ(tree, data, Q = Q, subs = subs, bf = bf, w = w, g = g, inv = inv, INV = INV, 
                ll.0 = ll.0, rate = rate, k = k)
            Q = res[[1]]
            eig = edQt(Q = Q, bf = bf)
            if (trace > 0) 
                cat("optimize rate matrix: ", ll, "-->", res[[2]], 
                  "\n")
            ll = res[[2]]
        }
        if(optInv) {
            res = optimInv(tree, data, inv = inv, INV = INV, Q = Q, 
                bf = bf, eig = eig, k = k, shape = shape, rate = rate)
            inv = res[[1]]
            w = rep(1/k, k)
            g = discrete.gamma(shape, k)
            w = (1 - inv) * w
            if (wMix > 0) 
                w <- (1 - wMix) * w
            g = g/(1 - inv)
            g <- g * rate
            ll.0 = as.matrix(INV %*% (bf * inv))
            if (wMix > 0) 
                ll.0 <- ll.0 + llMix
            if (trace > 0) 
                cat("optimize invariant sites: ", ll, "-->", res[[2]], "\n")
            ll = res[[2]]
        }
        if(optGamma) {
            res = optimGamma(tree, data, shape = shape, k = k, 
                inv = inv, INV = INV, Q = Q, bf = bf, eig = eig, 
                ll.0 = ll.0, rate = rate)
            shape = res[[1]]
            w = rep(1/k, k)
            g = discrete.gamma(shape, k)
            if (inv > 0) {
                w = (1 - inv) * w
                g = g/(1 - inv)
            }
            if (wMix > 0) 
                w <- (1 - wMix) * w
            g <- g * rate
            if (trace > 0) 
                cat("optimize shape parameter: ", ll, "-->", 
                  res[[2]], "\n")
            ll = res[[2]]
        }
        if(optRate) {
            res = optimRate(tree, data, rate = rate, inv = inv, 
                INV = INV, Q = Q, bf = bf, eig = eig, k = k, 
                shape = shape, w = w, ll.0 = ll.0)
            if (res[[2]] > ll)rate = res[[1]]
            g = discrete.gamma(shape, k)
            w = rep(1/k, k)
            if (inv > 0) {
                w = (1 - inv) * w
                g = g/(1 - inv)
            }
            if (wMix > 0) 
                w <- (1 - wMix) * w
            g <- g * rate
            if (trace > 0) 
                cat("optimize rate: ", ll, "-->", res[[2]], "\n")
            ll = res[[2]]
        }
        if (optEdge) {  
           res <- optimEdge(tree, data, eig=eig, Q=Q, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0,
                 control = pml.control(epsilon = 1e-08, maxit = 5, trace=trace - 1)) 
           if (trace > 0) 
              cat("optimize edge weights: ", ll, "-->", res[[2]], "\n")
           if (res[[2]] > ll){  
              ll <- res[[2]]
              tree <- res[[1]]
           }
        }
        if(optRooted){
	        object = update(object, tree=tree, bf=bf, Q=Q, shape=shape, inv=inv) 
	        res <- optimRooted(object, control = pml.control(epsilon = 1e-07, maxit = 25, trace = trace-1))
	        if (res[[2]] > ll){  
                ll <- res[[2]]
                tree <- res[[1]]
            }     
	    }
        if(optNni) {
            swap = 0
            iter = 1
            while (iter < 4) {
                tmp <- pml.nni(tree, data, w, g, eig, bf, ll.0, ll, ...) 
                swap = swap + tmp$swap
                res <- optimEdge(tmp$tree, data, eig=eig, Q=Q, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0,
                   control = pml.control(epsilon = 1e-08, maxit = 3, trace=0))
                tree <- res[[1]]
                if (trace > 0) 
                  cat("optimize topology: ", ll, "-->", res[[2]], "\n")
                ll = res[[2]]
                iter = iter + 1
                if (tmp$swap == 0) {
                  iter = 4
                }
            }
            if (trace > 0) 
                cat(swap, "\n")
            if (swap > 0) 
                rounds = 1
            if (swap == 0) 
                optNni = FALSE
        }
        rounds = rounds + 1
        if(rounds > control$maxit) opti <- FALSE
        if (( ll1 - ll ) / ll  < control$eps) #abs(ll1 - ll)
            opti <- FALSE
        ll1 = ll
    }  
    if(type=="CODON"){
        object$dnds = dnds
        object$tstv = tstv
    }
    object <- update(object, tree = tree, data = data, bf = bf, Q = Q, inv = inv, shape = shape, k = k, 
        rate = rate)

# change call
    extras = pairlist(bf = bf, Q = Q, inv = inv, shape = shape, rate = rate)[c(optBf, optQ, optInv, optGamma, optRate)]
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    object$call = call   
    object
}


preptmp = function(child, parent, eig, bf){
    res = vector("list", length(child))
    tmp = (t(eig[[3]]) * bf)
    for(i in 1:length(child)){
        res[[i]] = (child[[i]] %*% eig[[2]]) *  ( parent[[i]] %*% tmp )
    }
    res 
}


fs <- function (old.el, eig, parent.dat, child.dat, weight, g=g, 
    w=w, bf=bf, ll.0=ll.0, evi, getA=TRUE, getB=TRUE) 
{
    if (old.el < 1e-8) old.el <- 1e-8
    lg = length(parent.dat)
    P <- getP(old.el, eig, g)
    nr = as.integer(length(weight))
    nc = as.integer(length(bf))
    eve = eig[[2]]
    dad <- .Call("getDAD", parent.dat, child.dat, P, nr, nc) 
    X <- .Call("getPrep", dad, child.dat, eig[[2]], evi, nr, nc) 
    .Call("FS4", eig, as.integer(length(bf)), as.double(old.el), 
            as.double(w), as.double(g), X, child.dat, dad, as.integer(length(w)), 
            as.integer(length(weight)), as.double(bf), as.double(weight), 
            as.double(ll.0), as.integer(getA), as.integer(getB))
}


# in optim.pml und optimEdge verwendbar (weniger memory) 
fs3 <- function (old.el, eig, parent.dat, child, weight, g=g, 
    w=w, bf=bf, ll.0=ll.0, contrast, contrast2, evi, ext=TRUE, getA=TRUE, getB=TRUE) # child.dat
{
    if (old.el < 1e-8) old.el <- 1e-8
    lg = length(parent.dat)
    P <- getP(old.el, eig, g)
    nr = as.integer(length(weight))
    nc = as.integer(length(bf))
    if(ext==FALSE){ 
       child.dat <- child
       eve = eig[[2]]
       dad <- .Call("getDAD", parent.dat, child.dat, P, nr, nc) 
       X <- .Call("getPrep", dad, child.dat, eig[[2]], evi, nr, nc) 
    }
    else {
        nco = as.integer(nrow(contrast))
        dad <- .Call("getDAD2", parent.dat, child, contrast, P, nr, nc, nco)
        child.dat <- vector("list", lg)
        for (i in 1:lg)child.dat[[i]] <- contrast[child, , drop=FALSE]
        X <- .Call("getPrep2", dad, child, contrast2, evi, nr, nc, nco)
    }
    .Call("FS4", eig, as.integer(length(bf)), as.double(old.el), 
            as.double(w), as.double(g), X, child.dat, dad, as.integer(length(w)), 
            as.integer(length(weight)), as.double(bf), as.double(weight), 
            as.double(ll.0), as.integer(getA), as.integer(getB))
}



prep <- function(child, parent, eig, bf, evi){
    res = vector("list", length(child))
    for(i in 1:length(parent)){
        res[[i]] = (child[[i]] %*% eig[[2]]) *  ( parent[[i]] %*% evi ) 
    }
    res 
}

prep1 <- function(child, parent, eve, bf, evi){
    res = vector("list", length(child))
    for(i in 1:length(parent)){
        res[[i]] = (child[[i]] %*% eve) *  ( parent[[i]] %*% evi ) 
    }
    res 
}

prep2 <- function(child, parent, eig, bf, contrast2, evi){
    res = vector("list", length(child))
    for(i in 1:length(parent)){
        res[[i]] = (contrast2[child, , drop=FALSE]) *  ( parent[[i]] %*% evi ) 
    }
    res 
}


fs2 <- function (old.el, eig, parent.dat, child, weight, g=g, 
    w=w, bf=bf, ll.0=ll.0, contrast, contrast2, evi, getA=TRUE, getB=TRUE) 
{
    if (old.el < 1e-8) old.el <- 1e-8
    lg = length(parent.dat)
    P <- getP(old.el, eig, g)
    nr = as.integer(length(weight))
    nc = as.integer(length(bf))
    nco = as.integer(nrow(contrast))
    dad <- .Call("getDAD2", parent.dat, child, contrast, P, nr, nc, nco)
    child.dat <- vector("list", lg)
    for (i in 1:lg)child.dat[[i]] <- contrast[child, , drop=FALSE]
    X <- .Call("getPrep2", dad, child, contrast2, evi, nr, nc, nco)

    .Call("FS4", eig, as.integer(length(bf)), as.double(old.el), 
            as.double(w), as.double(g), X, child.dat, dad, as.integer(length(w)), 
            as.integer(length(weight)), as.double(bf), as.double(weight), 
            as.double(ll.0), as.integer(getA), as.integer(getB))
}


optimEdge <- function (tree, data, eig=eig, Q=Q, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0,
              control = pml.control(epsilon = 1e-08, maxit = 10, trace=0), ...) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") 
        tree <- reorderPruning(tree)
    el <- tree$edge.length
    tree$edge.length[el < 0] <- 1e-08
    oldtree = tree
    dat <- NULL    
    old.ll <- pml6(tree, data, bf, Q, eig, ll.0, w, g)
    start.ll = old.ll
    contrast <- attr(data, "contrast")
    contrast2 <- contrast %*% eig[[2]] 
    evi = (t(eig[[3]]) * bf)
    weight <- attr(data, "weight")
    eps = 1
    iter = 0
    child = tree$edge[, 2]
    parent = tree$edge[, 1]
    nTips = min(parent) - 1
    n = length(tree$edge.length)
    k = length(w)
    data = subset(data, tree$tip)  
    lt = length(tree$tip)
    for(i in 1:k)dat[i, 1:lt]=data    
    child.dat = vector("list", k)
    while (eps > control$eps && iter < control$maxit) {
        for (j in n:1) {
            parent.dat = dat[, parent[j]]
            old.el = tree$edge.length[j] 
            if (child[j] > nTips){ 
               newEL <- fs(old.el, eig, parent.dat, dat[, child[j]], weight, 
                    g = g, w = w, bf = bf, ll.0 = ll.0, evi, getA=TRUE, getB=TRUE)
            }
            else{
                newEL <- fs2(old.el, eig, parent.dat, data[[child[j] ]], weight, 
                    g = g, w = w, bf = bf, ll.0 = ll.0, contrast, contrast2, evi, getA=TRUE, getB=FALSE)
                }          
            el[j] = newEL[[1]]
            dat[, parent[j]] = newEL[[2]]
            if (child[j] > nTips) {
                dat[, child[j]] = newEL[[3]]
            }
        }
        tree$edge.length = el
        iter = iter + 1
        dat <- NULL
        newll <- pml6(tree, data, bf, Q, eig, ll.0, w, g)
        eps = ( old.ll - newll ) / newll
        if(eps<0) return(list(oldtree, old.ll))
        oldtree = tree
        old.ll = newll
    }
    if(control$trace>0) cat(start.ll, " -> ", newll, "\n")
    list(tree, newll)
}


optimRooted <- function(fit, control = pml.control(epsilon = 1e-08, maxit = 25, trace = 1)){
    tree = fit$tree
    g = fit$g
    w = fit$w
    eig = fit$eig
    bf = fit$bf
    ll.0 = fit$ ll.0
   
    if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") 
        tree <- ape:::reorder.phylo(tree, "pruningwise") # reorderPruning???
    if(!is.rooted(tree))stop("tree must be rooted")
  
    weight = attr(fit$data , "weight")
    optRoot <- function(t, dat, el1, el2, g, w, eig, bf, ll.0, weight){
        P = getP(c(el1+t,el2-t), eig, g)
        res=vector("list",length(g))
        for(i in 1:length(g)){
             res[[i]]=(dat[[i,1]] %*% P[[i,1]])
             for(j in 2:(dim(dat)[2])){
                 res[[i]] = res[[i]] * (dat[[i,j]] %*% P[[i,j]])
                 }
             }
        result <- ll.0
        for(i in 1:length(g))result <- result +  res[[i]]  %*% (w * bf)
        loglik <- sum(weight * log(result))
        loglik
    }
    optRoot1 <- function(t, dat, el, g, w, eig, bf, ll.0, weight, logLik=TRUE){
        P = getP(c(el+t), eig, g)
        res=vector("list",length(g))
        for(i in 1:length(g)){
             res[[i]]=(dat[[i,1]] %*% P[[i,1]])
             for(j in 2:(dim(dat)[2])){
                 res[[i]] = res[[i]] * (dat[[i,j]] %*% P[[i,j]])
                 }
             }
        if(!logLik)return(res)
        result <- ll.0  
        for(i in 1:length(g))result <- result +  res[[i]]  %*% (w * bf)
        loglik <- sum(weight * log(result))
        loglik
    }
# use pml2 instead of pml5
    scaleEdges <- function(t,fit){
        fit$tree$edge.length = fit$tree$edge.length*t
        pml5(fit)
    }

    child = tree$edge[, 2]
    parent = tree$edge[, 1]
    nTips = min(parent) - 1
 
    ll= fit$logLik
    eps=10
    iter = 1
    while(eps>control$eps && iter < control$maxit){
        t <- optimize(f=scaleEdges, interval=c(0.25,4),fit, maximum=TRUE)
        fit$tree$edge.length = fit$tree$edge.length*t[[1]]
        tree = fit$tree
        el = tree$edge.length
        .dat=NULL 
        dat = rnodes(fit) # replace
        dat2 = .dat
        for(i in 1:length(parent)){
            if(child[i]>nTips){
                dad = child[i]
                kids = which(parent==dad)
                children = child[kids]  
                kidsEl = el[kids] 
                minEl = min(kidsEl) 
                kidsEl = kidsEl - minEl
                maxEl = minEl + el[i] # el[dad]
                t <- optimize(f=optRoot,interval=c(0,maxEl),dat=cbind(dat2[,children, drop=FALSE],dat[,dad, drop=FALSE]),el1=kidsEl,el2=maxEl, g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, weight=weight, maximum=TRUE)
                el[kids] = kidsEl+t[[1]]
                dat2[,dad] = optRoot1(0,dat=dat2[,children, drop=FALSE],el=el[kids], g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, weight=weight, logLik=FALSE)
                el[i] = maxEl-t[[1]]
                tree$edge.length = el
                fit$tree=tree
                }
            }
        kids = which(parent==(nTips+1))
        children = child[kids]  
        kidsEl = el[kids] 
        minEl = min(kidsEl) 
        kidsEl = kidsEl - minEl
        t <- optimize(f=optRoot1,interval=c(0,3),dat=dat2[,children, drop=FALSE],el=kidsEl, g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, weight=weight, maximum=TRUE)
        el[kids] = kidsEl+t[[1]]
        tree$edge.length = el
        fit$tree=tree
        eps = (ll - t[[2]]) / t[[2]]
        ll=t[[2]]
        iter = iter+1
        }
    list(tree=tree, logLik=ll, c(eps=eps, iter=iter))
}


pml5 <- function (object, ...) 
{
    tree = object$tree
    Q = object$Q
    bf = object$bf
    eig = object$eig
    w = object$w
    g = object$g
    data = object$data
    lll <- object$ll.0
    weight = attr(data, "weight")
    m = 1
    p = length(g)
    q = length(tree$edge[, 1]) + 1
    resll = vector("list", p)
#    resll <- matrix(0, nr, p)
    dat = vector("list", q * p)
    attr(dat, "dim") = c(p, q)
    asdf <- NULL
    while (m <= p) {
        resll[[m]] = ll(data, tree, bf = bf, g = g[m], Q = Q, 
            eig = eig, assign.dat = TRUE, ...)
        dat[m, ] <- asdf
        m = m + 1
    }
    attr(dat, "dimnames") = list(NULL, attr(asdf, "names"))
    for (i in 1:p) lll = lll + resll[[i]] * w[i]
    siteLik <- lll
    siteLik <- log(siteLik)
    ll0 = sum(weight * siteLik)
    assign("dat", dat, envir = parent.frame(n = 1))
    ll0
}


pml6 <- function (tree, data, bf, Q, eig, ll.0, w, g, ...) 
{
    lll=ll.0
    weight = attr(data, "weight")
    m = 1
    p = length(g)
    q = length(tree$edge[, 1]) + 1
    resll = vector("list", p)
#    resll <- matrix(0, nr, p)
    dat = vector("list", q * p)
    attr(dat, "dim") = c(p, q)
    asdf <- NULL
    while (m <= p) {
        resll[[m]] = ll(data, tree, bf = bf, g = g[m], Q = Q, 
            eig = eig, assign.dat = TRUE, ...)
        dat[m, ] <- asdf
        m = m + 1
    }
    attr(dat, "dimnames") = list(NULL, attr(asdf, "names"))
    for (i in 1:p) lll = lll + resll[[i]] * w[i]
    siteLik <- lll 
    siteLik <- log(siteLik)
    ll0 = sum(weight * siteLik)
    assign("dat", dat, envir = parent.frame(n = 1))
    ll0
}


#
# pmlPart + pmlCluster
#
optimPartQ <- function (object, Q = c(1, 1, 1, 1, 1, 1), ...) 
{
    l = length(Q)
    Q = Q[-l]
    Q = sqrt(Q)
    fn = function(Q, object, ...) {
        result <- 0
        Q = c(Q^2, 1)
        n <- length(object)
        for (i in 1:n) result <- result + update(object[[i]], Q = Q, ...)$logLik
        result
    }
    res = optim(par = Q, fn = fn, gr = NULL, method = "L-BFGS-B", 
        lower = 0, upper = Inf, control = list(fnscale = -1, 
            maxit = 25), object = object, ...)
    res[[1]] = c(res[[1]]^2, 1)
    res
}


optimPartQGeneral <- function (object, Q = c(1, 1, 1, 1, 1, 1), subs=rep(1,length(Q)), ...) 
{
    m = length(Q)
    n = max(subs)
    ab = numeric(n)
    for(i in 1:n) ab[i]=log(Q[which(subs==i)[1]])
    fn = function(ab, object, m, n, subs, ...) {
        Q = numeric(m)
        for(i in 1:n)Q[subs==i] = ab[i]
        Q = exp(Q)
        result = 0
        for (i in 1:length(object)) result <- result + update(object[[i]], Q = Q, ...)$logLik
        result
    }
    res = optim(par = ab, fn = fn, gr = NULL, method = "L-BFGS-B", 
        lower = -Inf, upper = Inf, control = list(fnscale = -1, 
            maxit = 25), object = object, m=m, n=n, subs=subs, ...)
    Q = rep(1, m)
    for(i in 1:n) Q[subs==i] = exp(res[[1]][i])
    res[[1]] = Q
    res
}


optimPartBf <- function (object, bf = c(0.25, 0.25, 0.25, 0.25), ...) 
{
    l = length(bf)
    nenner = 1/bf[l]
    lbf = log(bf * nenner)
    lbf = lbf[-l]
    fn = function(lbf, object, ...) {
        result <- 0
        bf = exp(c(lbf, 0))
        bf = bf/sum(bf)
        n <- length(object)
        for (i in 1:n) result <- result + update(object[[i]], 
            bf = bf, ...)$logLik
        result
    }
    res = optim(par = lbf, fn = fn, gr = NULL, method = "Nelder-Mead", 
        control = list(fnscale = -1, maxit = 500), object, ...)
    print(res[[2]])
    bf = exp(c(res[[1]], 0))
    bf = bf/sum(bf)
}


optimPartInv <- function (object, inv = 0.01, ...) 
{
    fn = function(inv, object, ...) {
        result <- 0
        n <- length(object)
        for (i in 1:n) result <- result + update(object[[i]], inv = inv, 
            ...)$logLik
        result
    }
    res = optimize(f = fn, interval = c(0, 1), lower = 0, upper = 1, 
        maximum = TRUE, tol = 1e-04, object, ...)
    print(res[[2]])
    res[[1]]
}


optimPartGamma <- function (object, shape = 1, ...) 
{
    fn = function(shape, object, ...) {
        result <- 0
        n <- length(object)
        for (i in 1:n) result <- result + update(object[[i]], shape = shape, 
            ...)$logLik
        result
    }    
    res = optimize(f = fn, interval = c(0, 100), lower = 0, upper = 100, 
        maximum = TRUE, tol = 0.01, object, ...)
    res
}


dltmp <- function (fit, i=1, transform=transform) # i = weights
{
    tree = fit$tree 
    data = getCols(fit$data, tree$tip) 
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorderPruning(tree)
    q = length(tree$tip.label)
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    m = length(edge) + 1  # max(edge)
    dat = vector(mode = "list", length = m)
    eig = fit$eig
    w = fit$w[i]
    g = fit$g[i]
    bf = fit$bf
    el <- tree$edge.length
    P <- getP(el, eig, g)
    nr <- as.integer(attr(data, "nr"))
    nc <- as.integer(attr(data, "nc"))
    node = as.integer(node - min(node))
    edge = as.integer(edge - 1)
    nTips = as.integer(length(tree$tip))
    mNodes = as.integer(max(node) + 1)
    contrast = attr(data, "contrast")
    nco = as.integer(dim(contrast)[1])
    dat[(q + 1):m] <- .Call("LogLik2", data, P, nr, nc, node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")
    result = dat[[q+1]] %*% (bf * w)

    parent <- tree$edge[, 1]
    child <- tree$edge[, 2]
    nTips = min(parent) - 1
    datp = vector("list", m)
    el = tree$edge.length 
    if (transform) dP = getdP(tree$edge.length, eig, g)
    else dP = getdP2(tree$edge.length, eig, g)
   
    datp[(nTips + 1)] = dat[(nTips + 1)]
    l = length(child)
    dl = matrix(0, nr, l)
    for (j in (m - 1):1) {
        # tips have factor format, internal edges are matrices
        if (child[j] > nTips){
             tmp2 = (datp[[parent[j]]]/(dat[[child[j]]] %*% P[[j]]))
             dl[, j] = (tmp2 * (dat[[child[j]]] %*% dP[[j]])) %*% (w * bf)
             datp[[child[j]]] = (tmp2 %*% P[[j]]) * dat[[child[j]]]  
             }
        else{
             tmp2 = (datp[[parent[j]]]/((contrast %*% P[[j]])[data[[child[j]]],] ))
             dl[, j] = (tmp2 * ((contrast %*% dP[[j]])[data[[child[j]]],]) ) %*% (w * bf)    
             }
    }
    dl
}


dl <- function(x, transform = TRUE){
  w = x$w 
  l=length(x$w)
  dl = dltmp(x, 1, transform)
  i=2
  while(i < (l+1)){
    dl = dl + dltmp(x, i, transform)
    i = i + 1
  } 
  dl
}


# add control and change edge
optimPartEdge <- function (object, ...) 
{
    tree <- object[[1]]$tree
    theta <- object[[1]]$tree$edge.length
    n <- length(object)
    l <- length(theta)
    nrv <- numeric(n)
    for (i in 1:n) nrv[i] = attr(object[[i]]$data, "nr")
    cnr <- cumsum(c(0, nrv))
    weight = numeric(sum(nrv))
    dl <- matrix(NA, sum(nrv), l)
    for (i in 1:n) weight[(cnr[i] + 1):cnr[i + 1]] = attr(object[[i]]$data, 
        "weight")
    ll0 = 0
    for (i in 1:n) ll0 = ll0 + object[[i]]$logLik
    eps = 1
    scalep =1
    k = 1
    while (eps > 0.001 & k<50) {
        if(scalep==1){
            for (i in 1:n) {
                lv = drop(exp(object[[i]]$site))
                dl[(cnr[i] + 1):cnr[i + 1], ] = dl(object[[i]], TRUE)/lv
            }
            sc = colSums(weight * dl)
            F = crossprod(dl * weight, dl) #+ diag(l)*1e-10
        }
        thetaNew = log(theta) + scalep * solve(F, sc)
        tree$edge.length = as.numeric(exp(thetaNew))
        for (i in 1:n) object[[i]] <- update(object[[i]], tree = tree)
        ll1 = 0
        for (i in 1:n) ll1 = ll1 + object[[i]]$logLik
        eps <- ll1 - ll0
        if (eps < 0) {
            scalep = scalep/2
            eps = 1
            thetaNew = log(theta)
            ll1 = ll0
        }
        else scalep = 1
        theta = exp(thetaNew)
        ll0 <- ll1
        k=k+1
    }
    object
}


makePart <- function(fit, weight=~index+genes){
    dat <- fit$data 
    if(class(weight)[1]=="formula")     
        weight <- xtabs(weight, data=attr(dat, "index"))
    fits <- NULL 
    for(i in 1:dim(weight)[2]){ 
       ind <- which(weight[,i] > 0)
       dat2 <- getRows(dat, ind)
       attr(dat2, "weight") <- weight[ind,i]
       fits[[i]] <- update(fit, data = dat2)
    }
    names(fits) = colnames(fits)
    fits
}


pmlPart <- function (formula, object, control=pml.control(epsilon=1e-8, maxit=10, trace=1), ...) 
{
    call <- match.call()
    form <- phangornParseFormula(formula)
    opt <- c("nni", "bf", "Q", "inv", "shape", "edge", "rate")
    optAll <- match(opt, form$left)
    optPart <- match(opt, form$right)
    AllNNI <- !is.na(optAll[1])
    AllBf <- !is.na(optAll[2])
    AllQ <- !is.na(optAll[3])
    AllInv <- !is.na(optAll[4])
    AllGamma <- !is.na(optAll[5])
    AllEdge <- !is.na(optAll[6])
    PartNni <- !is.na(optPart[1])
    PartBf <- !is.na(optPart[2])
    PartQ <- !is.na(optPart[3])
    PartInv <- !is.na(optPart[4])
    PartGamma <- !is.na(optPart[5])
    PartEdge <- !is.na(optPart[6])
    PartRate <- !is.na(optPart[7])
 
    if(class(object)=="pml") fits <- makePart(object, ...)   
    if(class(object)=="pmlPart") fits <- object$fits
    if(class(object)=="list") fits <- object

    trace = control$trace
    epsilon = control$epsilon
    maxit = control$maxit

    p <- length(fits)
    m = 1
    logLik = 0
    for (i in 1:p) logLik = logLik + fits[[i]]$log
    eps = 10
    while (eps > epsilon & m < maxit) {
        loli = 0
        if(any(c(PartNni, PartBf, PartInv, PartQ, PartGamma, PartEdge, PartRate))){
            for (i in 1:p) {
                fits[[i]] = optim.pml(fits[[i]], PartNni, PartBf, 
                    PartQ, PartInv, PartGamma, PartEdge, PartRate, 
                    control = pml.control(maxit = 3, epsilon = 1e-8, trace-1))
            }
        } 
        if (AllQ) {
            Q = fits[[1]]$Q
            subs = c(1:(length(Q)-1), 0)
            newQ <- optimPartQGeneral(fits, Q=Q, subs=subs)
            for (i in 1:p) fits[[i]] <- update(fits[[i]], Q = newQ[[1]])
        }
        if (AllBf) {
             bf = fits[[1]]$bf
            newBf <- optimPartBf(fits, bf=bf)
            for (i in 1:p) fits[[i]] <- update(fits[[i]], bf = newBf)
        }
        if (AllInv) {
            inv = fits[[1]]$inv
            newInv <- optimPartInv(fits, inv=inv)
            for (i in 1:p) fits[[i]] <- update(fits[[i]], inv = newInv)
        }
        if (AllGamma) {
            shape = fits[[1]]$shape
            newGamma <- optimPartGamma(fits, shape=shape)[[1]]
            for (i in 1:p) fits[[i]] <- update(fits[[i]], shape = newGamma)
        }
        if (AllNNI){
            fits <- optimPartNNI(fits,AllEdge)
            if(trace>0) cat(attr(fits,"swap"), " NNI operations performed")
        }
        if (AllEdge) 
            fits <- optimPartEdge(fits)
        if (PartRate){
            tree = fits[[1]]$tree
            rate=numeric(p)
            wp =numeric(p) 
            for(i in 1:p){
                wp[i]=sum(fits[[i]]$weight)
                rate[i] <- fits[[i]]$rate
                }          
            ratemult = sum(wp) / sum(wp*rate)
            tree$edge.length = tree$edge.length/ratemult  
            for(i in 1:p)fits[[i]] = update(fits[[i]], tree=tree, rate=rate[i]*ratemult)   
        }
        loli <- 0
        for (i in 1:p) loli <- loli + fits[[i]]$log
        eps = (logLik - loli)/loli
        if(trace>0) cat("loglik:", logLik, "-->", loli, "\n")
        logLik <- loli
        m = m + 1
    }
    
    df <- matrix(1, 6 ,2)
    colnames(df) <- c("#df", "group")
    rownames(df) <- c("Edge", "Shape", "Inv", "Bf", "Q", "Rate")
    df[1,1] <- length(fits[[1]]$tree$edge.length)
    df[2,1] <- fits[[1]]$k > 1
    df[3,1] <- fits[[1]]$inv > 0
    df[4,1] <- length(unique(fits[[1]]$bf)) - 1
    df[5,1] <- length(unique(fits[[1]]$Q)) - 1
    df[6,1] <- 0 # rates 
    if(PartEdge) df[1,2] = p
    if(PartGamma) df[2,2] = p
    if(PartInv) df[3,2] = p
    if(PartBf) df[4,2] = p
    if(PartQ) df[5,2] = p
    if(PartRate) df[6,1] = p-1     
    attr(logLik, "df") = sum(df[,1]*df[,2])
    object <- list(logLik = logLik, fits = fits, call = call, df=df)
    class(object) <- "pmlPart" 
    object
}


#
# Distance Matrix methods
#
bip <- function (x) 
{
    if (is.null(attr(x, "order")) || attr(x, "order") == "cladewise") 
        x = reorderPruning(x)
    nTips = length(x$tip)
    res = vector("list", max(x$edge))
    res[1:nTips]=1:nTips
    tmp = bipart(x)
    res[attr(tmp, "nodes")] = tmp
    res
}

bip2 <- function (x) 
{
    if(is.null(attr(x,"order")) || attr(x, "order")=="cladewise") x = reorderPruning(x)
    nNode = x$Nnode
    nTips = length(x$tip)
    parent <- as.integer(x$edge[, 1])
    child <- as.integer(x$edge[, 2])
    res = vector("list", max(x$edge))
    p = parent[1]
    tmp = NULL
    for(i in 1:nTips) res[[i]]=i
    for (i in 1:length(parent)) {
        pi = parent[i]
        ci = child[i]
        if (pi == p) {
            if (ci < (nTips + 1)) 
                tmp = cisort(tmp, ci)
            else tmp = cisort(tmp, res[[ci]])
        }
        else {
            res[[p]] = (tmp)
            if (ci < (nTips + 1)) 
                tmp = ci
            else tmp = res[[ci]]
            p = pi
        }
    }
    res[[p]] = (tmp)
    res
}

# as.Matrix, sparse = TRUE, 
designTree <- function(tree, method="unrooted", sparse=FALSE, ...){
    if (!is.na(pmatch(method, "all"))) 
        method <- "unrooted"
    METHOD <- c("unrooted", "rooted")
    method <- pmatch(method, METHOD)
    if (is.na(method)) stop("invalid method")
    if (method == -1) stop("ambiguous method")
    if(!is.rooted(tree) & method==2) stop("tree has to be rooted")  
    if(method==1){ X <- designUnrooted(tree,...)
        if(sparse) X = Matrix(X)  
    }
    if(method==2) X <- designUltra(tree, sparse=sparse,...)
    X
}


designUnrooted = function(tree,order=NULL){
    if(is.rooted(tree))tree = unroot(tree)
    p=bipartition(tree)
    if(!is.null(order)) p=p[,order]
    n = dim(p)[1]
    m = dim(p)[2]
    res = matrix(0,(m-1)*m/2, n)
    k = 1
    for(i in 1:(m-1)){
        for(j in (i+1):m){
            res[k,] = p[,i]!=p[,j]
            k=k+1
        }
    }
    colnames(res) = paste(tree$edge[,1],tree$edge[,2],sep="<->")
    res
    }

    
designUltra <- function (tree, sparse=FALSE) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") 
        tree = reorderPruning(tree)
    leri = allChildren(tree)
    bp = bip(tree)
    n = length(tree$tip)
    l = tree$Nnode   
    nodes = integer(l)
    k = 1L
    u=numeric( n * (n - 1)/2)
    v=numeric( n * (n - 1)/2)
    m = 1L
    for (i in 1:length(leri)) {
        if (!is.null(leri[[i]])) {
            ind =  getIndex(bp[[leri[[i]][1] ]], bp[[leri[[i]][2] ]], n) 
            li = length(ind)
            v[m: (m+li-1)]=k
            u[m: (m+li-1)]=ind   
            nodes[k]=i
            m = m+li
            k = k + 1L
        }
    }
    if(sparse) X = sparseMatrix(i=u,j=v, x=2L)
    else{
        X = matrix(0L, n * (n - 1)/2, l)              
        X[cbind(u,v)]=2L
    }
    colnames(X) = nodes
    attr(X, "nodes") = nodes
    X
}
nnls.tree <- function(dm, tree, rooted=FALSE, trace=1){
    if(is.rooted(tree) & rooted==FALSE){
        tree = unroot(tree)
        warning("tree was rooted, I unrooted the tree!")
    }
    tree = reorderPruning(tree)
    dm = as.matrix(dm)
    k = dim(dm)[1]
    labels = tree$tip
    dm = dm[labels,labels]
    y = dm[lower.tri(dm)]
#computing the design matrix from the tree   
    if(rooted) X = designUltra(tree) 
    else X = designUnrooted(tree)
# na.action
    if(any(is.na(y))){
        ind = which(is.na(y))
        X = X[-ind,,drop=FALSE]
        y= y[-ind]
    }
# LS solution 
    fit = lm.fit(X,y)
    betahat = fit$coefficients
    if(rooted){
        bhat = numeric(max(tree$edge))
        bhat[as.integer(colnames(X))] = betahat
        betahat = bhat[tree$edge[,1]] - bhat[tree$edge[,2]]
    }
    if(!any(betahat<0)){
         RSS = sum(fit$residuals^2)
         if(trace)print(paste("RSS:", RSS))
         attr(tree, "RSS") = RSS
         tree$edge.length[] = betahat
         return(tree)
    }
# non-negative LS
    n = dim(X)[2]
    Dmat <- crossprod(X) # cross-product computations
    dvec <- crossprod(X, y)
    if(rooted){
        l = nrow(tree$edge)
        Amat = matrix(0, n, l)
        ind = match(tree$edge[,1], colnames(X))
        Amat[cbind(ind, 1:l)] = 1
        ind = match(tree$edge[,2], colnames(X))
        Amat[cbind(ind, 1:l)] = -1  
    }
    else Amat <- diag(n)
    betahat <- quadprog::solve.QP(Dmat,dvec,Amat)$sol # quadratic programing solving
    RSS = sum((y-(X%*%betahat))^2) 

    if(rooted){
        bhat = numeric(max(tree$edge))
        bhat[as.integer(colnames(X))] = betahat
        betahat = bhat[tree$edge[,1]] - bhat[tree$edge[,2]]
    }
    tree$edge.length[] = betahat
    if(trace)print(paste("RSS:", RSS))
    attr(tree, "RSS") = RSS
    tree
}


designSplits <- function (x, splits = "all", ...) 
{
    if (!is.na(pmatch(splits, "all"))) 
        splits <- "all"
    SPLITS <- c("all", "star") #,"caterpillar")
    splits <- pmatch(splits, SPLITS)
    if (is.na(splits)) stop("invalid splits method")
    if (splits == -1) stop("ambiguous splits method")  
    if(splits==1) X <-  designAll(x)
    if(splits==2) X <-  designStar(x)
    return(X)
}


designAll <- function(n){
    Y = matrix(0L, n*(n-1)/2, n)
    k = 1
    for(i in 1:(n-1)){
    for(j in (i+1):n){
          Y[k,c(i,j)]=1L
          k=k+1L
        }
    }
    m <- n-1L
    X <- matrix(0L, m+1, 2^m)
    for(i in 1:m)
    X[i, ] <- rep(rep(c(0L,1L), each=2^(i-1)),2^(m-i))
    X <- X[,-1]
    list(X=(Y%*%X)%%2,Splits=t(X))
}


designStar = function(n){
    res=NULL
    for(i in 1:(n-1)) res = rbind(res,cbind(matrix(0,(n-i),i-1),1,diag(n-i)))
    res
}


cisort <- function(x,y){
    k = length(x)
    l=length(y)
    .C("cisort",as.integer(x),as.integer(y),k,l, integer(k+l), DUP=FALSE)[[5]]
}


bipart = function(obj){
    if (is.null(attr(obj, "order")) || attr(obj, "order") == "cladewise") 
        obj <- reorderPruning(obj)
    maxP  = max(obj$edge)
    nTips = length(obj$tip)
    res <- .Call("bipart", as.integer(obj$edge[,1]) , as.integer(obj$edge[,2]), as.integer(nTips), as.integer(maxP), as.integer(obj$Nnode))
    attr(res, "nodes") = unique(obj$edge[,1])
    res    
}


bipartition <- function (tree) 
{
    if(is.rooted(tree))tree <- unroot(tree)
    if(is.null(attr(tree,"order")) || attr(tree, "order")=="cladewise") tree <- reorderPruning(tree)
    bp <- bipart(tree)
    nTips = length(tree$tip)
    l = length(bp)
    m = length(bp[[l]])
    k = length(tree$edge[, 1])
    result = matrix(0L, l, m)
    res = matrix(0L, k, m)
    for (i in 1:l) result[i, bp[[i]]] = 1L
    result = result[-l, ,drop=FALSE]
    for (i in 1:nTips) res[(tree$edge[, 2] == i), i] = 1L     
#    res[tree$edge[, 2] > nTips, ] = result
    res[ match(unique(tree$edge[,1]),tree$edge[,2])[-l], ] = result
    colnames(res) = tree$tip.label
    rownames(res) = tree$edge[,2]
    res[res[, 1] == 1, ] = 1L - res[res[, 1] == 1, ]
    res
}


#
# some generic tree functions
#
Ancestors <- function (x, node, type = c("all", "parent")) 
{
    parents <- x$edge[, 1]
    child <- x$edge[, 2]
    pvector <- numeric(max(x$edge)) # parents
    pvector[child] <- parents    
    type <- match.arg(type)
    if (type == "parent") 
        return(pvector[node])
    res <- numeric(0)
    repeat {
        anc <- pvector[node]
        if (anc == 0) break
        res <- c(res, anc)
        node <- anc
    }
    res
}


allChildren <- function(x){
   parent = x$edge[,1]
   children = x$edge[,2]
   res = vector("list", max(x$edge))
   for(i in 1:length(parent)) res[[parent[i]]] = c(res[[parent[i]]], children[i])
   res
}


allChildrenNew <- function(x){
   if (is.null(attr(x, "order")) || attr(x, "order") == "cladewise") 
        x <- reorderPruning(x)
   parent = unique(x$edge[,1])
   tab = tabulate(x$edge[,1])[parent]
   children = x$edge[,2]
   .Call("allChildren", as.integer(children), as.integer(parent), as.integer(tab), as.integer(max(x$edge)))
}


Children <- function(x, node){
   if(length(node)==1)return(x$edge[x$edge[,1]==node,2])
   allChildren(x)[node]
}


Descendants = function(x, node, type=c("tips","children","all")){
    type <- match.arg(type)
    if(type=="children") return(Children(x, node))
    if(type=="tips") return(bip(x)[node])
    ch = allChildren(x) # out of the loop
    desc = function(x, node, type){
        isInternal = logical(max(x$edge))
        isInternal[ unique(x$edge[,1]) ] =TRUE       
        if(!isInternal[node])return(node)   
        res = NULL
        while(length(node)>0){
            tmp = unlist(ch[node])
            res = c(res, tmp)
            node = tmp[isInternal[tmp]]
        }
    if(type=="tips") return(res[!isInternal[res]])
    res
    }
    if(length(node)>1) return(lapply(node, desc, x=x, type=type))
    desc(x,node, type)
}


Siblings = function (x, node, include.self = FALSE) 
{
    l = length(node)
    if(l==1){
        v <- Children(x, Ancestors(x, node, "parent"))
        if (!include.self) 
            v <- v[v != node]
        return(v)
    }
    else{    
        parents <- x$edge[, 1]
        child <- x$edge[, 2]
        pvector <- numeric(max(x$edge)) # parents
        pvector[child] <- parents
        root <- as.integer(parents[!match(parents, child, 0)][1])
        res = vector("list", l)
        ch = allChildren(x)
        k = 1
        for(i in node){
            if(i != root){
                tmp <- ch[[ pvector[i] ]]
                res[[k]] = tmp[tmp != i]
            } 
            k=k+1    
        }     
    }
    res
}


getIndex = function(left, right, n){
    if(n<max(left) | n<max(right)) stop("Error")  
    left = as.integer(left)
    right = as.integer(right)
    ll = length(left)
    lr = length(right)
    .C("giveIndex", left, right, ll, lr, as.integer(n), integer(ll*lr))[[6]]+1
}



pmlCluster.fit <- function (formula, fit, weight, p = 4, part = NULL, control=pml.control(epsilon=1e-8, maxit=10, trace=1), ...) 
{
    call <- match.call()
    form <- phangornParseFormula(formula)
    opt <- c("nni", "bf", "Q", "inv", "shape", "edge", "rate")
    optAll <- match(opt, form$left)
    optPart <- match(opt, form$right)
    AllNNI <- !is.na(optAll[1])
    AllBf <- !is.na(optAll[2])
    AllQ <- !is.na(optAll[3])
    AllInv <- !is.na(optAll[4])
    AllGamma <- !is.na(optAll[5])
    AllEdge <- !is.na(optAll[6])
    PartNni <- !is.na(optPart[1])
    PartBf <- !is.na(optPart[2])
    PartQ <- !is.na(optPart[3])
    PartInv <- !is.na(optPart[4])
    PartGamma <- !is.na(optPart[5])
    PartEdge <- !is.na(optPart[6])
    PartRate <- !is.na(optPart[7])
    nrw <- dim(weight)[1]
    ncw <- dim(weight)[2]
    if (is.null(part)){ 
        part = rep(1:p, length=ncw)
        part = sample(part)
        }
    Part = part
    Gtrees = vector("list", p)
    dat <- fit$data
    attr(fit$orig.data, "index") <- attr(dat, "index") <- NULL
    for (i in 1:p) Gtrees[[i]] = fit$tree
    fits = vector("list", p)
    for (i in 1:p) fits[[i]] = fit
    trace = control$trace
    eps = 0
    m = 1
    logLik = fit$log
    trees = list()
    weights = matrix(0, nrw, p)
    lls = matrix(0, nrw, p)
    loli = fit$log
    oldpart = part
    eps2 = 1
    iter = 0
    swap = 1
    while (eps < ncw || abs(eps2) > control$eps) {
        df2 = 0
        
        if(any(c(PartNni, PartBf, PartInv, PartQ, PartGamma, PartEdge, PartRate))){
            for (i in 1:p) {
                weights[, i] = rowSums(weight[, which(part == i), 
                    drop = FALSE])
                ind <- which(weights[, i] > 0)
                dat2 <- getRows(dat, ind)
                attr(dat2, "weight") <- weights[ind, i]
                fits[[i]] <- update(fits[[i]], data = dat2)
                fits[[i]] = optim.pml(fits[[i]], PartNni, PartBf, 
                    PartQ, PartInv, PartGamma, PartEdge, PartRate, 
                    control = pml.control(epsilon = 1e-8, maxit = 3, trace-1))
                lls[, i] = update(fits[[i]], data = dat)$site
                Gtrees[[i]] = fits[[i]]$tree
            }
        }
        if (AllQ) {
            Q = fits[[1]]$Q
            subs = c(1:(length(Q)-1), 0)
            newQ <- optimPartQGeneral(fits, Q=Q, subs=subs)[[1]]
            for (i in 1:p) fits[[i]] <- update(fits[[i]], Q = newQ)
            df2 = df2 + length(unique(newQ)) - 1
        }
        if (AllBf) {
	        bf = fits[[1]]$bf
            newBf <- optimPartBf(fits, bf=bf)
            for (i in 1:p) fits[[i]] <- update(fits[[i]], bf = newBf)
            df2 = df2 + length(unique(newBf)) - 1
        }
        if (AllInv) {
            inv = fits[[1]]$inv
            newInv <- optimPartInv(fits, inv=inv)
            for (i in 1:p) fits[[i]] <- update(fits[[i]], inv = newInv) #there was an Error
            df2 = df2 + 1
        }
        if (AllGamma) {
            shape = fits[[1]]$shape
            newGamma <- optimPartGamma(fits, shape=shape)[[1]]        
            for (i in 1:p) fits[[i]] <- update(fits[[i]], shape = newGamma)
            df2 = df2 + 1
        }
        if (AllNNI) {
            fits <- optimPartNNI(fits, AllEdge)
            if(trace>0)cat(attr(fits, "swap"), " NNI operations performed")
            swap <- attr(fits, "swap")
        }
        if (AllEdge) {
            fits <- optimPartEdge(fits)
            df2 = df2 + length(fits[[1]]$tree$edge.length)
        }
        if (PartRate) {
            tree = fits[[1]]$tree
            rate = numeric(p)
            wp = numeric(p)
            for (i in 1:p) {
                wp[i] = sum(fits[[i]]$weight)
                rate[i] <- fits[[i]]$rate
            }
            ratemult = sum(wp)/sum(wp * rate)
            tree$edge.length = tree$edge.length/ratemult
            for (i in 1:p) fits[[i]] = update(fits[[i]], tree = tree, 
                rate = rate[i] * ratemult)
        }
        for (i in 1:p) lls[, i] = update(fits[[i]], data = dat)$site
        trees[[m]] = Gtrees
        LL = t(weight) %*% lls       
# choose partitions which change        
        tmp =(LL[cbind(1:ncw,part)] - apply(LL, 1, max))/colSums(weight)
        fixi = numeric(p)
        for(i in 1:p){
            tmpi = which(part == i)
            fixi[i] = tmpi[which.max(tmp[tmpi])]     
            }
        oldpart = part
# restrict the number of elements changing groups 
# If more than 25% would change, only the 25% with the highest increase per site change       
        if( sum(tmp==0)/length(tmp) < .75){
           medtmp = quantile(tmp, .25)
           medind = which(tmp<=medtmp)
           part[medind] = apply(LL[medind,], 1, which.max)
           }
        else part = apply(LL, 1, which.max)
# force groups to have at least one member
        part[fixi] = 1:p
        Part = cbind(Part, part)
        eps = sum(diag(table(part, oldpart)))
        eps2 = loli
        loli = sum(apply(LL, 1, max))
        eps2 = (eps2 - loli)/loli
        logLik = c(logLik, loli)
        if(trace>0) print(loli)
        Part = cbind(Part, part)
        df2 = df2 + df2
        if (eps == ncw & swap == 0) 
            AllNNI = FALSE
        m = m + 1
        if (eps == ncw) 
            iter = iter + 1
        if (iter == 3) 
            break
    }
    df <- matrix(1, 6, 2)
    colnames(df) <- c("#df", "group")
    rownames(df) <- c("Edge", "Shape", "Inv", "Bf", "Q", "Rate")
    df[1, 1] <- length(fits[[1]]$tree$edge.length)
    df[2, 1] <- fits[[1]]$k - 1
    df[3, 1] <- fits[[1]]$inv > 0
    df[4, 1] <- length(unique(fits[[1]]$bf)) - 1
    df[5, 1] <- length(unique(fits[[1]]$Q)) - 1
    df[6, 1] <- 0
    if (PartEdge) 
        df[1, 2] = p
    if (PartGamma) 
        df[2, 2] = p
    if (PartInv) 
        df[3, 2] = p
    if (PartBf) 
        df[4, 2] = p
    if (PartQ) 
        df[5, 2] = p
    if (PartRate) 
        df[6, 1] = p - 1
    attr(logLik, "df") = sum(df[, 1] * df[, 2])
    res = list(logLik = logLik, Partition = Part, trees = trees) # intermediate results
    result <- list(logLik = loli, fits = fits, Partition = part, df = df, res = res, call = call)
    class(result) <- c("pmlPart")
    result
}


pmlCluster <- function (formula, fit, weight, p = 1:5, part = NULL, nrep = 10, control = pml.control(epsilon = 1e-08,
   maxit = 10, trace = 1), ...)
{
   call <- match.call()
   form <- phangornParseFormula(formula)
   if(any(p==1)){
       opt2 <- c("nni", "bf", "Q", "inv", "shape", "edge")
       tmp1 <- opt2 %in% form$left
       tmp1 <- tmp1 | (opt2 %in% form$right)
       fit <- optim.pml(fit, tmp1[1], tmp1[2], tmp1[3], tmp1[4],
       tmp1[5], tmp1[6])
   }

   p=p[p!=1]
   if(length(p)==0)return(fit)
   n = sum(weight)
   k=2

   BIC = matrix(0, length(p)+1, nrep)
   BIC[1,] = AIC(fit, k = log(n))
   LL = matrix(NA, length(p)+1, nrep)
   LL[1,] = logLik(fit)

   P = array(dim=c(length(p)+1, nrep, dim(weight)[2]))
   tmpBIC = Inf
   choice = c(1,1) 
   for(j in p){
       tmp=NULL
       for(i in 1:nrep){
           tmp = pmlCluster.fit(formula, fit, weight, p=j, part=part, control=control,...)
           P[k,i,] = tmp$Partition
           BIC[k,i] = AIC(tmp, k = log(n))
           LL[k,i] = logLik(tmp)
           if(BIC[k,i]<tmpBIC){
                tmpBIC = BIC[k,i]
                result = tmp
                choice = c(k,i) 
           }
       }
       k=k+1
   }      

   p = c(1,p)
   result$choice = choice 
   result$BIC = BIC
   result$AllPartitions = P
   result$AllLL = LL
   result$p = p 
   class(result) = c("pmlCluster", "pmlPart")
   result
}


plot.pmlCluster <- function(x, which = c(1L:3L), caption = list("BIC", "log-likelihood", "Partitions"), ...){
   show <- rep(FALSE, 3)
   show[which] <- TRUE
   choice = x$choice
   if(show[1]){
       X <- x$AllPartitions[choice[1],,]
       d <- dim(X)
       ind = order(X[choice[2],])
       im  = matrix(0, d[2], d[2])
       for(j in 1:d[1]){for(i in 1:d[2]) im[i,] = im[i,] + (X[j,] == X[j,i]) }
       image(im[ind, ind], ...)
   }

   if(show[1])matplot(x$p, x$BIC, ylab="BIC", xlab="number of clusters")
   if(show[1])matplot(x$p, x$AllLL, ylab="log-likelihood", xlab="number of clusters")
}


readAArate <- function(file){
    tmp <- read.table(system.file(file.path("extdata", file), package = "phangorn"), col.names = 1:20, fill=TRUE)
    Q <- tmp[1:19,1:19]
    names <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l", "k", "m", "f", "p", "s", "t", "w",  "y", "v")
    Q <- as.numeric(Q[lower.tri(Q,TRUE)])
    bf <- as.numeric(as.character(unlist(tmp[20,])))
    names(bf) <- names
    list(Q=Q, bf=bf)
}

#.LG <- readAArate("lg.dat")
#.WAG <- readAArate("wag.dat")
#.Dayhoff <- readAArate("dayhoff-dcmut.dat")
#.JTT <- readAArate("jtt-dcmut.dat")
#.cpREV <- readAArate("cpREV.dat")
#.mtmam <- readAArate("mtmam.dat")
#.mtArt <- readAArate("mtArt.dat")
# save(.LG,.WAG,.Dayhoff,.JTT,.cpREV,.mtmam,.mtArt, file = "sysdata2.rda")


getModelAA <- function(model, bf=TRUE, Q=TRUE){
    model <- match.arg(eval(model), c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24"))
    tmp = get(paste(".", model, sep=""), environment(pml))
    if(Q) assign("Q", tmp$Q, envir=parent.frame())
    if(bf) assign("bf", tmp$bf, envir=parent.frame())
}


pml <- function (tree, data, bf = NULL, Q = NULL, inv = 0, k = 1, shape = 1, 
    rate = 1, model=NULL, ...) 
{
    call <- match.call()
    extras <- match.call(expand.dots = FALSE)$...
    pmla <- c("wMix", "llMix") 
    existing <- match(pmla, names(extras))
    wMix <- ifelse(is.na(existing[1]), 0, eval(extras[[existing[1]]], parent.frame()) )  
    llMix <- ifelse(is.na(existing[2]), 0, eval(extras[[existing[2]]], parent.frame()) )
  
    if(class(tree)!="phylo") stop("tree must be of class phylo") 
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorderPruning(tree)
    if (class(data)[1] != "phyDat") stop("data must be of class phyDat")
    if(is.null(tree$edge.length)) stop("tree must have edge weights") 
    if(any(is.na(match(tree$tip, attr(data, "names"))))) stop("tip labels are not in data")  
    levels <- attr(data, "levels")
    weight <- attr(data, "weight")
    nr <- attr(data, "nr")
    type <- attr(data,"type")
    if(type=="AA" & !is.null(model)){
        model <- match.arg(model, c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24"))
        getModelAA(model, bf=is.null(bf), Q=is.null(Q)) 
    }  
    if(type=="CODON") Q <- as.numeric(.syn > 0)
    if (is.null(bf)) 
        bf <- rep(1/length(levels), length(levels))
    if (is.null(Q)) 
        Q <- rep(1, length(levels) * (length(levels) - 1)/2)
    m <- 1
    eig <- edQt(bf = bf, Q = Q)
    w <- rep(1/k, k)
    if (inv > 0) 
        w <- (1 - inv) * w
    if (wMix > 0) 
        w <- wMix * w  
    g <- discrete.gamma(shape, k)
    if (inv > 0) 
        g <- g/(1 - inv)
    g <- rate * g
    INV <- Matrix(lli(data, tree), sparse=TRUE)
    ll.0 <- as.matrix(INV %*% (bf * inv))
    if(wMix>0) ll.0 <- ll.0 + llMix
    resll <- matrix(0, nr, k)
    while (m <= k) {
        resll[,m] = ll(data, tree, bf = bf, g = g[m], Q = Q, eig = eig, assign.dat = FALSE, ...)
        m = m + 1
    }
    lll <- resll %*% w
    siteLik <- lll + ll.0
    siteLik <- log(siteLik)
    loglik = sum(weight * siteLik)
    if(type=="CODON"){ 
        df <- length(tree$edge.length) + (k>1) + (inv > 0) + length(unique(bf)) - 1 
        }
    else df = length(tree$edge.length) + (k>1) + (inv > 0) + length(unique(bf)) - 1 + length(unique(Q)) - 1
    result = list(logLik = loglik, inv = inv, k = k, shape = shape,
        Q = Q, bf = bf, rate = rate, siteLik = siteLik, weight = weight, 
        g = g, w = w, eig = eig, data = data, model=model, INV = INV, 
        ll.0 = ll.0, tree = tree, lv = resll, call = call, df=df, wMix=wMix, llMix=llMix)
    if(type=="CODON"){
        result$dnds <- 1
        result$tstv <- 1
    }
    class(result) = "pml"
    result
}


print.pml = function(x,...){
    cat("\n loglikelihood:", x$logLik, "\n")
    w <- x$weight
    w <- w[w>0]    
    type <- attr(x$data, "type")
    levels <- attr(x$data, "levels")
    nc <- attr(x$data, "nc")
    ll0 = sum(w*log(w/sum(w)))
    cat("\nunconstrained loglikelihood:", ll0, "\n")
    if(x$inv > 0)cat("Proportion of invariant sites:",x$inv,"\n")
    if(x$k >1){
        cat("Discrete gamma model\n")
        cat("Number of rate categories:",x$k,"\n")        
        cat("Shape parameter:",x$shape,"\n")
        }
    if(type=="AA") cat("Rate matrix:",x$model, "\n")    
    if(type=="DNA"){
        cat("\nRate matrix:\n")    
        QM = matrix(0, nc, nc, dimnames = list(levels,levels))    
        QM[lower.tri(QM)] = x$Q    
        QM = QM+t(QM)
        print(QM)
        cat("\nBase frequencies:  \n")
        bf = x$bf
        names(bf) = levels 
        cat(bf, "\n")
    }
    if(type=="AA") cat("Rate matrix:",x$model, "\n")
    if(type=="CODON") {
         cat("dn/ds:",x$dnds, "\n")
         cat("ts/tv:",x$tstv, "\n") 
    }
    if(type=="USER" & length(x$bf)<11){         
        cat("\nRate matrix:\n")    
        QM = matrix(0, nc, nc, dimnames = list(levels,levels))    
        QM[lower.tri(QM)] = x$Q    
        QM = QM+t(QM)
        print(QM)
        cat("\nBase frequencies:  \n")
        bf = x$bf
        names(bf) = levels 
        cat(bf, "\n")
    }        
}


optEdgeMulti <- function (object, control = pml.control(epsilon = 1e-8, maxit = 10, trace=1), ...) 
{
    tree <- object$tree
    theta <- object$tree$edge.length
    weight <- attr(object$data, "weight")
    ll0 = object$logLik
    eps = 1
    iter = 0
    iter2 = 0
    scale = 1
    # l = length(theta)
    while (abs(eps) > control$eps && iter < control$maxit) {
        dl = score(object)
        thetaNew = log(theta) + scale * solve(dl[[2]], dl[[1]]) #+ diag(l)*1e-10
        newtheta = exp(thetaNew)
        tree$edge.length = as.numeric(newtheta)
        object <- update(object, tree = tree)
        ll1 = object$logLik 
        eps <- ( ll0 - ll1 ) / ll1 
        if(eps < 0){
             newtheta = theta
             scale = scale / 2
             tree$edge.length = as.numeric(theta)  
             ll1 = ll0  
             iter2 <- iter2+1             
        }
        else{
            scale=1
            iter2 = 0
        }  
        theta = newtheta 
        if(iter2==0 && control$trace>0) cat("loglik: ",ll1,"\n")
        ll0 <- ll1
        if(iter2==10)iter2=0  
        if(iter2==0)iter <- iter+1
    }
    object <- update(object, tree = tree) 
    object
}


# add data for internal use parent.frame(n) for higher nestings 
update.pmlNew <- function (object, ..., evaluate = TRUE){
    call <- object$call
    if (is.null(call)) 
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if (evaluate) 
        eval(call, object, parent.frame())
    else call
}


# call is now updated
update.pml <- function (object, ...) 
{
    extras <- match.call(expand.dots = FALSE)$...
    pmla <- c("tree", "data", "bf", "Q", "inv", "k", "shape", 
        "rate", "model", "wMix", "llMix", "...") 
    names(extras) <- pmla[pmatch(names(extras), pmla[-length(pmla)])]
    call = object$call
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }    
    existing <- match(pmla, names(extras))
    
    updateRates <- FALSE
    if (is.na(existing[1])) tree <- object$tree
    else tree <- eval(extras[[existing[1]]], parent.frame())
    if(is.null(attr(tree,"order")) || attr(tree,"order")=="cladewise")tree <- reorderPruning(tree)
    if (is.na(existing[2])){
        data <- object$data
        INV <- object$INV
        }
    else{ 
        data <- eval(extras[[existing[2]]], parent.frame())
        ll.0 <- numeric(attr(data,"nr"))
        INV <- Matrix(lli(data, tree), sparse=TRUE)
    }
    nr <- attr(data, "nr")  
    if (is.na(existing[3])) bf <- object$bf
    else bf <- eval(extras[[existing[3]]], parent.frame())
    if (is.na(existing[4])) Q <- object$Q
    else Q <- eval(extras[[existing[4]]], parent.frame())
#    model <- object$model
    type <- attr(object$data, "type")
    model<-NULL
    if (type == "AA") {
        if(!is.na(existing[9]) ){
        model <- match.arg(eval(extras[[existing[9]]], parent.frame()), c("WAG", "JTT", 
            "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", 
            "mtREV24"))
        getModelAA(model, bf = is.na(existing[3]), Q = is.na(existing[4]))
        } 
#        else model <- object$model
    }
   
    if(is.na(existing[5])) inv <- object$inv
    else{
        inv <- eval(extras[[existing[5]]], parent.frame())
        updateRates <- TRUE
    }
    if(is.na(existing[6])) k <- object$k
    else{
        k <- eval(extras[[existing[6]]], parent.frame())
        updateRates <- TRUE
    }
    if(is.na(existing[7])) shape <- object$shape
    else{
        shape <- eval(extras[[existing[7]]], parent.frame())
        updateRates <- TRUE
    }
    rate <- ifelse(is.na(existing[8]), object$rate, eval(extras[[existing[8]]], parent.frame()))
    wMix <- ifelse(is.na(existing[10]), object$wMix, eval(extras[[existing[10]]], parent.frame()))
    if(is.na(existing[11])) llMix <- object$llMix
    else llMix <- eval(extras[[existing[11]]], parent.frame())
    levels <- attr(data, "levels")
    weight <- attr(data, "weight")
    eig <- edQt(bf = bf, Q = Q)
    g <- discrete.gamma(shape, k)
    g <- rate * g 
    if (inv > 0) g <- g/(1 - inv)
    ll.0 <- as.matrix(INV %*% (bf * inv))
    if(wMix>0) ll.0 <- ll.0 + llMix
    w = rep(1/k, k)
    if (inv > 0) 
        w <- (1 - inv) * w
    if (wMix > 0) 
        w <- wMix * w                  
    m <- 1
    resll <- matrix(0, nr, k)
    
    while (m <= k) {
        resll[,m] = ll(data, tree, bf, g[m], Q, eig, assign.dat = FALSE)
        m = m + 1
    }
    lll = resll %*% w
    siteLik = log(lll + ll.0)
    loglik = sum(weight * siteLik)
    if(type=="CODON"){ 
        df <- length(tree$edge.length) + (k>1) + (inv > 0) + length(unique(bf)) - 1 + as.numeric(object$dnds != 1)    + as.numeric(object$tstv != 1)
        } 
    else df <- length(tree$edge.length) + (k>1) + (inv > 0) + length(unique(bf)) - 1 + length(unique(Q)) - 1
     result = list(logLik = loglik, inv = inv, k = k, shape = shape,
        Q = Q, bf = bf, rate = rate, siteLik = siteLik, weight = weight, 
        g = g, w = w, eig = eig, data = data, model=model, INV = INV, 
        ll.0 = ll.0, tree = tree, lv = resll, call = call, df=df, wMix=wMix, llMix=llMix)
    if(type=="CODON"){
        result$dnds <- object$dnds
        result$tstv <- object$tstv
    }

    class(result) = "pml"
    result
}


optimMixQ <- function(object, Q=c(1, 1, 1, 1, 1, 1), omega,...){
    l = length(Q)
    Q = Q[-l]
    Q = sqrt(Q)
    fn = function(Q, object, omega,...) {
        Q = c(Q^2, 1)
        weight <- object[[1]]$weight
        n <- length(omega)
        p <- length(weight)
        result <- numeric(p)
        for(i in 1:n)result <- result + as.numeric(update(object[[i]], Q=Q, ...)$lv) * omega[i]
        result <- sum(weight %*% log(result))
        result 
    }
    res = optim(par=Q, fn=fn, gr=NULL, method="L-BFGS-B", lower=0, 
            upper=Inf, control=list(fnscale = -1, maxit=25), 
            object=object, omega=omega,...)
    res[[1]] = c(res[[1]]^2, 1)
    res
}


optimMixBf <- function(object, bf=c(.25,.25,.25,.25), omega,...){
    l = length(bf)
    nenner = 1/bf[l]
    lbf = log(bf * nenner)
    lbf = lbf[-l]
    fn = function(lbf, object, omega,...) {
    bf = exp(c(lbf,0))
    bf = bf/sum(bf)
    weight <- object[[1]]$weight
        n <- length(omega)
        p <- length(weight)
        result <- numeric(p)
        for(i in 1:n)result <- result + as.numeric(update(object[[i]], bf=bf, ...)$lv) * omega[i]
        result <- sum(weight %*% log(result))
        result 
    }
    res = optim(par=lbf, fn=fn, gr=NULL, method="Nelder-Mead", 
        control=list(fnscale=-1, maxit=500), object, omega=omega,...)
    print(res[[2]])
    bf = exp(c(res[[1]],0))
    bf = bf/sum(bf)
}


optimMixInv <- function(object, inv=0.01, omega,...){
    fn = function(inv, object, omega,...) {
        n <- length(omega)
        weight <- object[[1]]$weight
        p <- length(weight)
        result <- numeric(p)
         for(i in 1:n)result <- result + as.numeric(update(object, inv=inv, ...)$lv) * omega[i]
        result <- sum(weight %*% log(result))
        result 
    }
    res = optimize(f=fn, interval = c(0,1), lower = 0, upper = 1, maximum = TRUE,
        tol = .0001, object, omega=omega,...)
    print(res[[2]]) 
    res[[1]]
}


pml2 <- function (tree, data, bf = rep(1/length(levels), length(levels)), 
    shape = 1, k = 1, Q = rep(1, length(levels) * (length(levels) - 1)/2), 
    levels = attr(data, "levels"), inv = 0, rate = 1, g = NULL, w = NULL, 
    eig = NULL, INV = NULL, ll.0 = NULL, llMix = NULL, wMix = 0, ...) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorderPruning(tree)
    if (class(data)[1] != "phyDat") 
        stop("data must be of class phyDat")
    weight = attr(data, "weight")
    l = length(bf)
    lll = matrix(0, length(weight), l)
    m = 1
    if (is.null(eig)) 
        eig = edQt(bf = bf, Q = Q)
    if (is.null(w)) {
        w = rep(1/k, k)
        if (inv > 0) 
            w <- (1 - inv) * w
        if (wMix > 0) 
            w <- (1 - wMix) * w           
    }
    if (is.null(g)) {
        g = discrete.gamma(shape, k)
        if (inv > 0) 
            g <- g/(1 - inv)
        g <- g * rate     
    } 
    if (is.null(INV)) 
        INV <- Matrix(lli(data, tree), sparse=TRUE)
#        INV = lli(data, tree)
    if (is.null(ll.0)) 
        ll.0 <- numeric(attr(data,"nr"))
    if(inv>0)
        ll.0 <- as.matrix(INV %*% (bf * inv)) 
#        ll.0 = INV %*% (inv * bf)               
    if (wMix > 0)
         ll.0 <- ll.0 + llMix           
    p = length(g)
    nr <- attr(data, "nr")   
    resll <- matrix(0, nr, k)
    while (m <= k) {
        resll[,m] = ll(data, tree, bf = bf, g = g[m], Q = Q, eig = eig, assign.dat = FALSE, ...)
        m = m + 1
    }
    lll <- resll %*% w         
    siteLik <- lll + ll.0
    if(wMix >0) siteLik <- siteLik * (1-wMix) + llMix
    siteLik = log(siteLik)
    sum(weight * siteLik)
}


optimMixRate <- function (fits, ll, weight, omega, rate=rep(1,length(fits))) 
{
    r <- length(fits)
    rate0 <- rate[-r]   

    fn<-function(rate, fits, ll, weight, omega){
        r <-  length(fits)
        rate <- c(rate, (1- sum(rate *omega[-r]))/omega[r])
        for (i in 1:r) fits[[i]]<-update(fits[[i]], rate = rate[i])
        for (i in 1:r) ll[, i] <- fits[[i]]$lv
        sum(weight*log(ll%*%omega)) 
    }
    ui=diag(r-1)
    ui <- rbind(-omega[-r], ui)
    ci <- c(-1, rep(0, r-1))
    res <- constrOptim(rate0, fn, grad=NULL, ui=ui, ci=ci, mu = 1e-04, control = list(fnscale=-1),
        method = "Nelder-Mead", outer.iterations = 50, outer.eps = 1e-05, fits=fits, ll=ll, weight=weight, omega=omega)
    rate <- res[[1]]
    res[[1]] <- c(rate, (1- sum(rate *omega[-r]))/omega[r])
    res
}


optW <- function (ll, weight, omega,...) 
{
    k = length(omega)
    nenner = 1/omega[1]
    eta = log(omega * nenner)
    eta = eta[-1]
    fn = function(eta, ll, weight) {
        eta = c(0,eta)
        p = exp(eta)/sum(exp(eta))
        res = sum(weight * log(ll %*% p)) 
        res
    }
    if(k==2)res = optimize(f =fn , interval =c(-3,3) , lower = -3, upper = 3, maximum = TRUE, tol = .Machine$double.eps^0.25, ll = ll, weight = weight) 
    else res = optim(eta, fn = fn, method = "L-BFGS-B", lower=-5, upper=5,control = list(fnscale = -1, 
        maxit=25), gr = NULL, ll = ll, weight = weight)

    p = exp(c(0,res[[1]]))
    p = p/sum(p)
    result = list(par = p, value = res[[2]])
    result
}


optimMixEdge <- function(object, omega, trace=1,...){
    tree <- object[[1]]$tree
    theta <- object[[1]]$tree$edge.length
    weight = as.numeric(attr(object[[1]]$data,"weight"))
    n <- length(omega)
    p <- length(weight)
    q <- length(theta)
    lv1 = numeric(p)
    for(i in 1:n) lv1 = lv1 + as.numeric(object[[i]]$lv) * omega[i]
    ll0 <- sum(weight * log(lv1))
    eps=1
    iter <- 0
    scalep <- 1
    if(trace>0) cat(ll0)
    while(abs(eps)>.001 & iter<10){
        dl <- matrix(0,p,q)
        for(i in 1:n)dl <- dl + dl(object[[i]],TRUE) * omega[i]
        dl <- dl/lv1
        sc = colSums(weight * dl)
        F = crossprod(dl * weight, dl)+diag(q)*1e-6
        blub <- TRUE
        iter2 <- 0
        while(blub & iter2<10){
        thetaNew = log(theta) + scalep * solve(F, sc)
        tree$edge.length = as.numeric(exp(thetaNew))
        for(i in 1:n)object[[i]] <- update(object[[i]],tree=tree)
        lv1 = numeric(p)
        for(i in 1:n) lv1 = lv1 + as.numeric(object[[i]]$lv)  * omega[i]
        ll1 <- sum(weight * log(lv1))
        eps <- ll1 - ll0     
        if (eps < 0) {
            scalep = scalep/2
            eps = 1
            thetaNew = log(theta)
            ll1 = ll0
            iter2 <- iter2+1
        }
        else{
             scalep = 1;
             theta = exp(thetaNew)  
             blub=FALSE  
            }     
        }             
        iter <- iter+1
        ll0 <- ll1
    }       
    tree$edge.length <- theta
    for(i in 1:n)object[[i]] <- update(object[[i]],tree=tree)
    if(trace>0) cat("->", ll1, "\n")
    object
}


pmlMix <- function (formula, fit, m = 2, omega = rep(1/m, m), control=pml.control(epsilon=1e-8, maxit=20, trace=1), ...) 
{
    call <- match.call()
    form <- phangornParseFormula(formula)
    opt <- c("nni", "bf", "Q", "inv", "shape", "edge", "rate")
    optAll <- match(opt, form$left)
    optPart <- match(opt, form$right)
    AllBf <- !is.na(optAll[2])
    AllQ <- !is.na(optAll[3])
    AllInv <- !is.na(optAll[4])
    AllGamma <- !is.na(optAll[5])
    AllEdge <- !is.na(optAll[6])
    MixNni <- !is.na(optPart[1])
    MixBf <- !is.na(optPart[2])
    MixQ <- !is.na(optPart[3])
    MixInv <- !is.na(optPart[4])
    MixGamma <- !is.na(optPart[5])
    MixEdge <- !is.na(optPart[6])
    MixRate <- !is.na(optPart[7])
    if (class(fit) == "list") 
        fits <- fit
    else {
        fits <- vector("list", m) 
        for (i in 1:m) fits[[i]] <- fit
    }
    dat <- fits[[1]]$data
    p = attr(dat, "nr")
    weight = attr(dat, "weight")
    r = m
    ll = matrix(0, p, r)
    for (i in 1:r) ll[, i] = fits[[i]]$lv

    for (i in 1:r){
         pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
         fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
    }

    if(MixRate) rate <- rep(1,r)

    llstart = sum(weight * log(ll %*% omega))
    llold <- llstart
    ll0 <- llstart
    ll3 <- llstart
    eps0 <- 1
    iter0 <- 0
    trace = control$trace
    while (eps0 > control$eps & iter0 < control$maxit) {  #while (eps0 > 1e-6 & iter0 < 20) {
        eps1 <- 100
        iter1 <- 0
        
        if (AllQ) {
            newQ <- optimMixQ(fits, Q = fits[[1]]$Q, 
                omega = omega)[[1]]
            for (i in 1:m) fits[[i]] <- update(fits[[i]], Q = newQ)
        }
        if (AllBf) {
            newBf <- optimMixBf(fits, bf = fits[[1]]$bf, 
                omega = omega)
            for (i in 1:m) fits[[i]] <- update(fits[[i]], bf = newBf)
        }
        if (AllInv) {
            newInv <- optimMixInv(fits, inv = fits[[1]]$inv, 
                omega = omega)
            for (i in 1:m) fits[[i]] <- update(fits[[i]], Inv = newInv)
        }
        if (AllEdge) 
            fits <- optimMixEdge(fits, omega, trace=trace-1)
        for (i in 1:r) ll[, i] <- fits[[i]]$lv

        while ( abs(eps1) > 0.001 & iter1 < 3) {
             if(MixRate){
                 rate <- optimMixRate(fits, ll, weight, omega, rate)[[1]]
                 for (i in 1:r) fits[[i]] <- update(fits[[i]], rate=rate[i]) 
                 for (i in 1:r) ll[, i] <- fits[[i]]$lv
            }
            for (i in 1:r){
                pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
                fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
            }

            for (i in 1:r) {
                pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
                fits[[i]] <- optim.pml(fits[[i]], MixNni, MixBf, MixQ, MixInv, MixGamma, 
                    MixEdge, optRate=FALSE, control = pml.control(epsilon = 1e-8, maxit = 3,
                    trace-1), llMix = pl0, wMix = omega[i])
                 ll[, i] <- fits[[i]]$lv 

            res = optW(ll, weight, omega)
               omega = res$p
            
            if(MixRate){
                blub <- sum(rate*omega)
                rate <- rate / blub 
                tree <- fits[[1]]$tree
                tree$edge.length <-   tree$edge.length*blub
                for (i in 1:r) fits[[i]]<-update(fits[[i]], tree=tree, rate = rate[i])
                for (i in 1:r) ll[, i] <- fits[[i]]$lv
             }
             for (i in 1:r){
                 pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
                 fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
             }
             
         }
         ll1 = sum(weight * log(ll %*% omega))
         res = optW(ll, weight, omega)
         omega = res$p
         if(MixRate){
                blub <- sum(rate*omega)
                rate <- rate / blub 
                tree <- fits[[1]]$tree
                tree$edge.length <-   tree$edge.length*blub
                for (i in 1:r) fits[[i]]<-update(fits[[i]], tree=tree, rate = rate[i])
                     if(trace>0) print(rate)
                     for (i in 1:r) ll[, i] <- fits[[i]]$lv
                }
         for (i in 1:r){
             pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
             fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
        }

        ll2 = sum(weight * log(ll %*% omega)) 
        eps1 = llold - ll2
        iter1 <- iter1 + 1
        llold = ll2
        }   

        ll1 <- sum(weight * log(ll %*% omega))
        eps0 <- (ll3 - ll1) / ll1
        ll3 <- ll1
        iter0 <- iter0 + 1
        if(trace>0) print(iter0)
    }
    parameter <- c(AllBf=AllBf, AllQ=AllQ, AllInv=AllInv, AllGamma=AllGamma, AllEdge=AllEdge, MixNni=MixNni, 
       MixBf=MixBf, MixQ=MixQ, MixInv=MixInv, MixGamma=MixGamma, MixEdge=MixEdge, MixRate=MixRate)
    
    df <- matrix(1, 6 ,2)
    colnames(df) <- c("#df", "group")
    rownames(df) <- c("Edge", "Shape", "Inv", "Bf", "Q", "Rate")
    df[1,1] <- length(fits[[1]]$tree$edge.length)
#    df[2,1] <- fits[[1]]$k - 1     
    df[2,1] <- fits[[1]]$k > 1
    df[3,1] <- fits[[1]]$inv > 0
    df[4,1] <- length(unique(fits[[1]]$bf)) - 1
    df[5,1] <- length(unique(fits[[1]]$Q)) - 1
    df[6,1] <- 0  
    if(MixEdge) df[1,2] = r
    if(MixGamma) df[2,2] = r
    if(MixInv) df[3,2] = r
    if(MixBf) df[4,2] = r
    if(MixQ) df[5,2] = r
    if(MixRate) df[6,1] = r-1     
    attr(logLik, "df") = sum(df[,1]*df[,2])
    converge <- c(iter=iter0, eps=eps0)
    result <- list(logLik = ll1, omega = omega, fits = fits, call = call, converge=converge, parameter=parameter, df=df)
    class(result) <- "pmlMix"
    result
}


print.pmlMix <- function(x,...){
    nc <- attr(x$fits[[1]]$data, "nc")
    nr <- attr(x$fits[[1]]$data, "nr")
    levels <- attr(x$fits[[1]]$data, "levels")
    r <- length(x$fits)   
    w <- x$fits[[1]]$weight
    w <- w[w>0] 
    type <- attr(x$fits[[1]]$data, "type")
    nc <- attr(x$fits[[1]]$data, "nc")
    ll0 = sum(w*log(w/sum(w)))

    
    bf <- matrix(0,r,nc)
    dimnames(bf) <- list(1:r, levels)
    Q <- matrix(0, r, nc*(nc-1)/2)
    dimnames(Q) <- list(1:r, NULL)

    rate <- numeric(r)
    inv <- x$fits[[1]]$inv
    shape <- numeric(r)

    for(i in 1:r){
        bf[i, ] <- x$fits[[i]]$bf
        Q[i, ] <- x$fits[[i]]$Q
        rate[i] <- x$fits[[i]]$rate
        shape[i] <- x$fits[[i]]$shape
    }
    cat("\nloglikelihood:", x$logLik, "\n")
    cat("\nunconstrained loglikelihood:", ll0, "\n") 
    cat("AIC: ", AIC(x), " BIC: ", AIC(x, k=log(nr)), "\n\n")
    cat("\nposterior:", x$omega ,"\n")   
    if(inv > 0)cat("Proportion of invariant sites:",inv,"\n")
    cat("\nRates:\n")
    cat(rate,"\n")
    cat("\nBase frequencies:  \n")
    print(bf)
    cat("\nRate matrix:\n")
    print(Q)
}


logLik.pmlMix <- function (object, ...) 
{
    res <- object$logLik
    attr(res, "df") <- sum(object$df[,1] * object$df[,2])
    class(res) <- "logLik"
    res
}
 

print.pmlPart <- function(x,...){
    df <- x$df
    nc <- attr(x$fits[[1]]$data, "nc")
    levels <- attr(x$fits[[1]]$data, "levels")
    r <- length(x$fits)   
    nc <- attr(x$fits[[1]]$data, "nc")
    nr <- attr(x$fits[[1]]$data, "nr")
    k <- x$fits[[1]]$k    

    lbf=x$df["Bf",2]
    bf <- matrix(0, lbf, nc)
    if(lbf>1)dimnames(bf) <- list(1:r, levels)
    lQ = x$df["Q",2]
    Q <- matrix(0, lQ, nc*(nc-1)/2)
    if(lQ>1)dimnames(Q) <- list(1:r, NULL)
    type <- attr(x$fits[[1]]$data, "type")
    
    loli <- numeric(r)
    rate <- numeric(r)
    shape <- numeric(r)
    sizes <- numeric(r)
    inv <- numeric(r)      
    for(i in 1:r){
        loli[i] <- x$fits[[i]]$logLik
        if(i <= lbf)bf[i, ] <- x$fits[[i]]$bf
        if(i <= lQ)Q[i, ] <- x$fits[[i]]$Q
        rate[i] <- x$fits[[i]]$rate
        shape[i] <- x$fits[[i]]$shape
        inv[i] <- x$fits[[i]]$inv
        sizes[i] <- sum(attr(x$fits[[i]]$data,"weight"))
    }
    cat("\nloglikelihood:", x$logLik, "\n")
    cat("\nloglikelihood of partitions:\n ", loli, "\n")
    cat("AIC: ", AIC(x), " BIC: ", AIC(x, k=log(sum(sizes))), "\n\n")    
    cat("Proportion of invariant sites:",inv,"\n")
    cat("\nRates:\n")
    cat(rate,"\n")
    if(k>1){
        cat("\nShape parameter:\n") 
        cat(shape,"\n")
    }
    if(type=="AA") cat("Rate matrix:",x$fits[[1]]$model, "\n")
    else{
        cat("\nBase frequencies:  \n")
        print(bf)
        cat("\nRate matrix:\n")
        print(Q)
    }
}


logLik.pmlPart <- function (object, ...) 
{
    res <- object$logLik
    attr(res, "df") <- sum(object$df[,1] * object$df[,2])
    class(res) <- "logLik"
    res
}


pmlPen <- function(object, lambda, ...){
    if(class(object)=="pmlPart") return(pmlPartPen(object, lambda,...))
    if(class(object)=="pmlMix") return(pmlMixPen(object, lambda,...))
    else stop("object has to be of class pmlPart or pmlMix")
    }
       
   
pmlPartPen <- function(object, lambda, control=pml.control(epsilon=1e-8, maxit=20, trace=1),...){
    fits <- object$fits
    
    m <- length(fits)
    K = -diag(length(fits[[1]]$tree$edge.length))
    Ktmp=K
    for(i in 1:(m-1))Ktmp = cbind(Ktmp,K)
    KM = Ktmp
    for(i in 1:(m-1))KM = rbind(KM,Ktmp)
    diag(KM) = m-1
    theta=NULL
    l = length(fits[[1]]$tree$edge.length)
    loglik = 0
    for(i in 1:m){
        theta = c(theta,fits[[i]]$tree$edge.length)
        loglik = loglik + fits[[i]]$logLik
    }
    print(loglik)
    pen = - 0.5 * lambda * t(theta)%*%KM%*%theta
    loglik = loglik - 0.5 * lambda * t(theta)%*%KM%*%theta 
    eps=1
    H  = matrix(0, m * l, m * l)
    iter=0
    trace = control$trace
    while( abs(eps)>control$eps & iter<control$maxit){
        theta=NULL
        sc = NULL
        for(i in 1:m){
            theta = c(theta,fits[[i]]$tree$edge.length)
            scoretmp = score(fits[[i]], TRUE)
            sc = c(sc,scoretmp$sc)
            H[(1:l)+l*(i-1), (1:l)+l*(i-1)] = scoretmp$F
        }
        sc = sc - lambda * KM%*% log(theta)
        thetanew = log(theta) +  solve(H + lambda*KM, sc)
        for(i in 1:m) fits[[i]]$tree$edge.length = exp(thetanew[(1:l)+(i-1)*l])
        for(i in 1:m) fits[[i]] = update.pml(fits[[i]], tree=fits[[i]]$tree)
        loglik1 = 0
        for(i in 1:m) loglik1 = loglik1 + fits[[i]]$logLik
        logLik <- loglik1
        if(trace>0)print(loglik1)
        loglik0 = loglik1
        pen = - 0.5 * lambda * t(theta)%*%KM%*%theta
        loglik1 = loglik1 - 0.5 * lambda * t(thetanew)%*%KM%*%thetanew
        eps =  (loglik - loglik1) / loglik1   
        loglik = loglik1
        theta = exp(thetanew)
        iter = iter+1
        if(trace>0)print(iter)
    }
    df = sum( diag(solve(H + lambda* KM, H)))
    
    object$df[1,1] = df
    object$df[1,2] = 1
    object$fits = fits
    object$logLik = loglik0
    attr(object$logLik, "df") = sum(object$df[,1]*object$df[,2])
    object$logLik.pen = loglik
    attr(object$logLik.pen, "df") = sum(object$df[,1]*object$df[,2])      
    object
}


pmlMixPen = function (object, lambda, optOmega=TRUE, control=pml.control(epsilon=1e-8, maxit=20, trace=1), ...) 
{
    fits <- object$fits
    m <- length(fits)
    K = -diag(length(fits[[1]]$tree$edge.length))
    tree <- fits[[1]]$tree
    Ktmp = K
    for (i in 1:(m - 1)) Ktmp = cbind(Ktmp, K)
    KM = Ktmp
    for (i in 1:(m - 1)) KM = rbind(KM, Ktmp)
    diag(KM) = m - 1
    theta = NULL
    l = length(fits[[1]]$tree$edge.length)
    omega <- object$omega
    dat <- fits[[1]]$data
    nr = attr(dat, "nr")
    weight = drop(attr(dat, "weight"))
    ll = matrix(0, nr, m)
    for (i in 1:m) ll[, i] = fits[[i]]$lv
    lv = drop(ll %*% omega)
    loglik = sum(weight * log(lv))
    for (i in 1:m) theta = c(theta, fits[[i]]$tree$edge.length)
    pen = - 0.5 * lambda * t(theta) %*% KM %*% theta
    loglik = loglik + pen
    print(loglik)    
    eps0 = 1 
    dl <- matrix(0, nr, m * l)
    iter0 = 0
    trace = control$trace 
    while (abs(eps0) > control$eps & iter0 < control$maxit) {
      eps = 1
      iter = 0      
      while (abs(eps) > 0.01 & iter < 5) {
        for (i in 1:m) {
            dl[, (1:l) + l * (i - 1)] <- dl(fits[[i]], TRUE) * 
                omega[i]
        }
        dl <- dl/lv
        sc = colSums(weight * dl) - lambda * KM %*% log(theta)
        H = crossprod(dl * weight, dl)
        thetanew = log(theta) + solve(H + lambda * KM, sc)
        for (i in 1:m) fits[[i]]$tree$edge.length = exp(thetanew[(1:l) + 
            (i - 1) * l])
        for (i in 1:m) {
            tree$edge.length = exp(thetanew[(1:l) + (i - 1) * l])
            fits[[i]] = update.pml(fits[[i]], tree = tree)
            ll[, i] = fits[[i]]$lv
        }
        lv = drop(ll %*% omega)
        loglik1 = sum(weight * log(lv))
        pen =  - 0.5 * lambda * t(thetanew) %*% KM %*% thetanew
        loglik1 = loglik1 + pen
        eps = abs(loglik1 - loglik)
        theta = exp(thetanew)
        loglik <- loglik1
        iter = iter + 1  
       }
       if(optOmega){
            res = optWPen(ll, weight, omega, pen)
            omega = res$p
            for (i in 1:m) {
                pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
                fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
                }
            } 
        lv = drop(ll %*% omega)
        loglik1 = sum(weight * log(lv))
        loglik0 =loglik1
        loglik1 = loglik1 - 0.5 * lambda * t(thetanew) %*% KM %*% thetanew
        eps0 = (loglik - loglik1) / loglik1
        theta = exp(thetanew)
        loglik <- loglik1
        iter0 = iter0 + 1
        if(trace>0) print(loglik)  
    }

    for (i in 1:m) {
        pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
        fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
    }
    df = sum(diag(solve(H + lambda * KM, H)))
    penalty <- list(lambda=lambda, K=KM, thetanew=thetanew, ploglik=loglik)
    object$omega = omega
    object$df[1, 1] = df
    object$df[1, 2] = 1
    object$fits = fits
    object$logLik = loglik0
    object$penalty = penalty
    object
}


optWPen = function (ll, weight, omega, pen, ...) 
{
    k = length(omega)
    nenner = 1/omega[1]
    eta = log(omega * nenner)
    eta = eta[-1]
    fn = function(eta, ll, weight, pen) {
        eta = c(0, eta)
        p = exp(eta)/sum(exp(eta))
        res = sum(weight * log(ll %*% p)) + pen
        res
    }
    if (k == 2) 
        res = optimize(f = fn, interval = c(-3, 3), lower = -3, 
            upper = 3, maximum = TRUE, tol = .Machine$double.eps^0.25, 
            ll = ll, weight = weight, pen = pen)
    else res = optim(eta, fn = fn, method = "L-BFGS-B", lower = -5, 
        upper = 5, control = list(fnscale = -1, maxit = 25), 
        gr = NULL, ll = ll, weight = weight, pen=pen)
    p = exp(c(0, res[[1]]))
    p = p/sum(p)
    result = list(par = p, value = res[[2]])
    result
} 


optNNI <- function(fit, INDEX){    
       tree = fit$tree
       ll.0 <- fit$ll.0
       loli <- fit$logLik
       bf = fit$bf
       eig = fit$eig
       k = fit$k
       w = fit$w
       g = fit$g
       rootEdges <- attr(INDEX, "root")
       .dat <- NULL
       parent = tree$edge[, 1]
       child = tree$edge[, 2]
       datp = rnodes(fit) # raus
       evector <- numeric(max(parent))
       evector[child] <- tree$edge.length
       m <- dim(INDEX)[1]
       k = min(parent)
       loglik = numeric(2 * m)
       edgeMatrix <- matrix(0, 2 * m, 5)
       for (i in 1:m) {
           ei = INDEX[i, ]
           el0 = evector[INDEX[i, ]]
           l = length(datp[, 1])
           weight = fit$weight
           datn = vector("list", 4 * l)
           attr(datn, "dim") = c(l, 4)
           datn <- .dat[, ei[1:4], drop = FALSE]
           if (!(ei[5] %in% rootEdges)) 
                datn[, 1] = datp[, ei[1], drop = FALSE]
           new1 <- optim.quartet(el0[c(1, 3, 2, 4, 5)], 
               eig, bf, datn[, c(1, 3, 2, 4), drop = FALSE], g, 
               w, weight, ll.0, llcomp = fit$log)
           new2 <- optim.quartet(el0[c(1, 4, 3, 2, 5)], 
               eig, bf, datn[, c(1, 4, 3, 2), drop = FALSE], g, 
               w, weight, ll.0, llcomp = fit$log)
           loglik[(2 * i) - 1] = new1[[2]]
           loglik[(2 * i)] = new2[[2]]
           edgeMatrix[(2 * i) - 1, ] = new1[[1]]
           edgeMatrix[(2 * i), ] = new2[[1]]
           }
       list(loglik=loglik, edges = edgeMatrix)
       }


optimPartNNI <- function (object, AllEdge=TRUE,...) 
{
    tree <- object[[1]]$tree
    INDEX <- indexNNI(tree)   
    l = length(object)
    loglik0 = 0
    for(i in 1:l)loglik0 = loglik0 + logLik(object[[i]])    
    
    l = length(object)
    TMP=vector("list", l)
    for(i in 1:l){
        TMP[[i]] = optNNI(object[[i]], INDEX)
        }
    loglik=TMP[[1]][[1]] 
    for(i in 2:l)loglik=loglik+TMP[[i]][[1]]

    swap <- 0
    candidates <- loglik > loglik0

    while (any(candidates)) {
        ind = which.max(loglik)
        loglik[ind] = -Inf
        if (ind%%2) 
            swap.edge = c(2, 3)
        else swap.edge = c(2, 4)
        tree2 <- changeEdge(tree, INDEX[(ind + 1)%/%2, swap.edge], 
            INDEX[(ind + 1)%/%2, ], TMP[[1]][[2]][ind, ])
        tmpll = 0                 
        for(i in 1:l){
            if(!AllEdge)tree2 <- changeEdge(object[[i]]$tree, INDEX[(ind + 1)%/%2, swap.edge], 
                INDEX[(ind + 1)%/%2, ], TMP[[i]][[2]][ind, ]) 
            tmpll <- tmpll + update(object[[i]], tree = tree2)$logLik
            }

        if (tmpll < loglik0) 
            candidates[ind] = FALSE
        if (tmpll > loglik0) {

            swap = swap + 1
            tree <- tree2
            indi <- which(rep(colSums(apply(INDEX, 1, match, 
                INDEX[(ind + 1)%/%2, ], nomatch = 0)) > 0, each = 2))
            candidates[indi] <- FALSE
            loglik[indi] <- -Inf

            for(i in 1:l){
                if(!AllEdge)tree2 <- changeEdge(object[[i]]$tree, INDEX[(ind + 1)%/%2, swap.edge], 
                    INDEX[(ind + 1)%/%2, ], TMP[[i]][[2]][ind, ]) 
                object[[i]] <- update(object[[i]], tree = tree2)
                }
            loglik0 = 0
            for(i in 1:l)loglik0 = loglik0 + logLik(object[[i]])    
            cat(loglik0, "\n")
        }
    }
    if(AllEdge)object <- optimPartEdge(object)
    attr(object,"swap") = swap
    object
}


#
# add codon models, change to phyDat statt 3* 
#
simSeq = function(tree, l=1000, Q=NULL, bf=NULL, rootseq=NULL, type = "DNA", model="USER",
    levels = NULL, rate=1, ancestral=FALSE){
    
    pt <- match.arg(type, c("DNA", "AA", "USER"))
    if (pt == "DNA") 
        levels <- c("a", "c", "g", "t")
    if (pt == "AA") 
        levels <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", 
        "l", "k", "m", "f", "p", "s", "t", "w", "y", "v")
    if (pt == "USER") 
        if(is.null(levels))stop("levels have to be supplied if type is USER")
     
    lbf = length(levels)
    
    if (type == "AA" & !is.null(model)) {
        model <- match.arg(model, c("USER", "WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24"))
        if(model!="USER")getModelAA(model, bf=is.null(bf), Q=is.null(Q))
    }

    if(is.null(bf)) bf = rep(1/lbf,lbf)
    if(is.null(Q)) Q = rep(1,lbf*(lbf-1)/2)
    if(is.matrix(Q)) Q=Q[lower.tri(Q)]
    eig = edQt(Q, bf)
    
    m = length(levels)    
    
    if(is.null(rootseq))rootseq = sample(levels, l, replace=TRUE, prob=bf)
    tree = ape:::reorder.phylo(tree) #
    edge = tree$edge
    nNodes = max(edge)
    res = matrix(NA, l, nNodes)
    parent <- as.integer(edge[, 1])
    child <- as.integer(edge[, 2])
    root <- as.integer(parent[!match(parent, child, 0)][1])  
    res[, root] = rootseq   
    tl = tree$edge.length
    for(i in 1:length(tl)){
        from = parent[i] 
        to = child[i]
        P = getP(tl[i], eig, rate)[[1]]
        for(j in 1:m){
            ind = res[,from]==levels[j]
            res[ind,to] = sample(levels, sum(ind), replace=TRUE, prob=P[,j])
        }
    }
    k = length(tree$tip)
    label = c(tree$tip, as.character((k+1):nNodes))
    colnames(res)=label 
    if(!ancestral)res = res[, tree$tip, drop=FALSE]
    if(pt=="DNA") return(phyDat.DNA(as.data.frame(res), return.index=TRUE))
    if(pt=="AA") return(phyDat.AA(as.data.frame(res), return.index=TRUE))
    if(pt=="USER") return(phyDat.default(as.data.frame(res), levels = levels, return.index=TRUE))
}        

      
SH.test <- function (..., B = 10000, data = NULL)
{
   fits <- list(...)
   p = 1
   if (class(fits[[1]]) == "pmlPart") {
       fits = fits[[1]]$fits
       p = length(fits)
   }
   k = length(fits)
   if (is.null(data))
       data = fits[[1]]$data
   res = NULL
   for (h in 1:p) {
       if (p > 1)
           data = fits[[h]]$data
       weight = attr(data, "weight")
       lw = length(weight)
       siteLik = matrix(0, lw, k)
       for (i in 1:k) siteLik[, i] = update(fits[[i]], data = data)$site
       ntree = k
       Lalpha <- drop(crossprod(siteLik, weight))
       Talpha <- max(Lalpha) - Lalpha
       M <- matrix(NA, k, B)
#        S <- matrix(NA, k, B)
       wvec <- rep(1L:lw, weight)
       size = length(wvec)
       for (i in 1:B) {
           boot = tabulate(sample(wvec, replace=TRUE), nbins=lw)
           M[, i] <- crossprod(siteLik, boot)
       }
       M <- M - rowMeans(M)
#        for (i in 1:B) for (j in 1:ntree) S[j, i] <- max(M[j, i] - M[, i])
       S = matrix(apply(M,2,min), k, B, byrow=TRUE)
       S = M - S
       count <- numeric(ntree)
       for (j in 1:ntree) count[j] <- sum(S[j, ] > Talpha[j])
       count <- count/B
       trees <- 1:k
       if (p == 1)
           res = cbind(trees, Lalpha, Talpha, count)
       else res = rbind(res, cbind(h, trees[-h], Lalpha[-h],
           Talpha[-h], count[-h]))
   }
   if (p == 1)
       colnames(res) <- c("Trees", "ln L", "Diff ln L", "p-value")
   else colnames(res) <- c("Partition", "Trees", "ln L", "Diff ln L",
       "p-value")
   res
}


#
# Bootstrap functions 
# multicore support
#
bootstrap.pml = function (x, bs = 100, trees = TRUE, multicore=FALSE,  ...) 
{
    data = x$data
    weight = attr(data, "weight")
    v = rep(1:length(weight), weight)
    BS = vector("list", bs)
    for (i in 1:bs) BS[[i]] = tabulate(sample(v, replace = TRUE), 
        length(weight))
    pmlPar <- function(weights, fit, trees = TRUE, ...) {
        data = fit$data
        ind <- which(weights > 0)
        data <- getRows(data, ind)
        attr(data, "weight") <- weights[ind]
        fit = update(fit, data = data)
        fit = optim.pml(fit, ...)
        if (trees) {
            tree = fit$tree
            return(tree)
        }
        attr(fit, "data") = NULL
        fit
    }
    eval.success <- FALSE
    if (!eval.success & multicore) {
        if (!require(parallel) || .Platform$GUI!="X11") {
            warning("package 'parallel' not found or GUI is used, 
            bootstrapping is performed in serial")
        } else {       
            res <- mclapply(BS, pmlPar, x, trees = trees, ...)
            eval.success <- TRUE
        } 
    }
    if (!eval.success) res <- lapply(BS, pmlPar, x, trees = trees, ...)
    if (trees) {
        class(res) = "multiPhylo"
    }
    res
}


bootstrap.phyDat <- function (x, FUN, bs = 100, multicore=FALSE, ...) 
{
    weight = attr(x, "weight")
    v = rep(1:length(weight), weight)
    BS = vector("list", bs)
    for(i in 1:bs)BS[[i]]=tabulate(sample(v, replace=TRUE),length(weight)) 
 
    fitPar <- function(weights, data, ...){     
         ind <- which(weights > 0)
         data <- getRows(data, ind)
         attr(data, "weight") <- weights[ind]
         FUN(data,...)        
    }
    eval.success <- FALSE
    if (!eval.success & multicore) {         
        if (!require(parallel) || .Platform$GUI!="X11") {
            warning("package 'parallel' not found or GUI is used, 
            bootstrapping is performed in serial")
        } else { 
            res <- mclapply(BS, fitPar, x, ...) 
            eval.success <- TRUE
        } 
    }
    if (!eval.success) res <- lapply(BS, fitPar, x, ...) 
    res 
}


matchEdges = function(tree1, tree2){
    bp1 = bip(tree1)
    bp2 = bip(tree2)
    l = length(tree1$tip)
    fn = function(x, y){
        if(x[1]==1)return(x)
        else return(y[-x])
        } 
    bp1[] = lapply(bp1, fn, 1:l)
    bp2[] = lapply(bp2, fn, 1:l)
    match(bp1, bp2)
}


plotBS = function(tree, BStrees, type="unrooted", bs.col="black", bs.adj=NULL, ...){
	if(type=='phylogram' | type=='cladogram'){
		if(!is.rooted(tree)) tree2 = midpoint(tree)
		plot(tree2, type=type, ...)	   
		}
    else plot(tree, type=type, ...)
    x = prop.clades(tree, BStrees)
    x = round( (x / length(BStrees)) * 100 )
    label = c(rep("", length(tree$tip)),x)
    ind <- get("last_plot.phylo", envir = .PlotPhyloEnv)$edge[,2]
    if(type=='phylogram' | type=='cladogram'){
	    root = getRoot(tree)
	    label = c(rep(0, length(tree$tip)),x)
	    label[root] = 0
	    ind2 = matchEdges(tree2,tree)
	    label = label[ind2]
	    ind = which(label > 0)
	    if(is.null(bs.adj))bs.adj = c(1, 1)
	    nodelabels(text=label[ind], node=ind, frame="none", col=bs.col, adj=bs.adj, ...)
	    }
    else{
	     if(is.null(bs.adj))bs.adj = c(0.5, 0.5)
	     edgelabels(label[ind], frame="none", col=bs.col, adj=bs.adj, ...)
     }
}


#
# splits format, networx, Matrix, lento plot 
#
as.splits <- function (x, ...){
    if(inherits(x, "splits")) return(x)
    UseMethod("as.splits")
}


as.Matrix <- function (x, ...){
    if (class(x) == "Matrix") return(x)
    UseMethod("as.Matrix")
}


as.matrix.splits <- function(x, zero.print = 0L, one.print=1L, ...){
   m = length(x)
   labels = attr(x, "labels")
   n = length(labels)    
   res = matrix(zero.print, m, n)
   for(i in 1:m)res[i,x[[i]]]=one.print
   dimnames(res) = list(names(x), labels)
   res
}


#
# needs Matrix package, Matrix(t(as.matrix))
# 

# include weights
as.Matrix.splits <- function(x){
    l = length(x)
    j = unlist(x)
    i = rep(1:l, sapply(x, length))
    sparseMatrix(i,j) # ,x = rep(1L, length(i))
}


print.splits <- function (x, maxp = getOption("max.print"), 
    zero.print = ".", one.print="|", ...)
{
    x.orig <- x
    cx <- as.matrix(x, zero.print = zero.print, one.print=one.print)
    print(cx, quote = FALSE, right = TRUE, max = maxp)
    invisible(x.orig)
}


"[.splits" = function(x, i){
   result = unclass(x)[i]
   if(!is.null(attr(x, "weights"))) attr(result, "weights") = attr(x, "weights")[i] 
   if(!is.null(attr(x, "data"))) attr(result, "data") = attr(x, "data")[i,, drop=FALSE] 
   attr(result, "labels") = attr(x, "labels")
   class(result) = c("splits", "prop.part")
   result
}


c.splits <- function (..., recursive=FALSE) 
{
    x <- list(...)
    n <- length(x)
    match.names <- function(a, b) {
        if (any(!(a %in% b))) 
            stop("names do not match previous names")
    }
    if (n == 1) 
        return(x[[1]])

    labels <- attr(x[[1]], "labels")
    for (i in 2:n) {
        match.names(labels, attr(x[[i]], "labels"))
    }
    res = structure(NextMethod("c"), class=c("splits", "prop.part"))
    attr(res, "labels") = labels
    attr(res, "weight") = as.vector(sapply(x, attr, "weight"))
    res
}


# computes splits from phylo
as.splits.phylo <- function(x, ...){
    result = bip(x)
    edge.weights = numeric(max(x$edge))
    edge.weights[x$edge[,2]] = x$edge.length
    attr(result, "weights") = edge.weights
    attr(result, "labels") <- x$tip
    class(result) = c('splits', 'prop.part')
    result 
}


# computes splits from multiPhylo object (e.g. bootstrap, MCMC etc.)
as.splits.multiPhylo <- function(x, ...){
    if(class(x)=="multiPhylo")x = .uncompressTipLabel(x)
    lx = length(x)
    if(class(x)=="multiPhylo")class(x)='list'  # prop.part allows not yet multiPhylo
    firstTip = x[[1]]$tip[1]
    x = lapply(x, root, firstTip) # old trick  
    splits <- prop.part(x)
    class(splits)='list'
    weights = attr(splits, 'number')
    lab = attr(splits,'labels')
    attr(splits,'labels') <- attr(splits, 'number') <- NULL
    l = length(lab)
    splitTips = vector('list', l)
    for(i in 1:l) splitTips[[i]] = i
    result = c(splitTips,splits)
    attr(result, "weights") = c(rep(lx, l), weights)
    attr(result, "labels") <- lab
    class(result) = c('splits', 'prop.part')
    result  
}


as.splits.prop.part <- function(x, ...){
    if(is.null(attr(x, "number"))) attr(x, "weights") = rep(1, length(x)) 
	else attr(x, "weights") = attr(x, "number")
    class(x) = c('splits', 'prop.part')	
    x
}


as.prop.part <- function (x, ...){
    if (class(x) == "prop.part") return(x)
    UseMethod("as.prop.part")
}


as.prop.part.splits <- function(x, ...){
    attr(x, "number") = attr(x, "weights")
    attr(x, "weights") = NULL
    class(x) = c('prop.part')	
    x
}


as.phylo.splits <- function (x, result = "phylo", ...) 
{
    result <- match.arg(result, c("phylo", "all"))
    labels = attr(x, "labels")
    weights = attr(x, "weights")
    nTips = length(labels)
    dm = as.matrix(compatible(x))
    rs = rowSums(dm)
    ind = which(rs == 0)
    if (any(rs > 0)) {
        tmp = which(rs > 0)
        candidates = tmp[order(rs[tmp])]
        for (i in candidates) {
            if (sum(dm[ind, i]) == 0) 
                ind = c(ind, i)
        }
    }
    splits = x[ind]
    weights = weights[ind]
    l = length(ind)
    res = matrix(0L, l, nTips)
    for (i in 1:l) res[i, splits[[i]]] = 1L
    dm2 = (crossprod(res * weights, 1 - res))
    dm2 = dm2 + t(dm2)
    dimnames(dm2) = list(labels, labels)
    tree <- di2multi(NJ(dm2), tol = 1e-08)
    tree <- reorderPruning(tree)
    if (result == "phylo") 
        return(tree)  
    tree = reroot(tree, Ancestors(tree, 1, "parent")) 
    list(tree = tree, index = tree$edge[, 2], split = as.splits(tree), rest = x[-ind])
}


# computes compatible splits
compatible <- function(obj){
    labels = attr(obj, "labels")
    if(!inherits(obj, "splits"))stop("obj needs to be of class splits")
    
    l = length(labels)
    n = length(obj)
    
    bp = matrix(0L, n, l)
    for(i in 1:n)bp[i,obj[[i]]] = 1L
    bp[bp[, 1] == 0L, ] = 1L - bp[bp[, 1] == 0L, ]
    k=1
    res = matrix(0L, n, n) 
            
    tmp1 = tcrossprod(bp) #sum(bp[i,]* bp[j,])
    tmp2 = tcrossprod(1L - bp) #sum((1L - bp[i,])*(1L - bp[j,]))
    tmp3 = tcrossprod(bp, 1L - bp) #sum(bp[i,]*(1L - bp[j,]))
    tmp4 = tcrossprod(1L - bp, bp) #sum((1L - bp[i,])*bp[j,]) 
    res[(tmp1 * tmp2 * tmp3 * tmp4)>0]=1L
    k = k+1
    
    res = res[lower.tri(res)]
    attr(res, "Size") <- n
    attr(res, "Diag") <- FALSE
    attr(res, "Upper") <- FALSE
    class(res) <- "dist"
    return(res)
    }

    
compatible2 <- function (obj1, obj2=NULL) 
{   
    if (!inherits(obj1, "splits")) 
        stop("obj needs to be of class splits")
    labels = attr(obj1, "labels")    
    l = length(labels)
    n = length(obj1)
# use Matrix instead
    bp1 = matrix(0L, n, l)
    for (i in 1:n) bp1[i, obj1[[i]]] = 1L
    bp1[bp1[, 1] == 0L, ] = 1L - bp1[bp1[, 1] == 0L, ] 
    if(!is.null(obj2)){
        m = length(obj2) 
        bp2 = matrix(0L, m, l)
        for (i in 1:m) bp2[i, obj2[[i]]] = 1L

        labels2 = attr(obj2, "labels")
        bp2 = bp2[, match(labels2, labels), drop=FALSE]

        bp2[bp2[, 1] == 0L, ] = 1L - bp2[bp2[, 1] == 0L, ]
    }
    else bp2 = bp1

    if(is.null(obj2)) res = matrix(0L, n, n)
    else res = matrix(0L, n, m)

    tmp1 = tcrossprod(bp1, bp2)
    tmp2 = tcrossprod(1L - bp1, 1L - bp2)
    tmp3 = tcrossprod(bp1, 1L - bp2)
    tmp4 = tcrossprod(1L - bp1, bp2)
    res[(tmp1 * tmp2 * tmp3 * tmp4) > 0] = 1L
    if(is.null(obj2)){
        res = res[lower.tri(res)]
        attr(res, "Size") <- n
        attr(res, "Diag") <- FALSE
        attr(res, "Upper") <- FALSE
        class(res) <- "dist"
    }
    return(res)
}


compatible3 <- function(x, y=NULL) 
{
    if (!inherits(x, "splits")) 
        stop("x needs to be of class splits")
    if(is.null(y)) y <- x
        
    if (!inherits(y, "splits")) 
        stop("y needs to be of class splits")
    xlabels = attr(x, "labels")
    ylabels = attr(y, "labels")
    if(identical(xlabels, ylabels)) labels = xlabels 
    else labels = intersect(xlabels, ylabels)
#    if(length(labels) maybe warning
    nx = length(x)
    ny = length(y)   
    bp1 = as.matrix(x)[,labels, drop=FALSE]
    bp2 = as.matrix(y)[,labels, drop=FALSE]
    rs1 = rowSums(bp1)
    rs2 = rowSums(bp2)
    res = matrix(0L, nx, ny)
    tmp1 = tcrossprod(bp1, bp2)
    res = matrix(0L, nx, ny)
    for(i in 1:nx){
        for(j in 1:ny){            
            if(tmp1[i, j]==rs1[i]) res[i,j] = 1
            if(tmp1[i, j]==rs2[j]) res[i,j] = 2
            if(tmp1[i, j]==rs1[i] & tmp1[i, j]==rs2[j])res[i,j] = 3
        }
    }      
    if(is.null(y)){
        res = res[lower.tri(res)]
        attr(res, "Size") <- length(x)
        attr(res, "Diag") <- FALSE
        attr(res, "Upper") <- FALSE
        class(res) <- "dist"
    }
    return(res)
}
    

#
# networx
#
addEdge <- function(network, desc, split){
    edge = network[[1]]$edge
    parent = edge[,1]
    child = edge[,2]

    l = length(attr(split, "labels"))
    if(split[[1]][1]==1) split[[1]] <- c(1:l)[-split[[1]]] # ==1
 
    index = network[[2]]
    ind = which(compatible2(split, desc[index]) == 1)
    if(is.null(ind)) return(network)
    add = TRUE

    while(add){
        tmp = ind
        for(i in ind){
             tmp2 = which(compatible2(desc[index][i], desc[index]) == 1)
             tmp = union(tmp, tmp2)
        }
        if(identical(ind, tmp)){add=FALSE}
        ind=tmp 
    }    

    oldNodes = unique(as.vector(edge[ind,]))
    newNodes = (max(parent)+1L) : (max(parent)+length(oldNodes))

    ind2 = index[-ind]
    edge2 = edge[-ind,, drop=FALSE] 

    for(i in 1:length(oldNodes)){
        ind3 = which(edge2[,1] == oldNodes[i])
#        ind3 = which(edge2 == oldNodes[i],arr=TRUE)[,1]
        for(j in ind3) {
             if(any( desc[[ ind2[j] ]] %in% split[[1]])){
                   edge2[j,edge2[j,]==oldNodes[i]] = newNodes[i]
             }
        } 
    } 
    edge[-ind,] = edge2

#alle Splits verdoppeln
    dSpl = edge[ind,]
    for(i in 1:length(oldNodes)) dSpl[dSpl==oldNodes[i]] = newNodes[i]
    edge = rbind(edge, dSpl, deparse.level = 0) # experimental: no labels
    network[[1]]$edge.length = c(network[[1]]$edge.length, network[[1]]$edge.length[ind])   
    index = c(index, index[ind])

#neu zu alt verbinden   
    edge = rbind(edge, cbind(oldNodes, newNodes), deparse.level = 0) #  rbind(edge, cbind(oldNodes, newNodes), deparse.level = 0)# experimental: no labels
    network[[1]]$edge.length = c(network[[1]]$edge.length, rep(attr(split,"weights"), length(oldNodes)))
    index = c(index, rep(max(index)+1, length(oldNodes)) )
    network[[2]] = index
    
    network[[1]]$edge = edge
    network[[1]]$Nnode = length(unique(edge[,1]))
    
    network[[3]]= c(desc, split)   #[length(sp)+1] = split
    network   
}


as.networx <- function (x, ...) 
{
    if (inherits(x, "networx")) 
        return(x)
    UseMethod("as.networx")
}


as.networx.splits <- function(x, ...){   
    tmp = as.phylo(x, "all")
    obj = tmp[1:2]
    splitsOld = tmp[[3]]
    splitsNew = tmp[[4]]  
    if(length(splitsNew)>0){  
        ord <- order(colSums(compatible2(c(splitsOld,splitsNew), splitsNew)))
        ord2 <- order(colSums(compatible2(splitsOld, splitsNew)))
        ord3 <- order(colSums(compatible2(splitsNew, splitsNew)))
        splitsNew = splitsNew[ord]

        for(i in 1:length(splitsNew)){
            tmp = addEdge(obj, splitsOld, splitsNew[i])
            splitsOld = tmp[[3]]
            obj = tmp[1:2]    
        } 
    }
    res = obj[[1]]    
    nTips = length(res$tip.label)
    res$Nnode = max(res$edge) - nTips # raus in addEdge
    res$split = obj[[2]]
    class(res) = c("networx", "phylo")
    res
}


consensusNet <- function(obj, prob=.3, ...){
    l = length(obj)
    spl = as.splits(obj)
    w = attr(spl, "weight")
    ind = (w/l) > prob 
    spl = spl[ind] 
    as.networx(spl)
}


reorder.networx <- function (x, order = "cladewise", ...) 
{
    order <- match.arg(order, c("cladewise", "pruningwise"))
    if (!is.null(attr(x, "order"))) 
        if (attr(x, "order") == order) 
            return(x)
    nb.node <- x$Nnode
    if (nb.node == 1) 
        return(x)
    nb.tip <- length(x$tip.label)
    nb.edge <- dim(x$edge)[1]
    neworder <- if (order == "cladewise") 
        .C("neworder_cladewise", as.integer(nb.tip), as.integer(x$edge[, 
            1]), as.integer(x$edge[, 2]), as.integer(nb.edge), 
            integer(nb.edge), PACKAGE = "ape")[[5]]
    else .C("neworder_pruningwise", as.integer(nb.tip), as.integer(nb.node), 
        as.integer(x$edge[, 1]), as.integer(x$edge[, 2]), as.integer(nb.edge), 
        integer(nb.edge), PACKAGE = "ape")[[6]]
    x$edge <- x$edge[neworder, ]
    if (!is.null(x$edge.length)) 
        x$edge.length <- x$edge.length[neworder]
    # 
    if (!is.null(x$split))x$split <- x$split[neworder]
    attr(x, "order") <- order
    x
}


coords <- function(obj, dim="3D"){
	if(is.null(attr(obj,"order")) || attr(obj, "order")=="pruningwise") obj <- reorder(obj)
    obj = reorder.networx(obj)

    l = length(obj$edge.length)
    ind1 = which(!duplicated(obj$split))

    n = max(obj$edge)
    adj = Matrix::spMatrix(n, n, i = obj$edge[,2], j = obj$edge[,1], x = rep(1, length(obj$edge.length)))
    g = graph.adjacency(adj, "undirected")

    g2 <- graph(t(obj$edge)-1, directed=FALSE)
    g2 <- set.edge.attribute(g, "weight", value=obj$edge.length)
    if(dim=="3D"){
        coord <- igraph::layout.kamada.kawai(g, dim=3)
        k = matrix(0, max(obj$split), 3)
        for(i in ind1){
            tmp = coord[obj$edge[i, 2],] - coord[obj$edge[i, 1],]
            k[obj$split[i], ] = kart2kugel(tmp[1], tmp[2], tmp[3])
        }
        k[obj$split[ind1],1] = obj$edge.length[ind1] 

        res = matrix(0, vcount(g), 3)
        for(i in 1:l){# unique(obj$split)
            j = obj$edge[i,1]
            m = obj$edge[i,2]
            p = obj$split[i]
            res[m,] = res[j,] + kugel2kart(k[p,1], k[p,2], k[p,3])     
        }            
    }
    else{
        coord <- layout.kamada.kawai(g, dim=2)
        k = matrix(0, max(obj$split), 2)
        for(i in ind1){
            tmp = coord[obj$edge[i, 2],] - coord[obj$edge[i, 1],]
            k[obj$split[i], ] = kart2kreis(tmp[1], tmp[2])
        }
        k[obj$split[ind1],1] = obj$edge.length[ind1] 
        res = matrix(0, vcount(g), 2)
        for(i in 1:l){# unique(obj$split)
            j = obj$edge[i,1]
            m = obj$edge[i,2]
            p = obj$split[i]
            res[m,] = res[j,] + kreis2kart(k[p,1], k[p,2])     
        }
    }  
    res  
}


kart2kugel <- function(x,y,z){
    r = sqrt(x*x+y*y+z*z)
    alpha = atan(sqrt(x*x+y*y) / z)
    if(z<0) alpha = alpha+pi
    beta = atan(y/x)
    if(x<0) beta = beta+pi 
    c(r,alpha,beta)
}

	
kart2kreis <- function(x,y){
    r = sqrt(x*x+y*y)
    alpha = atan(y/x) 
    if(x<0) alpha = alpha+pi
    c(r,alpha)
}	
	

kreis2kart <- function(r,alpha){
	c(r*cos(alpha), r*sin(alpha))
}


kugel2kart <- function(r,alpha,beta){
    x = r * sin(alpha) * cos(beta) 
    y = r * sin(alpha) * sin(beta) 
    z = r * cos(alpha)
    c(x,y,z)
}


edgeLabels <- function(xx,yy,zz=NULL, edge){
        XX <- (xx[edge[, 1]] + xx[edge[, 2]])/2
        YY <- (yy[edge[, 1]] + yy[edge[, 2]])/2
        if(!is.null(zz)){
	        ZZ <- (zz[edge[, 1]] + zz[edge[, 2]])/2
	        return(cbind(XX, YY, ZZ))
        }  
        cbind(XX, YY)  
}


plot.networx = function(x, type="3D", show.tip.label=TRUE, show.edge.label=FALSE, show.nodes=FALSE, tip.color = "blue", edge.color="grey", edge.width = 3, 
    font = 3, cex = 1, ...){
    type = match.arg(type, c("3D", "2D")) 
    n = max(x$edge)
    tip = rep(NA, n)
    tips = x$tip.label
    tip[1:length(tips)] = tips
    
    x = reorder(x)
    
    adj = spMatrix(n, n, i = x$edge[,2], j = x$edge[,1], x = rep(1, length(x$edge.length)))
    g = graph.adjacency(adj, "undirected")
    
    plot.success <- FALSE
    if (!plot.success & type=="3D") {
        if (!require(rgl)) {
            warning("package 'rgl' not found, can only plot in 2D")
        } else {       
             coord <- coords(x, dim="3D")
             plotRGL(coord, x, show.tip.label=show.tip.label, show.edge.label=show.edge.label, show.nodes=show.nodes, tip.color = tip.color,
             edge.color=edge.color, edge.width = edge.width, font = font, cex = cex)
             plot.success <- TRUE
        } 
    }
    if (!plot.success){
	    coord <- coords(x, dim="2D")
	    plot2D(coord, x, show.tip.label=show.tip.label, show.edge.label=show.edge.label, tip.color = tip.color, edge.color=edge.color, 
	    edge.width = edge.width, font = font, cex = cex, add=FALSE)
	    }    
}

    
plotRGL <- function(coords, net, show.tip.label=TRUE, show.edge.label=FALSE, show.nodes=FALSE, tip.color = "blue", edge.color="grey", 
    edge.width = 3, font = 3, cex = par("cex"), ...){
    edge = net$edge
  
    x = coords[,1]
    y = coords[,2]
    z = coords[,3]
     
    nTips = length(net$tip.label)

    for(i in 1:dim(edge)[1]){
        segments3d(x[edge[i,]],y[edge[i,]],z[edge[i,]], col=edge.color, lwd=edge.width)
    }
    radius=0
    if(show.nodes){
        radius = sqrt((max(x)-min(x))^2 + (max(y)-min(y))^2 + (max(z)-min(z))^2) / 200    
        spheres3d(x[1:nTips], y[1:nTips],z[1:nTips], radius=2*radius, color="cyan")
        spheres3d(x[-c(1:nTips)], y[-c(1:nTips)],z[-c(1:nTips)], radius=radius, color="magenta")
    }
    if(show.tip.label){
      rgl.texts(x[1:nTips]+2.05*radius,y[1:nTips],z[1:nTips],net$tip.label, color=tip.color, cex=cex, font=font)
    }
    if(show.edge.label){
	    ec = edgeLabels(x, y, z, edge)
	    rgl.texts(ec[,1], ec[,2], ec[,3], net$split, color=tip.color, cex=cex, font=font)     
    } 
}


plot2D <- function(coords, net, show.tip.label=TRUE, show.edge.label=FALSE, tip.color = "blue", edge.color="grey", edge.width = 3, 
    font = 3, cex = par("cex"), add=FALSE, ...){
   edge = net$edge
   label = net$tip.label
   xx = coords[,1]
   yy = coords[,2]
   nTips = length(label)

   cex=1
   
   xlim <- range(xx)
   ylim <- range(yy)
     
   if(show.tip.label){
       offset <- max(nchar(label)) * 0.018 * cex * diff(xlim)
       xlim = c(xlim[1]-offset, xlim[2]+offset)
       ylim = c(ylim[1]-0.03 * cex * diff(ylim), ylim[2]+0.03 * cex * diff(ylim))
   }
   if(!add){ 
       plot.new() 
       plot.window(xlim, ylim, asp=1)
   }
   cladogram.plot(edge, xx, yy, edge.color, edge.width, 1)
   if(show.tip.label){
        ind=match(1:nTips, edge[,2])
        pos = rep(4, nTips)
        XX <- xx[edge[ind, 1]] - xx[edge[ind, 2]]
        pos[XX>0] = 2
        YY <- yy[edge[ind, 1]] - yy[edge[ind, 2]]
        pos2 <- rep(3, nTips)
        pos2[YY>0] = 1
        pos[abs(YY)>abs(XX)] <- pos2[abs(YY)>abs(XX)] 	
        text(xx[1:nTips], yy[1:nTips], labels=label, pos=pos, col=tip.color, cex=cex, font=font)
    }
    if(show.edge.label){
	    ec = edgeLabels(xx,yy, edge=edge)
	    text(ec[,1], ec[,2], labels=net$split, col=tip.color, cex=cex, font=font)     
	    } 
}   
   
    
lento <- function (obj, xlim = NULL, ylim = NULL, main = "Lento plot", 
    sub = NULL, xlab = NULL, ylab = NULL, bipart=TRUE, trivial=FALSE, ...) 
{
    if (class(obj) == "phylo") 
        obj = as.splits(obj)
    if (class(obj) == "multiPhylo") 
        obj = as.splits(obj)    
    labels = attr(obj, "labels") 
    l = length(labels)
    if(!trivial){
        triv = sapply(obj, length)
        ind = logical(length(obj)) 
        ind[(triv >1) & (triv < (l-1))] = TRUE
        obj = obj[ind]
        }
    CM = compatible(obj)
    support = attr(obj, "weights")
    if (is.null(support)) 
        support = rep(1, length(obj))
    conflict = -as.matrix(CM) %*% support
    n = length(support)
    if (is.null(ylim)) {
        eps = (max(support) - min(conflict)) * 0.05
        ylim = c(min(conflict) - eps, max(support) + eps)
    }
    if (is.null(xlim)) {
        xlim = c(0, n + 1)
    }

    ord = order(support, decreasing = TRUE)
    support = support[ord]
    conflict = conflict[ord]
    plot.new()
    plot.window(xlim, ylim)
    title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...)
    segments(0:(n - 1), support, y1 = conflict, ...)
    segments(1:n, support, y1 = conflict, ...)
    segments(0:(n - 1), support, x1 = 1:n, ...)
    segments(0:(n - 1), conflict, x1 = 1:n, ...)
    abline(h = 0)
    axis(2, ...)
    aty = diff(ylim)/(l+1)
    at = min(ylim) + (1:l) * aty
    if(bipart){
        Y = rep(at, n)
        X = rep((1:n)-.5, each=l)
        Circles = matrix(1, l, n)
        for(i in 1:n) Circles[obj[[ord[i]]],i] = 19   
#    axis(4, labels=labels, at=at)
        text(x=n+.1,y=at, labels, pos=4, ...) 
        points(X,Y,pch = as.numeric(Circles), col = rgb(0,0,0,.5), ...)
        }
    invisible(cbind(support, conflict))
    }

    
write.splits = function (x, file = "", zero.print = ".", one.print = "|", print.labels = TRUE, ...) 
{
    labels = attr(x, "labels")
    x.orig <- x
    cx <- as.matrix(x, zero.print = zero.print, one.print = one.print)
    w = FALSE
    if (!is.null(attr(x, "names"))) {
        nam = TRUE
        vnames = format(attr(x, "names"))
    }
    nam = FALSE
    if (!is.null(attr(x, "weight"))) {
        w = TRUE
        weight = format(attr(x, "weight"))
    }
    d = FALSE
    if (!is.null(attr(x, "data"))) {
        d = TRUE
        data = attr(x, "data")
    }
    if(print.labels){for(i in 1:length(labels)) cat(labels[i], "\n", file = file, append = TRUE)}
    if (w) 
        cat("weight", "\t", file = file, append = TRUE)
    if (d) 
        cat(paste(colnames(data), "\t"), file = file, append = TRUE)
    cat("\n", file = file, append = TRUE) #"Matrix", 
    for (i in 1:length(x)) {
        if (nam) 
            cat(vnames[i], "\t", file = file, append = TRUE)
        if (d) 
            cat(paste(data[i, ], "\t"), file = file, append = TRUE)
        if (w) 
            cat(weight[i], "\t", file = file)
        cat("\n", paste(cx[i, ], collapse = ""),"\n",  file = file, append = TRUE)
    }
}
 
    
write.nexus.splits <- function (obj, file = "", weights=NULL) 
{
    if(is.null(weights))weight <- attr(obj, "weights")
    taxa.labels <- attr(obj, "labels")
    ntaxa = length(taxa.labels)
    nsplits = length(obj)
    
    if (is.null(weight)) 
        weight = numeric(nsplits) + 100
    cat("#NEXUS\n\n", file = file)
    cat("[Splits block for Spectronet]\n", file = file, append = TRUE)
    cat("[generated by phangorn:\n", file = file, append = TRUE)
    cat(format(citation("phangorn"), "text"), "]\n\n",
       file = file, append = TRUE)
    cat(paste("BEGIN TAXA;\n\tDIMENSIONS NTAX=", ntaxa, ";\n", 
        sep = ""), file = file, append = TRUE)
    cat("\tTAXLABELS", paste(taxa.labels, sep = " "), ";\nEND;\n\n", 
        file = file, append = TRUE)
    cat(paste("BEGIN ST_SPLITS;\n\tDIMENSIONS NSPLITS=", nsplits, 
        ";\n", sep = ""), file = file, append = TRUE)
    cat("\tFORMAT LABELS WEIGHTS;\n\tMATRIX\n", file = file, 
        append = TRUE)
    for (i in 1:nsplits) cat("\t\t", i, weight[i], paste(obj[[i]]), 
        ",\n", file = file, append = TRUE)
    cat("\t;\nEND;\n", file = file, append = TRUE)
}


#
# ancestral sequences ML
#
ll2 = function (dat1, tree, g = 1, eig, ...) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorderPruning(tree)
    q = length(tree$tip.label)
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    m = length(edge) + 1
    dat = vector(mode = "list", length = m)
    el <- tree$edge.length
    P <- getP(el, eig, g)
    nr <- as.integer(attr(dat1, "nr"))
    nc <- as.integer(attr(dat1, "nc"))
    node = as.integer(node - min(node))
    edge = as.integer(edge - 1)
    nTips = as.integer(length(tree$tip))
    mNodes = as.integer(max(node) + 1)
    contrast = attr(dat1, "contrast")
    nco = as.integer(dim(contrast)[1])
    res <- .Call("LogLik2", dat1[tree$tip.label], P, nr, nc, 
        node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")[[1]]
    res
}


ancestral.pml <- function (object, type=c("ml", "bayes")) 
{
    call <- match.call()
    type <- match.arg(type)
    pt <- match.arg(type, c("ml", "bayes"))   
    tree = object$tree 
    
    INV <- object$INV
    inv <- object$inv
    
    data = getCols(object$data, tree$tip) 
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorderPruning(tree)
    q = length(tree$tip.label)
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    m = length(edge) + 1  # max(edge)
    w = object$w
    g = object$g
    l = length(w)    
    nr <- attr(data, "nr")
    nc <- attr(data, "nc")
    dat = vector(mode = "list", length = m*l)
    result = vector(mode = "list", length = m)
    dim(dat) <- c(l,m)
    
    x = attributes(data)
    label = as.character(1:m)
    nam = tree$tip.label
    label[1:length(nam)] = nam
    x[["names"]] = label
  
    
    tmp = length(data)

    result = new2old.phyDat(data) 

    eig = object$eig

    bf = object$bf
    el <- tree$edge.length
    P <- getP(el, eig, g)
    nr <- as.integer(attr(data, "nr"))
    nc <- as.integer(attr(data, "nc"))
    node = as.integer(node - min(node))
    edge = as.integer(edge - 1)
    nTips = as.integer(length(tree$tip))
    mNodes = as.integer(max(node) + 1)
    contrast = attr(data, "contrast")
    nco = as.integer(dim(contrast)[1])
    for(i in 1:l)dat[i,(q + 1):m] <- .Call("LogLik2", data, P[i,], nr, nc, node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")

    parent <- tree$edge[, 1]
    child <- tree$edge[, 2]
    nTips = min(parent) - 1
   
    for(i in 1:l){     
        for (j in (m - 1):1) {
            if (child[j] > nTips){
                tmp2 = (dat[[i, parent[j]]]/(dat[[i,child[j]]] %*% P[[i,j]]))
                dat[[i, child[j]]] = (tmp2 %*% P[[i,j]]) * dat[[i, child[j]]]  
            }
        }
    }
    for (j in unique(parent)) {
        tmp <- matrix(0, nr, nc)
        if(inv>0) tmp = as.matrix(INV) * inv
        for(i in 1:l){  
            tmp = tmp + w[i] * dat[[i, j]]                                 
        }
        if (pt == "bayes") tmp = tmp * rep(bf, each=nr)
        tmp = tmp / rowSums(tmp)
        result[[j]] = tmp
    } 
    attributes(result) = x
    attr(result, "call") <- call
    result 
}


fast.tree  = function(tree, node){
   parent = c(node, Ancestors(tree, node))
   children = Descendants(tree, parent, 'children')
   l = sapply(children, length)
   edge = cbind(rep(parent, l), unlist(children))
   obj = list(edge=edge, Nnode=sum(l>0), tip.label=as.character(edge[is.na(match(edge[,2], edge[,1])),2]))
   class(obj) = 'phylo'
   obj
}

# schneller ???
fast.tree2  = function(tree, node){
   parent = c(node, Ancestors(tree, node))
   edge = tree$edge 
   ind = match(edge[,1], parent)
   edge=edge[which(!is.na(ind)),] 
   obj = list(edge=edge, Nnode=length(parent), tip.label=as.character(edge[is.na(match(edge[,2], edge[,1])),2]))
   class(obj) = 'phylo'
   obj
}

###########################################################################################


getClans = function (tree) 
{
	if (is.rooted(tree)) 
        tree = unroot(tree)
    bp = bip(tree)
    nTips = length(tree$tip)
    root = nTips + 1
    bp[root] = NULL
    X = matrix(0, length(bp) - nTips, nTips)
    k = 1
    nl = NULL
    if (!is.null(tree$node.label)) {
        nl = c(rep("-1", nTips), rep("-1", nTips), tree$node.label[-1], 
            tree$node.label[-1])
    }
    if(root<=length(bp)){
        for (i in root:length(bp)) {
           X[k, bp[[i]]] = 1
           k = k + 1
        }
    }
    res <- rbind(diag(nTips), 1 - diag(nTips), X, 1 - X)
    colnames(res) <- tree$tip
    if (!is.null(nl)) 
        rownames(res) = nl
    res
}


getSlices <- function(tree){
    nTips = length(tree$tip)
    clans = getClans(tree)
    m = dim(clans)[1]
    X = tcrossprod(clans)
    z = rowSums(clans)
    Z1 = matrix(z,m,m)
    Z2 = t(Z1)
    Z = matrix(0,m,m)
    Z[Z1<=Z2] = Z1[Z1<=Z2]
    Z[Z2<Z1] = Z2[Z2<Z1]

    diag(X)=0
    X[upper.tri(X)] = 0
    X[X==1] = 0
    X[X==Z] = 0 
    index = which(X>0,arr.ind=TRUE)
    l = dim(index)[1]
    nSlices = 2 * nTips^2 -10 * nTips + 12
    result = matrix(0, nSlices, nTips)
    strClan = do.call("paste", c(as.data.frame(clans), sep = ""))
    k=1
    for(i in 1:l){
        tmp1 = as.numeric((clans[index[i,1],] + clans[index[i,2],])==2)
        tmp = paste(tmp1,sep="",collapse="")
        if(is.na(match(tmp,strClan))){ 
            result[k,] = tmp1
            k=k+1
        }
    }
    if(k<nSlices) result = result[1:(k-1),] 
    colnames(result) <- tree$tip
    result   
}


getClips = function (tree, all = TRUE) 
{
    if (any(is.na(tree$edge.length))) 
        return(NULL)
    dm = cophenetic(tree)
    tips = tree$tip
    nTips = length(tips)
    res = numeric(nTips)
    result = NULL
    for (i in 1:nTips) {
        ord = order(dm[i, ])
        for (j in 2:(nTips - 1)) {
            ind = ord[1:j]
            
            if(i>min(ind) ) break()
            within = max(dm[ind, ind])
            between = min(dm[ind, -ind])
            if (within < between) {
                res = numeric(nTips)
                res[ind] = 1L
                result = rbind(result, res)
            }
        }
    }
    dimnames(result) = list(NULL, tips)  
    if (all) 
        return(result)
    ind = which.max(rowSums(result))
    result[ind, ]
}


shannon <- function (x, norm=TRUE) 
{
    p = x/sum(x)
    ShD = -sum(p * log10(p))
    if(norm){
        if (sum(x) == 1) return(0)
        ShD = ShD/log10(sum(x))
    }
    ShD
}


shannon2 <- function (x, norm=TRUE) 
{
    p = x/sum(x)
    ShD = -sum(p * log(p))
    if(norm){
        if (sum(x) == 1) return(0)
        ShD = ShD/log(sum(x))
    }
    ShD
}


getE = function (tree, x, clans = NULL, norm = TRUE) 
{
    if (is.rooted(tree)) 
        tree = unroot(tree)
    if (is.null(clans)) 
        clans = getClans(tree)
    labels = tree$tip.label
    x = x[labels]
    result = rep(NA, 12)
    names(result) = c("E* tree", "# natives", "# intruder", "# unknown", 
        "E* clan", "# intruder", "# unknown", "E* slice", "# intruder", 
        "# unknown", "bs 1", "bs 2")
    result[2] = sum(x == 1)
    result[3] = sum(x == 2)
    result[4] = sum(x == 3)
    if (result[2] == 0 || result[3] == 0) {
        if (result[2] > 1) 
            return(list(result, labels))
        else return(list(result, integer(0)))
    }
    LHG = E_Intruder(clans, x)
    d = dim(LHG)[1]
    if (d == 1) {
        result[1] = 0
        if (!is.null(tree$node.label)) 
            result[11] = as.numeric(rownames(LHG))
        return(list(result, labels[LHG == 0]))
    }
    intr = drop(LHG %*% as.numeric(x == 2))
    result[1] = shannon2(intr, norm = norm)
    o <- order(intr, decreasing = TRUE)
    if (!is.null(tree$node.label)) 
        result[11:12] = as.numeric(rownames(LHG)[o[c(1, 2)]])
    ind = which(LHG[o[1], ] == 1)
    result[6] = sum(x[-ind] == 2)
    result[7] = sum(x[-ind] == 3)


    if (length(x[-ind]) < 4) 
        return(list(result, NULL))
    result[5] = shannon2(intr[-o[1]], norm = norm)
    ind2 = c(which(LHG[o[1], ] == 1), which(LHG[o[2], ] == 1))
  
    spl = structure(list(which(colSums(LHG)==0)), labels=labels, weights=1)
    class(spl)="splits"

    if (d == 2) {
         return(list(result, spl))
    } 
    result[9] = sum(x[-ind2] == 2)
    result[10] = sum(x[-ind2] == 3)
    if (length(x[-ind2]) < 4){ 
         return(list(result, spl))
    } 
    result[8] = shannon2(intr[-c(o[1], o[2])], norm = norm)
    return(list(result, spl))
}


E_Intruder <- function (clans, x) 
{
    cp = drop(clans %*% as.numeric(x == 1))
    ci = drop(clans %*% as.numeric(x == 2))
    homo = which(cp == 0 & ci > 0)
    l = length(homo)
    if (l > 0) {
        HG = clans[homo, , drop = FALSE]
        lhg = rep(TRUE, l)
        rsh = rowSums(HG)
        Z = tcrossprod(HG)>0
        Z = Z * rsh
        zmax = apply(Z,2,max)
        lhg = !(zmax > rsh)  
        LHG = HG[lhg, , drop = FALSE]
        return(LHG)
    }
    return(NULL)
}


E_Intruder_2 <- function (clans, x, native=NULL) 
{     
    contr = attr(x, "contr")
    d = dim(contr)
    if(d[1]>d[2])contr[(d[2]+1):d[1],]=0
    cp = clans %*% contr[as.numeric(x),]
    
	homo = which(rowSums(cp > 0) == 1)
		
    l = length(homo)
    if (l > 0) {
        HG = clans[homo, , drop = FALSE]
        lhg = rep(TRUE, l)
        rsh = rowSums(HG)
        Z = tcrossprod(HG)>0
        Z = Z * rsh
        zmax = apply(Z,2,max)
        lhg = !(zmax > rsh)  
        LHG = HG[lhg, , drop = FALSE]
        return(LHG)
    }
    return(NULL)
}


getDiv <- function(tree, x, native=NULL){
    clans = getClans(tree)
    labels = tree$tip.label
    x = subset(x, labels)
    LHG = E_Intruder_2(clans, subset(x,,1))
    if(!is.null(native)){
	    ll = match(native, attr(x, "allLevels"))
	    ind = (as.numeric(x) %in% ll)
	    }    	    
	if(!is.null(native)){    
	    rs = rowSums(clans)
	    intr = clans %*% ind    
	    clans = clans[intr==0,]
	    d = which.max(rs[intr==0])
	    tree2 = drop.tip(tree, tip=labels[which(clans[d, ]==1)])
    } 
    else tree2=NULL
	list(c(shannon(rowSums(LHG)),      
    summary(factor(attr(x, "allLevels"))[as.numeric(subset(x,,1))]), parsimony(tree, x)), tree2 )     
}


getDiversity <- function (tree, x, norm = TRUE, var.names = NULL, labels="new") 
{
    k = 1
    if(class(tree) == "multiPhylo") 
        k = length(tree)
    l = attr(x, "nr")
    tmp = matrix(0, k * l, 12)

    tnam = 1
    if (class(tree) == "multiPhylo") {
        tnam = names(tree)
        if (is.null(tnam)) 
            tnam = 1:length(tree)
    }
    if(is.null(var.names)) var.names = 1:l
    PM = data.frame("t1", "a", stringsAsFactors = FALSE)
    colnames(PM) = c("Tree", "Var")
    PM = PM[FALSE,] 
    PM[1 :(k*l), ] = NA 
    perfect = names(x)
    L = vector("list",k*l)
    m = 1
    o = 1
    ok= 0
    for (i in 1:k) {
        if (class(tree) == "multiPhylo") 
            tmptree = tree[[i]]
        else tmptree = tree
        if (is.rooted(tmptree)) 
            tmptree = unroot(tmptree)
        clans = getClans(tmptree) 
        for (j in 1:l) {          
            TMP = getE(tmptree, getRows(x, j), clans, norm = norm)
            tmp[m, ] = TMP[[1]]
            L[[m]] = TMP[[2]] # if class =splits else NULL
            PM[m, 1] = tnam[i]
            PM[m, 2] = var.names[j]
            m = m + 1   
        }
    }

    tnam = rep(tnam, each = l)
    dnam = var.names
    dnam = rep(dnam, k)
    pscore = as.numeric(sankoff(tree, x, site = "site"))
    res = data.frame(tnam, dnam, tmp, pscore)
    if(labels=="old")names(res) = c("tree", "variable", "E tree", "# natives", 
        "# intruder", "# unknown", "E clan", "# intruder", "# unknown", 
        "E slice", "# intruder", "# unknown", "bs 1", "bs 2", "p-score")
    else{
        names(res) = c("tree", "variable", "E clan", "# natives", 
            "# intruder", "# unknown", "E slice", "# intruder", "# unknown", 
            "E melange", "# intruder", "# unknown", "bs 1", "bs 2", "p-score")    
        warning("The variable names have changed")       
    }    
    attr(res, "Perfect") = L
    class(res) = c("clanistics", "data.frame")
    res
}


summary.clanistics <- function(object, ...){
    res <- matrix(FALSE, nrow(object), 5)
    res[,1] = object[,4]>0 & object[,"p-score"]==0 # "natives"
    res[,2] = object[,5]>0 & object[,"p-score"]==0 # "intruder"
    res[,3] = object[,"p-score"]==1
    res[,4] = ( (object[,"p-score"]==2) & (object[,7]==0) & (!is.na(object[,7])) ) | 
              ( (object[,"p-score"]==2) & (object[,4]==2) & (is.na(object[,7])) )  
    res[,5] = object[,"p-score"]>=2 & (object[,7]>0) & (!is.na(object[,7]))
    res[] = as.numeric(res)
    tmp = data.frame(factor(object[,"variable"]), res)	
    colnames(tmp) = c("Variable", "Natives_only", "Intruder_only", "Clan", "Slice", "Melange")
#        colnames(res) = c("Natives only", "Intruder only", "Clan", "Melange")
    class(tmp) <- c("summary.clanistics", "data.frame")
    tmp
    }
	

print.summary.clanistics <- function(x, ...){
    print(aggregate(x[,-1], list(Variable=x[,1]), sum), ...)
}


compareSplits <- function(res, nam1, nam2){
    wide <- reshape(res[, c("tree", "E tree", "variable")], v.names="E tree", idvar="tree", timevar="variable", direction="wide")
    wideI <- reshape(res[, c("tree", "# natives", "variable")], v.names="# natives", idvar="tree", timevar="variable", direction="wide")
    for(i in 2:dim(wide)[2])colnames(wide)[i] = strsplit(colnames(wide)[i],"E tree.")[[1]][2]
    for(i in 2:dim(wide)[2])colnames(wideI)[i] = strsplit(colnames(wideI)[i],"# natives.")[[1]][2]
	ntrees = wide[,1]
	splits = attr(res, "Perfect")
	dat = attr(attr(res, "Perfect"), "data")
    res = matrix(NA, length(ntrees), length(nam1)*length(nam2))
    for(m in 1:length(trees)){
        k=1
        trnam=ntrees[m]
        for(i in nam1){
            E1 = wide[m, i]
            for(j in nam2){
                E2 = wide[m, j] 
                if(!is.na(E1) & !is.na(E2)){
                    if(E1 == E2){ # if(E1 == 0 & E2 == 0){
	                if( (wideI[m, i] >0) & (wideI[m, j]) >0){
                            ind1 = which(dat[,1]==trnam & dat[,2]==i)
                            sp1 = splits[[ind1]]                     
                            ind2 = which(dat[,1]==trnam & dat[,2]==j) 
                            sp2 = splits[[ind2]]
                            if(length(ind1)>0 & length(ind2)>0 )res[m, k] = drop(compatible3(sp1, sp2))
                        }
                    }
                }
            k=k+1 
            }
        }
    }    
    res
}


diversity <- function(tree, X){  
# from kknn
    contr.dummy <- function (n, contrasts = TRUE) 
    {
        if (length(n) <= 1) {
            if (is.numeric(n) && length(n) == 1 && n > 1) 
                levels <- 1:n
            else stop("contrasts are not defined for 0 degrees of freedom")
        }
        else levels <- n
        lenglev <- length(levels)
       cont <- array(0, c(lenglev, lenglev), list(levels, levels))
       cont[col(cont) == row(cont)] <- 1
       cont
    }


    l = dim(X)[2]
    m <- ifelse(class(tree)=="multiPhylo", length(tree), 1)
    
    contr = as.list(rep("contr.dummy", l))
    names(contr) = names(X)
    tmp = model.matrix(~.-1, X, contrast=contr)
    tmp1 <- phyDat.default(tmp, levels=c(1,0), compress = FALSE)
    attr(tmp1, "varnames")  = colnames(tmp)
    fd = sankoff(tree,tmp1,site = "site") 
    fd = matrix(fd, ncol=m) 

    if(m>1){
         if(is.null(names(tree))) tnames <- paste("tree", 1:m, sep=".")
         else tnames <- names(tree)  
    }
    else tnames = "tree"
    dimnames(fd) = list(colnames(tmp), tnames)
    res = stack(data.frame(fd))

    if(m>1)nt = rep(sapply(tree, function(x)length(x$tip)), each=dim(fd)[1])    
    else nt = rep(length(tree$tip), each=dim(fd)[1]) 
    if(m>1)res2 = as.vector(sapply(tree, function(x,y)colSums(y[x$tip,,drop=FALSE]) , y=tmp))
    else res2 = colSums(tmp[tree$tip,,drop=FALSE])
    result <- data.frame(tree = res[,2], variable=rep(colnames(tmp),m), pscore=res[,1], ntips=nt, natives=res2)
    result
}
 








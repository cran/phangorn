
# logarithic version of this!! l*(log(2) - log(theta))+
# in C programmieren; nonlog raus
coalBranch = function(coaltimes, tau, k, theta=1){
  tt = diff.default(sort(c(tau[1], coaltimes, tau[2])))
  #  if(length(tt)==1)return(1)
  #  if(k==1)return(1)
  if(length(tt)==1)return(0)
  if(k==1)return(0)
  res = 0 
  #  res=1
  l = length(t)  
  n=k
  if(l>0){
    for(i in 1:l){
      #      res = res * 2/theta * exp(-n*(n-1)*tt[i] / theta)
      res = res + log(2) - log(theta) + (-n*(n-1)*tt[i] / theta)
      n = n - 1L
    }
    # res = res + l*log(2/theta)    
  }
  #  if(n>1L) res = res * exp((-n*(n-1) / theta) * tt[l+1])
  if(n>1L) res = res + ((-n*(n-1) / theta) * tt[l+1])   
  res
}


MSC <- function(gTree, sTree, X, theta=1){
  
  snh <- nodeHeight(sTree) 
  l=length(sTree$tip)  
  tmp = matrix(1:l,l,1,dimnames=list(sTree$tip, NULL))
  Y=phyDat(tmp, "USER", 1:l) 
  sst = ancstat(sTree, Y)
  
  edge = sTree$edge
  desc = Descendants(sTree, 1:max(sTree$edge), type="all")
  
  # start a loop over the trees   
  # should be a list
  
  gnh <- nodeHeight(gTree)  # scaling on these  
  gst <- ancstat(gTree, X)    
  
  br = comp(gst, sst) 
  ind1 = logical(length(br))
  ind1[ (length(gTree$tip)+1): length(br) ] = TRUE
  
  # loop for trees, triplets (edge length optimisation) & quartets (NNI moves)  
  
  l = nrow(edge)
  ct = vector("list", l+1L)
  root=getRoot(sTree)
  
  for(i in 1:l){
    ind2= br %in% desc[[edge[i,2]]] # vielleicht ausserhalb der Schleife
    ind3=(gnh>=snh[edge[i,2]] & gnh<snh[edge[i,1]])
    ind = which(ind1 & ind2 & ind3)
    ct[[edge[i,2]]] = gnh[ind]    
    # nr of coalescents
  }
  ct[[root]] = gnh[gnh>max(snh)]   
  
  nrc = sapply(ct, length)
  nsp = length(sTree$tip)
  
  
  k = integer(l+1)
  for(i in 1:nsp)k[i] = sum(br[!ind1]==i)
  
  for(i in 1:l){
    ei = edge[i,1]
    k[ei] = k[ei]+k[edge[i,2]] - nrc[edge[i,2]]  
  }
  res = numeric(l+1L)
  
  nhm = rep(Inf, length(snh))
  for(i in 1:l) nhm[edge[i,2]] = snh[edge[i,1]]
  for(i in 1:(l+1L))res[i] = coalBranch(ct[[i]], c(snh[i], nhm[i]), k[i], theta)  
  sum(res)
  
}

# lapply(trees, MSC, sp, X)

# gives index of edges which are a certain height
indEdges <- function(tree, x){
  h <- nodeHeight(tree)
  ind <- which( h[tree$edge[,1]]>x & h[tree$edge[,2]]<x )
  tree$edge[ind,2]
}



ancstat = function(phy, x){
  contrast= attr(x, "contrast")
  storage.mode(contrast) = "integer"
  phy=reorder(phy, "postorder")
  res=matrix(0L, max(phy$edge), ncol(contrast))
  colnames(res) = attr(x,"levels")
  nTips=length(phy$tip.label)
  pa=phy$edge[,1]
  ch=phy$edge[,2]
  res[1:nTips, ] = contrast[as.numeric(x)[match(phy$tip, names(x))],, drop=FALSE]
  for(i in 1:length(pa)){
    res[pa[i],] = res[pa[i],] | res[ch[i],]    
  }
  res
}

plotA = function (tree, data, i = 1, col = NULL, ...) 
{
  if(class(data)=="phyDAt"){
    y = subset(data, , i)
    nc = attr(data, "nc")
    y = matrix(unlist(y[]), ncol = nc, byrow = TRUE)
    l = dim(y)[1]
    dat = matrix(0, l, nc)
    for (i in 1:l) dat[i, ] = y[[i]]
    levels = attr(data, "levels")
  }
  else{
    dat = data
    y = data
    nc = ncol(y)
    levels=colnames(data) 
  } 
  
  args <- list(...)
  CEX <- if ("cex" %in% names(args)) 
    args$cex
  else par("cex")
  xrad <- CEX * diff(par("usr")[1:2])/50
  
  plot(tree, label.offset = 1.1 * xrad, plot = FALSE, ...)
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  XX <- lastPP$xx
  YY <- lastPP$yy
  xrad <- CEX * diff(lastPP$x.lim * 1.1)/50
  par(new = TRUE)
  plot(tree, label.offset = 1.1 * xrad, plot = TRUE, ...)
  if (is.null(col)) 
    col = rainbow(nc)
  if (length(col) != nc) 
    warning("Length of color vector differs from number of levels!")
  BOTHlabels(pie = y, XX = XX, YY = YY, adj = c(0.5, 0.5), 
             frame = "rect", pch = NULL, sel = 1:length(XX), thermo = NULL, 
             piecol = col, col = "black", bg = "lightblue", horiz = FALSE, 
             width = NULL, height = NULL)
  legend("bottomright", levels, text.col = col)
}


# old version
#STAR <- function(trees, speciesTree, X){
#  trees = unclass(trees)
#  States = lapply(trees, ancstat, X)
#  NH = lapply(trees, nodeHeight)
#  SST = ancstat(speciesTree, X)
#  Y=matrix(Inf, length(NH), nrow(SST)) 
#  for(i in 1:length(NH)){
#    ind = comp(States[[i]],SST) 
#    for(j in 1:length(ind))Y[i, ind[j]] = min(Y[i, ind[j]], NH[[i]][j])
#  }
#  STH = apply(Y, 2, min)
#  speciesTree$edge.length = STH[speciesTree$edge[,1]] - STH[speciesTree$edge[,2]]
#  speciesTree
#}

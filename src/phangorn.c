/* 
 * phangorn.c
 *
 * (c) 2008-2014  Klaus Schliep (klaus.schliep@gmail.com)
 * 
 * 
 * This code may be distributed under the GNU GPL
 *
 */

# define USE_RINTERNALS

#include <Rmath.h>
#include <math.h>
#include <R.h> 
#include <R_ext/Lapack.h>
#include <Rinternals.h>



// off-diagonal
#define DINDEX(i, j) n*(i - 1) - i * (i - 1)/2 + j - i - 1
// with diagonal (+i), R index (+1)
#define DINDEX2(i, j) n*(i - 1) - i * (i - 1)/2 + j - 1

// index likelihood pml
// need to define nr, nc, nTips, nNodes k
#define LINDEX(i) (i-nTips) * (nr*nc) //+ k * nTips * (nr * nc)
#define LINDEX2(i, k) (i-nTips) * (nr*nc) + k * nTips * (nr * nc)
#define LINDEX3(i, k) (i-*nTips-1L) * (*nr* *nc) + k * *nTips * (*nr * *nc)

// index sankoff
#define SINDEX(i) i * (nr*nc) 

/* from coalescentMCMC
void get_single_index_integer(int *x, int *val, int *index)
{
    int i = 0, v = *val;
	while (x[i] != v) i++;
	*index = i + 1;
}

void get_two_index_integer(int *x, int *val, int *index)
{
	int i1 = 0, i2, v = *val;
	while (x[i1] != v) i1++;
	i2 = i1 + 1;
	while (x[i2] != v) i2++;
	index[0] = i1 + 1;
	index[1] = i2 + 1;
}
*/


void nodeH(int *edge, int *node, double *el, int *l,  double *res){
    int ei, i;
    for (i=*l-1L; i>=0; i--) {
        ei = edge[i] - 1L;
        res[ei] = res[node[i]-1L] + el[ei];
    }
}


SEXP rowMax(SEXP sdat, SEXP sn, SEXP sk){
    int i, h, n=INTEGER(sn)[0], k=INTEGER(sk)[0];  
    double x, *res, *dat;
    SEXP result;
    PROTECT(result = allocVector(REALSXP, n));
    res = REAL(result);
    PROTECT(sdat = coerceVector(sdat, REALSXP));
    dat = REAL(sdat);
    for(i = 0; i < n; i++){
        x = dat[i];
        for(h = 1; h< k; h++) {if(dat[i + h*n] > x) x=dat[i + h*n];}
        res[i] = x;               
        }
    UNPROTECT(2);
    return(result);        
}
    

static R_INLINE void getP00(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double tmp, res;
    for(i = 0; i < m; i++){
        tmp = exp(eva[i] * w * el);
        for(j=0; j<m; j++) evi[i + j*m] *= tmp;   
    }    
    for(i = 0; i < m; i++){    
        for(j = 0; j < m; j++){
            res = 0.0;    
            for(h = 0; h < m; h++) res += ev[i + h*m] * evi[h + j*m];
            result[i+j*m] = res;    
        }
    }
}


 
static R_INLINE void getPP(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double tmp[m];
    for(i = 0; i < m; i++) tmp[i] = exp(eva[i] * w * el);
    for(i = 0; i < m; i++){    
        for(j = 0; j < m; j++){
            result[i+j*m] = 0;
            for(h = 0; h < m; h++) result[i+j*m] += ev[i + h*m] * tmp[h] * evi[h + j*m];                
        }
    }
}


void getdP(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double tmp[m], res;
    for(i = 0; i < m; i++) tmp[i] = (eva[i] * w  * el) * exp(eva[i] * w * el);
    for(i = 0; i < m; i++){    
        for(j = 0; j < m; j++){
            res = 0.0;    
            for(h = 0; h < m; h++)    res += ev[i + h*m] * tmp[h] * evi[h + j*m];
            result[i+j*m] = res;    
        }
    }
}


void getdP2(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double tmp[m], res;
    for(i = 0; i < m; i++) tmp[i] = (eva[i] * w) * exp(eva[i] * w * el);
    for(i = 0; i < m; i++){    
        for(j = 0; j < m; j++){
            res = 0.0;    
            for(h = 0; h < m; h++)    res += ev[i + h*m] * tmp[h] * evi[h + j*m];
            result[i+j*m] = res;    
        }
    }
}


void getd2P(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double tmp[m], res;
    for(i = 0; i < m; i++) tmp[i] = (eva[i] * w * el) * (eva[i] * w * el) * exp(eva[i] * w * el);
    for(i = 0; i < m; i++){    
        for(j = 0; j < m; j++){
            res = 0.0;    
            for(h = 0; h < m; h++)    res += ev[i + h*m] * tmp[h] * evi[h + j*m];
            result[i+j*m] = res;    
        }
    }
}


void getd2P2(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double tmp[m], res;
    for(i = 0; i < m; i++) tmp[i] = (eva[i] * w) * (eva[i] * w) * exp(eva[i] * w * el);
    for(i = 0; i < m; i++){    
        for(j = 0; j < m; j++){
            res = 0.0;    
            for(h = 0; h < m; h++)    res += ev[i + h*m] * tmp[h] * evi[h + j*m];
            result[i+j*m] = res;    
        }
    }
}




SEXP getPM2(SEXP eig, SEXP nc, SEXP el, SEXP w){
    R_len_t i, j, nel, nw;
    int m=INTEGER(nc)[0], l=0;
    double *ws=REAL(w);
    double *edgelen=REAL(el);
    double *eva, *eve, *evei;
    SEXP P, RESULT;
    nel = length(el);
    nw = length(w);
    if(!isNewList(eig)) error("'eig' must be a list");    
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));
    PROTECT(RESULT = allocVector(VECSXP, nel*nw));    
    for(j=0; j<nel; j++){ 
        for(i=0; i<nw; i++){
            PROTECT(P = allocMatrix(REALSXP, m, m));
            getPP(eva, eve, evei, m, edgelen[j], ws[i], REAL(P));
            SET_VECTOR_ELT(RESULT, l, P);
            UNPROTECT(1);
            l++;
        }
    }
    UNPROTECT(1);//RESULT
    return(RESULT);
} 


SEXP getdPM(SEXP eig, SEXP nc, SEXP el, SEXP w){
    R_len_t i, j, nel, nw;
    int m=INTEGER(nc)[0], l=0;
    double *ws=REAL(w);
    double *edgelen=REAL(el);
    double *eva, *eve, *evei;
    SEXP P, RESULT;
    nel = length(el);
    nw = length(w);
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));
    PROTECT(RESULT = allocVector(VECSXP, nel*nw));    
    double *p;
    if(!isNewList(eig)) error("'dlist' must be a list");    
    for(j=0; j<nel; j++){
        for(i=0; i<nw; i++){
            PROTECT(P = allocMatrix(REALSXP, m, m));
            p = REAL(P);
            getdP(eva, eve, evei, m, edgelen[j], ws[i], p);
            SET_VECTOR_ELT(RESULT, l, P);
            UNPROTECT(1); //P
            l++;
        }
    }
    UNPROTECT(1);//RESULT
    return(RESULT);
} 


SEXP getdPM2(SEXP eig, SEXP nc, SEXP el, SEXP w){
    R_len_t i, j, nel, nw;
    int m=INTEGER(nc)[0], l=0;
    double *ws=REAL(w);
    double *edgelen=REAL(el);
    double *eva, *eve, *evei;
    SEXP P, RESULT;
    nel = length(el);
    nw = length(w);
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));
    PROTECT(RESULT = allocVector(VECSXP, nel*nw));    
    double *p;
    if(!isNewList(eig)) error("'dlist' must be a list");    
    for(j=0; j<nel; j++){
        for(i=0; i<nw; i++){
            PROTECT(P = allocMatrix(REALSXP, m, m));
            p = REAL(P);
            getdP2(eva, eve, evei, m, edgelen[j], ws[i], p);
            SET_VECTOR_ELT(RESULT, l, P);
            UNPROTECT(1); //P
            l++;
        }
    }
    UNPROTECT(1); //RESULT
    return(RESULT);
} 


SEXP getd2PM(SEXP eig, SEXP nc, SEXP el, SEXP w){
    R_len_t i, j, nel, nw;
    int m=INTEGER(nc)[0], l=0;
    double *ws=REAL(w);
    double *edgelen=REAL(el);
    double *eva, *eve, *evei;
    SEXP P, RESULT;
    nel = length(el);
    nw = length(w);
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));
    PROTECT(RESULT = allocVector(VECSXP, nel*nw));    
    double *p;
    if(!isNewList(eig)) error("'dlist' must be a list");    
    for(j=0; j<nel; j++){
        for(i=0; i<nw; i++){
            PROTECT(P = allocMatrix(REALSXP, m, m));
            p = REAL(P);
            getd2P(eva, eve, evei, m, edgelen[j], ws[i], p);
            SET_VECTOR_ELT(RESULT, l, P);
            UNPROTECT(1); //P
            l++;
        }
    }
    UNPROTECT(1); //RESULT
    return(RESULT);
} 


SEXP getd2PM2(SEXP eig, SEXP nc, SEXP el, SEXP w){
    R_len_t i, j, nel, nw;
    int m=INTEGER(nc)[0], l=0;
    double *ws=REAL(w);
    double *edgelen=REAL(el);
    double *eva, *eve, *evei;
    SEXP P, RESULT;
    nel = length(el);
    nw = length(w);
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));
    PROTECT(RESULT = allocVector(VECSXP, nel*nw));    
    double *p;
    if(!isNewList(eig)) error("'dlist' must be a list");    
    for(j=0; j<nel; j++){
        for(i=0; i<nw; i++){
            PROTECT(P = allocMatrix(REALSXP, m, m));
            p = REAL(P);
            getd2P2(eva, eve, evei, m, edgelen[j], ws[i], p);
            SET_VECTOR_ELT(RESULT, l, P);
            UNPROTECT(1); //P
            l++;
        }
    }
    UNPROTECT(1); //RESULT
    return(RESULT);
} 

    
static R_INLINE void emult(double *x, double *y, int n){
    for(int i=0; i<n; i++) x[i]*=y[i];
}



void tabulate(int *x, int *n, int *nbin, int *ans){
    int i, tmp;
    for (i=0; i < *nbin; i++) ans[i]=0L; 
    for (i=0; i < *n; i++) {
        tmp = x[i];
        if( (tmp>0) & (tmp<(*nbin+1L)) )   
        ans[tmp-1L] ++;
    }
}


void reorder(int *from, int *to, int *n, int *sumNode,  int *neworder, int *root){ 
    int i, j, sum=0, k, Nnode, ind, *ord, *csum, *tips, *stack, z=0;  // l, 
    double *parent;
    int m=sumNode[0];
    parent = (double *) R_alloc((*n), sizeof(double));
    tips = (int *) R_alloc(m, sizeof(int));
    ord = (int *) R_alloc((*n), sizeof(int));
    csum = (int *) R_alloc( (m+1), sizeof(int));
    stack = (int *) R_alloc(m, sizeof(int));
    for(j=0;j<(*n);j++) parent[j] = (double)from[j];
   
    for(j=0;j<(*n);j++) ord[j] = j;
    for(j=0;j<m;j++) tips[j] = 0;
        
    rsort_with_index(parent, ord, *n);
    tabulate(from, n, sumNode, tips);
    csum[0]=0;
    for(i=0;i<(*sumNode);i++){
        sum+=tips[i];                 
        csum[i+1] = sum;                        
    }      
    k = (*n)-1;
    Nnode = 0;
    stack[0] = *root;
    
    while(z > -1){
        j=stack[z];          
        if(tips[j]>0){   
            for(i=csum[j];i<csum[j+1];i++){
                ind = ord[i];                     
                neworder[k] = ind + 1;        
                stack[z] = to[ind]-1;
                k -=1;
                z++;                            
            }         
            Nnode += 1; 
            }
        z--;       
    }                
    root[0]=Nnode;     
}



SEXP AllChildren(SEXP children, SEXP parent, SEXP M){
    int i, j, k, l=0L, m=INTEGER(M)[0], *tab, p;   
    R_len_t n=length(parent); 
    SEXP RESULT, TMP;
    tab = (int*)R_alloc(m, sizeof(int));
    for(i=0; i<m; i++)tab[i]=0L;
    
    j=0;
    p = INTEGER(parent)[0];
    for(i=0; i<n; i++){
        if(INTEGER(parent)[i]!=p){
            p = INTEGER(parent)[i]; 
            j=j+1;
        } 
        tab[j] += 1L;
    }
    PROTECT(RESULT = allocVector(VECSXP, m));

    i=0L;    
    while(l<n){    
        k=tab[i];        
        PROTECT(TMP = allocVector(INTSXP, k));  
        p = INTEGER(parent)[l]-1;
        for(j=0; j<k; j++){
            INTEGER(TMP)[j] = INTEGER(children)[l];
            l++;
        } 
        SET_VECTOR_ELT(RESULT, p, TMP);
        UNPROTECT(1);
        i++;
    }
    UNPROTECT(1);
    return(RESULT);
}



/*
library(phangorn)
tree =  rtree(10)

allDesc = function(x, node){
  x = reorder(x, "postorder")
  parent = x$edge[, 1]
  children = x$edge[, 2]
  .Call("AllDesc", as.integer(children), as.integer(parent), as.integer(max(parent)), as.integer(node)) 
}

allDesc(tree, 14)
 */
 
SEXP AllDesc(SEXP child, SEXP parent, SEXP M, SEXP NODE){
    int i, m=INTEGER(M)[0]+1, *tab, *res, p, node=INTEGER(NODE)[0];   
    R_len_t n=length(parent); 
    SEXP RESULT;
    tab = (int*)R_alloc( m, sizeof(int));
    for(i=0; i<m; i++)tab[i]=0L;
    tab[node] = 1L;
    PROTECT(RESULT = allocVector(INTSXP, n));
        res = INTEGER(RESULT);
    for(i=0; i<n; i++)res[i]=0L;

   for(i=n-1L; i>=0L; i--){
        p = INTEGER(parent)[i];
        if(tab[p]==1L){
            res[i] = 1L;
            tab[INTEGER(child)[i]] = 1L; 
        } 
    }
    
    UNPROTECT(1);
    return(RESULT);
}


// copy from ape 3.0-5
void neworder_cladewise(int *n, int *edge1, int *edge2,
			int *N, int *neworder)
/* n: nb of tips, N: nb of edges */
{
    int i, j, k, node, *done, dn, *node_back, eb;
    done = &dn;
    node_back = &eb;

    /* done: indicates whether an edge has been collected
       node_back: the series of node from the root to `node'
       node: the current node */

    done = (int*)R_alloc(*N, sizeof(int));
    node_back = (int*)R_alloc(*N + 2 - *n, sizeof(int));
    memset(done, 0, *N * sizeof(int));

    j = k = 0;
    node = *n + 1;
    while (j < *N) {
        for (i = 0; i < *N; i++) {
  	    if (done[i] || edge1[i] != node) continue;
	    neworder[j] = i + 1;
	    j++;
	    done[i] = 1;
	    if (edge2[i] > *n) {
	        node_back[k] = node;
		k++;
		node = edge2[i];
		/* if found a new node, reset the loop */
		i = -1;
	    }
	}
	/* if arrived at the end of `edge', go down one node */
	k--;
	node = node_back[k];
    }
}



// combine two sorted vectors
void crsort(double *x, double *y, int *a, int *b, double *res){
   double xi, yi;
   int i, j, k;    
   i=0;
   j=0;
   k=0;
   xi=x[0];
   yi=y[0];  
   while(k<((*a)+(*b))){
      if(i<(*a)){
          if( (xi<yi) | (j==((*b))) ){  //-1L
              res[k]=xi;      
              i++;
              if(i<(*a))xi=x[i];   
              k++;     
          }
          else{
              j++;
              res[k]=yi;
              if(j<(*b))yi=y[j];  
              k++;
          }
        }
        else{
              j++;
              res[k]=yi;
              if(j<(*b))yi=y[j];  
              k++;
          }
    }  
}    


void cisort(int *x, int *y, int *a, int *b, int *res){
   int xi, yi;
   int i, j, k;    
   i=0;
   j=0;
   k=0;
   xi=x[0];
   yi=y[0];  
   while(k<((*a)+(*b))){
      if(i<(*a)){
          if( (xi<yi) | (j==((*b))) ){  //-1L
              res[k]=xi;      
              i++;
              if(i<(*a))xi=x[i];   
              k++;     
          }
          else{
              j++;
              res[k]=yi;
              if(j<(*b))yi=y[j];  
              k++;
          }
        }
        else{
              j++;
              res[k]=yi;
              if(j<(*b))yi=y[j];  
              k++;
          }
    }
}    

void cisort2(int *x, int *y, int a, int b, int *res){
   int xi, yi;
   int i, j, k;    
   i=0;
   j=0;
   k=0;
   xi=x[0];
   yi=y[0];  
   while(k<((a)+(b))){
      if(i<(a)){
          if( (xi<yi) | (j==b) ){  //-1L
              res[k]=xi;      
              i++;
              if(i<(a))xi=x[i];   
              k++;     
          }
          else{
              j++;
              res[k]=yi;
              if(j<(b))yi=y[j];  
              k++;
          }
        }
        else{
              j++;
              res[k]=yi;
              if(j<(b))yi=y[j];  
              k++;
          }
    }
}    



SEXP C_bipart(SEXP parent, SEXP child, SEXP nTips, SEXP maxP){ //, SEXP Nnode){
   int eins=1L, i, j, k, l=length(child), *tmp, *tmp2, *lch, *kl, pi, ci, p, nt=INTEGER(nTips)[0], mp=INTEGER(maxP)[0], ltmp; 
   SEXP ans, ktmp;
   int nnode=1L;
   for(i=1; i<l; i++){
       if(INTEGER(parent)[i-1L] != INTEGER(parent)[i])nnode+=1L;
   }
   tmp = (int *) R_alloc(mp, sizeof(int));
   tmp2 = (int *) R_alloc(mp, sizeof(int));
   lch = (int *) R_alloc(mp+1L, sizeof(int));
   kl = (int *) R_alloc(mp+1L, sizeof(int));
// Nnode  
   PROTECT(ans = allocVector(VECSXP, nnode)); //INTEGER(Nnode)[0]));  
   p=INTEGER(parent)[0];
   k=0L;
   kl[p]=0;
   lch[p]=1;
   tmp[0] = INTEGER(child)[0]; 
   ltmp=1L;
   for(i=1; i<l; i++){ 
        pi = INTEGER(parent)[i]; 
        ci = INTEGER(child)[i];
        if(pi==p){
             if(ci < (nt+1L)){
                 cisort(&ci, tmp, &eins, &ltmp, tmp2);            
                 ltmp += 1L;
                 for(j=0; j<ltmp; j++) tmp[j] = tmp2[j];
             }
             else{
                 cisort(INTEGER(VECTOR_ELT(ans, kl[ci])), tmp, &(lch[ci]), &ltmp, tmp2);                       
                 ltmp += lch[ci]; 
                 for(j=0; j<ltmp; j++) tmp[j] = tmp2[j];                                
             }
             kl[pi]=k; 
             lch[pi] = ltmp;
        }  
        else{
            PROTECT(ktmp = allocVector(INTSXP, ltmp));
            for(j=0; j<ltmp; j++)INTEGER(ktmp)[j] = tmp2[j];
// k???           
            SET_VECTOR_ELT(ans, k, ktmp); 
            UNPROTECT(1); // ktmp

            if(ci < (nt+1)){ 
                 tmp[0]=ci;
                 ltmp=1L; 
            } 
            else{ 
                ltmp=lch[ci];
                for(j=0; j<ltmp; j++)tmp[j] = INTEGER(VECTOR_ELT(ans, kl[ci]))[j];
            }
            k += 1L; 
            p = pi;
        }
   }
// k ??   
   PROTECT(ktmp = allocVector(INTSXP, ltmp));// mp+1L
   for(j=0; j<ltmp; j++)INTEGER(ktmp)[j] = tmp2[j];
   SET_VECTOR_ELT(ans, k, ktmp);
   UNPROTECT(2);
   return(ans);  
}





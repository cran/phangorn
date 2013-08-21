/* 
 * phangorn.c
 *
 * (c) 2008-2013  Klaus Schliep (klaus.schliep@gmail.com)
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
// #include "phangorn.h" 


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

// from coalescentMCMC
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
// 


// char *transa = "N", *transb = "N";
// double one = 1.0, zero = 0.0;
// int ONE = 1L;
//const double ScaleEPS = 1.0/4294967296.0; 
//const double ScaleMAX = 4294967296.0;
//const double ScaleLOG =  (double)log(ScaleEPS);


// static double *LL, *ROOT, *XX;

// static int *SC;

// int nr=0L;
// static double *weight;


// void ll_free(){
//    free(LL);
//    free(SC);
//}

// type of fitch depending on nc e.g. int, long generic C++,  int *data, int *m, 
//void ll_init(int *nr, int *nTips, int *nc, int *k)
//{
//    LL = (double *) calloc(*nr * *nc * *k * *nTips, sizeof(double));
//}


//void optEdge_free(){
//    free(ROOT);
//    free(XX);
//}


//void optEdge_init(int *nr, int *nc, int *k)
//{
//    ROOT = (double *) calloc(*nr * *nc * *k, sizeof(double));
//    XX = (double *) calloc(*nr * *nc * *k, sizeof(double));
//}



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
     
/*
static R_INLINE void matprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z)
{
    F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one, x, &nrx, y, &nry, &zero, z, &nrx);
}


 
SEXP matpro(SEXP X, SEXP Y, SEXP nrx, SEXP ncx, SEXP nry, SEXP ncy) {   
    int nrmx = INTEGER(nrx)[0], ncmx = INTEGER(ncx)[0], nrmy = INTEGER(nry)[0], ncmy = INTEGER(ncy)[0];   
    SEXP ans;
    PROTECT(ans = allocMatrix(REALSXP, nrmx, ncmy));   
    matprod(REAL(X), nrmx, ncmx, REAL(Y), nrmy, ncmy, REAL(ans));
    UNPROTECT(1);
    return(ans);
}

*/

// INLINE raus
// loop unrolling??
// res = tmp * eva
// res = res %*% * evi ?? mit BLAS?? 
// Ziel 100 
// tmp = (double *) R_alloc((*nc) *(*nrs), sizeof(double)); 


/*
static R_INLINE void getP(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double tmp[m], res;
    for(i = 0; i < m; i++) tmp[i] = exp(eva[i] * w * el);
    for(i = 0; i < m; i++){    
        for(j = 0; j < m; j++){
            res = 0.0;    
            for(h = 0; h < m; h++) res += ev[i + h*m] * tmp[h] * evi[h + j*m];
            result[i+j*m] = res;    
        }
    }
}
*/


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

/*
void getPmix(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double tmp[m], res[m*m];
    for(i = 0; i < m; i++) tmp[i] = exp(eva[i] * w * el);
    for(j = 0; j < m; j++){
        for(h = 0; h < m; h++)    res[h + j*m] = evi[h + j*m] * tmp[h];
    }
    matprod(ev, m, m, res, m, m, result);
}
*/

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


//(dad  * (child %*% P))     



static R_INLINE void emult(double *x, double *y, int n){
    for(int i=0; i<n; i++) x[i]*=y[i];
}

/*
//(dad  * (child %*% P))     
SEXP getM4(SEXP dad, SEXP child, SEXP P, SEXP nr, SEXP nc){
    R_len_t i, n=length(P);
    int ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0], j;
    SEXP TMP, RESULT;
    double *tmp; 
    PROTECT(RESULT = allocVector(VECSXP, n));
    for(i=0; i<n; i++){
        PROTECT(TMP = allocMatrix(REALSXP, nrx, ncx));
        tmp = REAL(TMP);
        for(j=0; j<(ncx * nrx); j++) tmp[j]=0.0;
        F77_CALL(dgemm)(transa, transb, &nrx, &ncx, &ncx, &one, REAL(VECTOR_ELT(child, i)), &nrx, REAL(VECTOR_ELT(P, i)), &ncx, &zero, tmp, &nrx);
        emult(tmp,  REAL(VECTOR_ELT(dad, i)), nrx * ncx);
        SET_VECTOR_ELT(RESULT, i, TMP);
        UNPROTECT(1); // TMP
        }
    UNPROTECT(1); //RESULT    
    return(RESULT);    
    }
*/




// fuer optimEdge und optimRooted (Ziel schnellere Konvergenz und robuster)

/*

void moveLL(double *LL, double *child, double *P, int *nr, int *nc, double *tmp){
    int j;
    F77_CALL(dgemm)(transa, transb, nr, nc, nc, &one, child, nr, P, nc, &zero, tmp, nr);
    for(j=0; j<(*nc * *nr); j++) LL[j]/=tmp[j];               
    F77_CALL(dgemm)(transa, transb, nr, nc, nc, &one, LL, nr, P, nc, &zero, tmp, nr);
    for(j=0; j<(*nc * *nr); j++) tmp[j] *= child[j];
} 

void moveLLNew(double *LL, double *child, double *P, int *nr, int *nc){
    double *tmp;
    int j;
    tmp = (double *) R_alloc(*nr * *nc, sizeof(double));
    F77_CALL(dgemm)(transa, transb, nr, nc, nc, &one, child, nr, P, nc, &zero, tmp, nr);
    for(j=0; j<(*nc * *nr); j++) LL[j]/=tmp[j];               
    F77_CALL(dgemm)(transa, transb, nr, nc, nc, &one, LL, nr, P, nc, &zero, tmp, nr);
    for(j=0; j<(*nc * *nr); j++) LL[j]=tmp[j] * child[j];
} 

void moveLL2(double *LL, double *child, double *P, int nr, int nc){
    double *tmp;
    int j;
    tmp = (double *) R_alloc(nr*nc, sizeof(double));
    F77_CALL(dgemm)(transa, transb, &nr, &nc, &nc, &one, child, &nr, P, &nc, &zero, tmp, &nr);
    for(j=0; j<(nc * nr); j++) LL[j]/=tmp[j];               
    F77_CALL(dgemm)(transa, transb, &nr, &nc, &nc, &one, LL, &nr, P, &nc, &zero, tmp, &nr);
    for(j=0; j<(nc * nr); j++) LL[j]=tmp[j] * child[j];
} 

*/


/*
in ml.c

void moveLLtmp(double *LL, double *child, double *P, int *nr, int *nc, double *tmp){
    int j;    
    F77_CALL(dgemm)(transa, transb, nr, nc, nc, &one, child, nr, P, nc, &zero, tmp, nr);
    for(j=0; j<(*nc * *nr); j++) LL[j]/=tmp[j];               
    F77_CALL(dgemm)(transa, transb, nr, nc, nc, &one, LL, nr, P, nc, &zero, tmp, nr);
    for(j=0; j<(*nc * *nr); j++) LL[j]=tmp[j] * child[j];
} 

// 
// #define LINDEX2(i, j) i * (nr*nc) + j * nTips * (nr * nc)
void moveLL0(int *loli, double *eva, double *eve, double *evi, double *el, double *g, int *nr, int *nc, int *k, int *nTips){
    double *tmp, *P;
    tmp = (double *) R_alloc(*nr * *nc, sizeof(double));
    P = (double *) R_alloc(*nc * *nc, sizeof(double));
    int j;
    for(j = 0; j < *k; j++){
        getP(eva, eve, evi, *nc, el[0], g[j], P);
        moveLLtmp(&ROOT[j * *nr * *nc], &LL[LINDEX3(*loli, j)], P, nr, nc, tmp);
    }
}

// function to get node likelihoods!  #define LINDEX3(i, k) (i-*nTips-1L) * (*nr* *nc) + k * *nTips * (*nr * *nc)
SEXP getLL(SEXP ax, SEXP bx, SEXP nrx, SEXP ncx, SEXP nTips){
    int i, nc = INTEGER(ncx)[0], nr = INTEGER(nrx)[0], ntips = INTEGER(nTips)[0],  a = INTEGER(ax)[0], b = INTEGER(bx)[0];
    SEXP RES;
    PROTECT(RES = allocMatrix(REALSXP, nr, nc));
    for(i=0; i<(nr*nc); i++) REAL(RES)[i] = LL[i + (a-ntips-1L) * (nr * nc) + b * ntips * (nr * nc)];
    return(RES);
}

*/


/*

    double *eva, *eve, *evi;  
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));


void prepFS(SEXP child, SEXP parent, SEXP eig, SEXP el, SEXP g, SEXP nrx, SEXP ncx, SEXP nTips){
    R_len_t i, n=length(g);
    int ch = INTEGER(child)[0], pa = INTEGER(parent)[0];
    int nc=INTEGER(ncx)[0], nr=INTEGER(nrx)[0]; //, j
    double *P, *tmp; 
    double *eva, *eve, *evi;  
    tmp = (double *) R_alloc(nr * nc, sizeof(double));
    P = (double *) R_alloc(nc * nc, sizeof(double)); 
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evi = REAL(VECTOR_ELT(eig, 2));      
    for(i=0; i<n; i++){
        getP(eva, eve, evi, nc, el[0], g[j], P);
        helpDAD(ROOT[i * nr * nc], &LL[LINDEX(ch, i)], P, nr, nc, &LL[LINDEX(pa, i)]);
        helpPrep(&LL[LINDEX(ch, i)], &LL[LINDEX(pa, i)], eve, evi, nr, nc, tmp, &XX[i*nr*nc])
        }
}        
  

void getPrepNew{
    int i;
    double *tmp, *X;
//    X = (double *) R_alloc(nr * nc * k, sizeof(double));
    tmp = (double *) R_alloc(nr * nc, sizeof(double));
    PROTECT(RESULT = allocVector(VECSXP, n));
    for(i=0; i<n; i++){
        helpPrep(&LL[LINDEX3(dad, i)], &LL[LINDEX3(child, i)], eve,  evi, nr, nc, tmp, X);
        } 
    }


void helpPrep(double *dad, double *child, double *eve, double *evi, int nr, int nc, double *tmp, double *res){
    F77_CALL(dgemm)(transa, transb, &nr, &nc, &nc, &one, child, &nr, eve, &nc, &zero, res, &nr);
    F77_CALL(dgemm)(transa, transb, &nr, &nc, &nc, &one, dad, &nr, evi, &nc, &zero, tmp, &nr);
    for(int j=0; j<(nc * nr); j++) res[j]*=tmp[j];               
} 

void helpDAD(double *dad, double *child, double *P, int nr, int nc, double *res){
    F77_CALL(dgemm)(transa, transb, &nr, &nc, &nc, &one, child, &nr, P, &nc, &zero, res, &nr);
    for(int j=0; j<(nc * nr); j++) res[j]=dad[j]/res[j];               
} 


SEXP getPrep(SEXP dad, SEXP child, SEXP eve, SEXP evi, SEXP nr, SEXP nc){
    R_len_t i, n=length(dad);
    int ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0]; //, j
    double *tmp;
    SEXP TMP, RESULT;
    tmp = (double *) R_alloc(nrx*ncx, sizeof(double));  
    PROTECT(RESULT = allocVector(VECSXP, n));
    for(i=0; i<n; i++){
        PROTECT(TMP = allocMatrix(REALSXP, nrx, ncx));
        helpPrep(REAL(VECTOR_ELT(dad, i)), REAL(VECTOR_ELT(child, i)), REAL(eve),  REAL(evi), nrx, ncx, tmp, REAL(TMP));
        SET_VECTOR_ELT(RESULT, i, TMP);
        UNPROTECT(1); // TMP
        }
    UNPROTECT(1); //RESULT    
    return(RESULT);    
    }


        
        P <- getP(old.el, eig, g)
        dad <- .Call("getDAD", LL, dat[, child[j]], P, nr, nc)
        X <- .Call("getPrep", dad, dat[, child[j]], eig[[2]], evi, nr, nc) 



getLL = function(data, i, j){
    nr = attr(data, "nr")
    nc = attr(data, nc)
    nTips = as.integer(length(data))
    .Call("getLL", as.integer(i), as.integer(j) , nr, nc, nTips)
}

SEXP moveLL(SEP loli, SEXP eig, SEXP el, SEXP g, SEXP nr, SEXP nc, SEXP k, SEXP nTips){
    double *eva, *eve, *evi; 
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evi = REAL(VECTOR_ELT(eig, 2));
    while(){
        moveLL0(int *loli, eva, eve, evi, double *el, double *g, int nr, int nc, int k, int nTips)
        loli = ;
    }
    return(NULL);
}

create index
i, j, nr, nc auch fuer partition models
und mixtures

 */

/*
SEXP MoveLL(SEXP LL, SEXP child, SEXP P, SEXP nr, SEXP nc){
    R_len_t i, n=length(LL), j;
    int ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0];
    double *tmp; 
    SEXP RES, TMP;
    PROTECT(RES = allocVector(VECSXP, n));
    for(i=0; i<n; i++){
        PROTECT(TMP = allocVector(REALSXP, nrx * ncx));
        tmp = REAL(TMP);
        for(j=0;ncx*nrx;j++)tmp[j] = REAL(VECTOR_ELT(LL, i))[j];
        moveLL(tmp, REAL(VECTOR_ELT(child, i)), REAL(VECTOR_ELT(P, i)), nrx, ncx);     
        SET_VECTOR_ELT(RES, i, TMP);
        UNPROTECT(1);  
    }
    UNPROTECT(1); 
    return(RES); 
}
*/



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


SEXP allChildren(SEXP children, SEXP parent, SEXP tab, SEXP m){
    int i, j, k, l=0L;   
    R_len_t n=length(parent); 
    SEXP RESULT, TMP;
    PROTECT(RESULT = allocVector(VECSXP, INTEGER(m)[0]));
    for(i=0; i<n; i++){ 
        k=INTEGER(tab)[i];
        PROTECT(TMP = allocVector(INTSXP, k));  
        for(j=0; j<k; j++){
            INTEGER(TMP)[j] = INTEGER(children)[l];
            l++;
        } 
        SET_VECTOR_ELT(RESULT, INTEGER(parent)[i]-1L, TMP);
        UNPROTECT(1); 
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



SEXP bipart(SEXP parent, SEXP child, SEXP nTips, SEXP maxP, SEXP Nnode){
   int eins=1L, i, j, k, l=length(child), *tmp, *tmp2, *lch, *kl, pi, ci, p, nt=INTEGER(nTips)[0], mp=INTEGER(maxP)[0], ltmp; 
   SEXP ans, ktmp;
   tmp = (int *) R_alloc(mp, sizeof(int));
   tmp2 = (int *) R_alloc(mp, sizeof(int));
   lch = (int *) R_alloc(mp+1L, sizeof(int));
   kl = (int *) R_alloc(mp+1L, sizeof(int));
 
   PROTECT(ans = allocVector(VECSXP, INTEGER(Nnode)[0]));  
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
   PROTECT(ktmp = allocVector(INTSXP, ltmp));// mp+1L
   for(j=0; j<ltmp; j++)INTEGER(ktmp)[j] = tmp2[j];
   SET_VECTOR_ELT(ans, k, ktmp);
   UNPROTECT(2);
   return(ans);  
}



// funktioniert leider noch nicht
/*
void lll2(SEXP dlist, double *eva, double *eve, double *evi, double *el, double g, int *nr, int *nc, int *node, int *edge, int nTips, double *contrast, int nco, int n, int *scaleTmp, double *bf, double *TMP, double *ans){
    int  ni, ei, j, i, rc; //    R_len_t i, n = length(node);
    double *rtmp, *P;

    ni = -1;
    rc = *nr * *nc;
    rtmp = (double *) R_alloc(*nr * *nc, sizeof(double));
    P = (double *) R_alloc(*nc * *nc, sizeof(double));

    for(j=0; j < *nr; j++) scaleTmp[j] = 0L;

    for(i = 0; i < n; i++) {
        getP(eva, eve, evi, *nc, el[i], g, P); 
        ei = edge[i]; 
        if(ni != node[i]){
            if(i>0)scaleMatrixNew(&ans[ni * rc], nr, nc, scaleTmp); // (ni-nTips)
            ni = node[i];
            if(ei < nTips) 
                matp(INTEGER(VECTOR_ELT(dlist, ei)), contrast, P, nr, nc, &nco, &ans[ni * rc]); 
            else 
                F77_CALL(dgemm)(transa, transb, nr, nc, nc, &one, &ans[(ei-nTips) * rc], nr, P, nc, &zero, &ans[ni * rc], nr);
        }
        else {
            if(ei < nTips) 
                matp(INTEGER(VECTOR_ELT(dlist, ei)), contrast, P, nr, nc, &nco, rtmp);
            else 
                F77_CALL(dgemm)(transa, transb, nr, nc, nc, &one, &ans[(ei-nTips) * rc], nr, P, nc, &zero, rtmp, nr);
            for(j=0; j < rc; j++) ans[ni * rc + j] *= rtmp[j];
        }            
    }
    scaleMatrixNew(&ans[ni * rc], nr, nc, scaleTmp);
    F77_CALL(dgemv)(transa, nr, nc, &one, &ans[ni * rc], nr, bf, &ONE, &zero, TMP, &ONE); 
}
*/

// 1955 Zeilen pmlPart in pml einfuegen
// Teile von pmlMix ermoeglichen
// Ziel 100 Zeilen weniger fuer 1.99-0 < 1700



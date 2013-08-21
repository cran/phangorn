/* 
 * ml.c
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


// index likelihood pml
// need to define nr, nc, nTips, nNodes k

#define LINDEX(i, k) (i - ntips - 1L) * (nr * nc) + k * ntips * (nr * nc)
#define LINDEX2(i, k) (i - *ntips - 1L) * (*nr* *nc) + k * *ntips * (*nr * *nc)


char *transa = "N", *transb = "N";
double one = 1.0, zero = 0.0;
int ONE = 1L;
const double ScaleEPS = 1.0/4294967296.0; 
const double ScaleMAX = 4294967296.0;


static double *LL, *ROOT, *XX;
static int *SCM;

void ll_free(){
    free(LL);
    free(SCM);
    free(ROOT);
//    free(XX);
}


void ll_init(int *nr, int *nTips, int *nc, int *k)
{
    LL = (double *) calloc(*nr * *nc * *k * *nTips, sizeof(double));
    ROOT = (double *) calloc(*nr * *nc * *k, sizeof(double));
//    XX = (double *) calloc(*nr * *nc * *k, sizeof(double)); 
    SCM = (int *) calloc(*nr * *k * 2L * *nTips, sizeof(int));  
}

/*
void optEdge_free(){
    free(ROOT);
    free(XX);
}

void optEdge_init(int *nr, int *nc, int *k)
{
    ROOT = (double *) calloc(*nr * *nc * *k, sizeof(double));
    XX = (double *) calloc(*nr * *nc * *k, sizeof(double));
}
*/

void matm(int *x, double *contrast, int *nr, int *nc, int *nco, double *result){
    int i, j;
    for(i = 0; i < (*nr); i++){ 
        for(j = 0; j < (*nc); j++) result[i + j*(*nr)] *= contrast[x[i] - 1L + j*(*nco)];  
    }
}


SEXP invSites(SEXP dlist, SEXP nr, SEXP nc, SEXP contrast, SEXP nco){
    R_len_t n = length(dlist);
    int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0], i, j;
    SEXP result;  
    PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
    double *res;    
    res = REAL(result);
    for(j=0; j < (nrx * ncx); j++) res[j] = 1.0;
    for(i=0; i < n; i++) matm(INTEGER(VECTOR_ELT(dlist, i)), REAL(contrast), INTEGER(nr), INTEGER(nc), INTEGER(nco), res);   
    UNPROTECT(1); // result 
    return(result);
}     


void scaleMatrixD(double *X, int *nr, int *nc, double *result){
    int i, j; 
    double tmp;
    for(i = 0; i < *nr; i++) {    
        tmp = 0.0; 
        for(j = 0; j < *nc; j++) tmp += X[i + j* *nr]; 
        result[i] +=log(tmp); 
        tmp = 1.0/tmp; 
        for(j = 0; j < *nc; j++) X[i + j* *nr] *=tmp;  
    } 
}

void scaleMatrix(double *X, int *nr, int *nc, int *result){
    int i, j; 
    double tmp;
    for(i = 0; i < *nr; i++) {    
        tmp = 0.0; 
        for(j = 0; j < *nc; j++) tmp += X[i + j* *nr];    
        while(tmp < ScaleEPS){
           for(j = 0; j < *nc; j++) X[i + j* *nr] *=ScaleMAX;
           result[i] +=1L;
           tmp *= ScaleMAX;
       }        
    } 
}


// contrast to full
void matp(int *x, double *contrast, double *P, int *nr, int *nc, int *nrs, double *result){
    int i, j;
    double *tmp; 
    tmp = (double *) R_alloc((*nc) *(*nrs), sizeof(double)); 
//    matprod(contrast, (*nrs), (*nc), P, (*nc), (*nc), tmp);  
    F77_CALL(dgemm)(transa, transb, nrs, nc, nc, &one, contrast, nrs, P, nc, &zero, tmp, nrs);
    for(i = 0; i < (*nr); i++){ 
        for(j = 0; j < (*nc); j++) result[i + j*(*nr)] = tmp[x[i] - 1L + j*(*nrs)];  
    }
}

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



SEXP getPM(SEXP eig, SEXP nc, SEXP el, SEXP w){
    R_len_t i, j, nel, nw, k;
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
            if(edgelen[j]==0.0){
                for(k=0; k<(m*m);k++)REAL(P)[k]=0.0;
                for(k=0; k<m; k++)REAL(P)[k+k*m]=1.0;
            }
            else getP(eva, eve, evei, m, edgelen[j], ws[i], REAL(P));
            SET_VECTOR_ELT(RESULT, l, P);
            UNPROTECT(1); 
            l++;
        }
    }
    UNPROTECT(1);//RESULT
    return(RESULT);
} 



// 
void lll(SEXP dlist, SEXP eig, double *el, double g, int *nr, int *nc, int *node, int *edge, int nTips, double *contrast, int nco, int n, int *scaleTmp, double *bf, double *TMP, double *ans){
    int  ni, ei, j, i, rc; //    R_len_t i, n = length(node);
    double *rtmp, *P;
    double *eva, *eve, *evei;
 
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));

    ni = -1;
    rc = *nr * *nc;
    rtmp = (double *) R_alloc(*nr * *nc, sizeof(double));
    P = (double *) R_alloc(*nc * *nc, sizeof(double));

    for(j=0; j < *nr; j++) scaleTmp[j] = 0L;

    for(i = 0; i < n; i++) {
        getP(eva, eve, evei, *nc, el[i], g, P); 
        ei = edge[i]; 
        if(ni != node[i]){
            if(ni>0)scaleMatrix(&ans[ni * rc], nr, nc, scaleTmp); // (ni-nTips)
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
    scaleMatrix(&ans[ni * rc], nr, nc, scaleTmp);
    F77_CALL(dgemv)(transa, nr, nc, &one, &ans[ni * rc], nr, bf, &ONE, &zero, TMP, &ONE);
}

void lll2(SEXP dlist, SEXP eig, double *el, double g, int *nr, int *nc, int *node, int *edge, int nTips, double *contrast, int nco, int n, double *scaleTmp, double *bf, double *TMP, double *ans){
    int  ni, ei, j, i, rc; //    R_len_t i, n = length(node);
    double *rtmp, *P;
    double *eva, *eve, *evei;
 
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));

    ni = -1;
    rc = *nr * *nc;
    rtmp = (double *) R_alloc(*nr * *nc, sizeof(double));
    P = (double *) R_alloc(*nc * *nc, sizeof(double));

    for(j=0; j < *nr; j++) scaleTmp[j] = 0L;

    for(i = 0; i < n; i++) {
        getP(eva, eve, evei, *nc, el[i], g, P); 
        ei = edge[i]; 
        if(ni != node[i]){
            if(ni>0)scaleMatrixD(&ans[ni * rc], nr, nc, scaleTmp); // (ni-nTips)
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
    scaleMatrixD(&ans[ni * rc], nr, nc, scaleTmp);
    F77_CALL(dgemv)(transa, nr, nc, &one, &ans[ni * rc], nr, bf, &ONE, &zero, TMP, &ONE);
}


SEXP PML2(SEXP dlist, SEXP EL, SEXP W, SEXP G, SEXP NR, SEXP NC, SEXP K, SEXP eig, SEXP bf, SEXP node, SEXP edge, SEXP NTips, SEXP root, SEXP nco, SEXP contrast, SEXP N){
    int nr=INTEGER(NR)[0], nc=INTEGER(NC)[0], k=INTEGER(K)[0], i, indLL; 
    int nTips = INTEGER(NTips)[0];
    double *g=REAL(G), *tmp;
    SEXP ans, TMP, SC;
    PROTECT(ans = allocVector(VECSXP, 2L));  // 3L 
    PROTECT(SC = allocMatrix(REALSXP, nr, k));
    PROTECT(TMP = allocMatrix(REALSXP, nr, k)); // changed
    tmp=REAL(TMP);
    for(i=0; i<(k*nr); i++)tmp[i]=0.0;
    indLL = nr * nc * nTips;
    for(i=0; i<k; i++){                  
        lll2(dlist, eig, REAL(EL), g[i], &nr, &nc, INTEGER(node), INTEGER(edge), nTips, REAL(contrast), INTEGER(nco)[0], INTEGER(N)[0], &(REAL(SC)[nr * i]), REAL(bf), &tmp[i*nr], &LL[indLL *i]);           
     }   
     SET_VECTOR_ELT(ans, 0L, TMP);
     SET_VECTOR_ELT(ans, 1L, SC); 
     UNPROTECT(3);
     return ans;     
}


// TODO 
// openMP pragma, requires pure C code, lll2 statt lll
// long vector support (maybe C++)
SEXP PML(SEXP dlist, SEXP EL, SEXP W, SEXP G, SEXP NR, SEXP NC, SEXP K, SEXP eig, SEXP bf, SEXP node, SEXP edge, SEXP NTips, SEXP root, SEXP nco, SEXP contrast, SEXP N){
    int nr=INTEGER(NR)[0], nc=INTEGER(NC)[0], k=INTEGER(K)[0], i, indLL; 
    int nTips = INTEGER(NTips)[0];
    double *g=REAL(G), *tmp;
    SEXP ans, TMP, SC;
    PROTECT(ans = allocVector(VECSXP, 2L));  // 3L 
    PROTECT(SC = allocMatrix(INTSXP, nr, k));
    PROTECT(TMP = allocMatrix(REALSXP, nr, k)); // changed
    tmp=REAL(TMP);
    for(i=0; i<(k*nr); i++)tmp[i]=0.0;
    indLL = nr * nc * nTips;
    for(i=0; i<k; i++){                  
        lll(dlist, eig, REAL(EL), g[i], &nr, &nc, INTEGER(node), INTEGER(edge), nTips, REAL(contrast), INTEGER(nco)[0], INTEGER(N)[0], &(INTEGER(SC)[nr * i]), REAL(bf), &tmp[i*nr], &LL[indLL *i]);           
     }   
     SET_VECTOR_ELT(ans, 0L, TMP);
     SET_VECTOR_ELT(ans, 1L, SC); 
     UNPROTECT(3);
     return ans;     
}


// etwas compakterer code weniger movement
SEXP PML0(SEXP dlist, SEXP EL, SEXP W, SEXP G, SEXP NR, SEXP NC, SEXP K, SEXP eig, SEXP bf, SEXP node, SEXP edge, SEXP NTips, SEXP root, SEXP nco, SEXP contrast, SEXP N){
    int nr=INTEGER(NR)[0], nc=INTEGER(NC)[0], k=INTEGER(K)[0], i, indLL; 
    int nTips = INTEGER(NTips)[0], *SC;
    double *g=REAL(G), *tmp, logScaleEPS;
    SEXP TMP;
    SC = (int *) R_alloc(nr * nc, sizeof(int));   
//    PROTECT(ans = allocVector(VECSXP, 2L));  // 3L 
//    PROTECT(SC = allocMatrix(INTSXP, nr, k));
    PROTECT(TMP = allocMatrix(REALSXP, nr, k)); // changed
    tmp=REAL(TMP);
    for(i=0; i<(k*nr); i++)tmp[i]=0.0;
    indLL = nr * nc * nTips;  
    for(i=0; i<k; i++){                  
        lll(dlist, eig, REAL(EL), g[i], &nr, &nc, INTEGER(node), INTEGER(edge), nTips, REAL(contrast), INTEGER(nco)[0], INTEGER(N)[0], &SC[nr * i], REAL(bf), &tmp[i*nr], &LL[indLL *i]);           
     } 
/*
    ScaleEPS = log(1.0/4294967296.0)
    resll = tmp[[2]] * ScaleEPS + log(tmp[[1]])
*/
    logScaleEPS = log(ScaleEPS);
    for(i=0; i<(k*nr); i++) tmp[i] = logScaleEPS * SC[i] + log(tmp[i]);     
//     SET_VECTOR_ELT(ans, 0L, TMP);
//     SET_VECTOR_ELT(ans, 1L, SC); 
     UNPROTECT(1);
     return TMP;     
}


// replace child with LL
void moveLLNew(double *LL, double *child, double *P, int *nr, int *nc, double *tmp){
    int j;
    F77_CALL(dgemm)(transa, transb, nr, nc, nc, &one, child, nr, P, nc, &zero, tmp, nr);
    for(j=0; j<(*nc * *nr); j++) LL[j]/=tmp[j]; // new child              
    F77_CALL(dgemm)(transa, transb, nr, nc, nc, &one, LL, nr, P, nc, &zero, tmp, nr);
    for(j=0; j<(*nc * *nr); j++) child[j] *= tmp[j];
} 

void moveLL2(int *loli, int *nloli, double *eva, double *eve, double *evi, double *el, double *g, int *nr, int *nc, int *k, int *ntips){
    double *tmp, *P;
    int j;
    tmp = (double *) R_alloc(*nr * *nc, sizeof(double));
    P = (double *) R_alloc(*nc * *nc, sizeof(double));
    for(j = 0; j < *k; j++){
        getP(eva, eve, evi, *nc, el[0], g[j], P);
        moveLLNew(&LL[LINDEX2(*loli, j)], &LL[LINDEX2(*nloli, j)], P, nr, nc, tmp);
    }
}

void moveLLtmp(double *LL, double *child, double *P, int *nr, int *nc, double *tmp){
    int j;    
    F77_CALL(dgemm)(transa, transb, nr, nc, nc, &one, child, nr, P, nc, &zero, tmp, nr);
    for(j=0; j<(*nc * *nr); j++) LL[j]/=tmp[j];               
    F77_CALL(dgemm)(transa, transb, nr, nc, nc, &one, LL, nr, P, nc, &zero, tmp, nr);
    for(j=0; j<(*nc * *nr); j++) LL[j]=tmp[j] * child[j];
} 


void moveLL0(int *loli, double *eva, double *eve, double *evi, double *el, double *g, int *nr, int *nc, int *k, int *ntips){
    double *tmp, *P;
    tmp = (double *) R_alloc(*nr * *nc, sizeof(double));
    P = (double *) R_alloc(*nc * *nc, sizeof(double));
    int j;
    for(j = 0; j < *k; j++){
        getP(eva, eve, evi, *nc, el[0], g[j], P);
        moveLLtmp(&ROOT[j * *nr * *nc], &LL[LINDEX2(*loli, j)], P, nr, nc, tmp);
    }
}

// old version
void moveLL(double *LL, double *child, double *P, int *nr, int *nc, double *tmp){
    int j;
    F77_CALL(dgemm)(transa, transb, nr, nc, nc, &one, child, nr, P, nc, &zero, tmp, nr);
    for(j=0; j<(*nc * *nr); j++) LL[j]/=tmp[j];               
    F77_CALL(dgemm)(transa, transb, nr, nc, nc, &one, LL, nr, P, nc, &zero, tmp, nr);
    for(j=0; j<(*nc * *nr); j++) tmp[j] *= child[j];
} 




void helpDAD3(double *dad, double *child, double *P, int *nr, int *nc, double *res){
    F77_CALL(dgemm)(transa, transb, nr, nc, nc, &one, child, nr, P, nc, &zero, res, nr);
    for(int j=0; j<(*nc * *nr); j++) dad[j]/=res[j];               
} 

void getDAD3(int *dad, int *child, double *eva, double *eve, double *evi, double *el, double *g, int *nr, int *nc, int *k, int *ntips){
    double *tmp, *P;
    int j;
    tmp = (double *) R_alloc(*nr * *nc, sizeof(double));
    P = (double *) R_alloc(*nc * *nc, sizeof(double));
    for(j = 0; j < *k; j++){
        getP(eva, eve, evi, *nc, el[0], g[j], P);
        helpDAD3(&LL[LINDEX2(*dad, j)], &LL[LINDEX2(*child, j)], P, nr, nc, tmp);
    }
}


// children
void helpDAD4(double *dad, int *child, double *contrast, double *P, int *nr, int *nc, int *nco, double *res){
    matp(child, contrast, P, nr, nc, nco, res); 
    for(int j=0; j<(*nc * *nr); j++) res[j]=dad[j]/res[j];               
} 


void getDAD4(int *dad, int *child, double *contrast, double *eva, double *eve, double *evi, double *el, double *g, int *nr, int *nc, int *nco, int *k, int *ntips){
    double *tmp, *P;
    int j;
    tmp = (double *) R_alloc(*nr * *nc, sizeof(double));
    P = (double *) R_alloc(*nc * *nc, sizeof(double));
    for(j = 0; j < *k; j++){
        getP(eva, eve, evi, *nc, el[0], g[j], P);
        helpDAD4(&LL[LINDEX2(*dad, j)], child, contrast, P, nr, nc, nco, tmp);
    }
}


// dad / child * P 
void helpDAD(double *dad, double *child, double *P, int nr, int nc, double *res){
    F77_CALL(dgemm)(transa, transb, &nr, &nc, &nc, &one, child, &nr, P, &nc, &zero, res, &nr);
    for(int j=0; j<(nc * nr); j++) res[j]=dad[j]/res[j];               
} 


SEXP getDAD(SEXP dad, SEXP child, SEXP P, SEXP nr, SEXP nc){
    R_len_t i, n=length(P);
    int ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0]; //, j
    SEXP TMP, RESULT;
    PROTECT(RESULT = allocVector(VECSXP, n));
    for(i=0; i<n; i++){
        PROTECT(TMP = allocMatrix(REALSXP, nrx, ncx));
        helpDAD(REAL(VECTOR_ELT(dad, i)), REAL(VECTOR_ELT(child, i)), REAL(VECTOR_ELT(P, i)), nrx, ncx, REAL(TMP));
        SET_VECTOR_ELT(RESULT, i, TMP);
        UNPROTECT(1); // TMP
        }
    UNPROTECT(1); //RESULT    
    return(RESULT);    
    }



void helpPrep(double *dad, double *child, double *eve, double *evi, int nr, int nc, double *tmp, double *res){
    F77_CALL(dgemm)(transa, transb, &nr, &nc, &nc, &one, child, &nr, eve, &nc, &zero, res, &nr);
    F77_CALL(dgemm)(transa, transb, &nr, &nc, &nc, &one, dad, &nr, evi, &nc, &zero, tmp, &nr);
    for(int j=0; j<(nc * nr); j++) res[j]*=tmp[j];               
} 



    
void prepFS(double *XX, int ch, int pa, double *eva, double *eve, double *evi, double el, double *g, int nr, int nc, int ntips, int k){    
    int i;
    double *P, *tmp; 
    tmp = (double *) R_alloc(nr * nc, sizeof(double));
    P = (double *) R_alloc(nc * nc, sizeof(double)); 
    for(i=0; i<k; i++){
        getP(eva, eve, evi, nc, el, g[i], P);
        helpDAD(&ROOT[i * nr * nc], &LL[LINDEX(ch, i)], P, nr, nc, &LL[LINDEX(pa, i)]);
        helpPrep(&LL[LINDEX(ch, i)], &LL[LINDEX(pa, i)], eve, evi, nr, nc, tmp, &XX[i*nr*nc]);
        }
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


void helpDAD2(double *dad, int *child, double *contrast, double *P, int nr, int nc, int nco, double *res){
    matp(child, contrast, P, &nr, &nc, &nco, res); 
    for(int j=0; j<(nc * nr); j++) res[j]=dad[j]/res[j];               
} 

SEXP getDAD2(SEXP dad, SEXP child, SEXP contrast, SEXP P, SEXP nr, SEXP nc, SEXP nco){
    R_len_t i, n=length(P);
    int ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0], nrs=INTEGER(nco)[0]; //, j
    SEXP TMP, RESULT;
    PROTECT(RESULT = allocVector(VECSXP, n));
    for(i=0; i<n; i++){
        PROTECT(TMP = allocMatrix(REALSXP, nrx, ncx));
        helpDAD2(REAL(VECTOR_ELT(dad, i)), INTEGER(child), REAL(contrast), REAL(VECTOR_ELT(P, i)), nrx, ncx, nrs, REAL(TMP));
        SET_VECTOR_ELT(RESULT, i, TMP);
        UNPROTECT(1); // TMP
        }
    UNPROTECT(1); //RESULT    
    return(RESULT);    
    }



void helpPrep2(double *dad, int *child, double *contrast, double *evi, int nr, int nc, int nrs, double *res){
    int i, j;
    F77_CALL(dgemm)(transa, transb, &nr, &nc, &nc, &one, dad, &nr, evi, &nc, &zero, res, &nr);
    for(i = 0; i < nr; i++){ 
        for(j = 0; j < nc; j++) res[i + j*nr] *= contrast[child[i] - 1L + j*nrs];  
    }                  
} 


SEXP getPrep2(SEXP dad, SEXP child, SEXP contrast, SEXP evi, SEXP nr, SEXP nc, SEXP nco){
    R_len_t i, n=length(dad);
    int ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0], ncs=INTEGER(nco)[0]; 
    SEXP TMP, RESULT; 
    PROTECT(RESULT = allocVector(VECSXP, n));
    for(i=0; i<n; i++){
        PROTECT(TMP = allocMatrix(REALSXP, nrx, ncx));
        helpPrep2(REAL(VECTOR_ELT(dad, i)), INTEGER(child), REAL(contrast),  REAL(evi), nrx, ncx, ncs, REAL(TMP));
        SET_VECTOR_ELT(RESULT, i, TMP);
        UNPROTECT(1); 
        }
    UNPROTECT(1);     
    return(RESULT);    
    }


void prepFSE(double *XX, int *ch, int pa, double *eva, double *eve, double *evi, double el, double *g, int nr, int nc, int ntips, int k, double *contrast, double *contrast2, int ncs){    
    int i;
    double *P; //, *tmp 
//    tmp = (double *) R_alloc(nr * nc, sizeof(double));
    P = (double *) R_alloc(nc * nc, sizeof(double)); 
    for(i=0; i<k; i++){
        getP(eva, eve, evi, nc, el, g[i], P);
        helpDAD2(&ROOT[i * nr * nc], ch, contrast, P, nr, nc, ncs, &LL[LINDEX(pa, i)]);
//        helpDAD(&ROOT[i * nr * nc], &LL[LINDEX(ch, i)], P, nr, nc, &LL[LINDEX(pa, i)]);
//        helpPrep(&LL[LINDEX(ch, i)], &LL[LINDEX(pa, i)], eve, evi, nr, nc, tmp, &XX[i*nr*nc]);
        helpPrep2(&LL[LINDEX(pa, i)], ch, contrast2,  evi, nr, nc, ncs, &XX[i*nr*nc]);
        }
}        



/*

library(phangorn)

getNodeLogLik = function(data, i, j){
    nr = attr(data, "nr")
    nc = attr(data, "nc")
    ntips = length(data)
    browser()
    .Call("getLL", as.integer(i), as.integer(j), as.integer(nr), as.integer(nc), as.integer(ntips))
}

getNodeLogLik = function(i=48, j=0, nr=1605, nc=4, ntips=47){
    .Call("getLL", as.integer(i), as.integer(j), as.integer(nr), as.integer(nc), as.integer(ntips))
}

getNodeLogLik2 = function(i=48, j=0, nr=1605, nc=4, ntips=47){
    .Call("getLL2", as.integer(i), as.integer(j), as.integer(nr), as.integer(nc), as.integer(ntips))
}

getRoot = function(i=0, nr=1605, nc=4){
    .Call("getROOT", as.integer(i), as.integer(nr), as.integer(nc))
}

getXX = function(i=0, nr=1605, nc=4){
    .Call("getXX", as.integer(i), as.integer(nr), as.integer(nc))
}

dyn.load("ml.so")
data(yeast)
getNodeLogLik(yeast, 10, 0)



*/

// function to get node likelihoods!  #define LINDEX3(i, k) (i-*nTips-1L) * (*nr* *nc) + k * *nTips * (*nr * *nc)
SEXP getLL2(SEXP ax, SEXP bx, SEXP nrx, SEXP ncx, SEXP nTips){
    int nc = INTEGER(ncx)[0], nr = INTEGER(nrx)[0], ntips = INTEGER(nTips)[0],  a = INTEGER(ax)[0], b = INTEGER(bx)[0];
//    int j, start=LINDEX(i,j);
    SEXP RES;
//    PROTECT(RES = allocMatrix(REALSXP, nr, nc));
  PROTECT(RES = allocVector(INTSXP, 2));
//    for(j=0; j<(nr*nc); j++) REAL(RES)[i] = LL[j + LINDEX(a, b)];
   INTEGER(RES)[0] = LINDEX(a, b);
   INTEGER(RES)[1] = (a-ntips-1) * (nr * nc) + b * ntips * (nr * nc);

    //(a-ntips-1L) * (nr * nc) + b * ntips * (nr * nc)
    UNPROTECT(1);
    return(RES);
}


// function to get node likelihoods!  #define LINDEX3(i, k) (i-*nTips-1L) * (*nr* *nc) + k * *nTips * (*nr * *nc)
SEXP getLL3(SEXP ax, SEXP bx, SEXP nrx, SEXP ncx, SEXP nTips){
    int *nc = INTEGER(ncx), *nr = INTEGER(nrx), *ntips = INTEGER(nTips),  a = INTEGER(ax)[0], b = INTEGER(bx)[0];
//    int j, start=LINDEX(i,j);
    SEXP RES;
//    PROTECT(RES = allocMatrix(REALSXP, nr, nc));
  PROTECT(RES = allocVector(INTSXP, 2));
//    for(j=0; j<(nr*nc); j++) REAL(RES)[i] = LL[j + LINDEX(a, b)];
   INTEGER(RES)[0] = LINDEX2(a, b);
   INTEGER(RES)[1] = 0;

    //(a-ntips-1L) * (nr * nc) + b * ntips * (nr * nc)
    UNPROTECT(1);
    return(RES);
}



SEXP getLL(SEXP ax, SEXP bx, SEXP nrx, SEXP ncx, SEXP nTips){
    int j, nc = INTEGER(ncx)[0], nr = INTEGER(nrx)[0], ntips = INTEGER(nTips)[0],  a = INTEGER(ax)[0], b = INTEGER(bx)[0];
    SEXP RES;
    PROTECT(RES = allocMatrix(REALSXP, nr, nc));
    for(j=0; j<(nr*nc); j++) REAL(RES)[j] = LL[j + LINDEX(a, b)];
    UNPROTECT(1);
    return(RES);
}

SEXP getROOT(SEXP ax, SEXP nrx, SEXP ncx){
    int j, nc = INTEGER(ncx)[0], nr = INTEGER(nrx)[0], a = INTEGER(ax)[0];
    SEXP RES;
    PROTECT(RES = allocMatrix(REALSXP, nr, nc));
    for(j=0; j<(nr*nc); j++) REAL(RES)[j] = ROOT[j + a * nr * nc];
    UNPROTECT(1);
    return(RES);
}

SEXP getXX(SEXP ax, SEXP nrx, SEXP ncx){
    int j, nc = INTEGER(ncx)[0], nr = INTEGER(nrx)[0], a = INTEGER(ax)[0];
    SEXP RES;
    PROTECT(RES = allocMatrix(REALSXP, nr, nc));
    for(j=0; j<(nr*nc); j++) REAL(RES)[j] = XX[j + a * nr * nc];
    UNPROTECT(1);
    return(RES);
}

// ch * (pa %*% P)  
void getMI(int child, int parent, double el, double *eva, double *eve, double *evi, double *g, int nr, int nc, int k, int ntips){
    int i, a, j;
    double *P;
    P = (double *) R_alloc(nc * nc, sizeof(double)); 
    for(i=0; i<k; i++){
        getP(eva, eve, evi, nc, el, g[i], P);
        a = i * nr * nc;
        F77_CALL(dgemm)(transa, transb, &nr, &nc, &nc, &one, &LL[LINDEX(parent, i)], &nr, P, &nc, &zero, &ROOT[a], &nr);       
        for(j=0; j<(nc * nr); j++) ROOT[a + j]*=LL[LINDEX(child, i) + j];
    }    
}
 
// (ch %*% P) * pa
void getME(int *child, int parent, double el, double *eva, double *eve, double *evi, double *g, int nr, 
   int nc, int k, int ntips, double *contrast, int nrs){
    int i, a, j;
    double *P;
    P = (double *) R_alloc(nc * nc, sizeof(double)); 
    for(i=0; i<k; i++){
        getP(eva, eve, evi, nc, el, g[i], P);
        a = i * nr * nc;
        matp(child, contrast, P, &nr, &nc, &nrs, &ROOT[a]);   
        for(j=0; j<(nc * nr); j++) ROOT[a + j]*=LL[LINDEX(parent, i) + j];
    }    
}


void NR55(double *eva, int nc, double el, double *w, double *g, SEXP X, int ld, int nr, double *f, double *res){
    int i, j, k; 
    double *tmp;  
    tmp = (double *) R_alloc(nc, sizeof(double));

    for(k=0; k<nr; k++) res[k] = 0.0;
        
    for(j=0;j<ld;j++){
        for(i=0; i<nc ;i++) tmp[i] = (eva[i] * g[j]  * el) * exp(eva[i] * g[j] * el);       
        F77_CALL(dgemv)(transa, &nr, &nc, &w[j], REAL(VECTOR_ELT(X, j)), &nr, tmp, &ONE, &one, res, &ONE); 
        }
    for(i=0; i<nr ;i++) res[i]/=f[i];                
} 


void NR555(double *eva, int nc, double el, double *w, double *g, SEXP X, int ld, int nr, double *f, double *res){
    int i, j, k; 
    double *tmp;  
    tmp = (double *) R_alloc(nc, sizeof(double));

    for(k=0; k<nr; k++) res[k] = 0.0;
        
    for(j=0;j<ld;j++){
        for(i=0; i<nc ;i++) tmp[i] = (eva[i] * g[j]) * exp(eva[i] * g[j] * el);       
        F77_CALL(dgemv)(transa, &nr, &nc, &w[j], REAL(VECTOR_ELT(X, j)), &nr, tmp, &ONE, &one, res, &ONE); 
        }
    for(i=0; i<nr ;i++) res[i]/=f[i];                
} 

 
void NR66(double *eva, int nc, double el, double *w, double *g, SEXP X, int ld, int nr, double *res){
    int i, j;   
    double *tmp; //*res,  *dF,
 
    tmp = (double *) R_alloc(nc, sizeof(double));
        
    for(j=0;j<ld;j++){
        for(i=0; i<nc ;i++) tmp[i] = exp(eva[i] * g[j] * el);      
        // alpha = w[j] 
        F77_CALL(dgemv)(transa, &nr, &nc, &w[j], REAL(VECTOR_ELT(X, j)), &nr, tmp, &ONE, &one, res, &ONE); 
        }               
} 



// spaeter raus
SEXP LogLik2(SEXP dlist, SEXP P, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP nTips, SEXP mNodes, SEXP contrast, SEXP nco){
    R_len_t i, n = length(node);
    int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0], nt=INTEGER(nTips)[0], mn=INTEGER(mNodes)[0];
    int  ni, ei, j, *edges=INTEGER(edge), *nodes=INTEGER(node);
    SEXP ans, result;  
    double *res, *rtmp;
    if(!isNewList(dlist)) error("'dlist' must be a list");
    ni = nodes[0];
    PROTECT(ans = allocVector(VECSXP, mn)); 
    PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
    res = REAL(result);
    rtmp = (double *) R_alloc(nrx*ncx, sizeof(double));
    for(j=0; j < (nrx * ncx); j++) res[j] = 1.0;
    for(i = 0; i < n; i++) {
        ei = edges[i]; 
        if(ni != nodes[i]){
            SET_VECTOR_ELT(ans, ni, result);
            UNPROTECT(1); //result
            PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
            res = REAL(result);
            ni = nodes[i];
            if(ei < nt) 
               matp(INTEGER(VECTOR_ELT(dlist, ei)), REAL(contrast), REAL(VECTOR_ELT(P, i)), INTEGER(nr), INTEGER(nc), INTEGER(nco), res);
            else 
               F77_CALL(dgemm)(transa, transb, &nrx, &ncx, &ncx, &one, REAL(VECTOR_ELT(ans, ei-nt)), &nrx, 
                   REAL(VECTOR_ELT(P, i)), &ncx, &zero, res, &nrx);
            }
        else {
            if(ei < nt) 
                matp(INTEGER(VECTOR_ELT(dlist, ei)), REAL(contrast), REAL(VECTOR_ELT(P, i)), INTEGER(nr), INTEGER(nc), INTEGER(nco), rtmp);
            else 
                F77_CALL(dgemm)(transa, transb, &nrx, &ncx, &ncx, &one, REAL(VECTOR_ELT(ans, ei-nt)), &nrx, 
                    REAL(VECTOR_ELT(P, i)), &ncx, &zero, rtmp, &nrx);
            for(j=0; j < (nrx*ncx); j++) res[j] *= rtmp[j];
        }
    }
    SET_VECTOR_ELT(ans, ni, result);  
    UNPROTECT(2); // result ans 
    return(ans);
}


static R_INLINE void matprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z)
{
    F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one, x, &nrx, y, &nry, &zero, z, &nrx);
}


SEXP getM3(SEXP dad, SEXP child, SEXP P, SEXP nr, SEXP nc){
    R_len_t i, n=length(P);
    int ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0], j;
    SEXP TMP, RESULT;
    double *tmp, *daddy;
    PROTECT(RESULT = allocVector(VECSXP, n));
    for(i=0; i<n; i++){
        PROTECT(TMP = allocMatrix(REALSXP, nrx, ncx));
        tmp = REAL(TMP);
        matprod(REAL(VECTOR_ELT(child, i)), nrx, ncx, REAL(VECTOR_ELT(P, i)), ncx, ncx, tmp);
        daddy = REAL(VECTOR_ELT(dad, i));
        for(j=0; j<(ncx * nrx); j++) tmp[j]*=daddy[j];
        SET_VECTOR_ELT(RESULT, i, TMP);
        UNPROTECT(1); // TMP
        }
    UNPROTECT(1); //RESULT    
    return(RESULT);    
    }
    
SEXP FS4(SEXP eig, SEXP nc, SEXP el, SEXP w, SEXP g, SEXP X, SEXP dad, SEXP child, SEXP ld, SEXP nr, 
         SEXP basefreq, SEXP weight, SEXP f0, SEXP retA, SEXP retB)
{
    SEXP RESULT, EL, P; 
    double *tmp, *f, *wgt=REAL(weight), edle, ledle, newedle, eps=10, *eva=REAL(VECTOR_ELT(eig,0)); 
    double ll, lll, delta=0.0, scalep = 1.0, *ws=REAL(w), *gs=REAL(g), l1=0.0, l0=0.0;
    int i, k=0, ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0];
    tmp = (double *) R_alloc(nrx, sizeof(double));
    f = (double *) R_alloc(nrx, sizeof(double));
    PROTECT(RESULT = allocVector(VECSXP, 4));
    edle = REAL(el)[0];    
        
    for(i=0; i<nrx; i++)f[i] = REAL(f0)[i];
    NR66(eva, ncx, edle, ws, gs, X, INTEGER(ld)[0], nrx, f); // ncx-1L !!!
    for(i=0; i<nrx ;i++) l0 += wgt[i] * log(f[i]);    

    while ( (eps > 1e-05) &&  (k < 5) ) {
        if(scalep>0.6){  
            NR55(eva, ncx-1L, edle, ws, gs, X, INTEGER(ld)[0], nrx, f, tmp);  
            ll=0.0;  
            lll=0.0;        
            for(i=0; i<nrx ;i++) ll+=wgt[i]*tmp[i];
            for(i=0; i<nrx ;i++) lll+=wgt[i]*tmp[i]*tmp[i];  
            delta = ((ll/lll) < 3) ? (ll/lll) : 3;
        } // end if        
        ledle = log(edle) + scalep * delta;
        newedle = exp(ledle);
// some error handling avoid too big small edges & too big steps
        if (newedle > 10.0) newedle = 10.0;
        if (newedle < 1e-8) newedle = edle/2; 
        if (newedle < 1e-8) newedle = 1e-8; // 1e-8 phyML      
  
        for(i=0; i<nrx; i++)f[i] = REAL(f0)[i]; 
        NR66(eva, ncx, newedle, ws, gs, X, INTEGER(ld)[0], nrx, f);
        l1 = 0.0;
        for(i=0; i<nrx ;i++) l1 += wgt[i] * log(f[i]); // + log
        eps = l1 - l0;
// some error handling              
        if (eps < 0 || ISNAN(eps)) {
            if (ISNAN(eps))eps = 0;
            else {
                scalep = scalep/2.0;
                eps = 1.0;
            }
            newedle = edle;
            l1 = l0;
        }
        else scalep = 1.0;
        edle=newedle;
        l0 = l1; 
        k ++;
    }   
    PROTECT(EL = ScalarReal(edle));
    PROTECT(P = getPM(eig, nc, EL, g));  
    SET_VECTOR_ELT(RESULT, 0, EL); 
    if(INTEGER(retA)[0]>0L)SET_VECTOR_ELT(RESULT, 1, getM3(child, dad, P, nr, nc)); 
    if(INTEGER(retB)[0]>0L)SET_VECTOR_ELT(RESULT, 2, getM3(dad, child, P, nr, nc)); 
// add variance ??
    SET_VECTOR_ELT(RESULT, 3, ScalarReal(l1)); 
    UNPROTECT(3);
    return (RESULT);
} 


SEXP FS5(SEXP eig, SEXP nc, SEXP el, SEXP w, SEXP g, SEXP X, SEXP ld, SEXP nr, SEXP basefreq, SEXP weight, SEXP f0)
{
    SEXP RESULT; // EL, P; 
    double *tmp, *f, *wgt=REAL(weight), edle, ledle, newedle, eps=10, *eva=REAL(VECTOR_ELT(eig,0)); 
    double ll, lll, delta=0.0, scalep = 1.0, *ws=REAL(w), *gs=REAL(g), l1=0.0, l0=0.0;
    int i, k=0, ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0];
    tmp = (double *) R_alloc(nrx, sizeof(double));
    f = (double *) R_alloc(nrx, sizeof(double));
    PROTECT(RESULT = allocVector(REALSXP, 3));
    edle = REAL(el)[0];    
        
    for(i=0; i<nrx; i++)f[i] = REAL(f0)[i];
    NR66(eva, ncx, edle, ws, gs, X, INTEGER(ld)[0], nrx, f); // ncx-1L !!!
    for(i=0; i<nrx ;i++) l0 += wgt[i] * log(f[i]);    

    while ( (eps > 1e-05) &&  (k < 5) ) {
        if(scalep>0.6){  
            NR55(eva, ncx-1L, edle, ws, gs, X, INTEGER(ld)[0], nrx, f, tmp);  
            ll=0.0;  
            lll=0.0;        
            for(i=0; i<nrx ;i++) ll+=wgt[i]*tmp[i];
            for(i=0; i<nrx ;i++) lll+=wgt[i]*tmp[i]*tmp[i];  
            delta = ((ll/lll) < 3) ? (ll/lll) : 3;
        } // end if        
        ledle = log(edle) + scalep * delta;
        newedle = exp(ledle);
// some error handling avoid too big small edges & too big steps
        if (newedle > 10.0) newedle = 10.0;
        if (newedle < 1e-8) newedle = edle/2; 
        if (newedle < 1e-8) newedle = 1e-8; // 1e-8 phyML      
  
        for(i=0; i<nrx; i++)f[i] = REAL(f0)[i]; 
        NR66(eva, ncx, newedle, ws, gs, X, INTEGER(ld)[0], nrx, f);
        l1 = 0.0;
        for(i=0; i<nrx ;i++) l1 += wgt[i] * log(f[i]); // + log
        eps = l1 - l0;
// some error handling              
        if (eps < 0 || ISNAN(eps)) {
            if (ISNAN(eps))eps = 0;
            else {
                scalep = scalep/2.0;
                eps = 1.0;
            }
            newedle = edle;
            l1 = l0;
        }
        else scalep = 1.0;
        edle=newedle;
        l0 = l1; 
        k ++;
    }   
// variance 
    NR555(eva, ncx-1L, edle, ws, gs, X, INTEGER(ld)[0], nrx, f, tmp);  
    lll=0.0;        
    for(i=0; i<nrx ;i++) lll+=wgt[i]*tmp[i]*tmp[i]; 
    REAL(RESULT)[0] = edle;
    REAL(RESULT)[1] = 1.0/ lll; // variance
    REAL(RESULT)[2] = l0;
//    SET_VECTOR_ELT(RESULT, 3, ScalarReal(l1)); 
    UNPROTECT(1);
    return (RESULT);
} 



# define USE_RINTERNALS

# define USE_RINTERNALS

#include <Rmath.h>
#include <math.h>
#include <R.h> 
#include <R_ext/Lapack.h>
#include <Rinternals.h>


/*

static double *data1;
static double *weight;



void sankoff_init(int *n, double *weights, int *nr)
{
    int i;
    data1 = (double *) calloc(*n, sizeof(double));
    weight = (double *) calloc(*nr, sizeof(double));
    for(i=0; i<*nr; i++)weight[i] = weights[i]; 
}


void sankoff_free(){
    free(data1);
    free(weight);
}

static R_INLINE void sankoff4(double *dat, int n, double *cost, int k, double *result){
    int i, j, h; 
    double tmp[k], x;
    for(i = 0; i < n; i++){
        for(j = 0; j < k; j++){
            for(h = 0; h< k; h++){tmp[h] = dat[i + h*n] + cost[h + j*k];}
            x = tmp[0];
            for(h = 1; h< k; h++) {if(tmp[h]<x) {x=tmp[h];}}
            result[i+j*n] += x;
        }                   
    }        
}    


void sankoffTips(int *x, double *tmp, int nr, int nc, int nrs, double *result){
    int i, j;
    for(i = 0; i < (nr); i++){ 
        for(j = 0; j < (nc); j++) result[i + j*(nr)] += tmp[x[i] - 1L + j*(nrs)];  
    }
}


// mNodes raus
// SEXP fuer verschiedene IO mit R
// ohne zuviele SEXP definieren, Ausnahme dlist
void SANKOFF4(SEXP dlist, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP mNodes, SEXP tips, SEXP contrast, SEXP nrs){
    R_len_t i, n = length(node), nt = length(tips);
    int nTips = (int)length(tips), nii;
    int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0], nrc = INTEGER(nrs)[0];
    int  ni, ei, j, *edges=INTEGER(edge), *nodes=INTEGER(node);
    int rc = ncx*nrx;
    double *cost, *tmp; 
    tmp = (double *) R_alloc(ncx * nrc, sizeof(double)); 
    cost = REAL(scost);
    sankoff4(REAL(contrast), nrc, cost, ncx, tmp); 
    if(!isNewList(dlist)) error("'dlist' must be a list");
    ni = 0L;     
    for(i = 0; i < n; i++) {
        ei = edges[i]; 
        nii = (nodes[i]-nTips) * rc; 
        if(ni == nodes[i]){            
            if(ei < nt) sankoffTips(INTEGER(VECTOR_ELT(dlist,ei)), tmp, nrx, ncx, nrc, &data1[nii]);
            else sankoff4(&data1[(ei-1L)*rc], nrx, cost, ncx, &data1[nii]);
            }
        else{          
            ni = nodes[i];
            for(j=0; j<rc; j++) data1[nii+j] = 0.0; 
            if(ei < nt) sankoffTips(INTEGER(VECTOR_ELT(dlist,ei)), tmp, nrx, ncx, nrc, &data1[nii]);
            else sankoff4(&data1[(ei-1L)*rc], nrx, cost, ncx, &data1[nii]); 
            }
    }
}

*/



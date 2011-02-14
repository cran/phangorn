/* 
 * phangorn.c
 *
 * (c) 2011  Klaus Schliep (klaus.schliep@gmail.com)
 * 
 * 
 * This code may be distributed under the GNU GPL
 *
 */



# define USE_RINTERNALS


#include <stdlib.h>
#include <string.h>
#include <Rmath.h>
#include <math.h>
#include <R.h> 

#include <R_ext/Utils.h>
#include <R_ext/Memory.h> 
#include <R_ext/Print.h>  
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>



#define DINDEX(i, j) n*(i - 1) - i * (i - 1)/2 + j - i - 1

char *transa = "N", *transb = "N";
double one = 1.0, zero = 0.0;
int ONE = 1L;


int give_index(int i, int j, int n)
{
    if (i > j) return(DINDEX(j, i));
    else return(DINDEX(i, j));
}
 
 
void giveIndex(int *left, int* right, int *ll, int *lr, int *n, int *res){
    int i, j, k;
    k=0;
    for (i = 0; i < *ll; i++){
        for (j = 0; j < *lr; j++){
             res[k] = give_index(left[i], right[j], *n);
             k++;
             }    
        }
    }





int bitcount(int x){ 
    int count;
    for (count=0; x != 0; x>>=1)
       if ( x & 01)
           count++;
    return count;
}



void bitCount(int *x, int *count){
    count[0]=bitcount(x[0]);
} 




void fhm(double *v, int *n){
    unsigned int level, i, j; 
    unsigned int start, step, num_splits;
    double vi, vj;
    num_splits = (1 << (*n));
    step = 1;
    for(level = 0; level < (*n); level++){
        start = 0;
        while(start < (num_splits-1)){
            for(i = start; i < (start + step); i++){
                j = i + step;
                vi = v[i];
                vj = v[j];
                v[i] = vi + vj;
                v[j] = vi - vj;
            }
            start = start + 2*step;    
        }
        step *= 2;        
    }
}


void distance_hadamard(double *d, int n) {
    unsigned int num_splits;
    unsigned int x, r, nr, p, b, e;
    unsigned int odd = 1;                // The inner while loop can only terminate with odd == 1 so we don't need to set it inside the for loop.
    double cost, best_cost;
        
    num_splits = (1 << (n - 1));
        
    for (x = 1; x < num_splits; ++x) {
        r = (x - 1) & x;                // r = x without LSB
        nr = (r - 1) & r;                // nr = r without LSB
            
        if (nr) {                        // If x contains 1 or 2 bits only, then it has already been computed as a pairwise distance.
            b = x - r;                    // b = LSB of x: the "fixed" taxon in the pair.
            best_cost = 1e20;
            e = 0;                        // e holds bits to the right of the current p.
                
            while (1) {
                p = r - nr;                // p = 2nd half of pair
                cost = d[nr + e] + d[p + b];
                if (cost < best_cost) best_cost = cost;
                    
                if (!nr && odd) break;    // Ensure we get the (LSB with reference taxon) pair when there are an odd number of taxa
                r = nr;
                e += p;
                nr = (r - 1) & r;        // nr = r without LSB
                odd ^= 1;
                }
                d[x] = best_cost;
            }
        }

        d[0] = 0.0;
    }
    
void pairwise_distances(double *dm, int n, int num_splits, double *d) {
    int k=0;
    unsigned int offset;
        for (int i = 0; i < (n-1); ++i) {
            for (int j = (i+1); j < n; ++j) {
// Calculate the offset within the array to put the next value
                offset = (1 << i);
                if (j < n - 1) {            // If i == n - 1 then this is a distance between the reference taxon and some other taxon.
                    offset += (1 << j);        // Note that "+" is safe since (1 << i) and (1 << j) are always different bits.
                }
                d[offset]=dm[k];
                k++;
            }
        }
    }


SEXP dist2spectra(SEXP dm, SEXP nx, SEXP ns) {   
    int n = INTEGER(nx)[0];
    int nsp = INTEGER(ns)[0];   
    double *res;
    SEXP result;
    PROTECT(result = allocVector(REALSXP, nsp));
    res = REAL(result);
    pairwise_distances(REAL(dm), n, nsp, res);
    distance_hadamard(res, n);
    UNPROTECT(1);
    return(result);
}



// speed up some code for NJ    
void out(double *d, double *r, int *n, int *k, int *l){
    int i, j; 
    double res, tmp;
    k[0]=1;
    l[0]=2;
    res = d[1] - r[0] - r[1];
    for(i = 0; i < (*n-1); i++){
        for(j = i+1; j < (*n); j++){
                tmp = d[i*(*n)+j] - r[i] - r[j];
                if(tmp<res){
                    k[0]=i+1;
                    l[0]=j+1;
                    res = tmp;
                    }
            }
                
        }        
    }


    
// hamming distance    
void distHamming(int *x, double *weight, int *nr, int *l, double *d){
    int i, j, k, m;
    k = 0L;
    for(i = 0; i< (*l-1L); i++){
        for(j = (i+1L); j < (*l); j++){
             for(m=0; m<(*nr); m++){
                 if(!(x[i*(*nr) + m] & x[j*(*nr) + m])) d[k] += weight[m]; 
                 } 
             k++;
        }
    }
}

    
    
     
SEXP rowMin(SEXP sdat, SEXP sn, SEXP sk){
    int i, h, n=INTEGER(sn)[0], k=INTEGER(sk)[0];  
    double x, *res, *dat;
    SEXP result;
    PROTECT(result = allocVector(REALSXP, n));
    res = REAL(result);
    PROTECT(sdat = coerceVector(sdat, REALSXP));
    dat = REAL(sdat);
    for(i = 0; i < n; i++){
        x = dat[i];
        for(h = 1; h< k; h++) {if(dat[i + h*n] < x) x=dat[i + h*n];}
        res[i] = x;               
        }
    UNPROTECT(2);
    return(result);        
}



void rowMin2(double *dat, int n,  int k, double *res){
    int i, h;  
    double x;
    for(i = 0; i < n; i++){
        x = dat[i];
        for(h = 1; h< k; h++) {if(dat[i + h*n] < x) x=dat[i + h*n];}
        res[i] = x;               
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
        

void fitch43(int *dat1, int *dat2, int *nr, int *pars, double *weight, double *w){
    int k;
    int tmp;
    for(k = 0; k < (*nr); k++){
        tmp = dat1[k] & dat2[k];
        if(tmp) dat1[k] = tmp;
        else {dat1[k] = dat1[k] | dat2[k];
            pars[k] += 1;
            (*w)+=weight[k];
        }
    } 
}


void fitch6(int *dat, int *nr, int *pars, int *node, int *edge, int *nl, double *weight, double *pvec, double *pscore) 
{   
    int i, ni, k;
    ni = 0;
    for (i=0; i< *nl; i++) {
	    if (ni == node[i]){
	         pvec[ni-1] += pvec[edge[i]-1];
	         fitch43(&dat[(ni-1) * (*nr)], &dat[(edge[i]-1L) * (*nr)], nr, pars, weight, &pvec[(ni-1L)]); //pvec[(ni-1L)]
        }                  
        else {
            ni = node[i];   
            pvec[(ni-1L)] += pvec[(edge[i]-1L)];  
            // memcpy       
            for(k = 0; k < (*nr); k++) dat[(ni-1)*(*nr) + k] = dat[(edge[i]-1)*(*nr) + k];                     
        }
    }
    pscore[0]=pvec[ni-1];
}


/*
ACCTRAN
*/

void fitchT(int *dat1, int *dat2, int *nr, double *pars, double *weight, double *w){
    int k;
    int tmp;
    for(k = 0; k < (*nr); k++){
        tmp = dat1[k] & dat2[k];
        if(tmp > 0L){
             dat1[k] = tmp;
             }
    } 
}


void fitchT3(int *dat1, int *dat2, int *nr, double *pars, double *weight, double *w){
    int k;
    int tmp;
    for(k = 0; k < (*nr); k++){
       tmp = dat1[k] & dat2[k];
       if(tmp==0L) {
             (*w)+=weight[k];
             pars[k] += 1;
             }
       if(tmp >0){
           if(tmp < dat2[k]){ 
              (*w)+= .5*weight[k];
              pars[k] += .5;
           }
       }

    } 
}


// return lower and upper bound for the number of changes 
// upper bound very conservative 
void countMPR(double *res, int *dat1, int *dat2, int *nr, double *weight, int *external){
    int k;
    int tmp, tmp2;
    for(k = 0; k < (*nr); k++){
        tmp = dat1[k] & dat2[k];

        if(tmp==0){
            res[0] += weight[k];
            res[1] += weight[k];
        }
        else{ 
            if( external[0]==0L){ 
                 if( bitcount(dat1[k] | dat2[k])>1L ) res[1] += weight[k]; // dat1[k] != dat2[k]
            }   
            else{ 
                 if( tmp  < dat2[k] ) res[1] += weight[k];
            }
        }
    } 
}


void ACCTRAN2(int *dat, int *nr, double *pars, int *node, int *edge, int *nl, double *weight, double *pvec, int *nTips) 
{   
    int i;
    for (i=0; i< *nl; i++) {       
        if(edge[i]>nTips[0]) fitchT(&dat[(edge[i]-1L) * (*nr)], &dat[(node[i]-1) * (*nr)], nr, pars, weight, &pvec[i]); 
        }
}


void ACCTRAN3(int *dat, int *nr, double *pars, int *node, int *edge, int *nl, double *weight, double *pvec, int *nTips) 
{   
    int i;
    for (i=0; i< *nr; i++)pars[i]=0.0;
    for(i=0; i< *nl; i++)pvec[i] = 0.0;
    for (i=0; i< *nl; i++) {               
        fitchT3(&dat[(edge[i]-1L) * (*nr)], &dat[(node[i]-1) * (*nr)], nr, pars, weight, &pvec[i]); 
    }            
}


void fitchN(int *dat1, int *dat2, int *nr){
    int k;
    int tmp;
    for(k = 0; k < (*nr); k++){
        tmp = dat1[k] & dat2[k];
        if(tmp) dat1[k] = tmp;
        else dat1[k] = dat1[k] | dat2[k];
        } 
}


void fitchN2(int *res, int *dat, int *node, int *edge, int *nr, int *nl) { 
    int i;
    for (i=0; i< *nl; i++) {
        fitchN(&res[(node[i]-1L) * (*nr)], &dat[(edge[i]-1L) * (*nr)], nr);              
    }
}


void fitchTriplet(int *res, int *dat1, int *dat2, int *dat3, int *nr) 
{   
    int i, ni, k;
    ni = 0;
    
    int *v1, *v2, *v3;
    v1 = (int *) R_alloc(*nr, sizeof(int));    
    v2 = (int *) R_alloc(*nr, sizeof(int));
    v3 = (int *) R_alloc(*nr, sizeof(int));

    for(k = 0; k < (*nr); k++) v1[k] = dat1[k];
    fitchN(v1, dat2, nr);
    fitchN(v1, dat3, nr);

    for(k = 0; k < (*nr); k++) v2[k] = dat1[k];
    fitchN(v2, dat3, nr);
    fitchN(v2, dat2, nr);

    for(k = 0; k < (*nr); k++) v3[k] = dat2[k];
    fitchN(v3, dat3, nr);
    fitchN(v3, dat1, nr);

    for(k = 0; k < (*nr); k++)res[k] = v1[k] & v2[k]; // &v3[k];  
    for(k = 0; k < (*nr); k++)res[k] = res[k] & v3[k];  
}


void fitchNACC2(int *root, int *dat, int *nr, double *pars, int *result, double *weight, double *pars1){
    int k;
    int tmp;
    for(k = 0; k < (*nr); k++){
//       result[k] = 0L;
       tmp = root[k] & dat[k];
       if(tmp==0L) {
             pars[0] += weight[k];
             pars1[k] += weight[k];
             }
       if(tmp >0){
           if(tmp < root[k]){ 
              pars[0] += .5*weight[k];
              pars1[k] += .5*weight[k];
              result[k] += 1L;
           }
       }
    }        
}


void fitchTripletACC4(int *root, int *dat1, int *dat2, int *dat3, int *nr, double *p1, double *p2, double *p3, double *weight, double *pars1, int *v1) 
{   
    int k;
       
    int tmp, a, b, c, t1, t2, t3;
    double d, f;
    for(k = 0; k < (*nr); k++){
        tmp = root[k];
        a = dat1[k] & dat2[k]; 
        b = dat1[k] & dat3[k];
        c = dat2[k] & dat3[k];
        if((a+b+c) == 0L){
           d = (2.0/3.0) * weight[k];
           p1[0] += d; 
           p2[0] += d;
           p3[0] += d;        
           pars1[k] += 2*weight[k]; 
           v1[k] = 2L; 
        }
        else{  
            f = 0.0;
            d = weight[k];
            t1 = 0.0;
            t2 = 0.0;
            t3 = 0.0;
            if( (dat1[k] & tmp)<tmp){ 
                t1 = d; 
                f+=1.0;
            }
            if( (dat2[k] & tmp)<tmp){ 
                t2 = d; 
                f+=1.0;
            }
            if( (dat3[k] & tmp)<tmp){ 
                t3 = d; 
                f+=1.0;
            }
            if(f>0.0){   
                pars1[k] += weight[k]; 
                p1[0] += t1/f; 
                p2[0] += t2/f;
                p3[0] += t3/f;
                v1[k] += 1L;
            }
        }
    }
}




/*
ACCTRAN

result <- .C("ACCTRAN", dat, as.integer(nr), as.integer(node), as.integer(edge), as.integer(length(edge)), as.double(weight), numeric(m))

result <- .C("fitch6", dat, as.integer(nr), as.integer(pars), as.integer(node), as.integer(edge), as.integer(length(edge)), as.double(weight), 
   numeric(m), as.double(0))
.Call("FITCH", data[, tree$tip.label], as.integer(nr), as.integer(node), as.integer(edge), as.integer(length(edge)), as.double(weight), as.integer(m), 
   as.integer(q))
*/


SEXP FITCH(SEXP dat, SEXP nrx, SEXP node, SEXP edge, SEXP l, SEXP weight, SEXP mx, SEXP q){   
    int *data, *nr=INTEGER(nrx), m=INTEGER(mx)[0], i, n=INTEGER(q)[0];   //=INTEGER(dat)
    double *pvtmp;  
    SEXP DAT, pars, pvec, pscore, RESULT;
    PROTECT(RESULT = allocVector(VECSXP, 4L));
    PROTECT(pars = allocVector(INTSXP, *nr));
    PROTECT(pscore = allocVector(REALSXP, 1L));
    PROTECT(DAT = allocMatrix(INTSXP, nr[0], m));
    PROTECT(pvec = allocVector(REALSXP, m));
    pvtmp = REAL(pvec);
    data = INTEGER(DAT);
    for(i=0; i<m; i++) pvtmp[i] = 0.0;
    for(i=0; i<*nr; i++) INTEGER(pars)[i] = 0L;
    REAL(pscore)[0]=0.0;
//    for(i=0; i<(*nr * m); i++)data[i] = 0L;
    for(i=0; i<(*nr * n); i++)data[i] = INTEGER(dat)[i];
    
    fitch6(data, nr, INTEGER(pars), INTEGER(node), INTEGER(edge), INTEGER(l), REAL(weight), pvtmp, REAL(pscore));
    
    SET_VECTOR_ELT(RESULT, 0, pscore);
    SET_VECTOR_ELT(RESULT, 1, pars);
    SET_VECTOR_ELT(RESULT, 2, DAT);
    SET_VECTOR_ELT(RESULT, 3, pvec);
    UNPROTECT(5);
    return(RESULT); 
}


// eine else Schliefe raus
void FN2(int *dat, int *res, int *nr, int *pars, int *node, int *edge, int *nl, int *pc, double *weight, double *tmpvec, double *pvec) { 
    int i, ni, k;
    ni = 0;
    for (i=0; i< *nl; i++) {
        if (ni == node[i]){
              if(pc[i]==0L){
	              pvec[(ni-1L)] += tmpvec[(edge[i]-1L)];
	              fitch43(&res[(ni-1L) * (*nr)], &dat[(edge[i]-1L) * (*nr)], nr, pars, weight, &pvec[(ni-1L)]);              
              }    
              else{ 
                      pvec[(ni-1L)] += pvec[(edge[i]-1L)];
	              fitch43(&res[(ni-1L) * (*nr)], &res[(edge[i]-1L) * (*nr)], nr, pars, weight, &pvec[(ni-1L)]);   //&pvec[(ni-1)] 
              }
        }                      
        else {
            ni = node[i];
            pvec[(ni-1L)] += tmpvec[(edge[i]-1L)];
            for(k = 0; k < (*nr); k++) res[(ni-1L)*(*nr) + k] = dat[(edge[i]-1L)*(*nr) + k];                 
        }
    }
}


SEXP FNALL(SEXP dat, SEXP nrx, SEXP node, SEXP edge,  SEXP node2, SEXP edge2, SEXP l, SEXP weight, SEXP mx, SEXP my, SEXP q, SEXP pc){   
    int *data=INTEGER(dat), *nr=INTEGER(nrx), m=INTEGER(mx)[0], i, n=INTEGER(q)[0], *pars2, *dat2;  
    double *pvtmp, *pvtmp2;  
    SEXP DAT, DAT2, pars, pvec, pvec2, pscore, RESULT;
  
    pars2 = (int *) R_alloc(*nr, sizeof(int));    

    PROTECT(RESULT = allocVector(VECSXP, 6L));
    PROTECT(pars = allocVector(INTSXP, *nr));
    PROTECT(pscore = allocVector(REALSXP, 1L));
    PROTECT(DAT = allocMatrix(INTSXP, nr[0], m));
    PROTECT(DAT2 = allocMatrix(INTSXP, nr[0], m));
    PROTECT(pvec = allocVector(REALSXP, m));
    PROTECT(pvec2 = allocVector(REALSXP, m)); 
    for(i=0; i<*nr; i++) INTEGER(pars)[i] = 0L;
    pvtmp = REAL(pvec);
    for(i=0; i<m; i++) pvtmp[i] = 0.0;
    pvtmp2 = REAL(pvec2);
    for(i=0; i<m; i++) pvtmp2[i] = 0.0;
    REAL(pscore)[0]=0.0;
    for(i=0; i<(*nr * n); i++)INTEGER(DAT)[i] = data[i];
    dat2 = INTEGER(DAT2); // without wired things happen
    for(i=0; i<(*nr * m); i++)dat2[i] = 0L;
    fitch6(INTEGER(DAT), nr, INTEGER(pars), INTEGER(node), INTEGER(edge), INTEGER(l), REAL(weight), pvtmp, REAL(pscore));   
    FN2(INTEGER(DAT), dat2, nr, pars2, INTEGER(node2), INTEGER(edge2), INTEGER(my), INTEGER(pc), REAL(weight), pvtmp, pvtmp2);

    SET_VECTOR_ELT(RESULT, 0, pscore);
    SET_VECTOR_ELT(RESULT, 1, pars);
    SET_VECTOR_ELT(RESULT, 2, DAT);
    SET_VECTOR_ELT(RESULT, 3, DAT2);
    SET_VECTOR_ELT(RESULT, 4, pvec);
    SET_VECTOR_ELT(RESULT, 5, pvec2);
    UNPROTECT(7);
    return(RESULT); 
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



SEXP sankoffQuartet(SEXP dat, SEXP sn, SEXP scost, SEXP sk){
    int j, n=INTEGER(sn)[0], k = INTEGER(sk)[0];  
    double *cost, *res, *rtmp;
    SEXP result;
    PROTECT(result = allocVector(REALSXP, n));
    rtmp = (double *) R_alloc(n*k, sizeof(double));
    res = (double *) R_alloc(n*k, sizeof(double));
    PROTECT(scost = coerceVector(scost, REALSXP));
    cost = REAL(scost);
    for(j=0; j<(n*k); j++) rtmp[j] = 0.0;
    for(j=0; j<(n*k); j++) res[j] = 0.0;   
    sankoff4(REAL(VECTOR_ELT(dat,0)), n, cost, k, rtmp);
    sankoff4(REAL(VECTOR_ELT(dat,1)), n, cost, k, rtmp);
    sankoff4(rtmp, n, cost, k, res);
    sankoff4(REAL(VECTOR_ELT(dat,2)), n, cost, k, res);
    sankoff4(REAL(VECTOR_ELT(dat,3)), n, cost, k, res);
    rowMin2(res, n, k, REAL(result));  //res, sn sk  
    UNPROTECT(2);    
    return(result);        
}    


SEXP sankoffTwin(SEXP dat1, SEXP dat2, SEXP sn, SEXP scost, SEXP sk){
    int j, n=INTEGER(sn)[0], k = INTEGER(sk)[0];  
    double *cost, *tmp;
    SEXP result;
    PROTECT(result = allocMatrix(REALSXP, n, k));
    tmp = REAL(result);
    PROTECT(scost = coerceVector(scost, REALSXP));
    cost = REAL(scost);
    for(j=0; j<(n*k); j++) tmp[j] = REAL(dat1)[j];
//    sankoff4(REAL(dat1), n, cost, k, tmp);
    sankoff4(REAL(dat2), n, cost, k, tmp);
    UNPROTECT(2);    
    return(result);        
} 


SEXP sankoff3(SEXP dlist, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP mNodes, SEXP tips){
    R_len_t i, n = length(node), nt = length(tips);
    int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0], mn=INTEGER(mNodes)[0];
    int  ni, ei, j, *edges=INTEGER(edge), *nodes=INTEGER(node);
    SEXP result, dlist2; //tmp, 
    double *res, *rtmp, *cost;
    cost = REAL(scost);
    if(!isNewList(dlist)) error("'dlist' must be a list");
    ni = nodes[0];
    PROTECT(dlist2 = allocVector(VECSXP, mn));
    PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
    res = REAL(result);
    rtmp = (double *) R_alloc(nrx*ncx, sizeof(double));
    for(i = 0; i < nt; i++) SET_VECTOR_ELT(dlist2, INTEGER(tips)[i], VECTOR_ELT(dlist, INTEGER(tips)[i]));
    for(j=0; j<(nrx * ncx); j++) res[j] = 0.0; 
 
    for(i = 0; i < n; i++) {
        ei = edges[i]; 
        if(ni == nodes[i]){
            sankoff4(REAL(VECTOR_ELT(dlist2,ei)), nrx, cost, ncx, res);
            }
        else{          
            SET_VECTOR_ELT(dlist2, ni, result);
            UNPROTECT(1); 
            PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
            res = REAL(result);
            for(j=0; j<(nrx * ncx); j++) res[j] = 0.0; 
            ni = nodes[i];
            sankoff4(REAL(VECTOR_ELT(dlist2,ei)), nrx, cost, ncx, res); 
            }
    }
    SET_VECTOR_ELT(dlist2, ni, result);    
//    rowMin2(REAL(result), nrx, ncx, );
//    rowMin(SEXP sdat, SEXP sn, SEXP sk)
    UNPROTECT(2); // result dlist2  return dlist2
    return(dlist2);
}


    
SEXP pNodes(SEXP data, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge){
    R_len_t n = length(node); //, nt = length(tips);, SEXP tips , SEXP mNodes , mn=INTEGER(mNodes)[0]
    int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0];
    int k, pj, i, j, start, *edges=INTEGER(edge), *nodes=INTEGER(node);
    SEXP result, dlist;  
    double *res, *tmp, *cost;
    cost = REAL(scost);
    pj = nodes[n-1L];
    start = n-1L;
    PROTECT(dlist = allocVector(VECSXP, length(data)));
    tmp = (double *) R_alloc(nrx*ncx, sizeof(double));    
    for(i=0; i<(nrx * ncx); i++) tmp[i] = 0.0;
    for(j=n-1L; j>=0; j--) {
        PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
        res = REAL(result);
        if (pj != nodes[j]) {
            for(i=0; i<(nrx * ncx); i++) tmp[i] = 0.0;
            sankoff4(REAL(VECTOR_ELT(dlist, nodes[j])), nrx, cost, ncx, tmp);
            for(i=0; i<(nrx * ncx); i++) res[i] = tmp[i] ;  
            pj = nodes[j];
            start = j;
        }
        else for(i=0; i<(nrx * ncx); i++) res[i] = tmp[i] ;
        k = start;
        while (k >= 0 && pj == nodes[k]) {
            if (k != j) 
                sankoff4(REAL(VECTOR_ELT(data, edges[k])), nrx, cost, ncx, res); 
                
            k--;
        }
        SET_VECTOR_ELT(dlist, edges[j], result);    
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return(dlist);
}


static R_INLINE void matprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z)
{
    char *transa = "N", *transb = "N";
    double one = 1.0, zero = 0.0;
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

/*
matprod = function(X, Y){
    k=dim(X)
    l=dim(Y)
    .Call("matpro",X,Y,k[1],k[2],l[1],l[2])
}
*/

 
static R_INLINE void getP(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double tmp[m], res;
    for(i = 0; i < m; i++) tmp[i] = exp(eva[i] * w * el);
    for(i = 0; i < m; i++){    
        for(j = 0; j < m; j++){
            res = 0.0;    
            for(h = 0; h < m; h++)    res += ev[i + h*m] * tmp[h] * evi[h + j*m];
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
            for(h = 0; h < m; h++)    result[i+j*m] += ev[i + h*m] * tmp[h] * evi[h + j*m];                
        }
    }
}


void getPmix(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double tmp[m], res[m*m];
    for(i = 0; i < m; i++) tmp[i] = exp(eva[i] * w * el);
    for(j = 0; j < m; j++){
        for(h = 0; h < m; h++)    res[h + j*m] = evi[h + j*m] * tmp[h];
    }
    matprod(ev, m, m, res, m, m, result);
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




void NR5(SEXP eig, int nc, double el, double *w, double *g, SEXP dad, SEXP child, int ld, int nr, double *bf, double *f, double *res){
    int i, j, k, eins=1L; 
    double *kid;  
    double *eva, *eve, *evei; 
    double *dP, *dtmp; //*res,  *dF,
//    char *transa = "N", *transb = "N"; //raus
//    double one = 1.0, zero = 0.0; //raus
    if(!isNewList(eig)) error("eig must be a list"); //nach oben
    if(!isNewList(dad)) error("dad must be a list"); //nach oben
    if(!isNewList(child)) error("child must be a list"); //nach oben   
    dP = (double *) R_alloc(nc*nc, sizeof(double));
    dtmp = (double *) R_alloc(nr*nc, sizeof(double));
        
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));

    for(k=0; k<nr; k++) res[k] = 0.0;
        
    for(j=0;j<ld;j++){
        kid = REAL(VECTOR_ELT(child, j)); // kid
        getdP(eva, eve, evei, nc, el, g[j], dP); //dP faellt weg
        F77_CALL(dgemm)(transa, transb, &nr, &nc, &nc, &one, REAL(VECTOR_ELT(dad, j)), &nr, dP, &nc, &zero, dtmp, &nr);
        for(i = 0; i < (nr*nc); i++) dtmp[i]*=kid[i]; //dtmp =  (dad %*% dP) * kid 
//  dgemv
        F77_CALL(dgemv)(transa, &nr, &nc, &w[j], dtmp, &nr, bf, &ONE, &one, res, &ONE); 
//        F77_CALL(dgemm)(transa, transb, &nr, &eins, &nc, &w[j], dtmp, &nr, bf, &nc, &one, res, &nr);
        }
    for(i=0; i<nr ;i++) res[i]/=f[i];                
} 


void NR6(SEXP eig, int nc, double el, double *w, double *g, SEXP dad, SEXP child, int ld, int nr, double *bf, double *res){
    int i, eins=1L, j; 
    double *kid;  
    double *eva, *eve, *evei; 
    double *dP, *dtmp; //*res,  *dF,
    if(!isNewList(eig)) error("'eig' must be a list");
    if(!isNewList(dad)) error("'dad' must be a list");
    if(!isNewList(child)) error("'child' must be a list");    
    dP = (double *) R_alloc(nc*nc, sizeof(double));
    dtmp = (double *) R_alloc(nr*nc, sizeof(double));
        
    eva = REAL(VECTOR_ELT(eig, 0)); // raus
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));

    for(j=0;j<ld;j++){
        kid = REAL(VECTOR_ELT(child, j)); // kid
        getP(eva, eve, evei, nc, el, g[j], dP); 
        F77_CALL(dgemm)(transa, transb, &nr, &nc, &nc, &one, REAL(VECTOR_ELT(dad, j)), &nr, dP, &nc, &zero, dtmp, &nr);
        for(i = 0; i < (nr*nc); i++) dtmp[i]*=kid[i]; //dtmp =  (dad %*% dP) * kid //BLAS???
// dgemv
        F77_CALL(dgemv)(transa, &nr, &nc, &w[j], dtmp, &nr, bf, &ONE, &one, res, &ONE); 
//        F77_CALL(dgemm)(transa, transb, &nr, &eins, &nc, &w[j], dtmp, &nr, bf, &nc, &one, res, &nr);
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
//    double *p;    
    for(j=0; j<nel; j++){ 
        for(i=0; i<nw; i++){
            PROTECT(P = allocMatrix(REALSXP, m, m));
//            p = REAL(P);
            if(edgelen[j]==0.0){
                for(k=0; k<(m*m);k++)REAL(P)[k]=0.0;
                for(k=0; k<m; k++)REAL(P)[k+k*m]=1.0;
                }
            else getP(eva, eve, evei, m, edgelen[j], ws[i], REAL(P));//p
            SET_VECTOR_ELT(RESULT, l, P);
            UNPROTECT(1); //P
            l++;
        }
    }
    UNPROTECT(1);//RESULT
    return(RESULT);
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
            getPP(eva, eve, evei, m, edgelen[j], ws[i], REAL(P));//p
            SET_VECTOR_ELT(RESULT, l, P);
            UNPROTECT(1); //P
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


// return variance estimate
SEXP FS3(SEXP eig, SEXP nc, SEXP el, SEXP w, SEXP g, SEXP dad, SEXP child, SEXP ld, SEXP nr, 
         SEXP basefreq, SEXP weight,SEXP f0, SEXP ll0, SEXP ff0, SEXP retA, SEXP retB)
{
    SEXP RESULT, EL, P; //, logLik; 
    double *tmp, *tmp2, *f, *wgt=REAL(weight), edle, ledle, newedle, eps=10, *bf=REAL(basefreq); //*df, 
    double ll, lll, delta=0.0, scalep = 1.0, *ws=REAL(w), *gs=REAL(g), l1=0.0, l0;
    int i, k=0, ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0];
    tmp = (double *) R_alloc(nrx, sizeof(double));
    tmp2 = (double *) R_alloc(nrx, sizeof(double));
    f = (double *) R_alloc(nrx, sizeof(double));
    PROTECT(RESULT = allocVector(VECSXP, 4));
    edle = REAL(el)[0];    
    l0 = REAL(ll0)[0];
    for(i=0; i<nrx ;i++) f[i]=REAL(ff0)[i];
// eps uebergeben
    while ( (eps > 1e-05) &&  (k < 5) ) {
        if(scalep>0.6){  
// NR5, NR6 zusammenfassen: P uebergeben / f nicht  eig raus ziehen und alle variablen die oefters definiert werden
// getdP - tmp 
// NR5/6 explicit einbauen 
            NR5(eig, ncx, edle, ws, gs, dad, child, INTEGER(ld)[0], nrx, bf, f, tmp);  
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
        if (newedle < 1e-6) newedle = edle/2; // new 
        if (newedle < 1e-6) newedle = 1e-6; // 1e-8 phyML      
  
        for(i=0; i<nrx; i++)f[i] = REAL(f0)[i]; //memcpy ???
// getP         
        NR6(eig, ncx, newedle, ws, gs, dad, child, INTEGER(ld)[0], nrx, bf, f);
        l1 = 0.0;
        for(i=0; i<nrx ;i++) l1 += wgt[i] * log(f[i]);
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
// add variance
    SET_VECTOR_ELT(RESULT, 3, ScalarReal(l1));
    UNPROTECT(3);
    return (RESULT);
} 


void reorder(int *from, int *to, int *n, int *sumNode,  int *neworder, int *root){ 
    int i, j, sum=0, k, l, Nnode, ind, *ord, *csum, *tips, *stack, z=0;  
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
    R_tabulate(from, n, sumNode, tips);
    csum[0]=0;
    for(i=0;i<(*sumNode);i++){
        sum+=tips[i];                 
        csum[i+1] = sum;                        
    }      
    k = (*n)-1;
    l = (*n)-1;
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

/*
static R_INLINE void matprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z)
{
    char *transa = "N", *transb = "N";
    double one = 1.0, zero = 0.0;
    F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one, x, &nrx, y, &nry, &zero, z, &nrx);
}
*/


// tmp extern erstellen??
void matp(int *x, double *contrast, double *P, int *nr, int *nc, int *nrs, double *result){
    int i, j;
    double *tmp; // vielleicht rausziehen etwas weniger Speicher 
    tmp = (double *) R_alloc((*nc) *(*nrs), sizeof(double)); 
//    F77_CALL(dgemm)(transa, transb, nrs, nc, nc, &one, contrast, nc, P, nc, &zero, tmp, nrs);
    matprod(contrast, (*nrs), (*nc), P, (*nc), (*nc), tmp);   
    for(i = 0; i < (*nr); i++){ 
        for(j = 0; j < (*nc); j++) result[i + j*(*nr)] = tmp[x[i] - 1L + j*(*nrs)];  
    }
}


// nco: nr of contrasts (ca. 5% schneller), durch LogLik4 ersetzen
// void LogLik
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

    char *transa = "N", *transb = "N"; // global??
    double one = 1.0, zero = 0.0;

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


void matm(int *x, double *contrast, int *nr, int *nc, int *nco, double *result){
    int i, j;
    for(i = 0; i < (*nr); i++){ 
        for(j = 0; j < (*nc); j++) result[i + j*(*nr)] *= contrast[x[i] - 1L + j*(*nco)];  
    }
}

// sparse Matrix
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


void scaleMatrix(double *X, int nr, int nc, double *result){
    int i, j; 
    double tmp;   
    for(i = 0; i < nr; i++) {    
        tmp = 0.0; 
        for(j = 0; j < nc; j++) {
           tmp += X[i + j*nr];
        }        
        for(j = 0; j < nc; j++) {
           X[i + j*nr] /= tmp;
        }
        result[i] += log(tmp);
    } 
}

void scaleMatrix2(double *X, int nr, int nc, double *result, double eps){
    int i, j; 
    double tmp;   
    for(i = 0; i < nr; i++) {    
        tmp = 0.0; 
        for(j = 0; j < nc; j++) {
           tmp += X[i + j*nr];
        }
        if(tmp < eps){          
            for(j = 0; j < nc; j++) {
                X[i + j*nr] /= tmp;
            }
            result[i] += log(tmp);
        }
    } 
}

void scaleMatrix3(double *X, int nr, int nc, double *weight, double *result, double erg, double eps){
    int i, j; 
    double tmp, ltmp;   
    for(i = 0; i < nr; i++) {    
        tmp = 0.0; 
        for(j = 0; j < nc; j++) {
           tmp += X[i + j*nr];
        }
        if(tmp < eps){          
            for(j = 0; j < nc; j++) {
                X[i + j*nr] /= tmp;
            }
            ltmp = log(tmp); 
            result[i] += ltmp;
            erg += (ltmp * weight[i]);
        }
    } 
}



SEXP LogLik4(SEXP dlist, SEXP P, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP nTips, SEXP mNodes, SEXP contrast, SEXP nco){
    R_len_t i, n = length(node);
    int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0], nt=INTEGER(nTips)[0], mn=INTEGER(mNodes)[0];
    int  ni, ei, j, *edges=INTEGER(edge), *nodes=INTEGER(node);
    SEXP erg, ans, result, scVec , scTmp;  
    double *res, *rtmp, *scaleTmp;
    if(!isNewList(dlist)) error("'dlist' must be a list");
    ni = nodes[0];
    PROTECT(erg = allocVector(VECSXP, 2L));
    PROTECT(ans = allocVector(VECSXP, mn)); 
    PROTECT(scVec = allocVector(VECSXP, mn));
    PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
    PROTECT(scTmp = allocVector(REALSXP, nrx));
    res = REAL(result);
    rtmp = (double *) R_alloc(nrx*ncx, sizeof(double));
    scaleTmp = REAL(scTmp);
// global
//    char *transa = "N", *transb = "N";
//    double one = 1.0, zero = 0.0;
    for(j=0; j < nrx; j++) scaleTmp[j] = 0.0;
    for(j=0; j < (nrx * ncx); j++) res[j] = 1.0;
    for(i = 0; i < n; i++) {
        ei = edges[i]; 
        if(ni != nodes[i]){
//          no scalMatrix here??
            scaleMatrix(res, nrx, ncx, scaleTmp);
            SET_VECTOR_ELT(ans, ni, result);
            SET_VECTOR_ELT(scVec, ni, scTmp);
            UNPROTECT(2); //result scTmp
            PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
            res = REAL(result);
            PROTECT(scTmp = allocVector(REALSXP, nrx));
            scaleTmp = REAL(scTmp);
            for(j=0; j < nrx; j++) scaleTmp[j] = 0.0;
            ni = nodes[i];
            if(ei < nt) 
                matp(INTEGER(VECTOR_ELT(dlist, ei)), REAL(contrast), REAL(VECTOR_ELT(P, i)), INTEGER(nr), INTEGER(nc), INTEGER(nco), res); // nr, nrx, INTEGER(nr)
            else{ 
                F77_CALL(dgemm)(transa, transb, &nrx, &ncx, &ncx, &one, REAL(VECTOR_ELT(ans, ei-nt)), &nrx, 
                       REAL(VECTOR_ELT(P, i)), &ncx, &zero, res, &nrx);
                for(j=0; j < nrx; j++) scaleTmp[j] += REAL(VECTOR_ELT(scVec, ei-nt))[j];
            }     
//          z=1L;      
        }
        else {
            if(ei < nt) 
                matp(INTEGER(VECTOR_ELT(dlist, ei)), REAL(contrast), REAL(VECTOR_ELT(P, i)), INTEGER(nr), INTEGER(nc), INTEGER(nco), rtmp);
            else{ 
                F77_CALL(dgemm)(transa, transb, &nrx, &ncx, &ncx, &one, REAL(VECTOR_ELT(ans, ei-nt)), &nrx, 
                    REAL(VECTOR_ELT(P, i)), &ncx, &zero, rtmp, &nrx);
                for(j=0; j < nrx; j++) scaleTmp[j] += REAL(VECTOR_ELT(scVec, ei-nt))[j];
            }
            for(j=0; j < (nrx*ncx); j++) res[j] *= rtmp[j];
//          z+=1L; // epsScale=50; statt+1 von vector
//          if(z>epsScale){
//              scaleMatrix(res, nrx, ncx, scaleTmp);   
//              z=1;     
//          }
        }            
    }
    scaleMatrix(res, nrx, ncx, scaleTmp);
    SET_VECTOR_ELT(ans, ni, result);
    SET_VECTOR_ELT(scVec, ni, scTmp);
    SET_VECTOR_ELT(erg, 0L, ans);
    SET_VECTOR_ELT(erg, 1L, scVec);
    UNPROTECT(5); // result ans 
    return(erg);
}


SEXP LogLik5(SEXP dlist, SEXP P, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP nTips, SEXP mNodes, SEXP contrast, SEXP nco, SEXP weight){
    R_len_t i, n = length(node);
    int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0], nt=INTEGER(nTips)[0], mn=INTEGER(mNodes)[0];
    int  ni, ei, j, *edges=INTEGER(edge), *nodes=INTEGER(node), *nChild, eps1;
    SEXP erg, ans, result, scVec , scTmp;  
    double *res, *rtmp, *scaleTmp, *tmpS, eps2;
    
    eps1 = 50L;
    eps2 = 1e-10;

    if(!isNewList(dlist)) error("'dlist' must be a list");
    ni = nodes[0];
    PROTECT(erg = allocVector(VECSXP, 2L));
    PROTECT(ans = allocVector(VECSXP, mn)); 
    PROTECT(scVec = allocVector(VECSXP, mn));
    PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
    PROTECT(scTmp = allocVector(REALSXP, nrx));
    res = REAL(result);
    tmpS = (double *) R_alloc(n, sizeof(double));
    nChild = (int *) R_alloc(n, sizeof(int));
    rtmp = (double *) R_alloc(nrx*ncx, sizeof(double));
    scaleTmp = REAL(scTmp);
    for(j=0; j < nrx; j++) scaleTmp[j] = 0.0;
    for(j=0; j < (nrx * ncx); j++) res[j] = 1.0;
    for(i = 0; i < n; i++) {
        ei = edges[i]; 
        if(ni != nodes[i]){
            scaleMatrix2(res, nrx, ncx, scaleTmp, 1e-10);
            SET_VECTOR_ELT(ans, ni, result);
            SET_VECTOR_ELT(scVec, ni, scTmp);
            UNPROTECT(2); //result scTmp
            PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
            res = REAL(result);
            PROTECT(scTmp = allocVector(REALSXP, nrx));
            scaleTmp = REAL(scTmp);
            for(j=0; j < nrx; j++) scaleTmp[j] = 0.0;
            ni = nodes[i];
            if(ei < nt){ 
                matp(INTEGER(VECTOR_ELT(dlist, ei)), REAL(contrast), REAL(VECTOR_ELT(P, i)), INTEGER(nr), INTEGER(nc), INTEGER(nco), res); 
            }
            else{ 
                F77_CALL(dgemm)(transa, transb, &nrx, &ncx, &ncx, &one, REAL(VECTOR_ELT(ans, ei-nt)), &nrx, 
                       REAL(VECTOR_ELT(P, i)), &ncx, &zero, res, &nrx);
                for(j=0; j < nrx; j++) scaleTmp[j] += REAL(VECTOR_ELT(scVec, ei-nt))[j];
            }  
//          z= 1L;
        }
        else {
            if(ei < nt){ 
                matp(INTEGER(VECTOR_ELT(dlist, ei)), REAL(contrast), REAL(VECTOR_ELT(P, i)), INTEGER(nr), INTEGER(nc), INTEGER(nco), rtmp);

            }
            else{ 
                F77_CALL(dgemm)(transa, transb, &nrx, &ncx, &ncx, &one, REAL(VECTOR_ELT(ans, ei-nt)), &nrx, 
                    REAL(VECTOR_ELT(P, i)), &ncx, &zero, rtmp, &nrx);
                for(j=0; j < nrx; j++) scaleTmp[j] += REAL(VECTOR_ELT(scVec, ei-nt))[j];
            }
            for(j=0; j < (nrx*ncx); j++) res[j] *= rtmp[j];
//          z+=1L; // epsScale=50; statt+1 von vector
//          if(z>epsScale){
//              scaleMatrix(res, nrx, ncx, scaleTmp);   
//              z=1;     
//          }
        }            
    }
    scaleMatrix(res, nrx, ncx, scaleTmp);
    SET_VECTOR_ELT(ans, ni, result);
    SET_VECTOR_ELT(scVec, ni, scTmp);
    SET_VECTOR_ELT(erg, 0L, ans);
    SET_VECTOR_ELT(erg, 1L, scVec);
    UNPROTECT(5); // result ans 
    return(erg);
}





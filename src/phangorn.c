/* 
 * phangorn.c
 *
 * (c) 2008  Klaus Schliep (k.p.schliep@massey.ac.nz)
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
    unsigned int odd = 1;				// The inner while loop can only terminate with odd == 1 so we don't need to set it inside the for loop.
    double cost, best_cost;
		
    num_splits = (1 << (n - 1));
		
    for (x = 1; x < num_splits; ++x) {
        r = (x - 1) & x;				// r = x without LSB
        nr = (r - 1) & r;				// nr = r without LSB
			
        if (nr) {						// If x contains 1 or 2 bits only, then it has already been computed as a pairwise distance.
            b = x - r;					// b = LSB of x: the "fixed" taxon in the pair.
            best_cost = 1e20;
            e = 0;						// e holds bits to the right of the current p.
				
            while (1) {
                p = r - nr;				// p = 2nd half of pair
                cost = d[nr + e] + d[p + b];
                if (cost < best_cost) best_cost = cost;
					
                if (!nr && odd) break;	// Ensure we get the (LSB with reference taxon) pair when there are an odd number of taxa
                r = nr;
                e += p;
                nr = (r - 1) & r;		// nr = r without LSB
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
				if (j < n - 1) {			// If i == n - 1 then this is a distance between the reference taxon and some other taxon.
					offset += (1 << j);		// Note that "+" is safe since (1 << i) and (1 << j) are always different bits.
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
	
	
	
	
void fitch(unsigned char *dat, int *n, int *m, int *i, int *j, unsigned char * result, int *pars){
	int k;
	unsigned char tmp;
	for(k = 0; k < (*n); k++){
		tmp = dat[*i + k*(*m)] & dat[*j + k*(*m)];
		if(tmp) result[k] = tmp;
		else {result[k] = dat[*i + k*(*m)] | dat[*j + k*(*m)];
			pars[k] += 1;
	    }
	} 
}		



void fitch2(int *dat, int *n, int *m, int i, int j, int *pars){
	int k;
	int tmp;
	for(k = 0; k < (*n); k++){
		tmp = dat[i + k*(*m)] & dat[j + k*(*m)];
		if(tmp) dat[i + k*(*m)] = tmp;
		else {dat[i + k*(*m)] = dat[i + k*(*m)] | dat[j + k*(*m)];
			pars[k] += 1;
	    }
	} 
}		


void fitch3(int *dat, int *n, int *m, int *pars, int *node, int *edge, int *nl) 
{   
	int i, ni, k;
    ni = 0;
    for (i=0; i< *nl; i++) {
        if (ni == node[i]) fitch2(dat, n, m, ni-1, edge[i]-1, pars);	        	        
        else {
	        ni = node[i];
	        for(k = 0; k < (*n); k++) dat[ni-1 + k*(*m)] = dat[edge[i]-1 + k*(*m)];                     
        }
    }
}
	
		
void sankoff(double *dat, int *n, double *cost, int *k, double *result){
	int i, j, h; 
	double tmp[*k], x;
	for(i = 0; i < (*n); i++){
		for(j = 0; j < (*k); j++){
		    for(h = 0; h< (*k); h++){tmp[h] = dat[i + h*(*n)] + cost[h + j*(*k)];}
		    x = tmp[0];
            for(h = 1; h< (*k); h++) {if(tmp[h]<x) {x=tmp[h];}}
		    result[i+j*(*n)] = x;
		    }	   			
		}		
}	


SEXP sankoff2(SEXP sdat, SEXP sn, SEXP scost, SEXP sk){
	int i, j, h, n=INTEGER(sn)[0], k = INTEGER(sk)[0];  
	double tmp[k], x, *cost, *res, *dat;
	SEXP result;
	PROTECT(result = allocMatrix(REALSXP, n, k));
	res = REAL(result);
	PROTECT(sdat = coerceVector(sdat, REALSXP));
	PROTECT(scost = coerceVector(scost, REALSXP));
	cost = REAL(scost);
	dat = REAL(sdat);
	for(i = 0; i < n; i++){
		for(j = 0; j < k; j++){
		    for(h = 0; h< k; h++) tmp[h] = dat[i + h*n] + cost[h + j*k];
		    x = tmp[0];
        for(h = 1; h< k; h++) {if(tmp[h]<x) x=tmp[h]; }
		res[i+j* n] = x;
	    }	   			
    }
    UNPROTECT(3);
    return(result);		
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
	SEXP result, tmp;
	PROTECT(result = allocMatrix(REALSXP, n, k));
	PROTECT(tmp = allocMatrix(REALSXP, n, k));
	res = REAL(result);
	PROTECT(scost = coerceVector(scost, REALSXP));
	cost = REAL(scost);
	rtmp = REAL(tmp);
    for(j=0; j<(n*k); j++) rtmp[j] = 0.0;
    for(j=0; j<(n*k); j++) res[j] = 0.0;   
    sankoff4(REAL(VECTOR_ELT(dat,0)), n, cost, k, rtmp);
    sankoff4(REAL(VECTOR_ELT(dat,1)), n, cost, k, rtmp);
    sankoff4(rtmp, n, cost, k, res);
    sankoff4(REAL(VECTOR_ELT(dat,2)), n, cost, k, res);
    sankoff4(REAL(VECTOR_ELT(dat,3)), n, cost, k, res);
    UNPROTECT(3);
    return(result);		
}	



SEXP sankoff3(SEXP dlist, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP mNodes, SEXP tips){
	R_len_t i, n = length(node), nt = length(tips);
	int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0], mn=INTEGER(mNodes)[0];
	int  ni, ei, j, *edges=INTEGER(edge), *nodes=INTEGER(node);
	SEXP result, dlist2; //tmp, 
	double *res, *rtmp, *cost;
	cost = REAL(scost);
	if(!isNewList(dlist)) error("‘dlist’ must be a list");
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
	UNPROTECT(2); // result dlist2  return dlist2
	return(dlist2);
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

 
static R_INLINE void getP(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
	int i, j, h;
	double tmp[m], res;
	for(i = 0; i < m; i++) tmp[i] = exp(eva[i] * w * el);
	for(i = 0; i < m; i++){	
		for(j = 0; j < m; j++){
			res = 0.0;	
			for(h = 0; h < m; h++)	res += ev[i + h*m] * tmp[h] * evi[h + j*m];
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
			for(h = 0; h < m; h++)	result[i+j*m] += ev[i + h*m] * tmp[h] * evi[h + j*m];				
		}
	}
}


void getPmix(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
	int i, j, h;
	double tmp[m], res[m*m];
	for(i = 0; i < m; i++) tmp[i] = exp(eva[i] * w * el);
	for(j = 0; j < m; j++){
	    for(h = 0; h < m; h++)	res[h + j*m] = evi[h + j*m] * tmp[h];
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
			for(h = 0; h < m; h++)	res += ev[i + h*m] * tmp[h] * evi[h + j*m];
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
			for(h = 0; h < m; h++)	res += ev[i + h*m] * tmp[h] * evi[h + j*m];
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
			for(h = 0; h < m; h++)	res += ev[i + h*m] * tmp[h] * evi[h + j*m];
			result[i+j*m] = res;	
		}
	}
}


SEXP LogLik(SEXP dlist, SEXP P, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP nTips, SEXP mNodes){
	R_len_t i, n = length(node);
	int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0], nt=INTEGER(nTips)[0], mn=INTEGER(mNodes)[0];
	int  ni, ei, j, *edges=INTEGER(edge), *nodes=INTEGER(node);
	SEXP ans, result;  
	double *res, *rtmp;
	if(!isNewList(dlist)) error("‘dlist’ must be a list");
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
			if(ei < nt) matprod(REAL(VECTOR_ELT(dlist, ei)), nrx, ncx, REAL(VECTOR_ELT(P, i)), ncx, ncx, res);
			else matprod(REAL(VECTOR_ELT(ans, ei-nt)), nrx, ncx, REAL(VECTOR_ELT(P, i)), ncx, ncx, res);
			}
		else {
			if(ei < nt) matprod(REAL(VECTOR_ELT(dlist, ei)), nrx, ncx, REAL(VECTOR_ELT(P, i)), ncx, ncx, rtmp);
			else matprod(REAL(VECTOR_ELT(ans, ei-nt)), nrx, ncx, REAL(VECTOR_ELT(P, i)), ncx, ncx, rtmp);
			for(j=0; j < (nrx*ncx); j++) res[j] *= rtmp[j];
			}			
	}
	SET_VECTOR_ELT(ans, ni, result);	
	UNPROTECT(2); // result ans 
	return(ans);
}



double *NR4(SEXP eig, int nc, double el, double *w, double *g, SEXP dad, SEXP child, int ld, int nr, double *bf, double *weight){

	int i, j, k; 
	double *kid;  
	double *eva, *eve, *evei; 
	double *dP, *dtmp, *res, *dF;
	if(!isNewList(eig)) error("‘eig’ must be a list");
	if(!isNewList(dad)) error("dad‘’ must be a list");
	if(!isNewList(child)) error("‘child’ must be a list");	
	dP = (double *) R_alloc(nc*nc, sizeof(double));
	dtmp = (double *) R_alloc(nr*nc, sizeof(double));
	dF = (double *) R_alloc(nr*nc, sizeof(double));
		
	eva = REAL(VECTOR_ELT(eig, 0));
	eve = REAL(VECTOR_ELT(eig, 1));
	evei = REAL(VECTOR_ELT(eig, 2));

	res = (double *) R_alloc(nr, sizeof(double));
    for(k=0; k<nr; k++) res[k] = 0.0;
	    
    for(j=0;j<ld;j++){
		kid = REAL(VECTOR_ELT(child, j)); // kid
		getdP(eva, eve, evei, nc, el, g[j], dP); //dP
    	matprod(REAL(VECTOR_ELT(dad, j)), nr, nc, dP, nc, nc, dtmp); // dad %*% dP
    	for(i = 0; i < (nr*nc); i++) dtmp[i]*=kid[i]; // (dad %*% dP) * kid
    	matprod(dtmp, nr, nc, bf, nc, 1L, dF);// ( (dad %*% dP) * kid ) %*% bf
    	for(k=0; k<nr; k++) res[k] += dF[k] * w[j]; // += ( (dad %*% dP) * kid ) %*% (bf * w)
    	}	        
    return(res);
} 


// test 
SEXP getPM(SEXP eig, SEXP nc, SEXP el, SEXP w){
	R_len_t i, j, nel, nw;
	int m=INTEGER(nc)[0], l=0;
	double *ws=REAL(w);
	double *edgelen=REAL(el);
	double *eva, *eve, *evei;
	SEXP P, RESULT;
	nel = length(el);
	nw = length(w);
	if(!isNewList(eig)) error("‘eig’ must be a list");	
	eva = REAL(VECTOR_ELT(eig, 0));
	eve = REAL(VECTOR_ELT(eig, 1));
	evei = REAL(VECTOR_ELT(eig, 2));
	PROTECT(RESULT = allocVector(VECSXP, nel*nw));	
//	double *p;	
    for(j=0; j<nel; j++){ 
        for(i=0; i<nw; i++){
	        PROTECT(P = allocMatrix(REALSXP, m, m));
//	        p = REAL(P);
            getP(eva, eve, evei, m, edgelen[j], ws[i], REAL(P));//p
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
	if(!isNewList(eig)) error("‘eig’ must be a list");	
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
	if(!isNewList(eig)) error("‘dlist’ must be a list");	
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
	if(!isNewList(eig)) error("‘dlist’ must be a list");	
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
	if(!isNewList(eig)) error("‘dlist’ must be a list");	
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


SEXP FS(SEXP eig, SEXP nc, SEXP el, SEXP w, SEXP g, SEXP dad, SEXP child, SEXP ld, SEXP nr, SEXP basefreq, SEXP weight,SEXP f0, SEXP ll0, SEXP ff0)
{
	R_len_t j, n=length(dad);	
	SEXP RESULT, EL, DAT, P; //, logLik; 
	double *df, *tmp, *tmp2, *f, *wgt=REAL(weight), edle, ledle, newedle, eps=10, *bf=REAL(basefreq);
	double ll, lll, delta, scalep = 1.0, *ws=REAL(w), *gs=REAL(g), l1=0.0, l0;
	int i, k=0, ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0];
	tmp = (double *) R_alloc(nrx, sizeof(double));
	tmp2 = (double *) R_alloc(nrx, sizeof(double));
	f = (double *) R_alloc(nrx, sizeof(double));
  	PROTECT(RESULT = allocVector(VECSXP, 4));
  	edle = REAL(el)[0];	
  	l0 = REAL(ll0)[0];
  	for(i=0; i<nrx ;i++) f[i]=REAL(ff0)[i];
  	
    while ( (eps > 1e-05) &&  (k < 5) ) {
	    df = NR4(eig, ncx, edle, ws, gs, dad, child, INTEGER(ld)[0], nrx, bf, wgt);   
	    for(i=0; i<nrx ;i++) tmp[i]=df[i]/f[i];
        ll=0.0;  
        lll=0.0;        
        for(i=0; i<nrx ;i++) ll+=wgt[i]*tmp[i];
        for(i=0; i<nrx ;i++) lll+=wgt[i]*tmp[i]*tmp[i];
  
        delta = ((ll/lll) < 3) ? (ll/lll) : 3;
        ledle = log(edle) + scalep * delta;
        newedle = exp(ledle);
        if (newedle > 10.0) newedle = 10.0;
        if (newedle < 1e-6) newedle = 1e-6;
        PROTECT(EL = ScalarReal(newedle)); 
        PROTECT(P = getPM(eig, nc, EL, g));
        PROTECT(DAT = getM3(child, dad, P, nr, nc));
  
        for(i=0; i<nrx; i++)f[i] = REAL(f0)[i];
        for(j=0; j<n; j++){
        	matprod(REAL(VECTOR_ELT(DAT, j)), nrx, ncx, bf, ncx, 1L, tmp2);
            for(i=0; i< nrx; i++)f[i] += (tmp2[i] * ws[j]);
        }
        l1 = 0.0;
        for(i=0; i<nrx ;i++) l1 += wgt[i] * log(f[i]);
        eps = l1 - l0;
// some error handling        
        if (eps < 0 || ISNAN(eps)) {
            if (ISNAN(eps))eps = 0;
            else {
                scalep = scalep/2;
                eps = 1;
            }
            newedle = edle;
            l1 = l0;
        }
        else scalep = 1;
        edle=newedle;
        l0 = l1; 
        k ++;
        UNPROTECT(3);// EL, DAT, P
    }   
    PROTECT(EL = ScalarReal(edle));
    PROTECT(P = getPM(eig, nc, EL, g));
	SET_VECTOR_ELT(RESULT, 0, EL);
	SET_VECTOR_ELT(RESULT, 1, getM3(child, dad, P, nr, nc));
	SET_VECTOR_ELT(RESULT, 2, getM3(dad, child, P, nr, nc));
	SET_VECTOR_ELT(RESULT, 3, ScalarReal(l1));
    UNPROTECT(3);
    return (RESULT);
}	



	

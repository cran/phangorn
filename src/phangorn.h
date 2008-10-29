#include <R.h>
#include <Rinternals.h>


void fhm(double *v, int *n);
void out(double *d, double *r, int *n, int *k, int *l);
void fitch3(int *dat, int *n, int *m, int *pars, int *node, int *edge, int *nl);
void sankoff4(double *dat, int n, double *cost, int k, double *result);


SEXP dist2spectra(SEXP dm, SEXP nx, SEXP ns);
SEXP sankoff2(SEXP sdat, SEXP sn, SEXP scost, SEXP sk);
SEXP sankoff3(SEXP dlist, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP mNodes, SEXP tips);
SEXP sankoffQuartet(SEXP dat, SEXP sn, SEXP scost, SEXP sk);
SEXP sankoffQuartet2(SEXP dat, SEXP sn, SEXP scost, SEXP sk);
SEXP rowMin(SEXP sdat, SEXP sn, SEXP sk);
SEXP rowMax(SEXP sdat, SEXP sn, SEXP sk);
SEXP matpro(SEXP X, SEXP Y, SEXP nrx, SEXP ncx, SEXP nry, SEXP ncy);
SEXP LogLik(SEXP dlist, SEXP P, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP nTips, SEXP mNodes);
SEXP getPM(SEXP eig, SEXP nc, SEXP el, SEXP w);
SEXP getdPM(SEXP eig, SEXP nc, SEXP el, SEXP w);
SEXP getdPM2(SEXP eig, SEXP nc, SEXP el, SEXP w);
SEXP getd2PM(SEXP eig, SEXP nc, SEXP el, SEXP w);
SEXP getM3(SEXP dad, SEXP child, SEXP P, SEXP nr, SEXP nc);
SEXP FS(SEXP eig, SEXP nc, SEXP el, SEXP w, SEXP g, SEXP dad, SEXP child, SEXP ld, SEXP nr, SEXP basefreq, SEXP weight,SEXP f0, SEXP ll0, SEXP ff0);

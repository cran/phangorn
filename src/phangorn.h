#include <R.h>
#include <Rinternals.h>


void fhm(double *v, int *n);
void out(double *d, double *r, int *n, int *k, int *l);
void fitch5(int *dat, int *nr, int *pars, int *node, int *edge, int *nl);
void FN(int *dat, int *res, int *nr, int *pars, int *node, int *edge, int *nl, int *pc);
void giveIndex(int *left, int* right, int *ll, int *lr, int *n, int *res);
void sankoff4(double *dat, int n, double *cost, int k, double *result);
void reorder(int *from, int *to, int *n, int *sumNode,  int *neworder, int *root);
void cisort(int *x, int *y, int *a, int *b, int *res);



SEXP dist2spectra(SEXP dm, SEXP nx, SEXP ns);
SEXP sankoff3(SEXP dlist, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP mNodes, SEXP tips);
SEXP pNodes(SEXP data, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge);
SEXP sankoffQuartet(SEXP dat, SEXP sn, SEXP scost, SEXP sk);
SEXP rowMin(SEXP sdat, SEXP sn, SEXP sk);
SEXP rowMax(SEXP sdat, SEXP sn, SEXP sk);
SEXP matpro(SEXP X, SEXP Y, SEXP nrx, SEXP ncx, SEXP nry, SEXP ncy);
SEXP LogLik2(SEXP dlist, SEXP P, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP nTips, SEXP mNodes, SEXP contrast, SEXP nco);
SEXP LogLik4(SEXP dlist, SEXP P, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP nTips, SEXP mNodes, SEXP contrast, SEXP nco);
SEXP getPM(SEXP eig, SEXP nc, SEXP el, SEXP w);
SEXP getdPM(SEXP eig, SEXP nc, SEXP el, SEXP w);
SEXP getdPM2(SEXP eig, SEXP nc, SEXP el, SEXP w);
SEXP getd2PM(SEXP eig, SEXP nc, SEXP el, SEXP w);
SEXP getd2PM2(SEXP eig, SEXP nc, SEXP el, SEXP w);
SEXP getM3(SEXP dad, SEXP child, SEXP P, SEXP nr, SEXP nc);
SEXP FS3(SEXP eig, SEXP nc, SEXP el, SEXP w, SEXP g, SEXP dad, SEXP child, SEXP ld, SEXP nr, SEXP basefreq, SEXP weight, SEXP f0, SEXP ll0, SEXP ff0, SEXP retA, SEXP retB);
SEXP invSites(SEXP dlist, SEXP nr, SEXP nc, SEXP contrast, SEXP nco);




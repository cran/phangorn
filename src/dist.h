#include <R.h>
#include <Rinternals.h>


void giveIndex(int *left, int* right, int *ll, int *lr, int *n, int *res);
void fhm(double *v, int *n);
void distHamming(int *x, double *weight, int *nr, int *l, double *d);
void out(double *d, double *r, int *n, int *k, int *l);

SEXP dist2spectra(SEXP dm, SEXP nx, SEXP ns);
SEXP PWI(SEXP LEFT, SEXP RIGHT, SEXP L, SEXP N, SEXP W, SEXP LI);


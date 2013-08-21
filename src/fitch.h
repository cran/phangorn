#include <R.h>
#include <Rinternals.h>


void fitch_free();
void fitch_init(int *data, int *m, int *n, double *weights, int *nr);
void fitch8(int *dat, int *nr, int *pars, int *node, int *edge, int *nl, double *weight, double *pvec, double *pscore);
void fitchQuartet(int *index, int *n, int *nr, double *psc1, double *psc2, double *weight, double *res);

SEXP FITCH(SEXP dat, SEXP nrx, SEXP node, SEXP edge, SEXP l, SEXP weight, SEXP mx, SEXP q);   
SEXP FITCH345(SEXP nrx, SEXP node, SEXP edge, SEXP l, SEXP mx, SEXP ps);
SEXP FITCHTRIP3(SEXP DAT3, SEXP nrx, SEXP edge, SEXP score, SEXP PS);
SEXP FNALL5(SEXP nrx, SEXP node, SEXP edge, SEXP l, SEXP mx, SEXP my, SEXP root);






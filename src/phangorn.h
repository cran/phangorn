#include <R.h>
#include <Rinternals.h>


void ll_init(int *nr, int *nTips, int *nc, int *k); 
void ll_free();

SEXP PML(SEXP dlist, SEXP EL, SEXP W, SEXP G, SEXP NR, SEXP NC, SEXP K, SEXP eig, SEXP bf, SEXP node, SEXP edge, SEXP NTips, SEXP root, SEXP nco, SEXP contrast, SEXP N);
SEXP rowMax(SEXP sdat, SEXP sn, SEXP sk);
SEXP rowMin(SEXP sdat, SEXP sn, SEXP sk);




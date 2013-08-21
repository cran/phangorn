#include <R.h>
#include <Rinternals.h>

void ll_free();
void ll_init(int *nr, int *nTips, int *nc, int *k);


SEXP invSites(SEXP dlist, SEXP nr, SEXP nc, SEXP contrast, SEXP nco);
SEXP getPM(SEXP eig, SEXP nc, SEXP el, SEXP w);
SEXP PML0(SEXP dlist, SEXP EL, SEXP W, SEXP G, SEXP NR, SEXP NC, SEXP K, SEXP eig, SEXP bf, SEXP node, SEXP edge, SEXP NTips, SEXP root, SEXP nco, SEXP contrast, SEXP N);
SEXP FS5(SEXP eig, SEXP nc, SEXP el, SEXP w, SEXP g, SEXP X, SEXP ld, SEXP nr, SEXP basefreq, SEXP weight, SEXP f0);

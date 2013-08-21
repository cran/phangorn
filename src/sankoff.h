#include <R.h>
#include <Rinternals.h>


SEXP sankoff3(SEXP dlist, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP mNodes, SEXP tips);
SEXP rowMin(SEXP sdat, SEXP sn, SEXP sk);
SEXP pNodes(SEXP data, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge);
SEXP sankoffQuartet(SEXP dat, SEXP sn, SEXP scost, SEXP sk);



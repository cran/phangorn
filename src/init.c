#include <R.h> 
#include <Rinternals.h>
#include "phangorn.h"
#include <R_ext/Rdynload.h>

R_CallMethodDef callMethods[] = {
{"dist2spectra", (DL_FUNC) &dist2spectra, 3},
{"sankoff3", (DL_FUNC) &sankoff3, 8},
{"pNodes", (DL_FUNC) &pNodes, 6},
{"sankoffQuartet", (DL_FUNC) &sankoffQuartet, 4},
{"rowMin", (DL_FUNC) &rowMin, 3},
{"rowMax", (DL_FUNC) &rowMax, 3},
{"matpro", (DL_FUNC) &matpro, 6},
{"LogLik2", (DL_FUNC) &LogLik2, 10},
{"LogLik4", (DL_FUNC) &LogLik4, 10},
{"getPM", (DL_FUNC) &getPM, 4},
{"getdPM", (DL_FUNC) &getdPM, 4},
{"getdPM2", (DL_FUNC) &getdPM2, 4},
{"getd2PM", (DL_FUNC) &getd2PM, 4},
{"getd2PM2", (DL_FUNC) &getd2PM, 4},
{"getM3", (DL_FUNC) &getM3, 5},
{"invSites", (DL_FUNC) &invSites, 5},
{NULL, NULL, 0}
};
     

R_CMethodDef cMethods[] = {
{"fhm", (DL_FUNC) &fhm, 2},
{"out", (DL_FUNC) &out, 5},
{"fitch6", (DL_FUNC) &fitch6, 9},
{"FN2", (DL_FUNC) &FN2, 11},
{"giveIndex", (DL_FUNC) &giveIndex, 6},
{"cisort", (DL_FUNC) &cisort, 5},
{"reorder", (DL_FUNC) &reorder, 6},
{NULL, NULL, 0}
};


void R_init_myLib(DllInfo *info)
{
R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
}

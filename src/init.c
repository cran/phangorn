#include "phangorn.h"
#include <R_ext/Rdynload.h>


#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#define CDEF(name, n)   {#name, (DL_FUNC) &name, n}


static R_CallMethodDef CallEntries[] = {
    CALLDEF(PML, 16),
    {NULL, NULL, 0}
};


static R_CMethodDef CEntries[] = {
    CDEF(ll_init, 4),
    CDEF(ll_free, 0),
    {NULL, NULL, 0}
};



void R_init_phangorn(DllInfo *info)
{
R_registerRoutines(info, CEntries, CallEntries, NULL, NULL);

#define RREGDEF(name)  R_RegisterCCallable("Matrix", #name, (DL_FUNC) name)
    RREGDEF(PML);
    RREGDEF(ll_init);
    RREGDEF(ll_free);
    RREGDEF(rowMin);
    RREGDEF(rowMax);
}

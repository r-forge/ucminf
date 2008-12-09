/*

  Diethelm Wuertz
  ETH Zurich
  www.rmetrics.org
*/


#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>


SEXP fn ;
SEXP gr ;
SEXP environment ;


/*
--------------------------------------------------------------------------------
*/


void cfunc(int *n, double x[], double value[])
{
    SEXP PAR ;
    int i ;
    PROTECT(PAR = findVarInFrame(environment, install(".x"))) ;
    for (i = 0; i < *n; i++) REAL(PAR)[i] = x[i] ;
    value[0] = asReal(eval(fn, environment)) ;
    UNPROTECT(1) ;
}

void cgrad(int *n, double x[], double grval[])
{
    SEXP PAR, OUT ;
    int i ;
    PROTECT(OUT = allocVector(REALSXP, *n)) ;
    PROTECT(PAR = findVarInFrame(environment, install(".x"))) ;
    for (i = 0; i < *n; i++) REAL(PAR)[i] = x[i] ;
    OUT = eval(gr, environment) ;
    for (i = 0; i < *n; i++) grval[i] = REAL(OUT)[i];
    UNPROTECT(2) ;
}


/*
--------------------------------------------------------------------------------
*/


void F77_SUB(cfunc)(int *n, double x[], double value[])
{ 
    cfunc(n, x, value) ;    
}  

void F77_SUB(cgrad)(int *n, double x[], double grval[])
{ 
    cgrad(n, x, grval) ;    
}  


/*
--------------------------------------------------------------------------------
*/


SEXP mfopt(SEXP fnstr, SEXP grstr, SEXP rho)
{
    
    SEXP N, PAR, DX, EPS, MAXFUN, W, IW, ICONTR, GRAD,  GRSTEP ;
   
    fn = fnstr ;
    gr = grstr ;
    environment = rho ;
    
    PROTECT(N     = findVarInFrame(rho, install(".n"))) ;
    PROTECT(PAR   = findVarInFrame(rho, install(".par"))) ;
    PROTECT(DX    = findVarInFrame(rho, install(".stepmax"))) ;
    PROTECT(EPS   = findVarInFrame(rho, install(".eps"))) ;
    PROTECT(MAXFUN= findVarInFrame(rho, install(".maxfun"))) ;
    PROTECT(W     = findVarInFrame(rho, install(".w"))) ;
    PROTECT(IW    = findVarInFrame(rho, install(".iw"))) ;
    PROTECT(ICONTR= findVarInFrame(rho, install(".icontr"))) ;
    PROTECT(GRAD  = findVarInFrame(rho, install(".grad"))) ;
    PROTECT(GRSTEP= findVarInFrame(rho, install(".grstep"))) ;
        
    /* 
    Call the Fortran Solver:
    */
    F77_CALL(ucminf)(INTEGER(N),REAL(PAR), REAL(DX), REAL(EPS), INTEGER(MAXFUN),
		     REAL(W), INTEGER(IW), INTEGER(ICONTR), 
		     INTEGER(GRAD),REAL(GRSTEP)) ;
    
    UNPROTECT(10) ;
    return R_NilValue;
}




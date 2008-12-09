ucminf = function(par, fn, gr = NULL, ..., control = list(), hessian=0) {
  con <- list(trace=0, grtol=1e-6, xtol=1e-12, stepmax=100, maxit=500,
              grad='forward',gradstep=c(1e-6,1e-8), hessian.init = NULL,
              method.args = NULL)
  con[(namc <- names(control))] <- control
  stopifnot(length(con$gradstep)==2)
  fnstr <- quote(fn(.x, ...))	
  grstr <- quote(gr(.x, ...))	
  rho = new.env(parent = environment())
  n <- length(par)
  eps <- c(con$grtol,con$xtol)
  if(!is.null(gr)) { grad <- 0 } else {
    grad <- ifelse(con$grad=='forward',1,2) #else central
  }
  iw <- n*ceiling(max(n+1,(n+11)/2)) + 10
  w <- rep(0,iw)
  trace <- con$trace>0
  icontr = 1+trace+2*!is.null(con$hessian.init)
  logicMat <- (matrix(-(1:n^2),n,n,byrow=TRUE)+matrix(1:n^2,n,n))<=0
  if(!is.null(con$hessian.init)){
    inv.hessian.init <- solve(con$hessian.init)
    w[(4*n+1):(4*n+n*(n+1)/2)] <-  inv.hessian.init[logicMat]
  }
  assign(".n",      as.integer(n)           , envir = rho) 
  assign(".x",      as.double(rep(0,n))     , envir = rho)
  assign(".par",    as.double(par)          , envir = rho)
  assign(".stepmax",as.double(con$stepmax)  , envir = rho)
  assign(".eps",    as.double(eps)          , envir = rho)
  assign(".maxfun", as.integer(con$maxit)   , envir = rho)
  assign(".w",      as.double(w)            , envir = rho)
  assign(".iw",     as.integer(iw)          , envir = rho)
  assign(".icontr", as.integer(icontr)      , envir = rho)
  assign(".grad",   as.integer(grad)        , envir = rho)
  assign(".grstep", as.double(con$gradstep) , envir = rho)
  #
  .Call("mfopt", fnstr, grstr, rho, PACKAGE = "ucminf")
  #
  W <- get(".w", envir = rho)
  icontr <- get(".icontr", envir = rho)
  ans = list(
    par = get(".par", envir = rho), 
    value = W[1],
    convergence = icontr,
    message = switch(as.character(icontr),
      '1' ='Stopped by small gradient (grtol).',
      '2' ='Stopped by small step (xtol).',
      '3' ='Stopped by iteration limit (maxit)',
      '4' ='Stopped by zero step from line search',
      '-2'="Computation did not start: length(par) = 0.",
      '-4'="Computation did not start: stepmax is too small.",
      '-5'="Computation did not start: grtol or xtol <= 0.",
      '-6'="Computation did not start: maxit <= 0.",
      '-7'="Computation did not start: given hessian not pos. definite.",
      '-8'="Computation did not start: iw too small."
      )
    )
  if(0<icontr) {
    if(hessian == 1) 
      ans$hessian <- hessian(fn, ans$par, method="Richardson",
                             method.args=con$method.args, ...)
    if(hessian == 2 | hessian == 3) {
      COV <- matrix(0,n,n)
      COV[logicMat] <- W[(4*n+1):(4*n+n*(n+1)/2)]
      COV <- t(COV)+COV-diag(diag(COV))
      ans$invhessian <- COV
    }
    if(hessian == 3)
      ans$hessian <- solve(COV)
    ans$info = c( maxgradient = W[2],
                  laststep    = W[3],
                  stepmax     =get(".stepmax", envir = rho),
                  neval       = get(".maxfun", envir = rho)
                )
  }
  nm <- names(par)
  if (!is.null(nm)) 
    names(ans$par) <- nm
  ans
}		




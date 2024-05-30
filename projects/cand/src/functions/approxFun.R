approxFun <- function(U_sol) {
  require(fields)
  obj <- list(x = U_sol$t, y = U_sol$x, z = t(U_sol$U))
  approxfun_2d <- function(t,x) {
    interp.surface(obj, cbind(t, x))
  }
  return(approxfun_2d)
}
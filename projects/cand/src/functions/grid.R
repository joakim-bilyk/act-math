sinh <- function(x) (exp(x)-exp(-x))/2
sinh_inv <- function(x) log(x+sqrt(x^2+1))
g_identity <- function(
    r0,r1
) {
  function(v) {
    r0 + (r1-r0)*v
  }
}
g_hyperbolic <- function(r0,r1,r_mid = NULL,c = (r1+r0)/2,d = (r1-r0)/50) {
  if (!is.null(r_mid)) {c <- r_mid}
  coeffs <- c(
    sinh_inv((r0-c)/d),sinh_inv((r1-c)/d)
  )
  function(v) {
    tmp <- c + d*sinh(coeffs[1]+v*(coeffs[2]-coeffs[1]))
    pmax(r0,pmin(r1,tmp))
  }
}
g_positive <- function(r0,r1,d=(r1-r0)/50) {
  coeff <- sinh_inv((r1-r0)/d)
  function(v) {
    tmp <- r0+d*sinh(v*coeff)
    pmax(r0,pmin(r1,tmp))
  }
}
g_negative <- function(r0,r1,d=(r1-r0)/50) {
  coeff <- sinh_inv((r0-r1)/d)
  function(v) {
    tmp <- r1+d*sinh((1-v)*coeff)
    pmax(r0,pmin(r1,tmp))
  }
}
g_comb <- function(
    r0,r_lower,r_upper,
    r1,v1,vlower = (1-v1)/2,
    d1=(r_lower-r0)/50,d2=(r1-r_upper)/50) {
  coeff <- c(
    sinh_inv((r0-r_lower)/d1),
    sinh_inv((r1-r_upper)/d2)
  )
  vs <- c(vlower,vlower+v1)
  function(v) {
    tmp <- (v<vs[1])*(r_lower+d1*sinh((1-v/vs[1])*coeff[1])) +
      (v>=vs[1])*(v<vs[2])*(r_lower + (r_upper-r_lower)*(v-vs[1])/v1)+
      (v>=vs[2])*(r_upper+d2*sinh((v-vs[2])/(1-vs[2])*coeff[2]))
    pmax(r0,pmin(r1,tmp))
  }
}
RungeKuttaStage <- function(mu,sigma,dsigma,stages=2,par = NULL) {
  if (stages == 1) {
    function(t,x,delta_t,delta_W) {
      x + delta_t*mu(t,x) + delta_W*sigma(t,x)
    }
  } else if (stages == 2) {
    function(t,x,dt,dW) {
      tmp <- numeric(4)
      tmp[1] <- sigma(t,x)
      tmp[2:3] <- c(
        mu(t,x)*dt, tmp[1]*dW
      )
      tmp[4] <- tmp[2] + tmp[3]
      x + 
        0.5 * tmp[2] + 0.5*tmp[3] +
        0.5 * mu(t + dt,x + tmp[4]) * dt +
        0.5 * sigma(t + dt,x + tmp[4]) * dW -
        0.5 * tmp[1] * dsigma(t,x) * dt
    }
  } else if (stages == 3) {
    if (is.null(par)) par <- 1
    function(t,x,dt,dW) {
      tmp <- numeric(4)
      tmp[1] <- sigma(t,x)
      tmp[2:3] <- c(
        mu(t,x)*dt, tmp[1]*dW
      )
      tmp[4] <- tmp[2] + tmp[3]
      x + 
        0.5 * tmp[2] + 0.5*tmp[3] +
        0.5 * mu(t + dt,x + tmp[4]) * dt +
        1/(2+6*par[1]**2) * sigma(t + dt,x + tmp[2]+par[1]*tmp[3]) * dW +
        3*par[1]**2/(2+6*par[1]**2) * sigma(t + dt,x + tmp[2]-tmp[3]/(3*par[1])) * dW +
        0.5 * tmp[1] * dsigma(t,x) * (dW**2-dt)
    }
  } else {
    stop("stages in 1 to 3 only defined")
  }
}
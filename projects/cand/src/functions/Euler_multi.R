Euler_multi <- function(mu,sigma) {
  function(t,x,dt,dW) {
    x <- x + mu(t,x)*dt + sigma(t,x) %*% dW
    as.numeric(x)
  }
}
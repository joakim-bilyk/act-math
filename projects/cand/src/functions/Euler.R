Euler <- function(mu,sigma) {
  function(t,x,delta_t,delta_W) {
    x + mu(t,x)*delta_t + sigma(t,x)*delta_W
  }
}
Milstein <- function(mu,sigma,dsigma) {
  function(t,x,delta_t,delta_W) {
    x + mu(t,x)*delta_t + sigma(t,x)*(delta_W +
                                        0.5*dsigma(t,x)*(delta_W**2 - delta_t))
  }
}
simIncr <- function(T,n) {
  Delta_W <- rnorm(n,mean = 0,sd = sqrt(T/n))
  W <- c(0,cumsum(Delta_W))
  return(list(W = W, t = T*0:n/n))
}
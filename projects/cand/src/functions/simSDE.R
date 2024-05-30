simSDEDiff <- function(x0,T,n,F,simWeiner = simBridge,weinerPath = NULL,d=1) {
  m <- length(x0)
  if (is.null(weinerPath)) {
    W <- sapply(1:d, function(i) simWeiner(T,n)$W)
  } else {
    W <- weinerPath %>% as.matrix()
  }
  n <- dim(W)[1]-1; d <- dim(W)[2]
  X <- vector("list", length = n +1)
  X[[1]] <- x0
  delta_t <- T/n
  for (i in 1:n) {
    X[[i+1]] <- as.numeric(F(delta_t*(i-1),X[[i]],delta_t,W[i+1,]-W[i,]))
  }
  if (d == 1) {
    X <- unlist(X)
  }
  return(list(X = X,t = T*0:n/n))
}
simSDEDiffs <- function(N,x0,T,n,F,simWeiner = simBridge,d=1) {
  sims <- vector("list",length = N)
  for (sim_n in 1:N) {
    sims[[sim_n]] <- simSDEDiff(x0,T,n,F,simWeiner=simWeiner,d= d)
  }
  return(sims)
}
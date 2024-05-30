RK_solve <- function(
    t0,y0,yhit,dy,yeps = 10^-3,teps = 10^-2, nmax = 1000,T=Inf,
    step = function(t,y) {
      min(teps, abs(yeps/dy(t,y)))
    }
) {
  # Stopping criteria
  if (y0>yhit) {
    continue <- function(t,y,n) {(y>yhit) && (n < nmax) && (t<T)}
  } else {
    continue <- function(t,y,n) {(y<yhit) && (n < nmax) && (t<T)}
  }
  # Setup
  n <- 1; t <- t0; y <- y0
  while (continue(t,y,n)) {
    t0 <- t; y0 <- y # Old values
    dt <- step(t0,y0)
    tmp <- RK_epoch(t0,y0,dt,dy)
    t <- tmp$t; y <- tmp$y # New values
    n <- n + 1
  }
  # Solve done
  # Set final value based on middle point yhit in [y0,y]
  correct <- abs((yhit-y0)/(y-y0)) # Percentage error
  tmp <- RK_epoch(t0,y0,dt*correct,dy) # take smaller step
  t <- tmp$t; y <- tmp$y # Final values
  if (n == nmax) {
    warning("Maximum iterations reached")
    return(list(t=T,y=y,n=n))
  }
  return(list(t=t,y=y,n=n))
}
simBridge <- function(T,n) {
  log2n <- ceiling(log2(n))
  n0 <- n; n <- 2**log2n
  W <- numeric(n+1) # Define vector
  W[n+1] <- rnorm(1,mean = 0,sd = sqrt(T))
  if (n>=2) {
    for (i in 1:log2n) {
      # Simulate
      Delta_W <- rnorm(n = 2**(i-1),mean=0,sd = sqrt(T/(2**(i+1))))
      
      # Calculate indices for left, right, and midpoint
      points <- 2^(log2n - i+1)*(0:(2^(i - 1))) + 1
      left <- points[-length(points)]
      right <- points[-1]
      mid <- left + 2**(log2n-i)
      
      # Insert midpoints
      W[mid] <- 0.5*W[left]+0.5*W[right]+Delta_W
    } 
  }
  return(list(W = W, t = T*0:n/n))
}
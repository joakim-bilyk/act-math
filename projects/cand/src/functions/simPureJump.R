simPureJump <- function(
    n,z0,t0,T,lambda, logProb = TRUE,samplertau = NULL,steps = 10^2
) {
  if (!is.function(lambda)) {
    lambda0 <- lambda
    lambda <- function(t) {return(lambda0)}
  }
  J <- dim(lambda(0))[1] # Number of states
  if (J==2) {
    # Change state if only two states
    samplerZ <- function(t,z) {(1:J)[-z]}
  } else {
    # Sample based on the transition rates at the transition time
    samplerZ <- function(t,z) {
      lambda_t <- lambda(t)
      probs <- -lambda_t[z,-z]/lambda_t[z,z]
      if (is.nan(sum(probs))) {
        z
      } else {
        sample(x = (1:J)[-z],size =1,prob = -lambda_t[z,-z]/lambda_t[z,z])
      }
    }
  }
  if (is.null(samplertau)){
    if (logProb) {
      dy <- function(tau,y) {lambda(t+tau)[z,z]}
      samplertau <- function(t,u,j,...) {
        # Solve H(t)<log u
        yeps <- -log(u)/steps; teps <- (T-t)/steps
        suppressWarnings(
          RK_solve(0,0,log(u),dy,yeps = yeps,teps = teps,nmax = 3*steps,T=T+teps)$t
        )
      }
    } else {
      dy <- function(tau,y) {lambda(t+tau)[z,z]*y}
      samplertau <- function(t,u,j,...) {
        # Solve H(t)<log u
        yeps <- (1-u)/steps; teps <- (T-t)/steps
        suppressWarnings(
          RK_solve(0,1,u,dy,yeps = yeps,teps = teps,nmax = 3*steps,T=T+teps)$t
        )
      }
    }
  }
  if (length(z0) == 1) {
    z0all <- z0
    inidist <- function() {z0all}
  } else {
    probs <- z0/sum(z0)
    inidist <- function() {sample(1:J,size = 1, prob = probs)}
  }
  paths <- list()
  for (i in 1:n) {
    t <- t0; tt <- NULL; zz <- NULL
    z0 <- inidist(); z <- z0
    while (t < T) {
      tt <- c(tt,t); zz <- c(zz,z) # Append
      U <- runif(1,0,1) # Draw Unif(0,1)
      t <- t + samplertau(t,U,z,T+t) # Find inverse
      z <- samplerZ(t,z)
    }
    tt <- c(tt,T); zz <- c(zz,rev(zz)[1])
    paths[[i]] <- list(times = tt,states = zz)
  }
  if (n== 1) {
    return(paths[1])
  } else {
    return(paths)
  }
}
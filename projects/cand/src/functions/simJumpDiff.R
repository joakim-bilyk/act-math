simJumpDiff <- function(
  N,
  z0,x0,
  t0,T,n,
  alpha,gamma, dgamma = NULL, ddgamma = NULL,beta,lambda,
  method = "RK",
  simWeiner = simBridge,
  saveH = FALSE,
  ...
) {
  # Same as pure jump
  if (!is.function(lambda)) {
    lambda0 <- lambda
    lambda <- function(t) {return(lambda0)}
  }
  J <- dim(lambda(t0,x0))[1] # Number of states
  m <- length(x0) # Number of diffusions
  d <- dim(gamma(t0,x0,z0))[2] # Number of Weiner
  if(is.null(d)) d <- 1 # if vector
  # Sampler for Z
  if (J==2) {
    # Change state if only two states
    samplerZ <- function(t,z,x) {(1:J)[-z]}
  } else {
    # Sample based on the transition rates at the transition time
    samplerZ <- function(t,z,x) {
      lambda_t <- lambda(t,x)
      probs <- -lambda_t[z,-z]/lambda_t[z,z]
      if (is.nan(sum(probs))) {
        z
      } else {
        sample(x = (1:J)[-z],size =1,prob = -lambda_t[z,-z]/lambda_t[z,z])
      }
    }
  }
  dH <- function(t,y) {
    -lambda(t,x)[z,z]
  }
  alpha_j <- function(t,x) alpha(t,x,z)
  gamma_j <- function(t,x) gamma(t,x,z)
  if (is.null(dgamma)) {
    dgamma_j <- NULL
  } else {
    dgamma_j <- function(t,x) dgamma(t,x,z)
  }
  if (is.null(ddgamma)) {
    ddgamma_j <- NULL
  } else {
    ddgamma_j <- function(t,x) ddgamma(t,x,z)
  }
  updator_x <- switch (method,
                       RK = RK_multi(mu = alpha_j,sigma = gamma_j,
                                     dsigma = dgamma_j, ddsigma = ddgamma_j,
                                     m = m, d= d,...),
                       Euler = Euler_multi(mu = alpha_j,sigma = gamma_j)
  )
  # Simulate
  sims <- vector("list", length = N)
  for (sim_n in 1:N) {
    # Simulate Weiner paths
    W <- sapply(1:d, function(i) simWeiner(T-t0,n)$W)
    # Reset number of time steps based on Wiener increments
    n <- dim(W)[1]-1
    tt <- seq(t0,T,length.out = n+1); dt <- tt[2]-tt[1];t <- t0
    x <- x0
    X <- vector("list", length = m)
    names(X) <- paste0("X",1:m)
    for (mi in 1:m) X[[mi]] <- numeric(n+1) # Allocate
    for (mi in 1:m) X[[mi]][1] <- x[mi] # Fill initial
    z <- z0; times_z <- t0; states_z <- z
    x_leftlim <- list();x_leftlim[[1]] <- x
    H <- numeric(n+1); H[1] <- 0
    if (saveH) { H_int <- numeric(n+1); H_int[1] <- dH(t0,0)}
    U <- runif(1); logU <- -log(U)
    logUs <- logU
    n_transitions <- 0
    i <- 1
    for (i in 1:n) {
      #print(alpha_j(t,x))
      #print(c(logU,H,x,W[i+1,]-W[i,]))
      # Compute H((t+dt)-)
      if (saveH) {H_int[i+1] <- dH(t,H[i])}
      H[i+1] <- RK_epoch(t0 = t,y0 = H[i],dt = dt,dy = dH)$y
      # Compute X((t+dt)-)
      x <- updator_x(t,x,dt,dW = W[i+1,]-W[i,])
      for (mi in 1:m) X[[mi]][i+1] <- x[mi] #save
      t <- tt[i+1]
      # If transition happened on (t,t+dt]
      if (H[i+1] > logU) {
        # Save H path
        n_transitions <- n_transitions + 1
        # Transition in z
        H[i+1] <- 0
        U <- runif(1); logU <- -log(U); logUs <- c(logUs,logU)
        z_pre <- z
        z <- samplerZ(t,z,x)
        times_z <- c(times_z,t); states_z <- c(states_z,z)
        x_leftlim[[n_transitions+1]] <- x
        if (saveH) {H_int[i+1] <- dH(t,H)}
        x <- x + beta(t,x,z_pre)[,z]
        for (mi in 1:m) X[[mi]][i+1] <- x[mi] #save
      }
      i <- i + 1
    }
    x_leftlim[[n_transitions+2]] <- x
    if (saveH) {
      sims[[sim_n]] <- list(
        Z = list(times = c(times_z,T),states = c(states_z,z)),
        X = X, times = tt, H = list(H = H, Intensity = H_int),
        x_left = x_leftlim
      )
    } else {
      sims[[sim_n]] <- list(
        Z = list(times = c(times_z,T),states = c(states_z,z)),
        X = X, times = tt, x_left = x_leftlim
      )
    }
  }
  return(sims)
}
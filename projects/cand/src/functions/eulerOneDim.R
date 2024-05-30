eulerOneDim <- function(
    t0,t1,
    x0 = NULL,x1 = NULL,
    K,J,
    alpha,beta,gamma,zeta,uT,
    method = "midpoint",
    grid = g_identity(x0,x1),
    u_boundary = NULL,
    u_r_boundary = NULL,
    u_rr_boundary = NULL,
    ...
) {
  tt <- seq(t0,t1,length.out = K + 1)
  dt <- (t1-t0)/K
  xx <- grid(0:J/J)
  U <- matrix(ncol = K + 1,nrow= J + 1)
  # Terminal condition
  U[,K+1] <- uT(xx)
  # Compute each point
  for (k in K:1) {
    # k is the known plane
    Ut <- switch(method,
                 explicit = eulerEpochOneDim_explicit(
                   U[,k+1],dt,tt[k+1],xx,
                   alpha,beta,
                   gamma,zeta,
                   u_r_boundary,
                   u_rr_boundary
                 ),
                 implicit = eulerEpochOneDim_midpoint(
                   U[,k+1],dt,tt[k+1],xx,
                   alpha,beta,
                   gamma,zeta,
                   u_r_boundary,
                   u_rr_boundary,
                   theta = 1
                 ),
                 midpoint = eulerEpochOneDim_midpoint(
                   U[,k+1],dt,tt[k+1],xx,
                   alpha,beta,
                   gamma,zeta,
                   u_r_boundary,
                   u_rr_boundary,
                   ...
                 )
    )
    # Apply boundary conditions on U
    if (!is.null(u_boundary)) Ut <- u_boundary(Ut)
    U[,k] <- Ut
  }
  return(list(U = U,t = tt, x = xx))
}
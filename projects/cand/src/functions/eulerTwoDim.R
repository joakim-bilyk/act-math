eulerTwoDim <- function(
  t0,t1,
  r0 = NULL,r1 = NULL,
  y0 = NULL,y1 = NULL,
  K,J,L,
  alpha,beta,
  gamma_r,gamma_y,
  zeta_r,zeta_y,zeta_ry,
  uT, method = "midpoint",
  u_boundary = NULL,
  u_r_boundary = NULL,
  u_rr_boundary = NULL,
  u_y_boundary = NULL,
  u_yy_boundary = NULL,
  grid_r = g_identity(r0,r1),
  grid_y = g_identity(y0,y1),
  ...
) {
  tt <- seq(t0,t1,length.out = K + 1)
  dt <- (t1-t0)/K
  rr <- grid_r(0:J/J)
  yy <- grid_y(0:L/L)
  U <- list()
  # Terminal condition
  U[[K+1]] <- outer(rr,yy, uT) %>% t()
  # Compute each point
  for (k in K:1) {
    Ut <- switch(method,
                 implicit = eulerEpochTwoDim(
                   U[[k+1]],dt,tt[k+1],rr,yy,
                   alpha,beta,
                   gamma_r,gamma_y,zeta_r,zeta_y,zeta_ry, theta = 1,
                   u_r_boundary = u_r_boundary,u_rr_boundary = u_rr_boundary,
                   u_y_boundary = u_y_boundary,u_yy_boundary = u_yy_boundary
                 ),
                 explicit = eulerEpochTwoDim(
                   U[[k+1]],dt,tt[k+1],rr,yy,
                   alpha,beta,
                   gamma_r,gamma_y,zeta_r,zeta_y,zeta_ry, theta = 0,
                   u_r_boundary = u_r_boundary,u_rr_boundary = u_rr_boundary,
                   u_y_boundary = u_y_boundary,u_yy_boundary = u_yy_boundary
                 ),
                 midpoint = eulerEpochTwoDim(
                   U[[k+1]],dt,tt[k+1],rr,yy,
                   alpha,beta,
                   gamma_r,gamma_y,zeta_r,zeta_y,zeta_ry,
                   u_r_boundary = u_r_boundary,u_rr_boundary = u_rr_boundary,
                   u_y_boundary = u_y_boundary,u_yy_boundary = u_yy_boundary,
                   ...
                 )
    )
    # Apply boundary conditions on U
    if (!is.null(u_boundary)) Ut <- u_boundary(Ut)
    U[[k]] <- Ut
    #print(paste0(Sys.time(),": Iteration ",K+1-k," of ",K))
  }
  return(list(U = U,t = tt, r = rr, y = yy))
}
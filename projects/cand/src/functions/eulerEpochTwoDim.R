eulerEpochTwoDim <- function(
  Ut,dt,t,r,y,
  alpha,beta,
  gamma_r,gamma_y,
  zeta_r,zeta_y,zeta_ry, theta = 0.5,
  u_r_boundary = NULL,
  u_rr_boundary = NULL,
  u_y_boundary = NULL,
  u_yy_boundary = NULL
) {
  # k is the known plane
  mtx1 <- eulerMatrixTwoDim(
    t,r,y,
    alpha,beta,
    gamma_r,gamma_y,zeta_r,zeta_y,zeta_ry,
    u_r_boundary,u_rr_boundary,
    u_y_boundary,u_yy_boundary
  )
  J <- dim(mtx1$p0)[2]-1; L <- dim(mtx1$p0)[1] - 1
  P1 <- eulerPMatrix(mtx1)
  mtx0 <- eulerMatrixTwoDim(
    t-dt,r,y,
    alpha,beta,
    gamma_r,gamma_y,zeta_r,zeta_y,zeta_ry,
    u_r_boundary,u_rr_boundary,
    u_y_boundary,u_yy_boundary
  )
  P0 <- eulerPMatrix(mtx0)
  I <- diag(1,ncol = dim(P0)[1],nrow = dim(P0)[1])
  mtx_invert <- solve(
    I + theta*dt*P0 + (1-theta)*dt*diag(as.numeric(t(mtx1$B)))
  )
  Uk_mtx <- I - (1-theta)*dt*P1 - theta*dt*diag(as.numeric(t(mtx0$B)))
  Uk <- as.numeric(t(Ut))
  Ut1 <- mtx_invert %*% (Uk_mtx %*% Uk -dt*(theta*as.numeric(t(mtx0$A)) +
                                              (1-theta)*as.numeric(t(mtx1$A))))
  Ut1 <- matrix(Ut1,ncol = J+1, byrow = TRUE)
  return(Ut1)
}
eulerEpochOneDim_explicit <- function(
  Ut,dt,t,x,
  alpha,beta,
  gamma,zeta,
  u_r_boundary,
  u_rr_boundary
) {
  matrices <- eulerMatrixOneDim(
    t,x,
    alpha,beta,gamma,zeta,
    u_r_boundary,u_rr_boundary
  )
  Ut1 <- (Ut- dt*matrices$P %*% Ut - dt*matrices$a) /
    (1+ dt*matrices$b)
  return(Ut1)
}
eulerEpochOneDim_midpoint <- function(
  Ut,dt,t,x,
  alpha,beta,
  gamma,zeta,
  u_r_boundary,
  u_rr_boundary,
  theta = 0.5
) {
  matrices1 <- eulerMatrixOneDim(
    t,x,
    alpha,beta,gamma,zeta,
    u_r_boundary,u_rr_boundary
  )
  matrices0 <- eulerMatrixOneDim(
    t-dt,x,
    alpha,beta,gamma,zeta,
    u_r_boundary,u_rr_boundary
  )
  p_dim <- dim(matrices0$P)[1]
  I <- diag(1,ncol = p_dim,nrow=p_dim)
  mtx_invert <- solve(
    I + dt*theta*matrices0$P + dt*(1-theta)*diag(matrices1$b,ncol = p_dim,nrow=p_dim)
  )
  Ut1 <- mtx_invert %*% (
    (I -dt*(1-theta)*matrices1$P - dt*theta*diag(matrices0$b,ncol = p_dim,nrow=p_dim)) %*% Ut -
      dt*theta*matrices0$a-dt*(1-theta)*matrices1$a)
  
  return(Ut1)
}
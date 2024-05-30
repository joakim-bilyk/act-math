eulerMatrixOneDim <- function(
    t,x,
    alpha,beta,gamma,zeta,
    u_r_boundary,u_rr_boundary
) {
  require(SparseM) # Matrix is sparse
  J <- length(x)-1
  dx <- x[2:(J+1)]- x[1:J]
  if (is.null(u_r_boundary)) {
    u_r_boundary_tmp <- c(0,0)
    I_r <- c(0,0)
  } else {
    u_r_boundary_tmp <- u_r_boundary(t)
    I_r <- ifelse(is.na(u_r_boundary_tmp),0,1)
    u_r_boundary_tmp <- ifelse(is.na(u_r_boundary_tmp),0,u_r_boundary_tmp)
  }
  if (is.null(u_rr_boundary)) {
    u_rr_boundary_tmp <- c(0,0)
    I_rr <- c(0,0)
  } else {
    u_rr_boundary_tmp <- u_rr_boundary(t)
    I_rr <- ifelse(is.na(u_rr_boundary_tmp),0,1)
    u_rr_boundary_tmp <- ifelse(is.na(u_rr_boundary_tmp),0,u_rr_boundary_tmp)
  }
  # Compute values
  alpha_tmp <- alpha(t,x)
  if (length(alpha_tmp) == 1) alpha_tmp <- rep(alpha_tmp,J+1)
  beta_tmp <- beta(t,x)
  if (length(beta_tmp) == 1) beta_tmp <- rep(beta_tmp,J+1)
  gamma_tmp <- gamma(t,x)
  if (length(gamma_tmp) == 1) gamma_tmp <- rep(gamma_tmp,J+1)
  zeta_tmp <- zeta(t,x)
  if (length(zeta_tmp) == 1) zeta_tmp <- rep(zeta_tmp,J+1)
  # Compute interior p's
  p0 <- ((dx[2:J]-dx[1:(J-1)])/(dx[1:(J-1)]*dx[2:J]))*gamma_tmp[2:J] - 
    (1/(dx[1:(J-1)]*dx[2:J]))*zeta_tmp[2:J]
  p1 <- (dx[1:(J-1)]/(dx[2:J]*(dx[2:J]+dx[1:(J-1)])))*gamma_tmp[2:J] +
    (1/(dx[2:J]*(dx[2:J]+dx[1:(J-1)])))*zeta_tmp[2:J]
  pm1 <- -(dx[2:J]/(dx[1:(J-1)]*(dx[2:J]+dx[1:(J-1)])))*gamma_tmp[2:J] +
    (1/(dx[1:(J-1)]*(dx[2:J]+dx[1:(J-1)])))*zeta_tmp[2:J]
  
  # Insert diagonal
  P <- diag(c(0,p0,0),ncol = J+1,nrow = J+1)
  # Insert lower diagonal
  P[-c(1,J+1),-c(J,J+1)] <- P[-c(1,J+1),-c(J,J+1)] + diag(pm1)
  # Insert upper diagonal
  P[-c(1,J+1),-c(1,2)] <- P[-c(1,J+1),-c(1,2)] + diag(p1)
  # Fix upper row
  P[1,1] <- -((dx[2]+2*dx[1])/(dx[1]*(dx[1]+dx[2])))*gamma_tmp[1]*(1-I_r[1]) +
    (1/(dx[1]*(dx[1]+dx[2])))*zeta_tmp[1]*(1-I_rr[1])
  P[1,2] <- ((dx[1]+dx[2])/(dx[1]*dx[2]))*gamma_tmp[1]*(1-I_r[1]) -
    1/(dx[2]*dx[1])*zeta_tmp[1]*(1-I_rr[1])
  P[1,3] <- -dx[1]/(dx[2]*(dx[1]+dx[2]))*gamma_tmp[1]*(1-I_r[1]) +
    (1/(dx[2]*(dx[1]+dx[2])))*zeta_tmp[1]*(1-I_rr[1])
  # Fix lower row
  P[J+1,J+1] <- ((dx[J-1]+2*dx[J])/(dx[J]*(dx[J]+dx[J-1])))*gamma_tmp[J+1]*(1-I_r[2]) +
    (1/(dx[J]*(dx[J]+dx[J-1])))*zeta_tmp[J+1]*(1-I_rr[2])
  P[J+1,J] <- -((dx[J]+dx[J-1])/(dx[J]*dx[J-1]))*gamma_tmp[J+1]*(1-I_r[2]) -
    1/(dx[J-1]*dx[J])*zeta_tmp[J+1]*(1-I_rr[2])
  P[J+1,J-1] <- dx[J]/(dx[J-1]*(dx[J]+dx[J-1]))*gamma_tmp[J+1]*(1-I_r[2]) +
    (1/(dx[J-1]*(dx[J]+dx[J-1])))*zeta_tmp[J+1]*(1-I_rr[2])
  
  a <- alpha_tmp
  a[1] <- a[1] +
    u_r_boundary_tmp[1]*I_r[1]*gamma_tmp[1] +
    u_rr_boundary_tmp[1]*I_rr[1]*zeta_tmp[1]
  a[J+1] <- a[J+1] +
    u_r_boundary_tmp[2]*I_r[2]*gamma_tmp[J+1] +
    u_rr_boundary_tmp[2]*I_rr[2]*zeta_tmp[J+1]
  b <- beta_tmp
  
  matrices <- list(P = P,a= a,b=b)
  return(matrices)
}
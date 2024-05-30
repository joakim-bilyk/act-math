eulerMatrixTwoDim <- function(
    t,r,y,
    alpha,beta,
    gamma_r,gamma_y,
    zeta_r,zeta_y,zeta_ry,
    u_r_boundary = NULL,
    u_rr_boundary = NULL,
    u_y_boundary = NULL,
    u_yy_boundary = NULL
) {
  J <- length(r)-1; L <- length(y) - 1
  dr <- r[2:(J+1)]- r[1:J]
  dy <- y[2:(L+1)]- y[1:L]
  # Compute boundaries
  if (is.null(u_r_boundary)) {
    u_r_boundary_tmp <- matrix(0,nrow = 2,ncol = L+1)
    I_r <- c(0,0)
  } else {
    u_r_boundary_tmp <- u_r_boundary(t,y)
    I_r <- rowSums(!is.na(u_r_boundary_tmp))/(L+1)
    u_r_boundary_tmp <- ifelse(is.na(u_r_boundary_tmp),0,u_r_boundary_tmp)
  }
  if (is.null(u_rr_boundary)) {
    u_rr_boundary_tmp <- matrix(0,nrow = 2,ncol = L+1)
    I_rr <- c(0,0)
  } else {
    u_rr_boundary_tmp <- u_rr_boundary(t,y)
    I_rr <- rowSums(!is.na(u_rr_boundary_tmp))/(L+1)
    u_rr_boundary_tmp <- ifelse(is.na(u_rr_boundary_tmp),0,u_rr_boundary_tmp)
  }
  
  if (is.null(u_y_boundary)) {
    u_y_boundary_tmp <- matrix(0,nrow = 2,ncol = J+1)
    I_y <- c(0,0)
  } else {
    u_y_boundary_tmp <- u_y_boundary(t,r)
    I_y <- rowSums(!is.na(u_y_boundary_tmp))/(J+1)
    u_y_boundary_tmp <- ifelse(is.na(u_y_boundary_tmp),0,u_y_boundary_tmp)
  }
  if (is.null(u_yy_boundary)) {
    u_yy_boundary_tmp <- matrix(0,nrow = 2,ncol = J+1)
    I_yy <- c(0,0)
  } else {
    u_yy_boundary_tmp <- u_yy_boundary(t,r)
    I_yy <- rowSums(!is.na(u_yy_boundary_tmp))/(J+1)
    u_yy_boundary_tmp <- ifelse(is.na(u_yy_boundary_tmp),0,u_yy_boundary_tmp)
  }
  # Compute values
  alpha_tmp <- outer(r,y, alpha, t = t) %>% t()
  beta_tmp <- outer(r,y, beta, t= t) %>% t()
  gamma_r_tmp <- outer(r,y, gamma_r, t= t) %>% t()
  gamma_y_tmp <- outer(r,y, gamma_y, t= t) %>% t()
  zeta_r_tmp <- outer(r,y, zeta_r, t= t) %>% t()
  zeta_y_tmp <- outer(r,y, zeta_y, t= t) %>% t()
  zeta_ry_tmp <- outer(r,y, zeta_ry, t= t) %>% t()
  # Compute inner p's
  delta_j <- dr[2:J]; delta_jm1 <- dr[1:(J-1)]
  eps_l <- dy[2:L]; eps_lm1 <- dy[1:(L-1)]
  p0 <- rbind(0,
    cbind(0,t(((delta_j-delta_jm1)*t(gamma_r_tmp[2:L,2:J])-t(zeta_r_tmp[2:L,2:J])) / 
                (delta_j*delta_jm1)) +
          ((eps_l-eps_lm1)*(gamma_y_tmp[2:L,2:J])-zeta_y_tmp[2:L,2:J]) / 
            (eps_l*eps_lm1),0),
    0)
  p_r <- list(
    rbind(0,
      cbind(0,t((-delta_j*t(gamma_r_tmp[2:L,2:J])+t(zeta_r_tmp[2:L,2:J])) / 
                  (delta_jm1*(delta_j+delta_jm1))),0),0),
    rbind(0,
      cbind(0,t((delta_jm1*t(gamma_r_tmp[2:L,2:J])+t(zeta_r_tmp[2:L,2:J])) / 
                  (delta_j*(delta_j+delta_jm1))),0),0)
  )
  p_y <- list(
    rbind(0,cbind(0,(-eps_l*gamma_y_tmp[2:L,2:J]+zeta_y_tmp[2:L,2:J]) / 
                    (eps_lm1*(eps_l+eps_lm1)),0),0),
    rbind(0,cbind(0,(eps_lm1*gamma_y_tmp[2:L,2:J]+zeta_y_tmp[2:L,2:J]) /
                    (eps_l*(eps_l+eps_lm1)),0),0)
  )
  p_ry <- rbind(0,cbind(0,(
    matrix(eps_l,ncol = 1) %*% matrix(delta_j,nrow = 1) +
      matrix(eps_l,ncol = 1) %*% matrix(delta_jm1,nrow = 1) +
      matrix(eps_lm1,ncol = 1) %*% matrix(delta_j,nrow = 1) +
      matrix(eps_lm1,ncol = 1) %*% matrix(delta_jm1,nrow = 1)
  )^(-1)*zeta_ry_tmp[2:L,2:J],0),0)
  p0 + p_r[[1]] + p_r[[2]] + p_y[[1]] + p_y[[2]] + p_ry # Test
  
  # Compute a and b
  B <- beta_tmp
  A <- alpha_tmp
  A[1,1:(J+1)] <- A[1,1:(J+1)] +
    I_y[1]*u_y_boundary_tmp[1,1:(J+1)]*gamma_y_tmp[1,1:(J+1)] +
    I_yy[1]*u_yy_boundary_tmp[1,1:(J+1)]*zeta_y_tmp[1,1:(J+1)]
  A[L+1,1:(J+1)] <- A[L+1,1:(J+1)] +
    I_y[2]*u_y_boundary_tmp[2,1:(J+1)]*gamma_y_tmp[L+1,1:(J+1)] +
    I_yy[2]*u_yy_boundary_tmp[2,1:(J+1)]*zeta_y_tmp[L+1,1:(J+1)]
  A[1:(L+1),1] <- A[1:(L+1),1] +
    I_r[1]*u_r_boundary_tmp[1,1:(L+1)]*gamma_r_tmp[1:(L+1),1] +
    I_rr[1]*u_rr_boundary_tmp[1,1:(L+1)]*zeta_r_tmp[1:(L+1),1]
  A[1:(L+1),J+1] <- A[1:(L+1),J+1] +
    I_r[2]*u_r_boundary_tmp[2,1:(L+1)]*gamma_r_tmp[1:(L+1),J+1] +
    I_rr[2]*u_rr_boundary_tmp[2,1:(L+1)]*zeta_r_tmp[1:(L+1),J+1]
  
  matrices <- list(
    p0 = p0,  p_r = p_r, p_y = p_y,p_ry =p_ry,
    B = B, A = A, coeff = list(gamma_r = gamma_r_tmp,
                               gamma_y = gamma_y_tmp,
                               zeta_r = zeta_r_tmp,
                               zeta_y = zeta_y_tmp,
                               zeta_ry = zeta_ry_tmp),
    dr = dr, dy = dy, I = list(I_r = I_r, I_rr = I_rr,
                               I_y = I_y, I_yy = I_yy)
  )
  return(matrices)
}
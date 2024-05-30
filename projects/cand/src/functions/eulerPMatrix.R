eulerPMatrix <- function(
  matrices
) {
  require(SparseM) # Matrix is sparse
  J <- dim(matrices$p0)[2]-1; L <- dim(matrices$p0)[1]-1
  
  # Insert diagonal
  # Keep in mind R inserts by column so the matrix is constructed
  # with r in the vertical axis and y in the horizontal axis
  # In other words, we have to transpose the matrix at the end
  P <- matrix(0,ncol = (J+1)*(L+1), nrow = (J+1)*(L+1))
  
  ## Inner points
  jj <- rep(2:J,L-1)
  ll <- sort(rep(2:L,J-1))
  P_index <- jj + (ll-1)*(J+1) # inner points
  # Below yields the diagonal position of the inner points
  # (P_index-1)*(J+1)*(L+1) + P_index
  P[(P_index-1)*(J+1)*(L+1) + P_index ] <- t(matrices$p0)[P_index]
  P[(P_index-1)*(J+1)*(L+1) + P_index + 1] <- t(matrices$p_r[[2]])[P_index] #        (j+1,l)
  P[(P_index-1)*(J+1)*(L+1) + P_index - 1] <- t(matrices$p_r[[1]])[P_index] #        (j-1,l)
  P[(P_index-1)*(J+1)*(L+1) + P_index + (J + 1)] <- t(matrices$p_y[[2]])[P_index] #  (j  ,l+1)
  P[(P_index-1)*(J+1)*(L+1) + P_index - (J + 1)] <- t(matrices$p_y[[1]])[P_index] #  (j  ,l-1)
  P[(P_index-1)*(J+1)*(L+1) + P_index + (J + 1) + 1] <- -t(matrices$p_ry)[P_index]#  (j+1,l+1)
  P[(P_index-1)*(J+1)*(L+1) + P_index - (J + 1) - 1] <- -t(matrices$p_ry)[P_index]#  (j-1,l-1)
  P[(P_index-1)*(J+1)*(L+1) + P_index + (J + 1) - 1] <- t(matrices$p_ry)[P_index]#   (j-1,l+1)
  P[(P_index-1)*(J+1)*(L+1) + P_index - (J + 1) + 1] <- t(matrices$p_ry)[P_index]#   (j+1,l-1)
  
  ## Corners
  # j = 0, l = 0
  j <- c(1,J+1,1,J+1)
  l <- c(1,1,L+1,L+1)
  sr <- c(1,-1,1,-1)
  sy <- c(1,1,-1,-1)
  for (it in 1:4) {
    P_index <- (l[it]-1)*(J+1) + j[it]
    
    g_r <- t(matrices$coeff$gamma_r)[P_index]
    g_y <- t(matrices$coeff$gamma_y)[P_index]
    z_r <- t(matrices$coeff$zeta_r)[P_index]
    z_y <- t(matrices$coeff$zeta_y)[P_index]
    z_ry <- t(matrices$coeff$zeta_ry)[P_index]
    if (sr[it] == 1) {
      delta0 <- matrices$dr[1]
      delta1 <- matrices$dr[2]
      g_r <- g_r * (1 - matrices$I$I_r[1])
      z_r <- z_r * (1 - matrices$I$I_rr[1])
    } else {
      delta0 <- matrices$dr[J]
      delta1 <- matrices$dr[J-1]
      g_r <- g_r * (1 - matrices$I$I_r[2])
      z_r <- z_r * (1 - matrices$I$I_rr[2])
    }
    if (sy[it] == 1) {
      eps0 <- matrices$dy[1]
      eps1 <- matrices$dy[2]
      g_y <- g_y * (1 - matrices$I$I_y[1])
      z_y <- z_y * (1 - matrices$I$I_yy[1])
    } else {
      eps0 <- matrices$dy[L]
      eps1 <- matrices$dy[L-1]
      g_y <- g_y * (1 - matrices$I$I_y[2])
      z_y <- z_y * (1 - matrices$I$I_yy[2])
    }
    
    (-sr[it]*(delta1+2*delta0)/(delta0*(delta0+delta1)))*g_r -
      sy[it]*(eps0/(eps1*(eps0+eps1)))*g_r +sr[it]*((delta0+delta1)/(delta0*delta1))*g_r
    
    # (j,l)
    P[(P_index-1)*(J+1)*(L+1) + P_index] <- (-sr[it]*(delta1+2*delta0)/(delta0*(delta0+delta1)))*g_r +
      (-sy[it]*(eps1+2*eps0)/(eps0*(eps0+eps1)))*g_y +
      (1/(delta0*(delta0+delta1)))*z_r +
      (1/(eps0*(eps0+eps1)))*z_y
    # (j,l+sy)
    P[(P_index-1)*(J+1)*(L+1) + P_index + sy[it]*(J + 1)] <- sy[it]*((eps0+eps1)/(eps0*eps1))*g_y -
      (1/(eps0*eps1))*z_y
    # (j,l+2sy)
    P[(P_index-1)*(J+1)*(L+1) + P_index + 2*sy[it]*(J + 1)] <- -sy[it]*(eps0/(eps1*(eps0+eps1)))*g_y +
      (1/(eps1*(eps0+eps1)))*z_y
    # (j + sr,l)
    P[(P_index-1)*(J+1)*(L+1) + P_index + sr[it]] <- sr[it]*((delta0+delta1)/(delta0*delta1))*g_r -
      (1/(delta0*delta1))*z_r
    # (j+2sr,l)
    P[(P_index-1)*(J+1)*(L+1) + P_index + 2*sr[it]] <- -sr[it]*(delta0/(delta1*(delta0+delta1)))*g_r +
      (1/(delta1*(delta0+delta1)))*z_r
    # Cross
    tmp <- 1/(sr[it]*sy[it]*delta1*eps1)
    P[(P_index-1)*(J+1)*(L+1) + P_index + sr[it] + sy[it]*(J + 1)] <- tmp * z_ry
    P[(P_index-1)*(J+1)*(L+1) + P_index + 2*sr[it] + sy[it]*(J + 1)] <- tmp * z_ry
    P[(P_index-1)*(J+1)*(L+1) + P_index + sr[it] + 2*sy[it]*(J + 1)] <- tmp * z_ry
    P[(P_index-1)*(J+1)*(L+1) + P_index + 2*sr[it] + 2*sy[it]*(J + 1)] <- tmp * z_ry
  }
  
  ## Edges
  sr <- c(0,0,1,-1)
  sy <- c(1,-1,0,0)
  ii <- list(2:J,
             (L)*(J+1)+2:J,
             (2:L-1)*(J+1)+1,
             (2:L)*(J+1))
  delta_j <- matrices$dr[2:J]; delta_jm1 <- matrices$dr[1:(J-1)]
  eps_l <- matrices$dy[2:L]; eps_lm1 <- matrices$dy[1:(L-1)]
  for (it in 1:4) {
    i <- ii[[it]]
    P_index <- (i-1)*(J+1)*(L+1) + i
    
    g_r <- t(matrices$coeff$gamma_r)[i]
    g_y <- t(matrices$coeff$gamma_y)[i]
    z_r <- t(matrices$coeff$zeta_r)[i]
    z_y <- t(matrices$coeff$zeta_y)[i]
    z_ry <- t(matrices$coeff$zeta_ry)[i]
    
    if (sr[it] == 1) {
      delta0 <- matrices$dr[1]
      delta1 <- matrices$dr[2]
      g_r <- g_r * (1 - matrices$I$I_r[1])
      z_r <- z_r * (1 - matrices$I$I_rr[1])
    } else {
      delta0 <- matrices$dr[J]
      delta1 <- matrices$dr[J-1]
      g_r <- g_r * (1 - matrices$I$I_r[2])
      z_r <- z_r * (1 - matrices$I$I_rr[2])
    }
    if (sy[it] == 1) {
      eps0 <- matrices$dy[1]
      eps1 <- matrices$dy[2]
      g_y <- g_y * (1 - matrices$I$I_y[1])
      z_y <- z_y * (1 - matrices$I$I_yy[1])
    } else {
      eps0 <- matrices$dy[L]
      eps1 <- matrices$dy[L-1]
      g_y <- g_y * (1 - matrices$I$I_y[2])
      z_y <- z_y * (1 - matrices$I$I_yy[2])
    }
    
    # (j,l)
    if (sr[it] == 0) {
      z_ry_tmp <- sy[it]*(1/((2*eps0+eps1)*(delta_j+delta_jm1)))*z_ry
      
      P[P_index] <- 
        # Central approximation r
        ((delta_j - delta_jm1)/(delta_j * delta_jm1))*g_r -
        (1/(delta_j * delta_jm1))*z_r +
        # Forward y
        (-sy[it]*(eps1+2*eps0)/(eps0*(eps0+eps1)))*g_y +
        (1/(eps0*(eps0+eps1)))*z_y
      # r directions
      P[P_index + 1] <- (delta_jm1*g_r+z_r)/(delta_j*(delta_j + delta_jm1)) -
        2*z_ry_tmp
      P[P_index - 1] <- (-delta_j*g_r+z_r)/(delta_jm1*(delta_j + delta_jm1)) +
        2*z_ry_tmp
      # y directions
      P[P_index + sy[it]*(J+1)] <-
        sy[it]*((eps0+eps1)/(eps0*eps1))*g_y -
        (1/(eps0*eps1))*z_y
      P[P_index + 2*sy[it]*(J+1)] <- 
        -sy[it]*(eps0/(eps1*(eps0+eps1)))*g_y +
        (1/(eps1*(eps0+eps1)))*z_y
      
      # Cross directions
      P[P_index + 1 + (J+1)*sy[it]] <- z_ry_tmp
      P[P_index - 1 + (J+1)*sy[it]] <- -z_ry_tmp
      P[P_index + 1 + 2*(J+1)*sy[it]] <- z_ry_tmp
      P[P_index - 1 + 2*(J+1)*sy[it]] <- -z_ry_tmp
    } else {
      
      z_ry_tmp <- sr[it]*(1/((2*delta0+delta1)*(eps_l+eps_lm1)))*z_ry
      
      P[P_index] <- 
        # Central approximation y
        ((eps_l - eps_lm1)/(eps_l * eps_lm1))*g_y -
        (1/(eps_l * eps_lm1))*z_y +
        # Forward r
        (-sr[it]*(delta1+2*delta0)/(delta0*(delta0+delta1)))*g_r +
        (1/(delta0*(delta0+delta1)))*z_r
      # r directions
      P[P_index + sr[it]] <- 
        sr[it]*((delta0+delta1)/(delta0*delta1))*g_r -
        (1/(delta0*delta1))*z_r
      P[P_index + 2*sr[it]] <- 
        -sr[it]*(delta0/(delta1*(delta0+delta1)))*g_r +
        (1/(delta1*(delta0+delta1)))*z_r
      # y directions
      P[P_index + (J+1)] <- (eps_lm1*g_y+z_y)/(eps_l*(eps_l + eps_lm1)) -
        2*z_ry_tmp
      P[P_index - (J+1)] <- (-eps_lm1*g_y+z_y)/(eps_l*(eps_l + eps_lm1)) +
        2*z_ry_tmp
      
      # Cross directions
      P[P_index + sr[it] + (J+1)] <- z_ry_tmp
      P[P_index + sr[it] - (J+1)] <- -z_ry_tmp
      P[P_index + 2*sr[it] + (J+1)] <- z_ry_tmp
      P[P_index + 2*sr[it] - (J+1)] <- -z_ry_tmp
    }
    
  }
  
  # Transpose
  P <- t(P)
  return(P)
}
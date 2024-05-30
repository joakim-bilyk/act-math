RK_multi <- function(mu,sigma,dsigma = NULL,ddsigma = NULL,m,d, c = 1/2) {
  I <- rep(1,d)
  function (t,x,dt,dW) {
    mut <- mu(t,x); sigmat <- sigma(t,x)
    R1 <- 0; R2 <- 0
    if (!is.null(dsigma)) {
      dsigmat <- dsigma(t,x)
      if (d == 1) {
        V <- -dt
      } else {
        V <- diag(rep(-1,d)); upper <- sample(c(-1,1),size = sum(upper.tri(V)),
                                              replace = TRUE)
        V[upper.tri(V)] <- upper; V <- t(V); V[upper.tri(V)] <- -upper
        V <- V*dt
      }
      tmp <- (dW %*% t(dW) + V)
      R1 <- calcR1(dt,dW,sigmat,dsigmat,V)
    }
    if (!is.null(ddsigma)) {
      ddsigmat <- ddsigma(t,x)
      if (d > 1) R2 <- calcR2(dt,dW,sigmat,dsigmat,ddsigmat)
    }
    tmp <- sigmat %*% dW*I
    mu_part <- 0.5 * (mut + mu(t + dt, x + mut*dt + sigmat %*% dW)) * dt
    sigma_part <- 0.5 * sigmat %*% dW +
      1/(2+6*c**2)*sigma(t+dt,x+mut*dt+c*tmp) %*% dW +
      3*c**2/(2+6*c**2)*sigma(t+dt,x+mut*dt-tmp/(3*c))%*% dW
    as.numeric(x + mu_part + sigma_part + R1 + R2)
  }
}
calcR1 <- function(dt,dW,sigmat,dsigmat,V) {
  d <- dim(sigmat)[2]
  m <- dim(sigmat)[1]
  tmp1 <- 0.5*t(sigmat) %*% dsigmat
  tmp2 <- dW %*% t(dW) + V
  R1 <- numeric(m)
  for (k in 1:m) R1[k] <- sum(tmp1[,1:d+(k-1)*d] * tmp2)
  return(R1)
}
calcR2 <- function(dt,dW,sigmat,dsigmat,ddsigmat) {
  d <- dim(sigmat)[2]
  m <- dim(sigmat)[1]
  tmp <- dt*dW/6
  diag_tmp <- kronecker(diag(d),sigmat)
  R2 <- numeric(m)
  diag_placement <- rep(1:d,d)^2+m^2*rep(0:(d-1),each = d)
  diag_0 <- (1:d)^2+(1:d-1)*d^2
  row_placement <- rep(1:d,each = d) + rep(0:(d-1),each = d)*d^2 + rep(0:(d-1),d)*d
  for (k in 1:m) {
    mtx <- t(sigmat) %*% ddsigmat[1:m+(k-1)*m,] %*% diag_tmp
    mtx <- mtx %*% diag(rep(tmp,each = m))
    mtx[diag_0] <- 0
    R2[k] <- sum(mtx[diag_placement]) - sum(mtx[row_placement])
  }
  return(R2)
}
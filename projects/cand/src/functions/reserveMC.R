reserveMC <- function(samples,b,B,BT,discount_m) {
  times <- samples[[1]]$times
  n_times <- length(times)
  T <- max(times); t <- min(times); dt <- (T-t)/(n_times-1)
  m <- length(samples[[1]]$X)
  N <- length(samples)
  J <- dim(B(0,rep(0,m)))[1]
  b_payment <- sapply(samples, function(L) {
    row_n <- 1
    # Prepare data
    b_tmp <- numeric(n_times)
    dat <- cbind(t = times,as.data.frame(L$X))
    j <- L$Z$states[1]; t0 <- L$Z$times[1]
    for (i in 1:(length(L$Z$states)-1)) {
      j <- L$Z$states[i]
      t0 <- L$Z$times[i]
      t1 <- L$Z$times[i+1]
      ns <- (1:n_times)[(times>=t0)&(times<t1)]
      b_tmp[ns] <- b(times[ns],dat[ns,2:(m+1)],j)*dt
    }
    b_tmp <- b_tmp * dat[,discount_m + 1]
    return(sum(b_tmp))
  })
  transition_payment <- sapply(samples,function(L) {
    if (length(L$Z$states) > 2) {
      mtx <- matrix(0,ncol = (J-1)*J,nrow = length(L$Z$states)-2)
      for (trans_n in 1:(length(L$Z$states) - 2)) {
        x_pre <- L$x_left[[trans_n+1]] #left-limit
        t_pre <- L$Z$times[trans_n+1]
        z_pre <- L$Z$states[trans_n]; z <- L$Z$states[trans_n+1]
        payment <- B(t_pre,x_pre)[z_pre,z]
        mtx[trans_n,(z_pre-1)*(J-1) + z-1] <- payment*x_pre[discount_m]
      }
      return(colSums(mtx))
    } else {
      return(rep(0,(J-1)*J))
    }
  }) %>% t()
  maturity_payment <- sapply(samples,function(L) {
    transitions <- length(L$x_left)-2
    BT(L$x_left[[transitions+2]],L$Z$states[transitions+2]) *
      L$x_left[[transitions+2]][discount_m]
  })
  discounted_payments <- cbind(
    b_payment,rowSums(transition_payment)
    ,maturity_payment,transition_payment) %>% as.data.frame()
  trans_names <- sapply(1:J, function(j) {
    sapply((1:J)[-j], function(k) {
      paste0("B_",j,"_",k)
    })
  }) %>% as.character()
  colnames(discounted_payments) <- c("b","B","B_T",trans_names)
  convergence_V <- discounted_payments
  for (i in 1:dim(convergence_V)[2]) { 
    convergence_V[,i] <- cumsum(convergence_V[,i])/1:N
  }
  convergence_V[,"Total"] <- rowSums(convergence_V[,1:3])
  return(
    list(
      reserve = convergence_V[N,"Total"],
      reserve_b = convergence_V[N,"b"],
      reserve_B = convergence_V[N,"B"],
      reserve_B_T = convergence_V[N,"B_T"],
      convergence_V = convergence_V
    )
  )
}
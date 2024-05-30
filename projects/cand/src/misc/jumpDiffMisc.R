plotJumpDiff <- function(sample) {
  require(latex2exp)
  m <- length(sample$X)
  transitions <- length(sample$Z$times)-2
  # Plot for H
  # Data frame for U
  if ("logU" %in% names(sample)) {
    tstart <- sapply(sample$H, function(h) h$times[1])
    tend <- sapply(sample$H, function(h) rev(h$times)[1])
    plotU <- data.frame(t = tstart, tend = tend,y = -log(sample$U))
    p1 <- ggplot() +
      # logU
      geom_segment(data = plotU, aes(x=t,xend = tend, y = y, yend = y),col="darkblue") +
      geom_point(data = plotU,aes(x= t, y=y),col="darkblue") +
      geom_point(data = plotU,aes(x= tend, y=y),shape =1,col="darkblue")
    plotH <- data.frame(t = tstart, tend = tend)
    for (i in 1:length(sample$H)) {
      H <- sample$H[[i]]
      plotH[i,c("Hstart","Hend")] <- c(H$H[1],rev(H$H)[1])
      plotH0 <- data.frame(t = H$times, y = H$H)
      plotH1 <- data.frame(t = H$times, y = H$Intensity)
      p1 <- p1 + geom_line(data = plotH0,aes(x=t,y=y),col="darkred")+ geom_line(data = plotH1,aes(x=t,y=y),col="darkgreen")
    }
    p1 <- p1 + geom_point(data = plotH,aes(x=t,y=Hstart),col="darkred") +
      geom_point(data = plotH,aes(x=tend,y=Hend),shape = 1,col="darkred") +
      labs(y=TeX(paste0("$H_t$ and $-\\log U_i$"))) + theme_custom()
  }
  
  # Plot of X
  X_df <- as.data.frame(sample$X)
  X_df[,"times"] <- sample$times
  tstart <- sample$Z$times[1:transitions]
  tend <- sample$Z$times[2:(transitions+1)]
  plotX <- data.frame(t = tstart, tend = tend)
  p2 <- ggplot()
  for (i in 1:(transitions+1)) {
    X <- X_df[(sample$times < tend[i]),]
    plotX0 <- melt(X,id.vars= "times")
    p2 <- p2 + geom_line(data = plotX0,aes(x=times,y=value,col=variable))
  }
  p2 <- p2 + geom_point(data = plotX,aes(x=t,y=Xstart,col=Var)) +
    geom_point(data = plotX,aes(x=tend,y=Xend,col=Var),shape = 1) +
    labs(y=TeX(paste0("$X_t$"))) + theme_custom()
  
  p2_all <- lapply( 1:m, function(l) {
    plotX <- data.frame(t = tstart, tend = tend)
    p <- ggplot()
    for (i in 1:length(sample$X_jump)) {
      X <- sample$X_jump[[i]] %>% as.data.frame()
      plotX[i,c("Xstart","Xend")] <- c(X[1,l+1],rev(X[,l+1])[1])
      plotX0 <- data.frame(t = X$times, y = X[,l+1])
      p <- p + geom_line(data = plotX0,aes(x=t,y=y),col="darkred")
    }
    p <- p + geom_point(data = plotX,aes(x=t,y=Xstart),col="darkred") +
      geom_point(data = plotX,aes(x=tend,y=Xend),shape = 1,col="darkred") +
      labs(y=TeX(paste0("$X_t$"))) + theme_custom()
  })
  
  # Plot of Z
  plotZ <- sample$Z %>% as.data.frame()
  plotZ[,"tend"] <- c(plotZ$times[-1],max(sample$X$t))
  p3 <- ggplot(plotZ,aes(x = times)) + geom_segment(aes(xend = tend, y =states,yend =states),col="darkred") +
    geom_point(aes(x=times,y=states),col="darkred") +
    geom_point(aes(x=tend,y=states),shape = 1,col="darkred") +
    labs(y=TeX(paste0("$Z_t$"))) + theme_custom()
  
  return(list(p1 = p1,p2 = p2,p3 = p3,p2all = p2_all))
}

# Plot transition probabilities
jumpDiffTransition <- function(samples,J) {
  require(AalenJohansen)
  # Make paths object
  paths <- lapply(samples, function(L) { L$Z})
  ts <- range(unlist(sapply(paths, function(L) L$times)))
  fit <- aalen_johansen(paths)
  # Hazards
  hazards <- matrix(unlist(lapply(fit$Lambda, FUN = function(L) rev(L[-seq(1,J^2,J+1)]))),
                    nrow=length(fit$t),byrow=TRUE) %>% as.data.frame() %>%
    mutate(t = fit$t)
  nm <- matrix(ncol = J,nrow = J)
  nm <- paste0(sort(rep(1:J,J))," -> ",rep(1:J,J))
  colnames(hazards)[-dim(hazards)[2]] <- nm[-seq(1,J^2,J+1)]
  rg <- range(fit$t)
  plotdf1 <- hazards %>% reshape2::melt(id.vars="t")
  p1 <- plotdf1 %>%
    ggplot() + geom_step(aes(x = t,y =value,col = variable)) +
    labs(col = "Transitions",y=TeX("$\\Lambda_{jk}(t)$")) + theme_custom() +
    xlim(ts[1],ts[2])
  
  # Transition probs
  probs <- matrix(unlist(fit$p),
                  nrow=length(fit$t),byrow=TRUE) %>% as.data.frame() %>%
    mutate(t = fit$t)
  colnames(probs)[1:J] <- paste0("j = ",1:J)
  plotdf2 <- probs %>% reshape2::melt(id.vars="t")
  p2 <- plotdf2 %>%
    ggplot() + geom_step(aes(x = t,y =value,col = variable)) +
    ylim(0,1) + labs(col = "State",y=TeX("$P(Z_t=j)$")) + theme_custom()+
    xlim(ts[1],ts[2])
  return(list(
    hazards = hazards, plots = list(p1 = p1,p2 = p2), probs =probs,fit = fit
  ))
}

# Diagnostics
jumpDiffDiag <- function(samples) {
  m <- length(samples[[1]]$X)-1
  N <- length(samples)
  meanplots <- lapply(1:m, function(l) {
    paths <- sapply(samples, function(L) {
      L$X[l+1]
    }) %>% as.data.frame()
    df <- data.frame(t = as.numeric(samples[[1]]$X$t))
    df[,"mean"] <- rowMeans(paths)
    df[,"variance"] <- matrixStats::rowVars(paths %>% as.matrix())
    df[,"lower"] <- df$mean - df$variance
    df[,"upper"] <- df$mean + df$variance
    colnames(paths) <- paste0("Sim",1:N)
    paths[,"t"] <- as.numeric(df$t)
    p1 <- ggplot(df,aes(x=t)) + geom_line(aes(y=mean)) + theme_classic()
    p2 <- ggplot(df,aes(x=t)) + geom_line(aes(y=variance)) + theme_classic()
    list(mean = p1,variance = p2,data = df,paths = paths)
  })
  return(list(
    means = meanplots
  ))
}

jumpDiffMean <- function(samples,J){
  require(AalenJohansen)
  # Make paths object
  paths <- lapply(samples, function(L) { L$Z})
  ts <- range(unlist(sapply(paths, function(L) L$times)))
  fit <- aalen_johansen(paths)
  # Hazards
  hazards <- matrix(unlist(lapply(fit$Lambda, FUN = function(L) rev(L[-seq(1,J^2,J+1)]))),
                    nrow=length(fit$t),byrow=TRUE) %>% as.data.frame() %>%
    mutate(t = fit$t)
  nm <- matrix(ncol = J,nrow = J)
  nm <- paste0(sort(rep(1:J,J))," -> ",rep(1:J,J))
  colnames(hazards)[-dim(hazards)[2]] <- nm[-seq(1,J^2,J+1)]
  
  # Transition probs
  probs <- matrix(unlist(fit$p),
                  nrow=length(fit$t),byrow=TRUE) %>% as.data.frame() %>%
    mutate(t = fit$t)
  colnames(probs)[1:J] <- paste0("j = ",1:J)
  
  #X mean
  m <- length(samples[[1]]$X)
  X <- vector("list",m)
  for (i in 1:m) {
    dat <- lapply(samples, function(L) {
      as.numeric(L$X[[i]])
    }) %>%as.data.frame() %>%as.matrix()
    mean <- rowMeans(dat)
    var <- matrixStats::rowVars(dat)
    X[[i]] <- list(mean = mean, var = var,times = samples[[1]]$times)
  }
  return(list(
    hazards = hazards, probs =probs, X = X
  ))
}

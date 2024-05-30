lambda_sin <- function() {
  x <- function(t) {
    sin(t)+1
  }
  lambda12 <- function(t) {x(t)}
  lambda21 <- function(t) {x(t+pi/2)}
  return(
    function(t) {
      lambda<- matrix(
        c(NA,lambda12(t),
          lambda21(t),NA),ncol=2,nrow=2,byrow = TRUE)
      diag(lambda) <- -rowSums(lambda,na.rm = TRUE)
      lambda
    }
  )
}

# Plots
diagnostic <- function(paths,lambda) {
  require(AalenJohansen)
  if (!is.function(lambda)) {
    lambda0 <- lambda
    lambda <- function(t) {return(lambda0)}
  }
  fit <- aalen_johansen(paths)
  J <- dim(lambda(1))[1]
  # Hazards
  hazards <- matrix(unlist(lapply(fit$Lambda, FUN = function(L) rev(L[-seq(1,J^2,J+1)]))),
                    nrow=length(fit$t),byrow=TRUE) %>% as.data.frame() %>%
    mutate(t = fit$t)
  nm <- matrix(ncol = J,nrow = J)
  nm <- paste0(sort(rep(1:J,J))," -> ",rep(1:J,J))
  colnames(hazards)[-dim(hazards)[2]] <- nm[-seq(1,J^2,J+1)]
  rg <- range(fit$t)
  hazards_true <- RK(rg[1],matrix(rep(0,J^2),ncol = J,nrow = J),rg[2],dt=(rg[2]-rg[1])/1000,
                     dy = function(t,y) {lambda(t)})
  plotdf1 <- hazards %>% reshape2::melt(id.vars="t")
  plotdf1b <- matrix(unlist(lapply(hazards_true$y, FUN = function(L) rev(L[-seq(1,J^2,J+1)]))),
                     nrow=length(hazards_true$y),byrow=TRUE) %>% as.data.frame() %>%
    mutate(t = hazards_true$t)
  colnames(plotdf1b)[-dim(plotdf1b)[2]] <- nm[-seq(1,J^2,J+1)]
  plotdf1b <- plotdf1b %>% reshape2::melt(id.vars="t")
  
  p1 <- plotdf1 %>%
    ggplot() + geom_step(aes(x = t,y =value,col = variable)) +
    geom_line(data = plotdf1b,aes(x = t,y =value,col = variable), linetype="dashed") +
    labs(col = "Transitions",y=TeX("$\\Lambda_{jk}(t)$")) + theme_custom()
  
  # Transition probs
  probs <- matrix(unlist(fit$p),
                  nrow=length(fit$t),byrow=TRUE) %>% as.data.frame() %>%
    mutate(t = fit$t)
  colnames(probs)[1:J] <- paste0("j = ",1:J)
  prd_int <- prodint(rg[1],rg[2],step_size = (rg[2]-rg[1])/1000,lambda)
  prd_int <- matrix(unlist(lapply(prd_int, FUN = function(L) fit$I0%*%L)),
                    nrow=length(prd_int),byrow=TRUE) %>% as.data.frame() %>%
    mutate(t = seq(rg[1],rg[2], by = (rg[2]-rg[1])/1000))
  colnames(prd_int)[1:J] <- paste0("j = ",1:J)
  plotdf2 <- probs %>% reshape2::melt(id.vars="t")
  plotdf2b <- prd_int %>% reshape2::melt(id.vars="t")
  p2 <- plotdf2 %>%
    ggplot() + geom_step(aes(x = t,y =value,col = variable)) +
    geom_line(data = plotdf2b,aes(x = t,y =value,col = variable), linetype="dashed")+
    ylim(0,1) + labs(col = "State",y=TeX("$P(Z_t=j)$")) + theme_custom()
  return(list(
    hazards = hazards, plots = list(p1 = p1,p2 = p2)
  ))
}
plotPureJump <- function(path) {
  path <- path %>% as.data.frame()
  path[,"tend"] <- c(path$times[-1],path$times[length(path$times)])
  ggplot(path) + geom_segment(aes(x = times, xend = tend,y = states, yend = states)) +
    geom_point(aes(x=times,y=states)) + geom_point(aes(x=tend,y=states), shape =21) +
    theme_custom() + labs(x="t",y=TeX("$Z_t$"))
}
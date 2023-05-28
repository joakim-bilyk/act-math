#### Appendix A.1 ####
M <- t(matrix(c(-3,2,1,
              3,-4,1,
              0,0,0),ncol=3,nrow=3))
lambda <- function(t){
  1/(1+0.5*t)
}
d_Lambda <- function(t,M,lambda){
  lambda(t)*M
}
jump_rate <- function(i,t,u){
  #the variable u is not used, it will be used later for time sojourned in the current state
  d_L_t <- d_Lambda(t,M,lambda)
  J <- dim(d_L_t)[1]
  vec <- (1:J==i)*1
  -vec%*%d_L_t%*%vec
}
mark_dist <- function(i,s,v){
  #the variable v is not used
  d_L_t <- d_Lambda(s,M,lambda)
  J <- dim(d_L_t)[1]
  vec <- (1:J==i)*1
  tmp <- (vec %*% d_L_t)* (1-vec)
  tmp / sum(tmp)
}

L <- 1000
set.seed(1)
R <- runif(L,0,10)
#Simulate paths
paths <- lapply(1:L,function(n) {
  sim_path(1,rates = jump_rate, dist = mark_dist,
         tn = R[n], bs= c(R[n]*3,R[n]*4,0))})
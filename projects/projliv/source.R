library(AalenJohansen)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(latex2exp)
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##### CODE FROM FURRER #####

# 1. Markov model with independent censoring
set.seed(2)

jump_rate <- function(i, t, u){
  if(i == 1){
    2 / (1+1/2*t)
  } else if(i == 2){
    3 / (1+1/2*t)
  } else{
    0
  }
}

mark_dist <- function(i, s, v){
  if(i == 1){
    c(0, 1/2, 1/2)
  } else if(i == 2){
    c(2/3, 0, 1/3)
  } else{
    0
  } 
}

lambda <- function(t){
  A <- matrix(c(2/(1+1/2*t)*mark_dist(1, t, 0), 3/(1+1/2*t)*mark_dist(2, t, 0), rep(0, 3)),
              nrow = 3, ncol = 3, byrow = TRUE)
  diag(A) <- -rowSums(A)
  A
}

n <- 1000
c <- runif(n, 0, 5)

sim <- list()
for(i in 1:n){
  sim[[i]] <- sim_path(sample(1:2, 1), rates = jump_rate, dists = mark_dist,
                       tn = c[i], bs = c(2*c[i], 3*c[i], 0))
}

sum(c == unlist(lapply(sim, FUN = function(z){tail(z$times, 1)}))) / n

fit <- aalen_johansen(sim)

v1 <- unlist(lapply(fit$Lambda, FUN = function(L) L[2,1]))
v0 <- fit$t
p <- unlist(lapply(fit$p, FUN = function(L) L[2]))
P <- unlist(lapply(prodint(0, 5, 0.01, lambda), FUN = function(L) (c(1/2, 1/2, 0) %*% L)[2]))

par(mfrow = c(1, 2))
par(mar = c(2.5, 2.5, 1.5, 1.5))

plot(v0, v1, type = "l", lty = 2, xlab = "", ylab = "", main = "Hazard")
lines(v0, 4*log(1+1/2*v0))
plot(v0, p, type = "l", lty = 2, xlab = "", ylab = "", main = "Probability")
lines(seq(0, 5, 0.01), P)

# 2. Markov model with independent censoring and covariates

#### END OF CODE FROM FURRER #####

rm(list = ls())

#### Appendix A.1 ####
M <- t(matrix(c(-3.5,2,1.5,
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

#Simulate paths
simulate_markov_inhomogenous <- function(L) {
  R <- runif(L,0,10)
  paths <- lapply(1:L,function(n) {
    sim_path(1,rates = jump_rate, dist = mark_dist,
             tn = R[n], bs= c(R[n]*3,R[n]*4,0))})
}
L <- 10000
set.seed(1)
paths <- simulate_markov_inhomogenous(L)

#Tests
jump_rate(2,20,1)
mark_dist(1,20,1)

#### Appendix A.2 ####
#### Appendix A.2.1 ####
#Convert path data to main_df
paths_to_df <- function(paths){
  L <- length(paths)
  times <- unlist(lapply(1:L, function(i){
    paths[[i]]$times
  }))
  states <- unlist(lapply(1:L, function(i){
    paths[[i]]$states
  }))
  obs <- unlist(lapply(1:L, function(i){
    rep(i,length(paths[[i]]$times))
  }))
  df <- data.frame(Id = obs, Start_Time = times, Start_State = states)
  #End time & end state
  df <- df %>%
    arrange(Start_Time) %>%
    group_by(Id) %>%
    mutate(End_Time = data.table::shift(Start_Time,-1),
           End_State = data.table::shift(Start_State,-1)) %>%
    ungroup() %>%
    replace(is.na(.), Inf) %>%
    arrange(Id, Start_Time) %>%
    mutate(End_State = ifelse(is.infinite(End_State),Start_State,End_State),
           Censored = ifelse(Start_State == End_State,
                             ifelse(is.finite(End_Time),TRUE,FALSE),FALSE)) %>%
    group_by(Id) %>%
    mutate(Censored = ifelse((!Censored) & is.infinite(End_Time) & (cumsum(Censored) >0),TRUE,Censored)) %>%
    ungroup()
  df <- df[,c("Id","Start_Time","End_Time","Start_State","End_State","Censored")]
  return(df)
}
main_df <- paths_to_df(paths)
#### Appendix A.2.2 ####
#Convert main_df to I
df_to_I <- function(df,num_states) {
  Init <- df %>%
    group_by(State = Start_State) %>%
    summarise(Time = 0,
              Change = sum(Start_Time == 0))
  I <- suppressMessages(df %>%
    filter(End_Time < Inf) %>%
    mutate(Count_from = TRUE) %>%
    bind_rows(filter(.,Censored == FALSE) %>% mutate(Count_from = FALSE)) %>%
    arrange(Id,End_Time) %>%
    mutate(State = ifelse(Count_from, Start_State,End_State)) %>%
    group_by(Time = End_Time, State) %>%
    summarise(Change = sum((!Count_from)*1)-sum(Count_from*1)) %>%
    ungroup() %>%
    bind_rows(Init) %>%
    arrange(Time) %>%
    group_by(State) %>%
    mutate(I_j = cumsum(Change)) %>%
    reshape2::dcast(Time ~ State, value.var = "I_j") %>%
    zoo::na.locf(na.rm = FALSE) %>%
    replace(is.na(.), 0))
  if ( sum(!(1:num_states %in% colnames(I))) >0) {
    I[,1:num_states[!(1:num_states %in% colnames(I))]] <- 0
  }
  I <- I[,c("Time",1:num_states)]
  I_left_limit <- I
  I_left_limit[2:(dim(I)[1]-1),2:dim(I)[2]] <- I[3:dim(I)[1],2:dim(I)[2]] 
  I_left_limit[dim(I)[1],2:dim(I)[2]] <- I[(dim(I)[1]-1),2:dim(I)[2]]
  return(list(I = I, I_left_limit = I_left_limit))
}
I_list <- df_to_I(main_df,3)
#### Appendix A.2.3 ####
#Convert main_df to N
df_to_N <- function(df, num_states) {
  N <- df %>%
    filter((Censored == FALSE) & (End_Time < Inf)) %>%
    arrange(End_Time) %>%
    group_by(Start_State, End_State) %>%
    mutate(Transitions = cumsum(End_State >= 0),
           #j_k = paste0(Start_State,",", End_State)
    ) %>%
    select(Time = End_Time,Start_State, End_State,Transitions) %>%
    bind_rows(filter(., Transitions == 1) %>% mutate(Time = 0, Transitions = 0)) %>%
    reshape2::dcast(Time ~ Start_State * End_State, value.var = "Transitions") %>%
    zoo::na.locf(na.rm = FALSE) %>%
    replace(is.na(.), 0)
  c <- paste0(unlist(lapply( 1:num_states, function(i) rep(i,num_states))),"_",rep(1:num_states,num_states))
  if ( sum(!(c %in% colnames(N))) >0) {
    N[,c[!(c %in% colnames(N))]] <- 0
  }
  N <- N[,c("Time",c)]
  for (i in 1:num_states) {
    target <- paste0(i,"_",(1:num_states)[!(1:num_states %in% i)])
    N[,paste0(i,"_",i)] <- -rowSums(N[,target])
  }
  delta_N <- N
  delta_N[2:dim(N)[1],2:dim(N)[2]] <- delta_N[2:dim(N)[1],2:dim(N)[2]] - delta_N[1:(dim(N)[1]-1),2:dim(N)[2]]
  return(list(N = N, delta_N = delta_N))
}
N_list <- df_to_N(main_df, 3)
#### Appendix A.2.4 ####
#Calculate Nelson-Aalen
N_I_to_NA <- function(I_list,N_list, num_states) {
  delta_N <- N_list$delta_N
  I_left_limit <- I_list$I_left_limit
  I_left_limit <- I_left_limit[I_left_limit$Time %in% delta_N$Time,]
  I_factor <- delta_N
  for (i in 1:num_states) {
    target <- paste0(i,"_",1:num_states)
    vec <- 1/I_left_limit[,colnames(I_left_limit) == i]
    vec <- ifelse(is.infinite(vec),0,vec)
    I_factor[,target] <- vec
  }
  delta_NelsonAalen <- delta_N
  delta_NelsonAalen[,2:dim(delta_N)[2]] <- delta_NelsonAalen[,2:dim(delta_N)[2]]*I_factor[,2:dim(delta_N)[2]]
  NelsonAalen <- delta_NelsonAalen
  NelsonAalen[,2:dim(delta_N)[2]] <- cumsum(delta_NelsonAalen[,2:dim(delta_N)[2]])
  return(list(NelsonAalen = NelsonAalen,delta_NelsonAalen=delta_NelsonAalen))
}
NelsonAalen_list <- N_I_to_NA(I_list,N_list,3)
#### Appendix A.2.5 ####
#Calculate product integral
NA_to_p <- function(I_list,N_list,NelsonAalen_list, num_states) {
  #This may be slow
  identity <- diag(as.numeric(1), ncol = num_states,nrow=num_states)
  Nelson <- NelsonAalen_list$delta_NelsonAalen
  Nelson <- lapply(1:dim(Nelson)[1], function(i) matrix(Nelson[i,2:dim(Nelson)[2]],nrow=num_states,ncol = num_states,byrow=TRUE))
  Prod_int <- list()
  Prod_int[[1]] <- identity
  for (i in 2:length(Nelson)) {
    Prod_int[[i]] <- Prod_int[[i-1]] %*% (identity + as.numeric(Nelson[[i]]))
  }
  P <- NelsonAalen_list$delta_NelsonAalen
  P[,2:dim(P)[2]] <- matrix(unlist(Prod_int),ncol = num_states**2, byrow = TRUE)
  return(P)
}
transitionprobs <- NA_to_p(I_list,N_list,NelsonAalen_list, 3)
#### Appendix A.2.6 ####
#Conditional probabilities
P_conditioned <- function(NelsonAalen_list,s,j,num_states) {
  Init <- (1:num_states==j)*1
  Nelson <- NelsonAalen_list$delta_NelsonAalen %>%
    filter(Time >= s)
  Nelson <- lapply(1:dim(Nelson)[1], function(i) matrix(Nelson[i,2:dim(Nelson)[2]],nrow=num_states,ncol = num_states,byrow=TRUE))
  Prod_int <- list()
  Prod_int[[1]] <- Init
  identity <- diag(as.numeric(1), ncol = num_states,nrow=num_states)
  for (i in 1:length(Nelson)) {
    Prod_int[[i+1]] <- Prod_int[[i]] %*% (identity + as.numeric(Nelson[[i]]))
  }
  p_con <- NelsonAalen_list$delta_NelsonAalen %>%
    filter(Time >= s) %>%
    bind_rows(filter(., row_number()==1) %>% mutate(Time = s)) %>%
    arrange(Time)
  p_con <- p_con[,1:(num_states + 1)]
  colnames(p_con) <- c("Time",1:num_states)
  p_con[,2:dim(p_con)[2]] <- matrix(unlist(Prod_int),ncol = num_states, byrow = TRUE)
  return(p_con)
}
#### Appendix A.2.7 ####
#Putting it all together with as-if variant
Estimate <- function(paths,num_states,s= NA, j = NA, as_if = FALSE,debug = TRUE) {
  start_time <- Sys.time()
  # 1. Start by converting to data frame
  if (is.na(s)) {
    s <- 0
    j <- 1
  }
  main_df_tmp <- paths_to_df(paths)
  if(debug) {
    print(paste0("Convert paths to main_df: ",round(as.numeric(Sys.time()-start_time,units = "secs"),digits=3)," seconds."))
    start_time <- Sys.time()
  }
  if (as_if) {
    tmp <- main_df_tmp %>%
      filter(Start_Time <= s) %>%
      group_by(Id) %>%
      mutate(max_Time = max(Start_Time)) %>%
      filter(max_Time == Start_Time) %>%
      filter(((is.finite(End_Time) | (Censored == FALSE)))) %>%
      filter(Start_State == j)
    main_df_tmp <- main_df_tmp %>% filter(Id %in% tmp$Id)
  }
  # 2. Calculate I and N
  I_list_tmp <- df_to_I(main_df_tmp,num_states)
  if (debug) {
    print(paste0("Convert main_df to I: ",round(as.numeric(Sys.time()-start_time,units = "secs"),digits=3)," seconds."))
    start_time <- Sys.time()
  }
  N_list_tmp <- df_to_N(main_df_tmp,num_states)
  if (debug) {
    print(paste0("Convert main_df to N: ",round(as.numeric(Sys.time()-start_time,units = "secs"),digits=3)," seconds."))
    start_time <- Sys.time()
  }
  # 3. Calculate Nelson-Aalen
  NelsonAalen_list_tmp <- N_I_to_NA(I_list_tmp,N_list_tmp,num_states)
  if (debug) {
    print(paste0("Convert I and N to Nelson-Aalen: ",round(as.numeric(Sys.time()-start_time,units = "secs"),digits=3)," seconds."))
    start_time <- Sys.time()
  }
  # 4. Calculate Aalen-Johansen
  P_tmp <- NA_to_p(I_list_tmp,N_list_tmp,NelsonAalen_list_tmp, 3)
  if (debug) {
    print(paste0("Convert Nelson-Aalen to transition probs: ",round(as.numeric(Sys.time()-start_time,units = "secs"),digits=3)," seconds."))
    start_time <- Sys.time()
  }
  # 5. Calculate probabilties
  p_con_tmp <- P_conditioned(NelsonAalen_list_tmp,s,j,num_states)
  if (debug) {
    print(paste0("Convert transition probs to conditioned probs: ",round(as.numeric(Sys.time()-start_time,units = "secs"),digits=3)," seconds."))
    start_time <- Sys.time()
  }
  return(list(I = I_list_tmp$I, I_left_limit = I_list_tmp$I_left_limit,
              N = N_list_tmp$N,delta_N = N_list_tmp$delta_N,
              NelsonAalen = NelsonAalen_list_tmp$NelsonAalen, delta_NelsonAalen= NelsonAalen_list_tmp$delta_NelsonAalen,
              P = P_tmp, p_con = p_con_tmp))
}
total <- Estimate(paths,3)
total <- Estimate(paths,3,s=1,j=1,as_if = TRUE)

#### Appendix A.3 ####
num_states <- 3
s <- 1
j <- 1
debug <- TRUE
pi <- 1
T <- 3
estimate_to_cashflow <- function(estimate,pi,T,n_steps=10000) {
  p_con <- estimate$p_con
  s <- min(p_con$Time)
  t_0 <- p_con[-dim(p_con)[1],1]
  t_1 <- p_con[-1,1]
  t <- p_con$Time
  p_con <- p_con[,2:dim(p_con)[2]]
  
  Nelson <- estimate$delta_NelsonAalen %>%
    filter(Time >= s)
  
  A1_increments <-  p_con[-1,1]*ifelse(t_1>T,t_1-pmax(T,t_0),0) #A_1(t)-A_1(s)=p(Z_t=1)*(t-s)
  A1 <- c(0,cumsum(A1_increments))
  A2_increments <- - pi * p_con[-1,1] * ifelse(t_0<T,pmin(T,t_1)-t_0,0) #A_2(t)-A_2(s)=-pi*p(Z_t=1)*(t-s)
  A2 <- c(0,cumsum(A2_increments))
  A3_increments <- p_con[-1,2] * (t_1-t_0)
  A3 <- c(0,cumsum(A3_increments))
  A4_increments <- p_con[-1,1] * Nelson[,"1_3"] + p_con[-1,2] * Nelson[,"2_3"]
  A4 <- c(0,cumsum(A4_increments))
  A <- data.frame(
    Time = t,
    A1 = A1, A2 = A2, A3 = A3, A4 = A4
  ) %>% mutate(A = A1 + A2 + A3 + A4)
  return(A)
}
estimate <- Estimate(paths,3)
cashflow <- estimate_to_cashflow(estimate,1,3)
#Runge-Kutta for true probabilities
runge_kutta <- function(f,a,b,y0,n) {
  y <- list()
  y[[1]] <- y0
  t <- list()
  t[[1]] <- a
  h <- (b-a)/n
  t0 <- a
  for (i in 1:n) {
    k1 <- f(t0,y0)
    k2 <- f(t0 + h/2, y0 + (h/2)*k1)
    k3 <- f(t0 + h/2, y0 + (h/2)*k2)
    k4 <- f(t0 + h, y0 + h*k3)
    y1 <- y0+ (h/6)*(k1+2*k2+2*k3+k4)
    y[[i+1]] <- y1
    t[[i+1]] <- t0+h
    t0 <- t0+h
    y0 <- y1
  }
  return(list(t=t,y=y))
}
n_steps <- 10000
plot_function <- function(paths,pi,T,num_states,s= 0, j = 1, debug = TRUE,n_steps=10000) {
  #Markov estimate
  estimate1 <- Estimate(paths,num_states,s= s, j = j, as_if = FALSE, debug = debug)
  start_time <- Sys.time()
  cashflow1 <- estimate_to_cashflow(estimate1,pi,T)
  if (debug) {
    print(paste0("Calculate cashflows: ",round(as.numeric(Sys.time()-start_time,units = "secs"),digits=3)," seconds."))
    start_time <- Sys.time()
  }
  #As-If-Markov estimate
  estimate2 <- Estimate(paths,num_states,s= s, j = j, as_if = TRUE, debug = debug)
  start_time <- Sys.time()
  cashflow2 <- estimate_to_cashflow(estimate2,pi,T)
  if (debug) {
    print(paste0("Calculate cashflows: ",round(as.numeric(Sys.time()-start_time,units = "secs"),digits=3)," seconds."))
    start_time <- Sys.time()
  }
  #Plots of probabilities
  plotdf <- estimate1$p_con
  colnames(plotdf)[2:4] <- paste0("j=",1:3,", Markov")
  plotdf <- plotdf%>% reshape2::melt(., id = "Time")
  plotdf2 <- estimate2$p_con
  colnames(plotdf2)[2:4] <- paste0("j=",1:3,", As-If")
  plotdf2 <- plotdf2%>% reshape2::melt(., id = "Time")
  #True values (approximated with 4th order runge kutta)
  derivative <- function(t,y) {
    y %*% d_Lambda(t,M,lambda)
  }
  y <- runge_kutta(derivative, s,10,(1:num_states == j)*1,n_steps)
  plotdf3 <- data.frame(Time = unlist(y$t),
                        matrix(unlist(y$y),ncol=3,byrow=TRUE))
  colnames(plotdf3)[2:4] <- paste0("j=",1:3,", True")
  plotdf3 <- plotdf3 %>% reshape2::melt(., id = "Time")
  asif_n <- sum(estimate2$I[1,])-estimate2$I[1,1]
  markov_n <- sum(estimate1$I[1,])-estimate1$I[1,1]
  p1 <- ggplot() +
    geom_step(data = plotdf ,mapping = aes(x=Time, y = value,col = variable)) + 
    geom_step(data = plotdf2,mapping = aes(x=Time, y = value,col = variable),linetype = "dashed") +
    geom_line(data = plotdf3,mapping = aes(x=Time, y = value,col = variable),linetype = "dotted",linewidth=1) +
    xlim(0,10) +
    theme_bw() +
    labs(title = TeX(paste0("Occupation probabilities under assumption $\\psi(Z_",s,")=",j,"$ (G)")),
         y =  TeX(paste0("$P(Z_t=j|G)$")),
         x = "t",
         subtitle = paste0("As-if estimate based on ",asif_n," observations,\nMarkov based on ",markov_n," observations")) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic")) +
    scale_color_manual(values=c("#FA8072", "#FF0000","#8B0000",
                                "#7B68EE", "#1E90FF","#00008B",
                                "#3CB371", "#32CD32","#006400"),
                       name = "Probability, Estimate",
                       labels = c(TeX("$P(Z_t=a|\\cdot)$, As-If"),
                                  TeX("$P(Z_t=a|\\cdot)$, Markov"),
                                  TeX("$P(Z_t=a|\\cdot)$, True"),
                                  TeX("$P(Z_t=b|\\cdot)$, As-If"),
                                  TeX("$P(Z_t=b|\\cdot)$, Markov"),
                                  TeX("$P(Z_t=b|\\cdot)$, True"),
                                  TeX("$P(Z_t=c|\\cdot)$, As-If"),
                                  TeX("$P(Z_t=c|\\cdot)$, Markov"),
                                  TeX("$P(Z_t=c|\\cdot)$, True")))
  #Plots of intensities
  plotdf1 <- estimate1$NelsonAalen %>% 
    select(Time,
           `Lambda(1,2), Markov`=`1_2`,
           `Lambda(1,3), Markov`=`1_3`,
           `Lambda(2,1), Markov`=`2_1`,
           `Lambda(2,3), Markov`=`2_3`) %>%
    filter(Time >= s)
  plotdf1[,2:dim(plotdf1)[2]] <- plotdf1[,2:dim(plotdf1)[2]] - as.data.frame(matrix(as.numeric(rep(plotdf1[1,2:dim(plotdf1)[2]],dim(plotdf1)[1])),ncol = dim(plotdf1)[2]-1,byrow = TRUE))
  plotdf1 <- plotdf1 %>%
    reshape2::melt(., id = "Time")
  plotdf2 <- estimate2$NelsonAalen %>% 
    select(Time,
           `Lambda(1,2), As-if`=`1_2`,
           `Lambda(1,3), As-if`=`1_3`,
           `Lambda(2,1), As-if`=`2_1`,
           `Lambda(2,3), As-if`=`2_3`) %>%
    filter(Time >= s)
  plotdf2[,2:dim(plotdf2)[2]] <- plotdf2[,2:dim(plotdf2)[2]] - as.data.frame(matrix(as.numeric(rep(plotdf2[1,2:dim(plotdf2)[2]],dim(plotdf2)[1])),ncol = dim(plotdf2)[2]-1,byrow = TRUE))
  plotdf2 <- plotdf2 %>%
    reshape2::melt(., id = "Time")
  times <- s+0:n_steps*(10-s)/n_steps
  plotdf3 <- data.frame(Time = times,
                        matrix(unlist(lapply(times, function(t) as.numeric(t(2*M*(log(1+0.5*t)-log(1+0.5*s)))))),ncol = num_states**2,byrow = TRUE))
  colnames(plotdf3) <- c("Time",unlist(lapply(1:num_states, function(i) paste0(i,"_",1:num_states))))
  plotdf3 <- plotdf3 %>% 
    select(Time,
           `Lambda(1,2)`=`1_2`,
           `Lambda(1,3)`=`1_3`,
           `Lambda(2,1)`=`2_1`,
           `Lambda(2,3)`=`2_3`) %>%
    reshape2::melt(., id = "Time")
  p2 <- ggplot() +
    geom_step(data = plotdf1 ,mapping = aes(x=Time, y = value,col = variable)) + 
    geom_step(data = plotdf2,mapping = aes(x=Time, y = value,col = variable),linetype = "dashed") +
    geom_line(data = plotdf3,mapping = aes(x=Time, y = value,col = variable),linetype = "dotted",linewidth=1) +
    xlim(0,10) +
    theme_bw() +
    labs(title = TeX(paste0("Cumulative transition rates under assumption $\\psi(Z_",s,")=",j,"$ (G)")),
         y =  TeX(paste0("$\\Lambda_{jk}(t|G)-\\Lambda_{jk}(s|G)$")),
         x = "t",
         subtitle = paste0("As-if estimate based on ",asif_n," observations,\nMarkov based on ",markov_n," observations")) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic")) +
    scale_color_manual(values=c("#FA8072", "#FF0000","#8B0000","#7B68EE",
                                "#1E90FF","#00008B","#3CB371", "#32CD32",
                                "#006400","#D2691E","#F4A460","#DEB887"),
                       name = "j to k, Estimate",
                       labels = c("j=a, k=b, True",
                                  "j=a, k=b, As-If",
                                  "j=a, k=b, Markov",
                                  "j=a, k=c, True",
                                  "j=a, k=c, As-If",
                                  "j=a, k=c, Markov",
                                  "j=b, k=a, True",
                                  "j=b, k=a, As-If",
                                  "j=b, k=a, Markov",
                                  "j=b, k=c, True",
                                  "j=b, k=c, As-If",
                                  "j=b, k=c, Markov"))
  #Plots of cashflows
  colnames(cashflow1)[2:dim(cashflow1)[2]] <- paste0(colnames(cashflow1)[2:dim(cashflow1)[2]],", Markov")
  plotdf1 <- cashflow1 %>% reshape2::melt(id = "Time")
  colnames(cashflow2)[2:dim(cashflow2)[2]] <- paste0(colnames(cashflow2)[2:dim(cashflow2)[2]],", As-If")
  plotdf2 <- cashflow2 %>% reshape2::melt(id = "Time")
  #Calculate true values
  cashflow3 <- cashflow1
  prob_to_1 <- approxfun(y$t,y=unlist(y$y)[1:length(y$t)*length(y$y[[1]])-2])
  prob_to_2 <- approxfun(y$t,y=unlist(y$y)[1:length(y$t)*length(y$y[[1]])-1])
  prob_to_3 <- approxfun(y$t,y=unlist(y$y)[1:length(y$t)*length(y$y[[1]])-0])
  A1_derivative <- function(t,y){
    (t>T)*prob_to_1(t)
  }
  A1 <- runge_kutta(A1_derivative,s,10,0,n_steps)
  A2_derivative <- function(t,y){
    -pi*(t<=T)*prob_to_1(t)
  }
  A2 <- runge_kutta(A2_derivative,s,10,0,n_steps)
  A3_derivative <- function(t,y){
    prob_to_2(t)
  }
  A3 <- runge_kutta(A3_derivative,s,10,0,n_steps)
  A4_derivative <- function(t,y){
    #prob_to_1(t-)=prob_to_1(t) same for prob_to_2
    prob_to_1(t)*d_Lambda(t,M,lambda)[1,3]+prob_to_2(t)*d_Lambda(t,M,lambda)[2,3]
  }
  A4 <- runge_kutta(A4_derivative,s,10,0,n_steps)
  cashflow3 <- data.frame(
    Time = unlist(A1$t),
    A1 = unlist(A1$y),
    A2 = unlist(A2$y),
    A3 = unlist(A3$y),
    A4 = unlist(A4$y)
  ) %>%
    mutate(A = A1 + A2 + A3 + A4)
  colnames(cashflow3) <- c("Time","A1_true","A2_true","A3_true","A4_true","A_true")
  plotdf3 <- cashflow3 %>% reshape2::melt(id = "Time")
  if (debug) {
    print(paste0("Calculate theoretical cashflows and more: ",round(as.numeric(Sys.time()-start_time,units = "secs"),digits=3)," seconds."))
    start_time <- Sys.time()
  }
  p3 <- ggplot() +
    geom_step(data = plotdf1 ,mapping = aes(x=Time, y = value,col = variable)) +
    geom_step(data = plotdf2,mapping = aes(x=Time, y = value,col = variable), linetype = "dashed") +
    geom_line(data = plotdf3,mapping = aes(x=Time, y = value,col = variable), linetype = "dotted", linewidth=1) +
    geom_vline(xintercept = T, col = "black",linetype = "dashed") +
    xlim(0,10) +
    theme_bw() +
    labs(title = TeX(paste0("Accumulated cash-flow under assumption $\\psi(Z_",s,")=",j,"$ (G)")),
         y =  TeX(paste0("A(t|G)-A(s|G)")),
         x = "t",
         subtitle = paste0("As-if estimate based on ",asif_n," observations,\nMarkov based on ",markov_n," observations")) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic")) +
    scale_color_manual(values=c("#FA8072", "#FF0000","#8B0000","#7B68EE",
                                "#1E90FF","#00008B","#3CB371", "#32CD32",
                                "#006400","#D2691E","#F4A460","#DEB887",
                                "#FF1493","#C71585","#FFB6C1"),
                       name = "Cash-flow, Estimate",
                       labels = c(TeX("Total ($A$), Markov"),
                                  TeX("Total ($A$), As-If"),
                                  TeX("Total ($A$), True"),
                                  TeX("Pension ($A_1$), Markov"),
                                  TeX("Pension ($A_1$), As-If"),
                                  TeX("Pension ($A_1$), True"),
                                  TeX("Premium ($A_2$), Markov"),
                                  TeX("Premium ($A_2$), As-If"),
                                  TeX("Premium ($A_2$), True"),
                                  TeX("Disabled ($A_3$), Markov"),
                                  TeX("Disabled ($A_3$), As-If"),
                                  TeX("Disabled ($A_3$), True"),
                                  TeX("Death ($A_4$), Markov"),
                                  TeX("Death ($A_4$), As-If"),
                                  TeX("Death ($A_4$), True")))
  return(list(p1= p1, p2 = p2, p3 = p3))
}
plot1 <- plot_function(paths,1,3,3,s= 0, j = 1)
plot2 <- plot_function(paths,1,3,3,s= 1, j = 1)
plot3 <- plot_function(paths,1,3,3,s= 3, j = 1)
plot4 <- plot_function(paths,1,3,3,s= 6, j = 1)
p1 <- ggarrange(plotlist = list(plot1$p1,plot2$p1,plot3$p1,plot4$p1),ncol = 2,nrow=2)
p2 <- ggarrange(plotlist = list(plot1$p2,plot2$p2,plot3$p2,plot4$p2),ncol = 2,nrow=2)
p3 <- ggarrange(plotlist = list(plot1$p3,plot2$p3,plot3$p3,plot4$p3),ncol = 2,nrow=2)
scaler <- 1.75
ggsave("plot1.png",p1,units = "px", width = 1920*scaler,height = 1080*scaler,scale = 1.5)
ggsave("plot2.png",p2,units = "px", width = 1920*scaler,height = 1080*scaler,scale = 1.5)
ggsave("plot3.png",p3,units = "px", width = 1920*scaler,height = 1080*scaler,scale = 1.5)

#### Appendix A.4 ####
L <- ceil(2*65**(2*0:100/100+1))
S <- 6
J <- 1
supremums <- function(L,S,J, debug = TRUE,n_steps = 10000) {
  
  results <- data.frame(L = NA,s=NA,j=NA,as_if = TRUE,`Lambda(1,2)`=NA,`Lambda(1,3)`=NA,`Lambda(2,1)`=NA,`Lambda(2,3)`=NA)
  L <- as.integer(L)
  counter <- 1
  set.seed(1)
  if (debug) {
    print(paste0(Sys.time()," Starting simulating ",max(L)," sample paths."))
  }
  paths <- simulate_markov_inhomogenous(max(L))
  if (debug) {
    print(paste0(Sys.time()," Done simulating ",max(L)," sample paths."))
  }
  for (l in L) {
    tmp_paths <- paths[1:l]
    for (s in S) {
      for (j in J) {
        #Generate paths
        estimate1 <- tryCatch(Estimate(tmp_paths,3, s= s, j=j,as_if = FALSE, debug = FALSE),
                              error = function(e) NA)
        estimate2 <- tryCatch(Estimate(tmp_paths,3,s= s, j=j,as_if = TRUE, debug = FALSE),
                              error = function(e) NA)
        if ((length(estimate1)>1) & (length(estimate2)>1)) {
          markov <- estimate1$NelsonAalen %>% 
            select(Time,
                   `Lambda(1,2)`=`1_2`,
                   `Lambda(1,3)`=`1_3`,
                   `Lambda(2,1)`=`2_1`,
                   `Lambda(2,3)`=`2_3`) %>%
            filter(Time >= s)
          markov[,2:dim(markov)[2]] <- markov[,2:dim(markov)[2]] -
            as.data.frame(matrix(as.numeric(rep(markov[1,2:dim(markov)[2]],dim(markov)[1])),ncol = dim(markov)[2]-1,byrow = TRUE))
          as_if_markov <- estimate2$NelsonAalen %>% 
            select(Time,
                   `Lambda(1,2)`=`1_2`,
                   `Lambda(1,3)`=`1_3`,
                   `Lambda(2,1)`=`2_1`,
                   `Lambda(2,3)`=`2_3`) %>%
            filter(Time >= s)
          as_if_markov[,2:dim(as_if_markov)[2]] <- as_if_markov[,2:dim(as_if_markov)[2]] -
            as.data.frame(matrix(as.numeric(rep(as_if_markov[1,2:dim(as_if_markov)[2]],dim(as_if_markov)[1])),ncol = dim(as_if_markov)[2]-1,byrow = TRUE))
          times <- sort(unique(c(s+0:n_steps*(10-s)/n_steps,as_if_markov$Time,markov$Time)),decreasing = FALSE)
          true_values <- data.frame(Time = times,
                                    matrix(unlist(lapply(times, function(t) as.numeric(t(2*M*(log(1+0.5*t)-log(1+0.5*s)))))),ncol = num_states**2,byrow = TRUE))
          colnames(true_values) <- c("Time",unlist(lapply(1:num_states, function(i) paste0(i,"_",1:num_states))))
          true_values <- true_values %>% 
            select(Time,
                   `Lambda(1,2),true`=`1_2`,
                   `Lambda(1,3),true`=`1_3`,
                   `Lambda(2,1),true`=`2_1`,
                   `Lambda(2,3),true`=`2_3`)
          markov <- merge(markov,true_values,all.x = TRUE)
          as_if_markov <- merge(as_if_markov,true_values,all.x = TRUE)
          
          results[counter,] <- c(
            l,s,j,FALSE,
            max(abs(markov$`Lambda(1,2)`-markov$`Lambda(1,2),true`)),
            max(abs(markov$`Lambda(1,3)`-markov$`Lambda(1,3),true`)),
            max(abs(markov$`Lambda(2,1)`-markov$`Lambda(2,1),true`)),
            max(abs(markov$`Lambda(2,3)`-markov$`Lambda(2,3),true`))
          )
          results[counter+1,] <- c(
            l,s,j,TRUE,
            max(abs(as_if_markov$`Lambda(1,2)`-as_if_markov$`Lambda(1,2),true`)),
            max(abs(as_if_markov$`Lambda(1,3)`-as_if_markov$`Lambda(1,3),true`)),
            max(abs(as_if_markov$`Lambda(2,1)`-as_if_markov$`Lambda(2,1),true`)),
            max(abs(as_if_markov$`Lambda(2,3)`-as_if_markov$`Lambda(2,3),true`))
          )
          counter <- counter + 2
          if (debug == TRUE) {
            print(paste0(Sys.time(),": Done with L=",l,", s=",s," and j=",j,"."))
          }
        }
      }
    }
    
  }
  return(results)
}
plot_function <- function(results, s, j) {
  
  plotdf <- results %>% filter((s==s) & (j==j)) %>% filter(as_if == TRUE) %>%
    select(-s,-j,-as_if) %>%
    melt(id = "L")
  plotdf2 <- results %>% filter((s==s) & (j==j)) %>% filter(as_if == FALSE) %>%
    select(-s,-j,-as_if) %>%
    melt(id = "L")
  
  ggplot() + geom_point(data = plotdf %>% filter(variable == "Lambda.1.3."), aes(x=L,y = value, col = variable))+
    geom_line(data = plotdf2 %>% filter(variable == "Lambda.1.3."), aes(x=L,y = value, col = variable)) +
    scale_x_continuous(trans = "log2") + scale_y_continuous(trans = "log2")
  
}

results1 <- supremums(L = 1:10*1000,S=0:6,J=1) #run-time ≈ 
results2 <- supremums(L = ceil(2*65**(2*0:25/25+1)),S=6,J=1,n_steps = 1000000) #run-time ≈ 11.5 min
results3 <- supremums(L = ceil(2*88**(2*0:25/25+1)),S=6,J=2) #run-time ≈ 

plot1 <- plot_function(results,6,1)
plot1
#Calcuæate P(Z_s=j)P(R>s)
s <- 6
c(1,0,0)%*%expm::expm(M*2*log(1+0.5*s))

##################################
###### Estimate Semi-Markov ######
##################################
#### Appendix A.5 ####
jump_rate <- function(i, t, u){
  if(i == 1){
    0.1 + 0.002*t
  } else if(i == 2){
    ifelse(u < 4, 0.29, 0.09) + 0.001*t
  } else{
    0
  }
}

mark_dist <- function(i, s, v){
  if(i == 1){
    c(0, 0.9, 0.1)
  } else if(i == 2){
    c(0, 0, 1)
  } else{
    0
  }
}
c <- runif(n, 10, 40)



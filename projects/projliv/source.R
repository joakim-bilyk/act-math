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
end_time <- Sys.time()
round(as.numeric(end_time-start_time,units = "secs"),digits = 3)
#Tests
jump_rate(2,20,1)
mark_dist(1,20,1)

#### Appendix A.2 ####
#Convert path data to main_df
paths_to_df <- function(paths){
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
#Calculate p(0,t) i.e. Aalen-Johansen
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
plot_function1 <- function(paths,num_states,s= 0, j = 1, debug = TRUE) {
  #Markov estimate
  total <- Estimate(paths,num_states,s= s, j = j, as_if = FALSE, debug = debug)
  #As-If-Markov estimate
  total2 <- Estimate(paths,num_states,s= s, j = j, as_if = TRUE, debug = debug)
  plotdf <- total$p_con
  colnames(plotdf)[2:4] <- paste0("j=",1:3,", Markov")
  plotdf <- plotdf%>% reshape2::melt(., id = "Time")
  plotdf2 <- total2$p_con
  colnames(plotdf2)[2:4] <- paste0("j=",1:3,", As-If")
  plotdf2 <- plotdf2%>% reshape2::melt(., id = "Time")
  times <- s+0:1000*(10-s)/1000
  plotdf3 <- data.frame(Time = times,
                        matrix(unlist(lapply(times, function(t) ((1:num_states == j)*1)%*%expm::expm(2*M*(log(1+0.5*t)-log(1+0.5*s))))),ncol=3,byrow=TRUE))
  colnames(plotdf3)[2:4] <- paste0("j=",1:3,", True")
  plotdf3 <- plotdf3 %>% reshape2::melt(., id = "Time")
  asif_n <- sum(total2$I[1,])-total2$I[1,1]
  markov_n <- sum(total$I[1,])-total$I[1,1]
  ggplot() +
    geom_step(data = plotdf ,mapping = aes(x=Time, y = value,col = variable)) + 
    geom_step(data = plotdf2,mapping = aes(x=Time, y = value,col = variable),linetype = "dashed") +
    geom_line(data = plotdf3,mapping = aes(x=Time, y = value,col = variable),linetype = "dotted",size=1) +
    theme_bw() +
    labs(title = TeX(paste0("Occupation probabilities under assumption $Z_",s,"=",j,"$")),
         y =  TeX(paste0("$P(Z_t=j|Z_",s,"=",j,")$")),
         subtitle = paste0("As-if estimate based on ",asif_n," observations,\nMarkov based on ",markov_n," observations")) +
    theme(legend.title = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic")) +
    scale_color_manual(values=c("#FA8072", "#FF0000","#8B0000",
                                "#7B68EE", "#1E90FF","#00008B",
                                "#3CB371", "#32CD32","#006400"))
}
plot1 <- plot_function1(paths,3,s= 0, j = 1)
plot2 <- plot_function1(paths,3,s= 1, j = 1)
plot3 <- plot_function1(paths,3,s= 3, j = 1)
plot4 <- plot_function1(paths,3,s= 6, j = 1)
ggarrange(plotlist = list(plot1,plot2,plot3,plot4),ncol = 2,nrow=2)
ggsave("plot1.png",units = "px", width = 1920,height = 1080,scale = 1.75)

#### Appendix A.4 ####
estimate_to_cashflow <- function(estimate,pi,T) {
  Nelson <- estimate$delta_NelsonAalen
  Nelson <- lapply(1:dim(Nelson)[1], function(i) matrix(Nelson[i,2:dim(Nelson)[2]],nrow=num_states,ncol = num_states,byrow=TRUE))
  p_con <- estimate$p_con
  times <- lapply(1:dim(p_con)[1], function(i) p_con[i,1])
  p_con <- lapply(1:dim(p_con)[1], function(i) p_con[i,2:dim(p_con)[2]])
  
  A <- list()
  A[[1]] <- rep(0,4)
  for (i in 2:length(p_con)) {
    A1 <- A[[i-1]][1] + p_con[[i]][1]*ifelse(times[[i]]>T,times[[i]]-max(T,times[[i-1]]),0)
    A2 <- A[[i-1]][2] - pi * p_con[[i]][1]*ifelse(times[[i-1]]<T,min(T,times[[i]])-times[[i-1]],0)
    A3 <- A[[i-1]][3] + p_con[[i]][2]*(times[[i]]-times[[i-1]])
    A4 <- A[[i-1]][4] + p_con[[i-1]][1]*Nelson[[i-1]][1,3] + p_con[[i-1]][2]*Nelson[[i-1]][2,3]
    A[[i]] <- c(A1, A2, A3, A4)
  }
  A <- cbind(unlist(times),matrix(unlist(A),ncol = 4,byrow=TRUE)) %>% as.data.frame()
  colnames(A) <- c("Time","A1","A2","A3","A4")
  A[,"A"] <- rowSums(A[,2:5])
  return(A)
}
estimate <- Estimate(paths,3)
cashflow <- estimate_to_cashflow(estimate,1,3)
plot_function2 <- function(paths,pi,T,num_states,s= 0, j = 1, debug = TRUE) {
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
  cashflow2 <- estimate_to_cashflow(estimate2,pi,T)
  if (debug) {
    print(paste0("Calculate cashflows: ",round(as.numeric(Sys.time()-start_time,units = "secs"),digits=3)," seconds."))
    start_time <- Sys.time()
  }
  plotdf1 <- cashflow1 %>% reshape2::melt(id = "Time")
  plotdf2 <- cashflow2 %>% reshape2::melt(id = "Time")
  asif_n <- sum(estimate2$I[1,])-estimate2$I[1,1]
  markov_n <- sum(estimate1$I[1,])-estimate1$I[1,1]
  #Calculate true values
  times <- s+0:1000*(10-s)/1000
  cashflow3 <- cashflow1
  for (i in 2:dim(cashflow3)[1]) {
    t_1 <- cashflow3$Time[i]
    t_0 <- cashflow3$Time[i-1]
    prod_int_0 <- expm::expm(2*M*(log(1+0.5*t_0)))
    delta_Prod_int <- expm::expm(2*M*(log(1+0.5*t_1)))-prod_int_0
    p_con_0 <- (1:num_states==j)%*%expm::expm(2*M*(log(1+0.5*t_0)-log(1+0.5*s)))
    p_con_1 <- (1:num_states==j)%*%expm::expm(2*M*(log(1+0.5*t_1)-log(1+0.5*s)))
    A1 <- cashflow3$A1[i-1] + p_con_1[1]*ifelse(t_1>T,t_1-max(T,t_0),0)
    A2 <- cashflow3$A2[i-1] - pi * p_con_1[1]*ifelse(t_0<T,min(T,t_1)-t_0,0)
    A3 <- cashflow3$A3[i-1] + p_con_1[2]*(t_1-t_0)
    A4 <- cashflow3$A4[i-1] + p_con_0[1]*delta_Prod_int[1,3] + p_con_0[2]*delta_Prod_int[2,3]
    cashflow3[i,2:5] <- c(A1,A2,A3,A4)
  }
  cashflow3[,6] <- rowSums(cashflow3[,2:5])
  colnames(cashflow3) <- c("Time","A1_true","A2_true","A3_true","A4_true","A_true")
  plotdf3 <- cashflow3 %>% reshape2::melt(id = "Time")
  if (debug) {
    print(paste0("Calculate theoretical cashflows and more: ",round(as.numeric(Sys.time()-start_time,units = "secs"),digits=3)," seconds."))
    start_time <- Sys.time()
  }
  ggplot() +
    geom_step(data = plotdf1 ,mapping = aes(x=Time, y = value,col = variable)) +
    geom_step(data = plotdf2,mapping = aes(x=Time, y = value,col = variable), linetype = "dashed") +
    geom_line(data = plotdf3,mapping = aes(x=Time, y = value,col = variable), linetype = "dotted", size=1) +
    geom_vline(xintercept = T, col = "black",linetype = "dashed") +
    theme_bw() +
    labs(title = TeX(paste0("Accumulated cash-flow under assumption $Z_",s,"=",j,"$")),
         y =  TeX(paste0("A(t)-A(s)")),
         subtitle = paste0("As-if estimate based on ",asif_n," observations,\nMarkov based on ",markov_n," observations")) +
    theme(legend.title = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic")) +
    scale_color_manual(values=c("#FA8072","#8B0000","#1E90FF","#7B68EE","#FF00FF","#800080",
                                "#32CD32","#006400","#F4A460","#DEB887"))
}
plot1 <- plot_function2(paths,1,3,3,s= 0, j = 1)
plot2 <- plot_function2(paths,1,3,3,s= 1, j = 1)
plot3 <- plot_function2(paths,1,3,3,s= 3, j = 1)
plot4 <- plot_function2(paths,1,3,3,s= 6, j = 1)
ggarrange(plotlist = list(plot1,plot2,plot3,plot4),ncol = 2,nrow=2)
ggsave("plot2.png",units = "px", width = 1920,height = 1080,scale = 1.75)
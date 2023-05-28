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
    print(paste0("Convert paths to main_df: ",round(Sys.time()-start_time,digits=2)," seconds."))
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
    print(paste0("Convert main_df to I: ",round(Sys.time()-start_time,digits=2)," seconds."))
    start_time <- Sys.time()
  }
  N_list_tmp <- df_to_N(main_df_tmp,num_states)
  if (debug) {
    print(paste0("Convert main_df to N: ",round(Sys.time()-start_time,digits=2)," seconds."))
    start_time <- Sys.time()
  }
  # 3. Calculate Nelson-Aalen
  NelsonAalen_list_tmp <- N_I_to_NA(I_list_tmp,N_list_tmp,num_states)
  if (debug) {
    print(paste0("Convert I and N to Nelson-Aalen: ",round(Sys.time()-start_time,digits=2)," seconds."))
    start_time <- Sys.time()
  }
  # 4. Calculate Aalen-Johansen
  P_tmp <- NA_to_p(I_list_tmp,N_list_tmp,NelsonAalen_list_tmp, 3)
  if (debug) {
    print(paste0("Convert Nelson-Aalen to transition probs: ",round(Sys.time()-start_time,digits=2)," seconds."))
    start_time <- Sys.time()
  }
  # 5. Calculate probabilties
  p_con_tmp <- P_conditioned(NelsonAalen_list_tmp,s,j,num_states)
  if (debug) {
    print(paste0("Convert transition probs to conditioned probs: ",round(Sys.time()-start_time,digits=2)," seconds."))
    start_time <- Sys.time()
  }
  return(list(I = I_list_tmp$I, I_left_limit = I_list_tmp$I_left_limit,
              N = N_list_tmp$N,delta_N = N_list_tmp$delta_N,
              NelsonAalen = NelsonAalen_list_tmp$NelsonAalen, delta_NelsonAalen= NelsonAalen_list_tmp$delta_NelsonAalen,
              P = P_tmp, p_con = p_con_tmp))
}
total <- Estimate(paths,3)
total <- Estimate(paths,3,s=1,j=1,as_if = TRUE)
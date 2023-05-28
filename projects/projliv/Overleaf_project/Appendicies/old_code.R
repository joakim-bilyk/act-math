
#### Appendix A.2 ####
markov <- function(paths,num_states) {
  require(dplyr)
  L <- length(paths)
  t <- unlist(lapply(1:L, function(i){
    paths[[i]]$times
  }))
  s <- unlist(lapply(1:L, function(i){
    paths[[i]]$states
  }))
  obs <- unlist(lapply(1:L, function(i){
    rep(i,length(paths[[i]]$times))
  }))
  df <- data.frame(id = obs, times = t, states = s)
  #Initiate N^jk, I^j,...
  N_0 <- matrix(rep(0,num_states**2),ncol = num_states, nrow=num_states)
  counts <- df %>% filter(times == 0) %>% group_by(states) %>%
    summarise(count = sum(states>0))
  I_0 <- unlist(lapply(1:num_states, function(i) {sum(counts$count[counts$states == i],na.rm = TRUE)}))
  Nelson_Aalen_0 <- N_0
  Aalen_Johansen_0 <- I_0/L
  results <- list(list(times = 0,N=N_0,I=I_0,Nelson_Aalen = Nelson_Aalen_0, Aalen_Johansen = Aalen_Johansen_0))
  #Go through all unique jump times and calculate N and I
  unique_times <- unique(df %>% arrange(times) %>% select(times))
  for (i in 2:dim(unique_times)[1]) {
    t_0 <- unique_times[i,1]
    #Find observations that jump
    tmp <- df %>%
      #Get history up until t_0
      filter(times <= t_0) %>%
      #Filter out non-jumping Z's
      filter(id %in% id[times == t_0]) %>%
      #Get the latest state and the new state
      group_by(id) %>%
      summarise(sojourn_since = max(times[times < t_0]),
                old_state = states[times == sojourn_since],
                new_state = states[times == t_0])
    #Update N^jk and I^j
    N_1 <- N_0
    I_1 <- I_0
    for (j in unique(tmp$old_state)) {
      for (k in unique(tmp$new_state)) {
        transitions <- sum((tmp$old_state == j)*(tmp$new_state == k),na.rm = TRUE)
        #Remember to through out censoring j->j
        N_1[j,k] <- N_1[j,k] + transitions*(j!=k)
        N_1[j,j] <- N_1[j,j] - transitions*(j!=k)
        I_1[j] <- I_1[j] - transitions
        I_1[k] <- I_1[k] + (transitions)*(j!=k)
      }
    }
    I_inverse <- 1/matrix(rep(I_0,num_states),ncol = 3)
    I_inverse[is.infinite(I_inverse)] <- 0
    Nelson_Aalen_1 <- Nelson_Aalen_0 + I_inverse*(N_1-N_0)
    Aalen_Johansen_1 <- Aalen_Johansen_0 %*% (diag(rep(1,num_states)) +Nelson_Aalen_1 -Nelson_Aalen_0)
    results[[i]] <-list(times = t_0,N=N_1,I=I_1,Nelson_Aalen = Nelson_Aalen_1, Aalen_Johansen = Aalen_Johansen_1)
    #Save new values
    N_0 <- N_1
    I_0 <- I_1
    Nelson_Aalen_0 <- Nelson_Aalen_1
    Aalen_Johansen_0 <- Aalen_Johansen_1
  }
  return(results)
}
results_markov <- markov(paths,3)

#Estimates
plotdf <- data.frame(
  times = rep(unlist(lapply(1:length(results_markov), function(i) results_markov[[i]]$times)),3),
  matrix(unlist(lapply(1:length(results_markov), function(i) results_markov[[i]]$Aalen_Johansen)),ncol=3,byrow=TRUE)
)
colnames(plotdf)[2:4] <- paste0("p_",1:3,"(t)")
plotdf <- plotdf%>% reshape2::melt(., id = "times")
#True values
times <- 0:1000/100
plotdf2 <- data.frame(times = times,
                      matrix(unlist(lapply(times, function(t) c(1,0,0)%*%expm::expm(2*M*log(1+0.5*t)))),ncol=3,byrow=TRUE))
colnames(plotdf2)[2:4] <- paste0("p_",1:3,"(t)")
plotdf2 <- plotdf2 %>% reshape2::melt(., id = "times")

ggplot() +
  geom_step(data = plotdf,mapping = aes(x=times, y = value,col = variable)) + 
  geom_line(data = plotdf2,mapping = aes(x=times, y = value,col = variable)) +
  theme_bw()

#### Appendix A.3 ####
markov_conditioned <- function(paths,num_states, j=1,s=0,asif = TRUE) {
  #We use our markov framework on the subset with Z_s=j.
  #Find observations
  times <- unlist(lapply(1:L, function(i){
    paths[[i]]$times
  }))
  states <- unlist(lapply(1:L, function(i){
    paths[[i]]$states
  }))
  obs <- unlist(lapply(1:L, function(i){
    rep(i,length(paths[[i]]$times))
  }))
  df <- data.frame(id = obs, times = times, states = states)
  tmp <- df %>%
    #Get history up until t_0
    filter(times <= s) %>%
    #Get the latest state and the new state
    group_by(id) %>%
    summarise(sojourn_since = max(times),
              Z_s = states[times == sojourn_since]) %>%
    #Filter only id with Z_s=j
    filter(Z_s == j)
  if (asif) {
    new_paths <- lapply(tmp$id, function(id) paths[[id]])
  } else {
    new_paths <- paths
  }
  #Get regular markov results
  results <- markov(new_paths, num_states)
  times <- unlist(lapply(1:length(results), function(i) results[[i]]$times))
  earlist_time_pre_s <- max(times[times<s])
  earlist_time_pre_s_id <- max(1:length(times)*(times == earlist_time_pre_s))
  #Recalculate aalen-Johansen
  Aalen_Johansen_0 <- results[[earlist_time_pre_s_id]]$I
  Aalen_Johansen_0 <- Aalen_Johansen_0/sum(Aalen_Johansen_0)
  results[[earlist_time_pre_s_id]]$Aalen_Johansen <- Aalen_Johansen_0
  #Forward
  for (i in (earlist_time_pre_s_id+1):length(times)) {
    #Recalculate
    Nelson_Aalen_0 <- results[[i-1]]$Nelson_Aalen
    Nelson_Aalen_1 <- results[[i]]$Nelson_Aalen
    Aalen_Johansen_1 <- Aalen_Johansen_0 %*% (diag(rep(1,num_states)) + Nelson_Aalen_1 - Nelson_Aalen_0)
    results[[i]]$Aalen_Johansen <- Aalen_Johansen_1
    #Save new values
    Aalen_Johansen_0 <- Aalen_Johansen_1
  }
  #Backward
  Aalen_Johansen_0 <- results[[earlist_time_pre_s_id]]$I
  Aalen_Johansen_0 <- Aalen_Johansen_0/sum(Aalen_Johansen_0)
  for (i in (earlist_time_pre_s_id-1):1) {
    #Recalculate
    Nelson_Aalen_0 <- results[[i]]$Nelson_Aalen
    Nelson_Aalen_1 <- results[[i+1]]$Nelson_Aalen
    Aalen_Johansen_1 <- Aalen_Johansen_0%*%(diag(1,nrow = num_states) -(Nelson_Aalen_1 - Nelson_Aalen_0))
    results[[i]]$Aalen_Johansen <- Aalen_Johansen_1
    #Save new values
    Aalen_Johansen_0 <- Aalen_Johansen_1
  }
  return(results)
}
results_asif <- markov_conditioned(paths,num_states = 3,j=2,s=1)

#Estimates
plotdf <- data.frame(
  times = rep(unlist(lapply(1:length(results_asif), function(i) results_asif[[i]]$times)),3),
  matrix(unlist(lapply(1:length(results_asif), function(i) results_asif[[i]]$Aalen_Johansen)),ncol=3,byrow=TRUE)
)
colnames(plotdf)[2:4] <- paste0("p_",1:3,"(t|Z_1=2)")
plotdf <- plotdf%>% reshape2::melt(., id = "times")
#True values
times <- 0:1000/100
plotdf2 <- data.frame(times = times,
                      matrix(unlist(lapply(times, function(t) c(1,0,0)%*%expm::expm(2*M*log(1+0.5*t)))),ncol=3,byrow=TRUE))
colnames(plotdf2)[2:4] <- paste0("p_",1:3,"(t|Z_1=2)")
plotdf2 <- plotdf2 %>% reshape2::melt(., id = "times")

ggplot() +
  geom_step(data = plotdf %>% filter(times >=1),mapping = aes(x=times, y = value,col = variable)) + 
  #geom_line(data = plotdf2 %>% filter(times >=1),mapping = aes(x=times, y = value,col = variable)) +
  theme_bw()

##################################
############# PART 2 #############
##################################
setwd("~/Documents/act-math/projects/skade2-2")
df_claim <- read.csv("bulkcargoclaim.csv")
df_policy <- read.csv("bulkcargopolicy.csv")
library(dplyr)
df_claim <- df_claim %>%
  mutate(
    Dedr = Dedr * 10^(-6), Clr = Clr * 10^(-6),
    Age = log(Age + 2), Brt = log(Brt), Dwt = log(Dwt), Value = log(Value), HP = log(HP),
    Code1 = factor(Code1, levels = c(1,2))
  ) %>% select(-Year)
df_policy <- df_policy %>%
  mutate(
    Dedr = Dedr * 10^(-6),
    Age = log(Age + 2), Brt = log(Brt), Dwt = log(Dwt), Value = log(Value), HP = log(HP),
    Code1 = factor(Code1, levels = c(1,2))
  ) %>% select(-Year)
summary(df_claim)
summary(df_policy)
beta <- beta_hat
beta <- c(1,rep(0,6))
phi <- 1
Z <- df_claim %>%
  mutate(
    Intercept = 1,
    Code1 = (Code1 == 2)*1
  ) %>%
  select(
    Intercept,Code1, Age, Brt, Dwt,Value,HP
  )
X <- df_claim$Clr
D <- df_claim$Dedr
P_G_loglikelihood <- function(beta,Z,phi,X,D) {
  beta_z <- as.numeric(as.matrix(Z) %*% beta)
  phi <- abs(phi)
  n <- dim(Z)[1]
  ll <- n*log(1+phi) + sum(
    ((phi+1)/phi) * log(exp(beta_z) + phi * D) - ((2*phi+1)/phi)*log(exp(beta_z) + phi* X)
  )
  return(ll)
}
P_G_loglikelihood(c(1,rep(0,6)),Z,1,X,D) 
# -280.7415

result <- optim(
  par = c(1,rep(0,6)), #(beta_0,...,beta_6)
  fn= function(x){P_G_loglikelihood(x[1:7],Z,100,X,D)},hessian =TRUE #phi = 100
)
beta_hat <- result$par
beta_hat
# 1.9894381  1.1660804  0.9396159 -1.6693281 -0.4314743 -1.6689529 -0.2078674
result$value
# 8.38396
exp(as.numeric(as.matrix(Z) %*% beta_hat))[1]

count_glm <- glm(
  Claims ~ (Code1 == 2) + Age + Brt + Dwt + Value + HP,
  data = df_policy,family = poisson(link = "log"),
  offset = log(Time)
)
count_glm$coefficients
# (Intercept) Code1 == 2TRUE            Age            Brt            Dwt          Value             HP 
# -4.3273710     -0.4248077      0.2485437     -0.2782506     -0.2958922      0.1344354      0.5553562

#Plotting
Beta_survival <- function(x,a,b) {
  1-zipfR::Ibeta(x,a,b)/zipfR::Cbeta(a,b)
}
Beta_survival(d/(mu/phi+d), 1 + 0:k, 1/phi+1-0:k)
exp_loss_less_ded_P <- function(d,mu,phi,k=1) {
  E <- zipfR::Cgamma(1+0:k)*zipfR::Cgamma(1+1/phi-0:k)*mu**0:k/(phi**0:k*zipfR::Cgamma(1+1/phi))
  sum(choose(k,0:k)*(-d)**(k-0:k)*E*Beta_survival(d/(mu/phi+d), 1 + 0:k, 1/phi+1-0:k))
}
mu_hat <- exp(as.numeric(as.matrix(Z) %*% beta_hat))
phi_hat <- 100
x <- 0:400/100
lambda <- count_glm[["fitted.values"]][1]
y <- unlist(lapply(x, function(d) exp_loss_less_ded_P(d,mean(df_claim$Clr),1)))
plot(x,y*lambda)

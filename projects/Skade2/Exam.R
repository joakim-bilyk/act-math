######################
####### SKADE2 #######
######################

######################
######## 2022 ########
######################

### OPGAVE 2A

X <- c(0.01,0.25,1.41,0.34,0.07,0.71,1.70,0.07,1.12,2.01)

log_likelihood_exp <- function(X,parameter) {
  lambda <- parameter[1]
  -lambda * sum(X)+length(X)*log(lambda)
}

log_likelihood_par <- function(X,parameter) {
  alpha <- parameter[1]
  theta <- parameter[2]
  length(X)*log(alpha)+length(X)*alpha*log(theta)-(1+alpha)*sum(log(1+X))
}

log_likelihood_exp(X,1)
log_likelihood_par(X,c(1,1))

MLE_exp <- function(X) {
  length(X) / sum(X)
}

MLE_par <- function(X) {
  length(X) / sum(log(1+X))
}

lambda_hat <- MLE_exp(X)
alpha_hat <- MLE_par(X)

AIC <- function(X,parameter,log_likelihood) {
  -2 * log_likelihood(X,parameter) + 2
}

#AIC EXP
AIC(
  X,
  lambda_hat,
  log_likelihood_exp
)
#AIC Par
AIC(
  X,
  c(alpha_hat,1),
  log_likelihood_par
)

# PROBLEM 2B
#We transform 
X_exp <- qnorm(exp(-lambda_hat*X))
X_exp
# [1]  2.2286165  0.5901558 -0.9950947  0.3655897  1.3594171 -0.2605618 -1.2285093  1.3594171 -0.7287895 -1.4519568
theta_hat <- 1
X_par <- qnorm(theta_hat**alpha_hat*(theta_hat+ X)**(-alpha_hat))
X_par
# [1]  2.0536296  0.3467102 -0.9638940  0.1304515  1.1340430 -0.4224335 -1.1122202  1.1340430 -0.7810565 -1.2444149

shapiro.test(X_exp) # p-value = 0.5614
shapiro.test(X_par) # p-value = 0.3557

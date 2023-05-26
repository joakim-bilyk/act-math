##################################
############# PART 1 #############
##################################
setwd("~/Documents/act-math/projects/skade2-1")
df <- read.csv("bulkcargoclaim.csv")
head(df)
dim(df)
eq1 <- function(X,D,theta) {
  s1 <- sum(1/(theta+X))
  s2 <- sum(1/(theta+D))
  s3 <- sum(log(theta+X)-log(theta+D))
  n <- length(X)
  total <- s1/(s2-s1)-n/s3
  return(total)
}
eq1(df$Clr,df$Dedr,100000)
uniroot(
  function(x){eq1(df$Clr,df$Dedr,x)},
  lower = 100000,
  upper = 1000000
)
theta <- round(uniroot(
  function(x){eq1(df$Clr,df$Dedr,x)},
  lower = 100000,
  upper = 1000000
)$root,digits = 1)
alpha <- length(df$Clr)/(sum(log(theta + df$Clr)-log(theta+df$Dedr)))
alpha
X <- df$Clr
D <- df$Dedr
mu <- 10
sigma_sqrd <- 1
lognormal_ll <- function(X,D,mu,sigma) {
  sigma_sqrd <- sigma**2
  n <- length(X)
  s1 <- n*log(sigma_sqrd)
  s2 <- n*log(sqrt(2*pi))
  s3 <- sum((log(X)-mu)^2)/(2*sigma_sqrd)
  s4 <- sum(log(pnorm(q = (log(D)-mu)/sqrt(sigma_sqrd),
                  mean = 0, sd = 1,lower.tail = FALSE)))
  total <- (s1+s2+s3+s4)
  return(total)
}
lognormal_ll(df$Clr,df$Dedr,10,1)

inital <- c(mean(log(df$Clr)),1)
result <- optim(
  par = inital,
  fn = function(x){lognormal_ll(df$Clr,df$Dedr,x[1],x[2])}
)
result$par
I_function <- function(theta,a,b,c) {
  require(zipfR)
  (theta**(1+a-b)*gamma(1+a)*gamma(b-a-1)/gamma(b))*zipfR::Ibeta(c/(theta + c),a+1,b-a-1)/beta(a+1,b-a-1)
}
exp_loss_less_ded_P <- function(d,alpha,theta) {
  prob <- (1+d/theta)**(-alpha)
  exp_value <- theta/(alpha-1)
  subtraction <- alpha*theta**(alpha)*I_function(theta,1,1+alpha,d)
  return(exp_value - subtraction - d*prob)
}
theta <- 437469.2
alpha <- 2.6236
exp_loss_less_ded_P(1,alpha,theta)

mu <- 12.56
sigma <- sqrt(0.3406)
exp_loss_less_ded_LN <- function(d,mu,sigma) {
  prop <- pnorm((log(d)-mu)/sigma,mean=0,sd=1,lower.tail = FALSE)
  start <- exp(mu+sigma**2/2)*pnorm((mu-log(d))/sigma+sigma,mean=0,sd=1,lower.tail = TRUE)
  return(start-d*prop)
}
exp_loss_less_ded_LN(1000000,mu,sigma)
plotdf <- data.frame(
  d = (0:20000)*100,
  P =  unlist(lapply((0:20000)*100, function(d){exp_loss_less_ded_P(d,alpha,theta)})),
  LN = unlist(lapply((0:20000)*100, function(d){exp_loss_less_ded_LN(d,mu,sigma)}))
)
library(ggplot2)
p <- ggplot(plotdf) + geom_line(aes(x=d,y=P), col = "red") +geom_line(aes(x=d,y=LN), col = "blue") +
  labs(x= "Deductible", y = "Expected loss") + theme_bw()
ggsave("plot1.png", plot = p,height = 5,width = 7.5)

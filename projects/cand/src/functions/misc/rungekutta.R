rungekutta <- function(t0,y0,t,steps=NULL,dy,eps=NULL) {
  if (is.null(eps)) {
    # Fixed step size
    dt <- (t-t0)/steps
    tt <- seq(t0,t, length.out = steps +1)
    yy <- numeric(steps+1); yy[1] <- y0
    # Solve differential equaiton using Runge-Kutta
    for (i in 1:steps) {
      # Update y by Runge-Kutta fourth order
      k1 <- dy(tt[i],yy[i])
      k2 <- dy(tt[i]+dt/2,yy[i]+dt*k1/2)
      k3 <- dy(tt[i]+dt/2,yy[i]+dt*k2/2)
      k4 <- dy(tt[i]+dt,yy[i]+dt*k3)
      yy[i+1] <- yy[i] + (dt/6)*(k1+2*k2+2*k3+k4)
    }
  } else {
    # Dynamics step size
    if (is.null(steps)) {
      dt_fun <- function() {abs(eps/dy(t0,y0))}
    } else {
      # Step is the minimum step length
      dt0 <- (t-t0)/steps
      dt_fun <- function() {min(dt0,abs(eps/dy(t0,y0)))}
    }
    tt <- t0; yy <- y0
    # Solve differential equaiton using Runge-Kutta
    while (t0 < t) {
      # Est dt by derivative
      dt <- dt_fun()
      # Update y by Runge-Kutta fourth order
      k1 <- dy(t0,y0)
      k2 <- dy(t0+dt/2,y0+dt*k1/2)
      k3 <- dy(t0+dt/2,y0+dt*k2/2)
      k4 <- dy(t0+dt,y0+dt*k3)
      y0 <- y0 + (dt/6)*(k1+2*k2+2*k3+k4)
      t0 <- t0 + dt
      tt <- c(tt,t0); yy <- c(yy,y0)
    }
  }
  list(t=tt,y=yy)
}
rungekutta_hittingtime <- function(t0,y0,yhit,steps=NULL,dy,eps=NULL,nmax =10^4) {
  if (y0 < yhit) {
    stop_criteria <- function(y) {y>yhit}
  } else {
    stop_criteria <- function(y) {y<yhit}
  }
  # Dynamics step size
  if (is.null(steps)) {
    dt_fun <- function() {abs(eps/dy(t0,y0))}
  } else {
    # Step is the minimum step length
    dt0 <- (t-t0)/steps
    dt_fun <- function() {min(dt0,abs(eps/dy(t0,y0)))}
  }
  # Solve differential equaiton using Runge-Kutta
  n <- 1
  while (!stop_criteria(y0) && (n < nmax)) {
    # Est dt by derivative
    dt <- dt_fun()
    # Update y by Runge-Kutta fourth order
    k1 <- dy(t0,y0)
    k2 <- dy(t0+dt/2,y0+dt*k1/2)
    k3 <- dy(t0+dt/2,y0+dt*k2/2)
    k4 <- dy(t0+dt,y0+dt*k3)
    y0 <- y0 + (dt/6)*(k1+2*k2+2*k3+k4)
    t0 <- t0 + dt
    n <- n + 1
  }
  list(t=t0,y=y0,n=n)
}
dy <- function(t,y) {y}
exp_test <- rungekutta(0,1,3,steps = 10,dy=dy) # By stepsize dt
exp_test <- rungekutta(0,1,3,eps = 10^-3,dy=dy) # By stepsize dy
exp_hit <- rungekutta_hittingtime(0,1,2,eps = 10^-2, dy = dy)
exp(exp_hit$t)-2 < 10^-2


dy_exp <- function(t,y) {y}
dy_log <- function(t,y) {y}
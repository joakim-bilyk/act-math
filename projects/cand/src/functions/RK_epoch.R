RK_epoch <- function(t0,y0,dt,dy) {
  # dy is a function of (t,y)
  # Update y by Runge-Kutta fourth order
  k1 <- dy(t0,y0)
  k2 <- dy(t0+dt/2,y0+dt*k1/2)
  k3 <- dy(t0+dt/2,y0+dt*k2/2)
  k4 <- dy(t0+dt,y0+dt*k3)
  y <- y0 + (dt/6)*(k1+2*k2+2*k3+k4)
  return(list(t=t0+dt,y=y))
}
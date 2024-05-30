RK <- function(t0,y0,t,dt,dy) {
  tt <- seq(t0,t,by = dt)
  yy <- list(); yy[[1]] <- y0
  for (i in 1:(length(tt)-1)) {
    yy[[i+1]] <- RK_epoch(tt[i],yy[[i]],dt,dy)$y
  }
  return(list(t = tt, y = yy))
}
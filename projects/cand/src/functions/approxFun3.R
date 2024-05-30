approxFun3 <- function(U_sol) {
  require(fields)
  t <- U_sol$t; r <- U_sol$r; y <- U_sol$y; u <- U_sol$U
  r_range <- range(r); y_range <- range(y)
  obj <- lapply(1:length(t), function(i) list(x = r, y = y, z = t(u[[i]])))
  weights <- function(s) {
    index <- min((1:length(t))[t>=s])
    index <- c(
      max(1,min(index-1,length(t))),
      max(1,min(index,length(t)))
    )
    t_range <- t[index]
    w <- rev(abs(s-t_range)/sum(abs(s-t_range)))
    w <- ifelse(is.nan(w),0.5,w)
    list(index = index, w = w, t_range = t_range)
  }
  approxfun_3d <- function(t,r,y) {
    w <- weights(t)
    tmp <- cbind(
      pmin(r_range[2],pmax(r_range[1],r)),
      pmin(y_range[2],pmax(y_range[1],y)))
    surf1 <- interp.surface(obj[[w$index[1]]], tmp)
    surf2 <- interp.surface(obj[[w$index[2]]], tmp)
    surf1*w$w[1]+surf2*w$w[2]
  }
  return(approxfun_3d)
}
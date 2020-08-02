#' calculating the log target value for the target function.
logtarget <- function(z,a0,h,b0,varc,lambda,v,w1,w2,Y,X){
  B <- balance(X,z)
  mplus <- sum(z == 1)
  mminus <- sum( z == -1)
  n <- dim(X)[1]
  p <- dim(X)[2]
  B1 <- cbind(rep(1,n),B)
  b0_ <- matrix(c(a0,b0),ncol = 1)
  V <- diag(c(h,varc))
  V_star_inv <- solve(V)+t(B1)%*% B1
  V_star <- solve(V_star_inv)
  u <- V_star %*% (solve(V)%*%b0_ + t(B1)%*% Y)
  b_star <- c(v*lambda+t(b0_)%*% V_star_inv %*% b0_ + t(Y) %*% Y - t(u) %*% V_star_inv %*% u)
  logtarget <- log(det(V_star_inv))*(-0.5) - 
    (v + n) / 2 * log(b_star/2) +
    mplus * log(w1) + log(w2) *mminus +  log(1- w1 - w2)*(p - mplus - mminus)
  return(list(logtarget=logtarget,post_beta=list(mean=u,sigma = b_star/(lambda + n)*V_star)))
}


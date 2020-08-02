#for probit regression
##This script is for the target distribution Y,Z|X 
##We first derive the function that's proportional to Y|X,Z, and multiply with the priors of Z
##a0,b0,v0,v1,w1,w2 are all hyperparameters
##Y_star is the unobserved continous outcome, n*1
##B is the balance index (calculated from X and Z), n*1
##X is the relative abudance matrix after 0 replacement, n*p
##Z is the indicator for balance group membership, p*1
##All vectors should be encoded as matrix 

##v4:simplify the target function ,fix a bug in likelihood (should be sigma2^n not sigma2^-1)




probit_logtarget <- function(z,a0,v0,b0,v1,w1,w2,Y_star,X){
  B <- balance(X,z)
  mplus <- sum(z == 1)
  mminus <- sum( z == -1)
  n <- dim(X)[1]
  p <- dim(X)[2]
  sig <- diag(c(v0,v1))
  B1 <- cbind(rep(1,n),B)
  b0_ <- matrix(c(a0,b0),ncol = 1)
  mu <- B1%*% b0_

  sig_star <- solve(B1 %*% sig %*% t(B1) +diag(rep(1,n)))
  logtarget <- 0.5* log(det(sig_star)) -
              0.5*t(Y_star-mu)%*%sig_star%*%(Y_star-mu) +
              mplus * log(w1) + log(w2) *mminus +  log(1- w1 - w2)*(p - mplus - mminus)
  sigma_tilde <- solve(solve(sig) + t(B1)%*%B1)
  mu_tilde <- sigma_tilde %*% (t(B1)%*% Y_star + solve(sig)%*%b0_)
  return(list(logtarget=logtarget,post_beta=list(mean=mu_tilde,sigma = sigma_tilde)))
}


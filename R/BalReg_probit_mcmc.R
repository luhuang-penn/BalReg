#This is the MCMC set up for probit regression
# initialize with z
# sample y*|z,y
# update z

#output: 

#depend on the following two functions:
#source('probit_logtarget_mat_v1.R')
#source('mh_update_v2.R')
#source('balance.R')

# it is possible that balances is all 0 when selected covs have the same value for every individual
# add this constriant
# X1 is the potential 




#' This is the MCMC set up for Balance Regression with binary outcome. The three major steps are:
#'    initialize with z
#'    sample y*|z,y
#'    update z
#'    
#' @param beta0  mean of a multivariate normal distribution which is the prior distribution for regression coefficients
#' @param V0 diagnol elements of the variance-covariance matrix of a multivariate normal distribution which is the prior distribution for regression coefficients
#' @param w1,w2 parameter of prior distribution for the balance configuration vector z
#' @param X compositional matrix
#' @param Y outcome vector
#' @param verbose boolean; indicating whether the program will output verbose progress report
#' @param iter total number of iteration to run
#' @param start starting value of Balance configuration vector z
#' @param filename name of the file to output verbose progress, default to log.txt



#' @expose
BalReg_probit_mcmc<-function(beta0,V0,w1,w2,Y,X,verbose,iter,start,filename = 'log.txt'){
  a0 <- beta0[1]
  b0 <- beta0[2]
  v0 <- V0[0]
  v1 <- V0[1]
  n<- dim(X)[1]
  p <- dim(X)[2]
  Z <- matrix(NA,nrow = iter,ncol = p)
  logtargetval <- rep(NA,iter)
  flip <- rep(0,iter)
  posterior_beta <- rep(list(NA),iter)
  ystar <- matrix(NA,nrow = iter,ncol = n)
  lower <- ifelse(Y==1,0,-Inf)
  upper <- ifelse(Y==0,0,Inf)
  Z[1,] <- start
  
  hidden <- sample_ystar(start,a0,v0,b0,v1,X,lower,upper)
  ltarget <- logtarget(Z[1,],a0,v0,b0,v1,w1,w2,c(hidden$ystar),X)
  ystar[1,] <- c(hidden$ystar)
  ystar_dist[[1]] <-hidden$ystar_dist
  logtargetval[1] <- ltarget$logtarget
  posterior_beta[[1]] <- ltarget$post_beta
  
  B_tmp <- balance(X,Z[1,])
  m <- .lm.fit(cbind(rep(1,n),B_tmp),c(hidden$ystar))
  beta_hat[1,] <- m$coefficients[2]
  
  
  for (i in 2:iter){
    if ( verbose && (i %%1e4 == 2)) {write.table(paste0("Performing Iteration # ", i, " at ", Sys.time()),filename,append=T)}
    
    update <- mh_update(Z[i-1,,drop=F],1/3,1/3,1/4,1/4,1/2,verbose,probit_logtarget,a0,v0,b0,v1,w1,w2,ystar[i-1,],X)
    if (enforce) {
    B_tmp <- balance(X,update$z)
    m <- .lm.fit(cbind(rep(1,n),B_tmp),ystar[i-1,])
    #beta<0 
    #flipping sign of z
      if (is.na(m$coefficients[2]) || (m$coefficients[2] < 0)) {
        update$z <- -1*update$z
        flip[i] <- flip[i] +1
      }
    }
  
    Z[i,] <- update$z
   hidden <- sample_ystar(Z[i,],a0,v0,b0,v1,X,lower,upper)
    ystar[i,] <- c(hidden$ystar)
    posterior_beta[[i]] <- update$post_beta
  }
  return(list(Z=Z,flip=flip,
              posterior_beta=posterior_beta,ystar=ystar))
}

sample_ystar<-function(z,a0,v0,b0,v1,X,lower,upper){
  n <- dim(X)[1]
  p <- dim(X)[2]
  B_tmp <- balance(X,z)
  B1 <- cbind(rep(1,n),B_tmp)
  b0_ <- matrix(c(a0,b0),ncol=1)
  sig <- diag(c(v0,v1))
  sig_star_inv <- diag(rep(1,n))+B1%*%sig%*%t(B1)
  m <- c(B1%*%b0_)
  ystar <- tmvtnorm::rtmvnorm(1, mean = m,
                              sigma = sig_star_inv,
                              lower=lower,
                              upper=upper,
                              algorithm="gibbs")
  return(list(ystar=ystar))
}

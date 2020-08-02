
#' This is the main function to perform MCMC sampling for Balance Regression in linear outcomes
#' 
#' @param beta0  mean of a multivariate normal distribution which is the prior distribution for regression coefficients
#' @param V0 diagnol elements of the variance-covariance matrix of a multivariate normal distribution which is the prior distribution for regression coefficients
#' @param lambda,v parameters of the inverse gamma distribution which is the prior distribution of model variance
#' @param w1,w2 parameter of prior distribution for the balance configuration vector z
#' @param X compositional matrix
#' @param Y outcome vector
#' @param verbose boolean; indicating whether the program will output verbose progress report
#' @param iter total number of iteration to run
#' @param start starting value of Balance configuration vector z
#' 
#' @output a liist of configuration vector z at each iteration and the parameters of posterior distribution of regression coefficients at each iteration
#' 
#' @export
BalReg_mcmc<-function(beta0,V0,lambda,v,w1,w2,Y,X,verbose,iter,start){
  a0 <- beta0[1]
  b0 <- beta0[2]
  h <- V0[0]
  c <- v0[1]
  n<- dim(X)[1]
  p <- dim(X)[2]
  Z <- matrix(NA,nrow = iter,ncol = p)
  beta_trace <- vector(mode = "list",length=iter)
  flip <- rep(0,iter)
  Z[1,] <- start
  
  for (i in 2:iter){
    if (verbose) {cat(sprintf("Performing Iteration # %s \n", i))}
    if ( verbose && (i %%1e4 == 2)) {write.table(paste0("Performing Iteration # ", i, " at ", Sys.time()),filename,append=T)}
    update <- mh_update(Z[i-1,,drop=F],1/3,1/3,1/4,1/4,1/2,verbose,a0,h,b0,c,lambda,v,w1,w2,Y,X)
    beta_trace[[i]] <- update$postbeta
      B_tmp <- balance(X,update$z)
      m <- .lm.fit(cbind(rep(1,n),B_tmp),Y)
      if ((m$coefficients[2] < 0)) {
        update$z <- -1*update$z
        flip[i] <- flip[i] + 1
        }
    Z[i,] <- update$z
    
    
  }
  return(list(Z=Z,beta_trace = beta_trace, flip=flip))
  
}

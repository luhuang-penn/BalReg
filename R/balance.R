#' This is the function to calculate balance; allowing NA values in z corresponds to undetermined group status
#' @param X is the relative abundance matrix
#' @param z is the group indicator vector; 
#' 

balance <- function(X,z){
  #replace NA with 0
  z[is.na(z)] <- 0
  bplus <- z == 1 
  mplus <- sum(bplus)
  bminus <- z == -1
  mminus <- sum(bminus)
  if (mplus == 0 && mminus == 0) {B<- rep(0,dim(X)[1])}
  if (mplus == 0 && mminus != 0){
    B <- apply(X,1,
               function(x){
                 x <- log(x) # take log
                 - sum(x[bminus]) / sqrt(mminus) 
               })
  }
  if (mminus == 0 && mplus != 0){
    B <- apply(X,1,
               function(x){
                 x <- log(x) # take log
                 sum(x[bplus]) / sqrt(mplus) 
               })
  } 
  if (mminus !=0 && mplus !=0){
    B <- apply(X,1,
               function(x){
                 x <- log(x) # take log
                 (sum(x[bplus]) / mplus - sum(x[bminus]) / mminus) *
                   sqrt(1 / (1/mplus + 1/mminus))
               })
  }
  
  
  return(B)
}

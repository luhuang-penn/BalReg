#' Metroplis Hastings Algorithm
#' proposal density:
#'  phi0 probability for changes between 0 and 1 (addition,deletion and switch)
#'  phi1 probability for changes between 0 and -1
#'  phi2 probaiblity for changes between 1 and -1
#'  pa, probability for addition
#'  pd,probability for deletion
#'  pw, probability for swtich

#return is z for next iteration; logtarget value; MH ratio; posterior for beta|z
#possible to add simulated annealing

mh_update <- function(zt,phi0,phi1,pa,pd,pw,verbose,...){
  znew <- zt
  phi2 <- 1 - phi0 - phi1
  mplus <- sum(zt == 1)
  mminus <- sum(zt == -1)
  m0 <- length(zt) - mplus - mminus
  logratio_proposal <- 0
  #type[1] between 0 and 1
  #type[2] between 0 and -1
  #type[3] between 1 and -1
  type <- rmultinom(1,1,c(phi0,phi1,phi2))
  if (m0 == 0 && (mplus == 1 && mminus == 1)){
    type <- c(0,0,1)
  } 
  while (type[1] && m0 == 0 && mplus == 1 && mminus > 1){
    type <- rmultinom(1,1,c(phi0,phi1,phi2))
  }
  while (type[2] && m0 == 0 && mplus > 1 && mminus == 1){
    type <- rmultinom(1,1,c(phi0,phi1,phi2))
  }
  #between 0 and 1
  #movetype == 1 : 0->1
  #movetype == 2 : 1->0
  #movetype == 3 : 0 <==>1 switch
  
  if (type[1]){
    movetype <- which(rmultinom(1,1,c(pa,pd,pw))==1)
    #when there's no 0 in zt
    if ( m0 == 0 && mplus > 1){
      movetype <- 2
    }
    #when there's only one 1 in zt and no zero
    while (movetype == 2  && mplus == 1 && m0 > 0){
      movetype <- which(rmultinom(1,1,c(pa,pd,pw))==1)
    }
    #movetype == 1 : 0->1
    if (movetype == 1){
      idx <- which(zt == 0)
      ad <-  idx[ceiling(runif(1,0,m0))]
      znew[ad] <- 1
      logratio_proposal <- log(pd/(mplus+1)) - 
        log(pa/m0)
      if (verbose) {
        cat(sprintf('  change 0 to 1 at position %s \n', ad))
      }
    }
    #movetype == 2 : 1->0
    if (movetype == 2){
      idx <- which(zt == 1)
      rm <-  idx[ceiling(runif(1,0,mplus))]
      znew[rm] <- 0
      logratio_proposal <- log(pa/(m0+1)) - 
        log(pd/mplus)
      if (verbose) {
        cat(sprintf('  change 1 to 0 at position %s \n', rm))
      }
    }
    #movetype == 3: 1<=>0 swtich
    #logratio_proposal does not change
    if (movetype == 3){
      idx1 <- which(zt == 0)
      idx2 <- which(zt == 1)
      sw1 <- idx1[ceiling(runif(1,0,m0))]
      sw2 <- idx2[ceiling(runif(1,0,mplus))]
      znew[sw1] <- 1
      znew[sw2] <- 0
      if (verbose) {
        cat(sprintf('  switch 0 and 1 bw pos %s and pos %s \n', sw1,sw2))
      }
    }
    
    
  }
  
  
  #between 0 and -1
  #movetype == 1 : 0->-1
  #movetype == 2 : -1->0
  #movetype == 3 : 0 <==>-1 switch
  if (type[2]){
    movetype <- which(rmultinom(1,1,c(pa,pd,pw))==1)
    #when there's no 0 in zt
    if ( m0 == 0 && mminus > 1){
      movetype <- 2
    }
    #when there's only one -1 in zt
    while (movetype == 2  && mminus == 1 && m0 > 0){
      movetype <- which(rmultinom(1,1,c(pa,pd,pw))==1)
    }
    #movetype == 1 : 0->-1
    if (movetype == 1){
      idx <- which(zt == 0)
      ad <-  idx[ceiling(runif(1,0,m0))]
      znew[ad] <- -1
      logratio_proposal <- log(pd/(mminus+1)) - 
        log(pa/m0)
      if (verbose) {
        cat(sprintf('  change 0 to -1 at position %s \n', ad))
      }
    }
    #movetype == 2 : -1->0
    if (movetype == 2){
      idx <- which(zt == -1)
      rm <-  idx[ceiling(runif(1,0,mminus))]
      znew[rm] <- 0
      logratio_proposal <- log(pa/(m0+1)) - 
        log(pd/mminus)
      if (verbose) {
        cat(sprintf('  change -1 to 0 at position %s \n', rm))
      }
    }
    #movetype == 3: -1<=>0 swtich
    #logratio_proposal does not change
    if (movetype == 3){
      idx1 <- which(zt == 0)
      idx2 <- which(zt == -1)
      sw1 <- idx1[ceiling(runif(1,0,m0))]
      sw2 <- idx2[ceiling(runif(1,0,mminus))]
      znew[sw1] <- -1
      znew[sw2] <- 0
      if (verbose) {
        cat(sprintf('  switch -1 and 0 bw pos %s and pos %s \n', sw1,sw2))
      }
    }
    
  }
  #between 1 and -1
  #movetype == 1 : 1 -> -1
  #movetype == 2 : -1 -> 1
  #movetype == 3 : 1 <==> -1 switch
  if (type[3]){
    movetype <- which(rmultinom(1,1,c(pa,pd,pw))==1)
 
    #when there's only one -1 in zt
    while (movetype == 2  && mminus == 1 && mplus > 1){
      movetype <- which(rmultinom(1,1,c(pa,pd,pw))==1)
    }
    #when there's only one 1 in zt
    while (movetype == 1  && mplus == 1 && mminus > 1){
      movetype <- which(rmultinom(1,1,c(pa,pd,pw))==1)
    }
    #when there's only one 1 and only one -1
    if (mplus == 1 && mminus ==1){
      movetype <- 3
    }
    #movetype == 1 : 1 -> -1
    if (movetype == 1){
      idx <- which(zt == 1)
      ad <-  idx[ceiling(runif(1,0,mplus))]
      znew[ad] <- -1
      logratio_proposal <- log(pd/(mminus+1)) - 
        log(pa/mplus)
      if (verbose) {
        cat(sprintf('  change 1 to -1 at pos %s \n', ad))
      }
    }
    #movetype == 2 : -1->1 pd
    if (movetype == 2){
      idx <- which(zt == -1)
      rm <-  idx[ceiling(runif(1,0,mminus))]
      znew[rm] <- 1
      logratio_proposal <- log(pa/(mplus+1)) - 
        log(pd/mminus)
      if (verbose) {
        cat(sprintf('  change -1 to 1 at pos %s \n', rm))
      }
    }
    #movetype == 3: 1<=> -1 swtich
    #logratio_proposal does not change
    if (movetype == 3){
      idx1 <- which(zt == 1)
      idx2 <- which(zt == -1)
      sw1 <- idx1[ceiling(runif(1,0,mplus))]
      sw2 <- idx2[ceiling(runif(1,0,mminus))]
      znew[sw1] <- -1
      znew[sw2] <- 1
      if (verbose) {
        cat(sprintf('  switch 1 and -1 bw pos %s and pos %s \n', sw1,sw2))
      }
    }
    
  }
  
  u <- runif(1,0,1)
  newltarget <- logtarget(znew,...)
  oldltarget <- logtarget(zt,...)
  logratio <- logratio_proposal + newltarget$logtarget - oldltarget$logtarget
  
  
  
  #print("New candicate of Z: ")
  #print(as.vector(znew),quote = F)
  #decide whether accept the new value of z 
 
  #logratio <- logratio_proposal + 
  #  logtarget(znew,a0,h,c,lambda,v,w1,w2,Y,X) - 
  #  logtarget(zt,a0,h,c,lambda,v,w1,w2,Y,X)
  
  if (logratio > 0 || exp(logratio) > u){
    if (verbose){cat(sprintf('Accept the new value \n'))}
    return(list(z=znew, post_beta = newltarget$post_beta))
  } 
  if (verbose) {cat(sprintf('Reject the new value \n'))}
  return(list(z=zt ,post_beta = oldltarget$post_beta))
}

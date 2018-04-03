###############################################################################
#                                                                        Fed 18
#        Fitting a state - space version of a CJS model to the RBT data 
#        Marginalized maximum likelihood version
#
#  Notes:
#  * 
#
#  To do: 
#  * Change the i and K indexing to match the other examples
#  * Ask Charles about the CH 000000001 --> discard ? 
#
###############################################################################
source("RBT_Functions.R", chdir = F)

tmpCH = collapse.ch(CH)[[1]]
sumFR = collapse.ch(CH)[[2]]

# Create vector with occasion of marking
get.first <- function(x) min(which(x != 0))
sumf <- apply(tmpCH, 1, get.first)

#
sumCH = tmpCH
sumCH[sumCH[,] == 0] = 2


NsumCH = nrow(sumCH)         # number of capture histories 
n.occasions = ncol(sumCH)    # number of sampling occasions


#-----------------------------------------------------------------------------#
# Maximum likelihood

rbt_nll <- function(par){
  s  <- par[1:(n.occasions - 1)]     # last s confounded with last p
  p  <- c(par[19:35], 1)
  tr <- array(0, dim = c(2, n.occasions - 1, 2))
  po <- array(0, dim = c(2, n.occasions - 1, 2))
  
  for(i in 1:n.occasions - 1){
    tr[1,i,1] <- s[i]
    tr[1,i,2] <- (1 - s[i])
    tr[2,i,2] <- 1
    
    po[1,i,1] <- p[i]
    po[1,i,2] <- 1 - p[i]
    po[2,i,2] <- 1
  }
  
  nll <- 0
  
  for(k in 1:NsumCH){
    pz <- matrix(0, nrow = n.occasions, ncol = 2)
    pz[sumf[k],1] <- 1
    
    for(i in sumf[k]:(n.occasions - 1)){
      if(sumf[k] == n.occasions) next             # ask Charles about this... need to fix the indexing
      pz[(i + 1),] <- (pz[i,] %*% tr[,i,]) * po[,i,sumCH[k,(i + 1)]]
    }
    
    nll <- sumFR[k] * (-1) * log(sum(pz[n.occasions,])) + nll
  }
  return(nll)
}

#-----------------------------------------------------------------------------#

m <- optim(rep(0.5, 35), rbt_nll, lower = rep(0.01, 35), upper = rep(0.99, 35),
           method = "L-BFGS")

m <- optim(rep(0.5, 35), rbt_nll, lower = rep(0.01, 35), upper = rep(0.99, 35),
           method = "L-BFGS", hessian = T)

sqrt(diag(solve(m$hessian))) #SE










#-----------------------------------------------------------------------------#
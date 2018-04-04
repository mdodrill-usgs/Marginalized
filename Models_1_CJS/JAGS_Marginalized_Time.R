###############################################################################
#                                                                        Feb 18
#        Fitting a state - space version of a CJS model to the RBT data 
#        Marginalized JAGS version - fixed time effects
#
#  Notes:
#  * 
#
#  To do: 
#  * Add parallel version to bottom?
#
###############################################################################
setwd(paste0(getwd(), '/Models_1_CJS'))
library(R2jags)

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

ones <- sumFR

pz = known.state.cjs(tmpCH)

#-----------------------------------------------------------------------------#
sink("JAGS_Marginalized_Time.jags")
cat("
model {
  
  for(k in 1:(n.occasions - 1)){

    s[k] ~ dunif(0,1)
    tr[1,k,1] <- s[k]
    tr[1,k,2] <- (1 - s[k])
    tr[2,k,1] <- 0
    tr[2,k,2] <- 1
    
    p[k] ~ dunif(0,1)
    pmat[1,k,1] <- p[k]
    pmat[1,k,2] <- (1 - p[k])
    pmat[2,k,1] <- 0
    pmat[2,k,2] <- 1
  }
  
  for(i in 1:NsumCH){
    pz[i,sumf[i],1] <- 1
    pz[i,sumf[i],2] <- 0
    
    for(k in (sumf[i] + 1):n.occasions){
      pz[i,k,1] <- pz[i,k-1,1] * tr[1,k-1,1] * pmat[1,k-1,sumCH[i,k]]
      pz[i,k,2] <- inprod(pz[i,k-1,], tr[,k-1,2]) * pmat[2,k-1,sumCH[i,k]]
    }
    
    like[i] <- sum(pz[i,n.occasions,])
    ones[i] ~ dbin(like[i], sumFR[i])
  }
}
    ", fill = TRUE)
sink() 
#-----------------------------------------------------------------------------#

# 
# JM.inits <- function(){list(s = runif((n.occasions - 1), 0, 1),
#                             p = runif((n.occasions - 1), 0, 1),
#                             z = cjs.init.z(CH, indf))}

JM.data <- list(NsumCH = NsumCH, n.occasions = n.occasions, sumCH = sumCH,
                sumf = sumf, sumFR = sumFR, ones = ones)

JM.par <- c('s', 'p')

ni <- 100
nt <- 1
nb <- 50


t1 <- proc.time()
JM.out <- jags(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_Time.jags",
               n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb)
t2 <- proc.time()

print(JM.out, digits = 3)

#-----------------------------------------------------------------------------#
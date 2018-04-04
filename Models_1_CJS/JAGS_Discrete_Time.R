###############################################################################
#                                                                        Feb 18
#        Fitting a state - space version of a CJS model to the RBT data 
#        Discrete JAGS version - fixed time effects
#
#  Notes:
#  * 
#
#  To do: 
#  * Add parallel version?
#
###############################################################################
setwd(paste0(getwd(), '/Models_1_CJS'))
library(R2jags)

source("RBT_Functions.R", chdir = F)

# format data for model fitting
indCH = CH

indCH[indCH[,] == 0] = 2

NindCH = nrow(indCH)         # number of capture histories 
n.occasions = ncol(indCH)    # number of sampling occasions

# Create vector with occasion of marking
get.first <- function(x) min(which(x != 0))
indf <- apply(CH, 1, get.first)

z = known.state.cjs(CH)

#-----------------------------------------------------------------------------#
sink("JAGS_Discrete_Time.jags")
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
  
  for(i in 1:NindCH){
    z[i, indf[i]] <- 1
    
    for(k in (indf[i] + 1):n.occasions){
      
      z[i,k] ~ dcat(tr[z[i, k-1], k-1, ])
      
      indCH[i,k] ~ dcat(pmat[z[i, k], k-1, ])
    }
  }
}
    
    ", fill = TRUE)
sink()    
#-----------------------------------------------------------------------------#


JD.data <- list(NindCH = NindCH, n.occasions = n.occasions, indCH = indCH, indf = indf, z = z)
JD.par <- c('s', 'p')

# JD.inits <- function(){list(s = runif((n.occasions - 1), 0, 1),
#                             p = runif((n.occasions - 1), 0, 1),
#                             z = cjs.init.z(CH, indf))}

ni <- 10
nt <- 1
nb <- 5

t1 <- proc.time()
JD.out <- jags(JD.data, inits = NULL, JD.par, "JAGS_Discrete_Time.jags",
               n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb)
t2 <- proc.time()

print(JD.out, digits = 3)

#-----------------------------------------------------------------------------#

trace_plots(JD.out, "s")

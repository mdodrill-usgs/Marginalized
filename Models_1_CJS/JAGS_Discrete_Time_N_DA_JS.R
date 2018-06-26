###############################################################################
#                                                                       June 18
#        Fitting a state - space version of a JS model to the RBT data 
#        Discrete JAGS version
#
#  Notes:
#  * Abundance (N) is computed using Data Augmentation (DA)
#  To do: 
#  * 
#
###############################################################################
setwd(paste0(getwd(), '/Models_1_CJS'))
# setwd("U:/marginalized_2/Models_1_CJS/")
library(R2jags)

source("RBT_Functions.R", chdir = F)

# format data for model fitting
indCH = CH

indCH[indCH[,] == 0] = 2

# NindCH = nrow(indCH)         # number of capture histories 
# n.occasions = ncol(indCH)    # number of sampling occasions

# Create vector with occasion of marking
get.first <- function(x) min(which(x != 0))
indf <- apply(CH, 1, get.first)

# DA parms 

# Add dummy occasion
CH.du <- cbind(rep(0, dim(CH)[1]), CH)

n.occasions = ncol(CH.du)    # number of sampling occasions

nz = 5000
CH.ms  = rbind(CH.du, array(0, dim = c(nz, n.occasions)))

# Recode CH matrix: a 0 is not allowed in WinBUGS!
CH.ms[CH.ms==0] <- 2                     # Not seen = 2, seen = 1

#-----------------------------------------------------------------------------#
sink("JAGS_Discrete_Time_working.jags")
cat("
model{
  
  for(t in 1:(n.occasions - 1)){
    s[t] ~ dunif(0,1)
    p[t] ~ dunif(0,1)
    gamma[t] ~ dunif(0, 1) # Prior for entry probabilities
  }
  
  for(i in 1:M){
    for(t in 1:(n.occasions - 1)){
      
      tr[1,i,t,1] <- 1 - gamma[t]
      tr[1,i,t,2] <- gamma[t]
      tr[1,i,t,3] <- 0
      tr[2,i,t,1] <- 0
      tr[2,i,t,2] <- s[t]
      tr[2,i,t,3] <- 1 - s[t]
      tr[3,i,t,1] <- 0
      tr[3,i,t,2] <- 0
      tr[3,i,t,3] <- 1
      
      pmat[1,i,t,1] <- 0
      pmat[1,i,t,2] <- 1
      pmat[2,i,t,1] <- p[t]
      pmat[2,i,t,2] <- 1 - p[t]
      pmat[3,i,t,1] <- 0
      pmat[3,i,t,2] <- 1
    }
  }
  
  for (i in 1:M){
    z[i,1] <- 1
    
    for(t in 2:n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(tr[z[i, t-1], i, t-1, ])
      
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(pmat[z[i,t], i, t-1,])
    }
  }
  
  # Calculate derived population parameters
  for (t in 1:(n.occasions-1)){
    qgamma[t] <- 1-gamma[t]
  }
  
  cprob[1] <- gamma[1]
  for (t in 2:(n.occasions-1)){
    cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
  } #t
  
  psi <- sum(cprob[])            # Inclusion probability
  for (t in 1:(n.occasions-1)){
    b[t] <- cprob[t] / psi      # Entry probability
  } #t
  
  for (i in 1:M){
    for (t in 2:n.occasions){
      al[i,t-1] <- equals(z[i,t], 2)
    } #t
    
    for (t in 1:(n.occasions-1)){
      d[i,t] <- equals(z[i,t]-al[i,t],0)
    } #t   
    
    alive[i] <- sum(al[i,])
  } #i
  
  for (t in 1:(n.occasions-1)){
    N[t] <- sum(al[,t])        # Actual population size
    B[t] <- sum(d[,t])         # Number of entries
  } #t
  
  for (i in 1:M){
    w[i] <- 1-equals(alive[i],0)
  } #i
  
  Nsuper <- sum(w[])            # Superpopulation size
}
    
    ", fill = TRUE)
sink()    
#-----------------------------------------------------------------------------#
JD.data <- list(y = CH.ms, n.occasions = dim(CH.ms)[2], M = dim(CH.ms)[1])

JD.par <- c('s', 'p', 'gamma', 'b', 'Nsuper', 'N', 'B')


inits <- function(){list(z = js.multistate.init(CH.du, nz))}

ni <- 20
nt <- 1
nb <- 10

# t1 <- proc.time()
# JD.out <- jags(JD.data, inits = inits, JD.par, "JAGS_Discrete_Time_working.jags",
#                n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb)
# t2 <- proc.time()

# print(JD.out, digits = 3)

#-----------------------------------------------------------------------------#
t1 <- proc.time()
JD.out <- jags.parallel(JD.data, inits = inits, JD.par, "JAGS_Discrete_Time_working.jags",
               n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb,
               export_obj_names = c("nb", "js.multistate.init", "CH.du", "nz",
                                    "ni", "nt"))
t2 <- proc.time()




# trace_plots(JD.out, "s")
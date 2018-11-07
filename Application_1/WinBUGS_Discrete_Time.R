###############################################################################
#                                                                        Nov 18
#  Fitting a multi-state version of a CJS model to the RBT data 
#  Discrete WinBUGS version - fixed time effects
#
#  Notes:
#  * 
#
#  To do: 
#  * 
#
###############################################################################
setwd(paste0(getwd(), '/Application_1'))
library(R2WinBUGS)

source("RBT_Functions.R", chdir = F)

# format data for model fitting
Y = CH

Y[Y[,] == 0] = 2

NY = nrow(Y)          # number of capture histories 
Nint = ncol(Y) - 1    # number of intervals

# Create vector with occasion of marking
get.first <- function(x) min(which(x != 0))
indf <- apply(CH, 1, get.first)

Z = known.state.cjs(CH)

#-----------------------------------------------------------------------------#
sink("WinBUGS_Discrete_Time.WinBUGS")
cat("
model{
  
  for(t in 1:Nint){ 
    phi[t] ~ dunif(0,1)
    p[t] ~ dunif(0,1)
    
    omega[1,1,t] <- phi[t]
    omega[1,2,t] <- 1 - phi[t]
    omega[2,1,t] <- 0
    omega[2,2,t] <- 1
    
    rho[1,1,t] <- p[t]
    rho[1,2,t] <- 1 - p[t]
    rho[2,1,t] <- 0
    rho[2,2,t] <- 1
  }
  
  for(i in 1:NY){
    Z[i,indf[i]] <- 1
    
    for(t in indf[i]:Nint){
      Z[i,(t + 1)] ~ dcat(omega[Z[i, t], , t])
      Y[i,(t + 1)] ~ dcat(rho[Z[i, (t + 1)], , t])
    }
  }
}
    ", fill = TRUE)
sink()    

#-----------------------------------------------------------------------------#
JD.data <- list(NY = NY, Nint = Nint, Y = Y, indf = indf, Z = Z)
JD.par <- c('phi', 'p')

JD.inits <- function(){list(phi = runif((Nint), .5, .5),
                            p = runif((Nint), .5, .5),
                            Z = cjs.init.z(CH, indf))}


ni <- 10
nt <- 1
nb <- 5

# t1 <- proc.time()
# JD.out <- R2WinBUGS::bugs(JD.data, inits = NULL, JD.par, "WinBUGS_Discrete_Time.WinBUGS",
JD.out <- R2WinBUGS::bugs(JD.data, inits = JD.inits, JD.par, "WinBUGS_Discrete_Time.WinBUGS",
               n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb, debug = TRUE)
# t2 <- proc.time()

print(JD.out, digits = 3)

#-----------------------------------------------------------------------------#
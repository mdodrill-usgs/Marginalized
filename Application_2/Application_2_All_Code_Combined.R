
# check the MCMC setting vs the paper...




###############################################################################
#                                                                     Spring 19
#  Dynamic two-species occupancy model with autologistic effects
#  (Barred owls and Northern Spotted owls)
#  JAGS Discrete Version
#
#  Notes:
#  * Need to set directory for data
#
###############################################################################
library(R2jags)

# Get data:
RF <-read.csv(".\\RF.csv", header = FALSE)
FOR <- read.csv(".\\FOR.csv", header = FALSE)
visited <- read.csv(".\\visited.csv", header = FALSE)

#-----------------------------------------------------------------------------#
# format data for model fitting:

# detection covariates: day, night, method1, method3, second half, and survey
# length effects on BO detection
pX <- array(NA, dim = c(158, 176, 7))
pX[,,1] <- matrix(unlist(visited), nrow = dim(visited)[1], ncol = dim(visited)[2])
pX[,,2] <- as.matrix(read.csv(".\\day.csv", header = FALSE))
pX[,,3] <- as.matrix(read.csv(".\\night.csv", header = FALSE))
pX[,,4] <- as.matrix(read.csv(".\\mtd1.csv", header = FALSE))
pX[,,5] <- as.matrix(read.csv(".\\mtd3.csv", header = FALSE))
pX[,1:88,6] <- 0
pX[,89:176,6] <- 1
pX[,,7] <- as.matrix(read.csv(".\\tott.csv", header = FALSE))

# If site is not visited, change all covariates to zero:
for(i in 1:158){
  for(j in 1:176){
    if(visited[i,j] == 0){pX[i,j,] <- 0}
  }
}

owls <- as.matrix(read.csv(".\\owls.csv", header = FALSE))
Y_BO <- ifelse(owls > 1, 1, 0)
Y_NSO <- ifelse(owls == 1 | owls == 3,1,0)
Y <- owls + 1

# Matrix 'BO' is whether or not barred owls were detected at each site at least
# once in each year
BO <- numeric()
for(i in 1:22){
  BO <- cbind(BO, apply(Y_BO[,(i * 8-7):(i * 8)], 1, sum))
}

1 -> BO[BO > 1]

# Matrix 'NSO' is whether or not spotted owls were detected at each site at
# least once in each year
NSO <- numeric()
for(i in 1:22){
  NSO <- cbind(NSO, apply(Y_NSO[,(i * 8-7):(i * 8)], 1, sum))
}

1 -> NSO[NSO > 1]

#-----------------------------------------------------------------------------#
# JAGS-Discrete
# For each species, two states: 0 - absent; 1 - present and two sets of process
# variables - colonization and extinction requires matrix of visited per
# site/year/8 visits with 0 when not visited and 1 when visited and y

sink("JAGS_Discrete_owls.txt")
cat("
model{
  
  psi0_BO ~ dunif(0,1)
  psi0_NSO ~ dunif(0,1)
  
  # time varying colonization
  for(t in 1:21){
    gamNSO_int[t] ~ dunif(-5,5)
  }
  
  gamNSO_b ~ dunif(-5,5) # Forest effect
  epsNSO_int ~ dunif(-5,5)
  epsNSO_b ~ dunif(-5,5) # BO effect
  gamBO_int ~ dunif(-5,5)
  epsBO_int ~ dunif(-5,5)
  
  gamBO_b[1] ~ dunif(-5,5)   # Riparian forest effect on colonization
  epsBO_b[1] ~ dunif(-5,5)   # Riparian forest effect on extinction
  gamBO_b[2] ~ dunif(-15,15) # Autologistic effecton colonization
  epsBO_b[2] ~ dunif(-15,15) # Autologistic effect on extinction
  gamBO_b[3] ~ dunif(-5,5)   # NSO effect on colonization
  epsBO_b[3] ~ dunif(-5,5)   # NSO effect on extinction
  
  # day, night, method1, method3 and BO effects on NSO detection
  for(i in 1:6){
    pNSO_b[i] ~ dunif(-5,5)
  }
  
  # day, night, method1, method3, second half, and survey length effects on BO detection
  for(i in 1:7){
    pBO_b[i] ~ dunif(-5,5)
  }
  
  # Initial occupancy at t=1
  for(k in 1:158){
    BO[k,1] ~ dbern(psi0_BO)
    NSO[k,1] ~ dbern(psi0_NSO)
    for(j in 1:8){
      pBO[k,j] <- pX[k,j,1] * ilogit(inprod(pBO_b,pX[k,j,]))
      pNSO[k,j] <- pX[k,j,1] * ilogit(inprod(pNSO_b[1:5],pX[k,j,1:5]) + pNSO_b[6] * BO[k,1])
      Y_BO[k,j] ~ dbern(BO[k,1] * pBO[k,j])
      Y_NSO[k,j] ~ dbern(NSO[k,1] * pNSO[k,j])
    }
  }
  
  # Colonization/extinction dynamics:
  for(t in 1:21){
    st_psiBO[t] <- (sum(BO[,t]) / 158 - .5) # autologistic BO effect - approximately centered
    for(k in 1:158){
      logit(gamBO[k,t]) <- gamBO_int + gamBO_b[1] * RF[k,t] + gamBO_b[2] * st_psiBO[t] + gamBO_b[3] * NSO[k,t]
      logit(epsBO[k,t]) <- epsBO_int + epsBO_b[1] * RF[k,t] + epsBO_b[2] * st_psiBO[t] + epsBO_b[3] * NSO[k,t]
      logit(gamNSO[k,t]) <- gamNSO_int[t] + gamNSO_b * FOR[k,t]
      logit(epsNSO[k,t]) <- epsNSO_int + epsNSO_b * BO[k,t]
      BO[k,(t+1)] ~ dbern(gamBO[k,t] * (1 - BO[k,t]) + (1 - epsBO[k,t]) * BO[k,t]) 
      NSO[k,(t+1)] ~ dbern(gamNSO[k,t] * (1 - NSO[k,t]) + (1 - epsNSO[k,t]) * NSO[k,t])
      # Use process model to update observation model	
      for(j in 1:8){
        pBO[k,(j + t * 8)] <- pX[k,(j + t * 8),1] * ilogit(inprod(pBO_b,pX[k,(j + t * 8),]))
        pNSO[k,(j + t * 8)] <- pX[k,(j + t * 8),1] * ilogit(inprod(pNSO_b[1:5],pX[k,(j + t * 8),1:5]) + pNSO_b[6] * BO[k,t+1])
        Y_BO[k,(j + t * 8)] ~ dbern(BO[k,(t + 1)] * pBO[k,(j + t * 8)])
        Y_NSO[k,(j + t * 8)] ~ dbern(NSO[k,(t + 1)] * pNSO[k,(j + t * 8)])
      }
    }
  }
}
    
    ", fill = TRUE)
sink()  

#-----------------------------------------------------------------------------#
owls_JD.data <- list(Y_BO = array(Y_BO, dim = c(dim(Y_BO)[1], dim(Y_BO)[2])),
                     Y_NSO = array(Y_NSO,dim = c(dim(Y_NSO)[1], dim(Y_NSO)[2])),
                     BO = array(BO, dim = c(dim(BO)[1], dim(BO)[2])),
                     NSO = array(NSO, dim = c(dim(NSO)[1], dim(NSO)[2])),
                     RF = array(as.vector(unlist(RF)), dim = c(dim(RF)[1], dim(RF)[2])),  
                     FOR = array(as.vector(unlist(FOR)), dim = c(dim(FOR)[1], dim(FOR)[2])),
                     pX = array(pX, dim = c(dim(pX)[1], dim(pX)[2], dim(pX)[3])))

owls_JD.par <- c('pBO_b', 'pNSO_b', 'psi0_BO', 'psi0_NSO', 'gamNSO_b', 'epsNSO_b',
                 'epsNSO_int', 'gamBO_int', 'epsBO_int', 'gamBO_b', 'epsBO_b', 'gamNSO_int')

owls_JD.out <- jags.parallel(owls_JD.data, inits = NULL, owls_JD.par,
                             "JAGS_Discrete_owls.txt", n.chains = 3,
                             n.iter = 2000, jags.seed = 1)

#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Dynamic two-species occupancy model with autologistic effects
#  (Barred owls and Northern Spotted owls)
#  JAGS Marginalized Version
#
#  Notes:
#  * Need to set directory for data
#
###############################################################################
library(R2jags)

# Get data:
RF <- read.csv(".\\RF.csv", header = FALSE)
FOR <- read.csv(".\\FOR.csv", header = FALSE)
visited <- read.csv(".\\visited.csv", header = FALSE)

#-----------------------------------------------------------------------------#
# format data for model fitting:

# JAGS-Marginalized - 4 states: 
# 1 - Neither present, 2 - NSO present, 3 - BO present, 4 - both species present
# requires matrix of visited per site/year/8 visits with 0 when not visited and
# 1 when visited, one is rep (1,158), and y

# detection covariates: day, night, method1, method3, second half, and survey
# length effects on BO detection
pX <- array(NA, dim = c(158, 176, 7))
pX[,,1] <- matrix(unlist(visited), nrow = dim(visited)[1], ncol = dim(visited)[2])
pX[,,2] <- as.matrix(read.csv(".\\day.csv", header = FALSE))
pX[,,3] <- as.matrix(read.csv(".\\night.csv", header = FALSE))
pX[,,4] <- as.matrix(read.csv(".\\mtd1.csv", header = FALSE))
pX[,,5] <- as.matrix(read.csv(".\\mtd3.csv", header = FALSE))
pX[,1:88,6] <- 0
pX[,89:176,6] <- 1
pX[,,7] <- as.matrix(read.csv(".\\tott.csv", header = FALSE))

# If site is not visited, change all covariates to zero:
for(i in 1:158){
  for(j in 1:176){
    if(visited[i,j] == 0){pX[i,j,] <- 0}
  }
}

owls <- as.matrix(read.csv(".\\owls.csv", header = FALSE))
Y_BO <- ifelse(owls > 1, 1, 0)
Y_NSO <- ifelse(owls == 1 | owls == 3,1,0)
Y <- owls + 1

#-----------------------------------------------------------------------------#
sink("JAGS_Marginalized_owls.txt")
cat("
model{
  
  # Unconditional occupancy estimates:
  psi0_BO ~ dunif(0,1)
  psi0_NSO ~ dunif(0,1)
  
  # time varying colonization
  for(t in 1:21){
    gamNSO_int[t] ~ dunif(-5,5)
  } 
  
  gamNSO_b ~ dunif(-5,5) # Forest effect
  epsNSO_int ~ dunif(-5,5)
  epsNSO_b ~ dunif(-5,5) # BO effect
  gamBO_int ~ dunif(-5,5)
  epsBO_int ~ dunif(-5,5)
  
  gamBO_b[1] ~ dunif(-5,5)   # Riparian forest effect on colonization
  epsBO_b[1] ~ dunif(-5,5)   # Riparian forest effect on extinction
  gamBO_b[2] ~ dunif(-15,15) # Autologistic effecton colonization
  epsBO_b[2] ~ dunif(-15,15) # Autologistic effect on extinction
  gamBO_b[3] ~ dunif(-5,5)   # NSO effect on colonization
  epsBO_b[3] ~ dunif(-5,5)   # NSO effect on extinction
  
  # day, night, method1, method3 and BO effects on NSO detection
  for(i in 1:6){
    pNSO_b[i] ~ dunif(-5,5)
  } 
  
  # day, night, method1, method3, second half, and survey length effects
  # on BO detection
  for(i in 1:7){
    pBO_b[i] ~ dunif(-5,5)
  } 
  
  # Observation model:
  for(k in 1:158){
    for(j in 1:176){
      pBO[k,j] <- pX[k,j,1] * ilogit(inprod(pBO_b, pX[k,j,]))
      pNSO_BO[k,j] <- pX[k,j,1] * ilogit(inprod(pNSO_b[1:5], pX[k,j,1:5]) + pNSO_b[6])
      pNSO_bo[k,j] <- pX[k,j,1] * ilogit(inprod(pNSO_b[1:5], pX[k,j,1:5]))
      p[k,j,1,1] <- 1 # prob that site in state 1 (no BO or NSO) is detected as state 1 
      p[k,j,1,2] <- 0
      p[k,j,1,3] <- 0
      p[k,j,1,4] <- 0
      p[k,j,2,1] <- 1 - pNSO_bo[k,j] # prob that site in state 2 (NSO only) is detected as state 1 (no NSO or BO)
      p[k,j,2,2] <- pNSO_bo[k,j] # prob that site in state 2 (NSO only) is detected as state 2 (NSO only)
      p[k,j,2,3] <- 0
      p[k,j,2,4] <- 0
      p[k,j,3,1] <- 1 - pBO[k,j] # prob that site in state 3 (BO only) is detected as state 1 (no NSO or BO)
      p[k,j,3,2] <- 0
      p[k,j,3,3] <- pBO[k,j] # prob that site in state 3 (BO only) is detected as state 3 (BO only)
      p[k,j,3,4] <- 0
      p[k,j,4,1] <- 1 - (pBO[k,j] + pNSO_BO[k,j] - pBO[k,j] * pNSO_BO[k,j]) # prob that site in state 4 (NSO & BO) is detected as state 1 (no NSO or BO)
      p[k,j,4,2] <- pNSO_BO[k,j] * (1 - pBO[k,j]) # prob that site in state 4 (NSO & BO) is detected as state 2 (NSO only)
      p[k,j,4,3] <- (1 - pNSO_BO[k,j]) * pBO[k,j] # prob that site in state 4 (NSO & BO) is detected as state 3 (BO only)
      p[k,j,4,4] <- pNSO_BO[k,j] * pBO[k,j] # prob that site in state 4 (NSO & BO) is detected as state 4 (NSO & BO)
    }
    for(t in 1:22){
      cp[k,(1+(t-1)*8),1] <- p[k,(1+(t-1)*8),1,Y[k,(1+(t-1)*8)]]
      cp[k,(1+(t-1)*8),2] <- p[k,(1+(t-1)*8),2,Y[k,(1+(t-1)*8)]]
      cp[k,(1+(t-1)*8),3] <- p[k,(1+(t-1)*8),3,Y[k,(1+(t-1)*8)]]
      cp[k,(1+(t-1)*8),4] <- p[k,(1+(t-1)*8),4,Y[k,(1+(t-1)*8)]]
      for (j in 2:8){
        cp[k,(j+(t-1)*8),1] <- cp[k,(j-1+(t-1)*8),1]*p[k,(j+(t-1)*8),1,Y[k,(j+(t-1)*8)]]
        cp[k,(j+(t-1)*8),2] <- cp[k,(j-1+(t-1)*8),2]*p[k,(j+(t-1)*8),2,Y[k,(j+(t-1)*8)]]
        cp[k,(j+(t-1)*8),3] <- cp[k,(j-1+(t-1)*8),3]*p[k,(j+(t-1)*8),3,Y[k,(j+(t-1)*8)]]
        cp[k,(j+(t-1)*8),4] <- cp[k,(j-1+(t-1)*8),4]*p[k,(j+(t-1)*8),4,Y[k,(j+(t-1)*8)]]
      }
    }
  }
  
  # Unconditional occupancy:
  unc_pz[1] <- (1-psi0_BO) * (1-psi0_NSO)
  unc_pz[2] <- (1-psi0_BO) * psi0_NSO
  unc_pz[3] <- psi0_BO * (1-psi0_NSO)
  unc_pz[4] <- psi0_BO * psi0_NSO
  
  # Use detection data at t=1 to update state probabilities:		
  for (k in 1:158){
    pz[k,1,1] <- unc_pz[1] * cp[k,8,1]
    pz[k,1,2] <- unc_pz[2] * cp[k,8,2]
    pz[k,1,3] <- unc_pz[3] * cp[k,8,3]
    pz[k,1,4] <- unc_pz[4] * cp[k,8,4]
    BO[k,1] <- sum(pz[k,1,3:4]) / sum(pz[k,1,1:4])
  }
  
  logit(epsNSO_BO) <- epsNSO_int + epsNSO_b
  logit(epsNSO_bo) <- epsNSO_int
  
  #Process model:
  for(t in 1:21){
    st_psiBO[t] <- (sum(BO[1:158,t]) / 158 - .5) #covariate for BO autologistic effect
    for(k in 1:158){
      logit(gamNSO[k,t]) <- gamNSO_int[t] + gamNSO_b * FOR[k,t]
      logit(gamBO_NSO[k,t]) <- gamBO_int + gamBO_b[1] * RF[k,t] + gamBO_b[2] * st_psiBO[t] + gamBO_b[3]
      logit(gamBO_nso[k,t]) <- gamBO_int + gamBO_b[1] * RF[k,t] + gamBO_b[2] * st_psiBO[t]
      logit(epsBO_NSO[k,t]) <- epsBO_int + epsBO_b[1] * RF[k,t] + epsBO_b[2] * st_psiBO[t] + epsBO_b[3]
      logit(epsBO_nso[k,t]) <- epsBO_int + epsBO_b[1] * RF[k,t] + epsBO_b[2] * st_psiBO[t]
      tr[k,t,1,1] <- (1 - gamNSO[k,t]) * (1 - gamBO_nso[k,t])
      tr[k,t,1,2] <- gamNSO[k,t] * (1 - gamBO_nso[k,t])
      tr[k,t,1,3] <- (1 - gamNSO[k,t]) * gamBO_nso[k,t]
      tr[k,t,1,4] <- gamNSO[k,t] * gamBO_nso[k,t]
      tr[k,t,2,1] <- epsNSO_bo * (1 - gamBO_NSO[k,t])
      tr[k,t,2,2] <- (1 - epsNSO_bo) * (1 - gamBO_NSO[k,t])
      tr[k,t,2,3] <- epsNSO_bo * gamBO_NSO[k,t]
      tr[k,t,2,4] <- (1 - epsNSO_bo) * gamBO_NSO[k,t]
      tr[k,t,3,1] <- (1 - gamNSO[k,t]) * epsBO_nso[k,t]
      tr[k,t,3,2] <- gamNSO[k,t] * epsBO_nso[k,t]
      tr[k,t,3,3] <- (1 - gamNSO[k,t]) * (1 - epsBO_nso[k,t])
      tr[k,t,3,4] <- gamNSO[k,t] * (1 - epsBO_nso[k,t])
      tr[k,t,4,1] <- epsNSO_BO * epsBO_NSO[k,t]
      tr[k,t,4,2] <- (1 - epsNSO_BO) * epsBO_NSO[k,t]
      tr[k,t,4,3] <- epsNSO_BO * (1 - epsBO_NSO[k,t])
      tr[k,t,4,4] <- (1 - epsNSO_BO) * (1 - epsBO_NSO[k,t])
      
      # Update process model using observation model:
      pz[k,(t+1),1] <- inprod(pz[k,t,], tr[k,t,,1]) * cp[k,(8*(t+1)),1]
      pz[k,(t+1),2] <- inprod(pz[k,t,], tr[k,t,,2]) * cp[k,(8*(t+1)),2]
      pz[k,(t+1),3] <- inprod(pz[k,t,], tr[k,t,,3]) * cp[k,(8*(t+1)),3]
      pz[k,(t+1),4] <- inprod(pz[k,t,], tr[k,t,,4]) * cp[k,(8*(t+1)),4]
      BO[k,(t+1)] <- sum(pz[k,(t+1),3:4]) / sum(pz[k,(t+1),1:4])
    }
  }
  
  for(k in 1:158){
    lik[k] <- sum(pz[k,22,])
    one[k] ~ dbern(lik[k])
  }
}
    ",fill=TRUE)
sink()    

#-----------------------------------------------------------------------------#
owls_JM.data <- list(Y = array(Y, dim = c(dim(Y)[1], dim(Y)[2])), 
                     RF = array(as.vector(unlist(RF)), dim = c(dim(RF)[1], dim(RF)[2])),
                     one = rep(1,158),
                     FOR = array(as.vector(unlist(FOR)), dim = c(dim(FOR)[1], dim(FOR)[2])),
                     pX = array(pX, dim = c(dim(pX)[1], dim(pX)[2], dim(pX)[3])))

owls_JM.par <- c('pBO_b', 'pNSO_b', 'psi0_BO', 'psi0_NSO', 'gamNSO_b',
                 'epsNSO_b', 'epsNSO_int', 'gamNSO_int', 'gamBO_int', 
                 'epsBO_int', 'gamBO_b', 'epsBO_b')

owls_JM.out <- jags.parallel(owls_JM.data, inits = NULL, owls_JM.par,
                             ".\\JAGS_Marginalized_owls.txt", n.chains = 3,
                             n.iter = 10)

#-----------------------------------------------------------------------------#
###############################################################################
#              Dynamic two-species occupancy model with autologistic effects (Barred owls and Northern Spotted owls)
#
#                               Stan Marginalized Version
#
#-----------------------------------------------------------------------------#

library(rstan)

#To run Stan in parallel:
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#Get data:
setwd("U:\\Marginalized likelihood\\BO_NSO")
RF<-read.csv(".\\RF.csv",header=FALSE)
FOR<-read.csv(".\\FOR.csv",header=FALSE)
visited<-read.csv(".\\visited.csv",header=FALSE)

#-----------------------------------------------------------------------------#
# format data for model fitting:

# Stan-M - 4 states: 1 - Neither present, 2 - NSO present, 3 - BO present, 4 - both species present
# requires list of visited sites & notvisited sites, matrix of visits per site/year with 0 visit for unvisited sites, and y, 

pX<-array(NA,dim=c(158,176,7))# detection covariates: day, night, method1, method3, second half, and survey length effects on BO detection
pX[,,1]<-matrix(unlist(visited),nrow=dim(visited)[1],ncol=dim(visited)[2])
pX[,,2]<-as.matrix(read.csv(".\\day.csv",header=FALSE))
pX[,,3]<-as.matrix(read.csv(".\\night.csv",header=FALSE))
pX[,,4]<-as.matrix(read.csv(".\\mtd1.csv",header=FALSE))
pX[,,5]<-as.matrix(read.csv(".\\mtd3.csv",header=FALSE))
pX[,1:88,6]<-0
pX[,89:176,6]<-1
pX[,,7]<-as.matrix(read.csv(".\\tott.csv",header=FALSE))
#If site is not visited, change all covariates to zero:
for(i in 1:158){for(j in 1:176){ if(visited[i,j]==0){pX[i,j,]<-0}}}  
owls<-as.matrix(read.csv(".\\owls.csv",header=FALSE))
Y_BO<-ifelse(owls>1,1,0)
Y_NSO<-ifelse(owls==1|owls==3,1,0)
Y<-owls+1

n_occ<-apply(visited[,1:8],1,sum)
start<-1:22*8-7
for(i in 2:length(start)){n_occ<-cbind(n_occ,apply(visited[,start[i]:(start[i]+7)],1,sum))}

#-----------------------------------------------------------------------------#

sink("Stan_Marginalized_owls.stan")
cat("
    
    data {
    int<lower=1,upper=4> Y[158,176];
    real RF [158,21];
    real FOR [158,21];
    matrix[176,7] pX[158];
    }
    
    parameters {
    real<lower=0,upper=1> psi0_BO; // initial barred owl (BO) occupancy
    real<lower=0,upper=1> psi0_NSO; // initial northern spotted owl (NSO) occupancy
    real<lower=-5,upper=5> gamNSO_int [21]; // time varying NSO colonization intercept
    real<lower=-5,upper=5> gamNSO_b; //Forest effect on NSO colonization
    real<lower=-5,upper=5> epsNSO_int; // constant NSO intercept on extinction
    real<lower=-5,upper=5> epsNSO_b; // BO effect on NSO extinction
    real<lower=-5,upper=5> gamBO_int; // constant BO intercept on colonization
    real<lower=-5,upper=5> gamBO_RF; //Riparian forest effect, 
    real<lower=-15, upper=15> gamBO_autolog;//autologistic effect, 
    real<lower=-5, upper=5> gamBO_xNSO; //and NSO effect on BO colonization
    real<lower=-5,upper=5> epsBO_int; // constant BO intercept on extinction
    real<lower=-5,upper=5> epsBO_RF; // Riparian forest effect, 
    real<lower=-15,upper=15> epsBO_autolog; //autologistic effect, 
    real<lower=-5,upper=5> epsBO_xNSO; //and NSO effect on BO extinction
    vector<lower=-5,upper=5>[6] pNSO_b; // intercept, day, night, method1, method3 and BO effects on NSO detection
    vector<lower=-5,upper=5>[7] pBO_b; // intercept, day, night, method1, method3, second half, and survey length effects on BO detection
    }
    
    transformed parameters {
    simplex [4] unc_pz; // initial probability of being in any state
    vector<lower=0,upper=1>[176] pBO;
    vector<lower=0,upper=1>[176] pNSO_BO;
    vector<lower=0,upper=1>[176] pNSO_bo;
    real<lower=0,upper=1> gamNSO;
    real<lower=0,upper=1> epsNSO_BO;
    real<lower=0,upper=1> epsNSO_bo;
    real<lower=0,upper=1> gamBO_NSO;
    real<lower=0,upper=1> gamBO_nso;
    real<lower=0,upper=1> epsBO_NSO;
    real<lower=0,upper=1> epsBO_nso;
    simplex [4] p [176,4];
    real<lower=0,upper=1> cp [158,22,4];//cumulative probability of observation history (8 occasions a year for 22 years for a total of 176 occasions)
    simplex [4] tr [4];
    real<lower=0,upper=1> pz [158,22,4];
    real<lower=0,upper=1> BO [158];
    real<lower=-0.5,upper=0.5> st_psiBO;
    vector[4] temp;
    
    unc_pz[1]=(1-psi0_BO)*(1-psi0_NSO);//unc_pz = probability of neither species being present in the study are at the start of the study.
    unc_pz[2]=(1-psi0_BO)*psi0_NSO;
    unc_pz[3]=psi0_BO*(1-psi0_NSO);
    unc_pz[4]=psi0_BO*psi0_NSO;
    
    for (k in 1:158){
    pBO=pX[k,,1].*inv_logit(pX[k,,]*pBO_b);
    pNSO_BO=pX[k,,1].*inv_logit(pX[k,,1:5]*pNSO_b[1:5]+pNSO_b[6]); // NSO detection probability if BO is present
    pNSO_bo=pX[k,,1].*inv_logit(pX[k,,1:5]*pNSO_b[1:5]);// NSO detection probability if BO is absent
    for (j in 1:176){
    p[j,1,1]=1; //if two species absent, then of course state of detection will be 1 (i.e., no species will be detected)
    p[j,1,2]=0;
    p[j,1,3]=0;
    p[j,1,4]=0;
    p[j,2,1]=1-pNSO_bo[j]; // if NSO present, probability of it not being detected
    p[j,2,2]=pNSO_bo[j]; // if NSO present, probability of it being detected
    p[j,2,3]=0;
    p[j,2,4]=0;
    p[j,3,1]=1-pBO[j];// if BO present, probability of it not being detected
    p[j,3,2]=0;
    p[j,3,3]=pBO[j];// if BO present, probability of it being detected
    p[j,3,4]=0;
    p[j,4,1]=1-(pBO[j]+pNSO_BO[j]-pBO[j]*pNSO_BO[j]);//I probably would have written this as : (1-pBO)*(1-pNSO_BO)*visited
    p[j,4,2]=pNSO_BO[j]*(1-pBO[j]); // if both species present, probability of only detecting BO
    p[j,4,3]=(1-pNSO_BO[j])*pBO[j];// if both species present, probability of only detecting NSO
    p[j,4,4]=pNSO_BO[j]*pBO[j];// probability of detecting both species
    }
    for (t in 1:22){
    //Note: 22 years and 8 occasions per year?  So, this is looking at the first occasion in each year (i.e., t= 1,9,17, etc)
    cp[k,t,1]=p[(1+(t-1)*8),1,Y[k,(1+(t-1)*8)]];//what is the probability that site is in state 1-4 given state is observed in state Y[k,,] 
    cp[k,t,2]=p[(1+(t-1)*8),2,Y[k,(1+(t-1)*8)]];
    cp[k,t,3]=p[(1+(t-1)*8),3,Y[k,(1+(t-1)*8)]];
    cp[k,t,4]=p[(1+(t-1)*8),4,Y[k,(1+(t-1)*8)]];
    for (j in 2:8){
    cp[k,t,1]=cp[k,t,1]*p[(j+(t-1)*8),1,Y[k,(j+(t-1)*8)]];//update probabilities given subsequent visits within the same year 
    cp[k,t,2]=cp[k,t,2]*p[(j+(t-1)*8),2,Y[k,(j+(t-1)*8)]];
    cp[k,t,3]=cp[k,t,3]*p[(j+(t-1)*8),3,Y[k,(j+(t-1)*8)]];
    cp[k,t,4]=cp[k,t,4]*p[(j+(t-1)*8),4,Y[k,(j+(t-1)*8)]];				}
    }
    }
    
    for (k in 1:158){
    pz[k,1,1]=unc_pz[1]*cp[k,1,1];//for first year, update probability of not being occupied by either species (kinda like a prior) by the observation history in the first year
    pz[k,1,2]=unc_pz[2]*cp[k,1,2];
    pz[k,1,3]=unc_pz[3]*cp[k,1,3];
    pz[k,1,4]=unc_pz[4]*cp[k,1,4];
    BO[k]=sum(pz[k,1,3:4])/sum(pz[k,1,]);//this is the probability of a BO being present on the site
    }
    
    epsNSO_BO=inv_logit(epsNSO_int+epsNSO_b);//probability of extinction of NSO given BO present
    epsNSO_bo=inv_logit(epsNSO_int);//probability of extinction of NSO given BO absent
    
    for (t in 1:21){
    st_psiBO=(sum(BO)/158-.5);//mean probability of BO presence (subtract 0.5 to center it)
    for (k in 1:158){
    gamNSO=inv_logit(gamNSO_int[t]+gamNSO_b*FOR[k,t]);//probability of colonization by NSO is a function of forest cover?
    gamBO_NSO=inv_logit(gamBO_int+gamBO_RF*RF[k,t]+gamBO_autolog*st_psiBO+gamBO_xNSO);//probability of BO colonization given NSO present (think this is a function of riparian forest, and autologistic effects)
    gamBO_nso=inv_logit(gamBO_int+gamBO_RF*RF[k,t]+gamBO_autolog*st_psiBO);
    epsBO_NSO=inv_logit(epsBO_int+epsBO_RF*RF[k,t]+epsBO_autolog*st_psiBO+epsBO_xNSO);//
    epsBO_nso=inv_logit(epsBO_int+epsBO_RF*RF[k,t]+epsBO_autolog*st_psiBO);
    tr[1,1]=(1-gamNSO)*(1-gamBO_nso);
    tr[1,2]=gamNSO*(1-gamBO_nso);
    tr[1,3]=(1-gamNSO)*gamBO_nso;
    tr[1,4]=gamNSO*gamBO_nso;
    tr[2,1]=epsNSO_bo*(1-gamBO_NSO);
    tr[2,2]=(1-epsNSO_bo)*(1-gamBO_NSO);
    tr[2,3]=epsNSO_bo*gamBO_NSO;
    tr[2,4]=(1-epsNSO_bo)*gamBO_NSO;
    tr[3,1]=(1-gamNSO)*epsBO_nso;
    tr[3,2]=gamNSO*epsBO_nso;
    tr[3,3]=(1-gamNSO)*(1-epsBO_nso);
    tr[3,4]=gamNSO*(1-epsBO_nso);
    tr[4,1]=epsNSO_BO*epsBO_NSO;
    tr[4,2]=(1-epsNSO_BO)*epsBO_NSO;
    tr[4,3]=epsNSO_BO*(1-epsBO_NSO);
    tr[4,4]=(1-epsNSO_BO)*(1-epsBO_NSO);
    for(j in 1:4)	{
    temp[j]=pz[k,t,j]*tr[j,1];}
    pz[k,(t+1),1]=sum(temp)*cp[k,(t+1),1];
    for(j in 1:4)	{
    temp[j]=pz[k,t,j]*tr[j,2];}
    pz[k,(t+1),2]=sum(temp)*cp[k,(t+1),2];
    for(j in 1:4)	{
    temp[j]=pz[k,t,j]*tr[j,3];}
    pz[k,(t+1),3]=sum(temp)*cp[k,(t+1),3];
    for(j in 1:4)	{
    temp[j]=pz[k,t,j]*tr[j,4];}
    pz[k,(t+1),4]=sum(temp)*cp[k,(t+1),4];
    BO[k]=sum(pz[k,(t+1),3:4])/sum(pz[k,(t+1),]);
    }
    }
    }
    
    
    model {
    for (k in 1:158) {  
    target +=log(sum(pz[k,22,])); }
    }
    
    
    ",fill=TRUE)
sink()  

#-----------------------------------------------------------------------------#
#Run Stan model:

owls_SM.data<-list (Y=array(Y,dim=c(dim(Y)[1],dim(Y)[2])),RF=array(as.vector(unlist(RF)),dim=c(dim(RF)[1],dim(RF)[2])),
                    FOR=array(as.vector(unlist(FOR)),dim=c(dim(FOR)[1],dim(FOR)[2])),pX=array(pX,dim=c(dim(pX)[1],dim(pX)[2],dim(pX)[3])))

owls_SM.par<-c('pBO_b','pNSO_b','psi0_BO','psi0_NSO','gamNSO_b','epsNSO_b','epsNSO_int','gamBO_int','epsBO_int','gamBO_RF','gamBO_autolog', 
               'gamBO_xNSO', 'epsBO_RF', 'epsBO_autolog','epsBO_xNSO', 'gamNSO_int')

owls_SM.out<- stan(".\\Stan_Marginalized_owls.stan",data=owls_SM.data,pars=owls_SM.par,chains = 1,iter=10,seed=1)


###############################################################################
###############################################################################
#              Dynamic two-species occupancy model with autologistic effects (Barred owls and Northern Spotted owls)
#
#                               Stan Marginalized Version with random effects 
#                                    and multi-level R^2 calculation
#
#-----------------------------------------------------------------------------#

library(rstan)

#To run Stan in parallel:
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


#Get data:
setwd("U:\\Marginalized likelihood\\BO_NSO")
RF<-read.csv(".\\RF.csv",header=FALSE)
FOR<-read.csv(".\\FOR.csv",header=FALSE)
visited<-read.csv(".\\visited.csv",header=FALSE)


#-----------------------------------------------------------------------------#
# format data for model fitting:

# Stan-MRE - 4 states: 1 - Neither present, 2 - NSO present, 3 - BO present, 4 - both species present
# requires list of visited sites & notvisited sites, matrix of visits per site/year with 0 visit for unvisited sites, and y, 

pX<-array(NA,dim=c(158,176,7))# detection covariates: day, night, method1, method3, second half, and survey length effects on BO detection
pX[,,1]<-matrix(unlist(visited),nrow=dim(visited)[1],ncol=dim(visited)[2])
pX[,,2]<-as.matrix(read.csv(".\\day.csv",header=FALSE))
pX[,,3]<-as.matrix(read.csv(".\\night.csv",header=FALSE))
pX[,,4]<-as.matrix(read.csv(".\\mtd1.csv",header=FALSE))
pX[,,5]<-as.matrix(read.csv(".\\mtd3.csv",header=FALSE))
pX[,1:88,6]<-0
pX[,89:176,6]<-1
pX[,,7]<-as.matrix(read.csv(".\\tott.csv",header=FALSE))
#If site is not visited, change all covariates to zero:
for(i in 1:158){for(j in 1:176){ if(visited[i,j]==0){pX[i,j,]<-0}}}  
owls<-as.matrix(read.csv(".\\owls.csv",header=FALSE))
Y_BO<-ifelse(owls>1,1,0)
Y_NSO<-ifelse(owls==1|owls==3,1,0)
Y<-owls+1

n_occ<-apply(visited[,1:8],1,sum)
start<-1:22*8-7
for(i in 2:length(start)){n_occ<-cbind(n_occ,apply(visited[,start[i]:(start[i]+7)],1,sum))}

#Mean values required for multi-level R^2 calculation:
mean_FOR_site<-apply(FOR,1,mean)
mean_FOR_time<-apply(FOR,2,mean)
mean_RF_site<-apply(RF,1,mean)
mean_RF_time<-apply(RF,2,mean)

#-----------------------------------------------------------------------------#

sink("Stan_Marginalized_owls_with_RE.stan")
cat("
    data {
    int<lower=1,upper=4> Y[158,176];
    real RF [158,21];
    real FOR [158,21];
    matrix[176,7] pX[158];
    vector[158] mean_FOR_site; // for multi-level R2 calculation
    vector[158] mean_RF_site; // for multi-level R2 calculation
    vector[21] mean_FOR_time; // for multi-level R2 calculation
    vector[21] mean_RF_time; // for multi-level R2 calculation
    
    }
    
    parameters {
    real<lower=0, upper=1> psi0_BO;
    real<lower=0, upper=1> psi0_NSO;
    
    real<lower=-5,upper=5> gamNSO_int; // time varying NSO colonization intercept
    real<lower=-5,upper=5> gamNSO_b; //Forest effect on NSO colonization
    vector[158] z_gamNSO_site;
    vector[21] z_gamNSO_time;
    real<lower=0, upper=5> sd_gamNSO_site;
    real<lower=0, upper=5> sd_gamNSO_time;
    
    real<lower=-5,upper=5> epsNSO_int; // constant NSO intercept on extinction
    real<lower=-5,upper=5> epsNSO_b; // BO effect on NSO extinction
    vector[158] z_epsNSO_site;
    vector[21] z_epsNSO_time;
    real<lower=0, upper=5> sd_epsNSO_site;
    real<lower=0, upper=5> sd_epsNSO_time;
    
    real<lower=-5,upper=5> gamBO_int; // constant BO intercept on colonization
    real<lower=-5, upper=5> gamBO_RF; //Riparian forest effect, 
    real<lower=-15, upper=15> gamBO_autolog; //autologistic effect
    real<lower=-5, upper=5> gamBO_xNSO; //NSO effect on BO colonization
    vector[158] z_gamBO_site;
    vector[21] z_gamBO_time;
    real<lower=0, upper=5> sd_gamBO_site;
    real<lower=0, upper=5> sd_gamBO_time;
    
    real<lower=-15,upper=15> epsBO_int; // constant BO intercept on extinction
    real<lower=-5, upper=5> epsBO_RF; //Riparian forest effect, 
    real<lower=-15, upper=15> epsBO_autolog; //autologistic effect
    real<lower=-5, upper=5> epsBO_xNSO; //NSO effect on BO extinction
    vector[158] z_epsBO_site;
    vector[21] z_epsBO_time;
    real<lower=0, upper=5> sd_epsBO_site;
    real<lower=0, upper=5> sd_epsBO_time;
    
    vector<lower=-5,upper=5>[6] pNSO_b; // intercept, day, night, method1, method3 and BO effects on NSO detection
    vector<lower=-5,upper=5>[7] pBO_b; // intercept, day, night, method1, method3, second half, and survey length effects on BO detection
    }
    
    transformed parameters {
    simplex [4] unc_pz; // initial probability of being in any state
    vector<lower=0,upper=1>[176] pBO;
    vector<lower=0,upper=1>[176] pNSO_BO;
    vector<lower=0,upper=1>[176] pNSO_bo;
    real<lower=0,upper=1> gamNSO;
    real<lower=0,upper=1> epsNSO_BO;
    real<lower=0,upper=1> epsNSO_bo;
    real<lower=0,upper=1> gamBO_NSO;
    real<lower=0,upper=1> gamBO_nso;
    real<lower=0,upper=1> epsBO_NSO;
    real<lower=0,upper=1> epsBO_nso;
    simplex [4] p [176,4];
    real<lower=0,upper=1> cp [158,22,4];//cumulative probability of observation history (8 occasions a year for 22 years for a total of 176 occasions)
    simplex [4] tr [4];
    real<lower=0,upper=1> pz [158,22,4];
    real<lower=0,upper=1> BO [158];
    real<lower=0,upper=1> NSO [158];
    vector<lower=-0.5,upper=0.5>[21] st_psiBO;
    vector<lower=0,upper=1>[21] psiNSO_time_mean;
    vector<lower=0,upper=1>[21] psiBO_time_mean;
    vector<lower=0, upper=1>[4] temp;
    
    //Unconditional occupancy probability (used to help determine occupancy state at t=1)
    unc_pz[1]=(1-psi0_BO)*(1-psi0_NSO);//unc_pz = probability of neither species being present in the study are at the start of the study.
    unc_pz[2]=(1-psi0_BO)*psi0_NSO;
    unc_pz[3]=psi0_BO *(1-psi0_NSO);
    unc_pz[4]=psi0_BO *psi0_NSO;
    
    //Observation model
    //22 years and 8 occasions per year
    for (k in 1:158){
    pBO=pX[k,,1].*inv_logit(pX[k,,]*pBO_b);
    pNSO_BO=pX[k,,1].*inv_logit(pX[k,,1:5]*pNSO_b[1:5]+pNSO_b[6]); // NSO detection probability if BO is present
    pNSO_bo=pX[k,,1].*inv_logit(pX[k,,1:5]*pNSO_b[1:5]);// NSO detection probability if BO is absent
    
    for (j in 1:176){
    p[j,1,1]=1; //if two species absent, then of course state of detection will be 1 (i.e., no species will be detected)
    p[j,1,2]=0;
    p[j,1,3]=0;
    p[j,1,4]=0;
    p[j,2,1]=1-pNSO_bo[j]; // if NSO present, probability of it not being detected
    p[j,2,2]=pNSO_bo[j]; // if NSO present, probability of it being detected
    p[j,2,3]=0;
    p[j,2,4]=0;
    p[j,3,1]=1-pBO[j];// if BO present, probability of it not being detected
    p[j,3,2]=0;
    p[j,3,3]=pBO[j];// if BO present, probability of it being detected
    p[j,3,4]=0;
    p[j,4,1]=1-(pBO[j]+pNSO_BO[j]-pBO[j]*pNSO_BO[j]);//I probably would have written this as : (1-pBO)*(1-pNSO_BO)*visited
    p[j,4,2]=pNSO_BO[j]*(1-pBO[j]); // if both species present, probability of only detecting BO
    p[j,4,3]=(1-pNSO_BO[j])*pBO[j];// if both species present, probability of only detecting NSO
    p[j,4,4]=pNSO_BO[j]*pBO[j];// probability of detecting both species
    }
    
    for (t in 1:22){
    
    cp[k,t,1]=p[(1+(t-1)*8),1,Y[k,(1+(t-1)*8)]];//what is the probability that site is in state 1-4 given state is observed in state Y[k,,] 
    cp[k,t,2]=p[(1+(t-1)*8),2,Y[k,(1+(t-1)*8)]];
    cp[k,t,3]=p[(1+(t-1)*8),3,Y[k,(1+(t-1)*8)]];
    cp[k,t,4]=p[(1+(t-1)*8),4,Y[k,(1+(t-1)*8)]];
    
    for (j in 2:8){
    cp[k,t,1]=cp[k,t,1]*p[(j+(t-1)*8),1,Y[k,(j+(t-1)*8)]];//update probabilities given subsequent visits within the same year 
    cp[k,t,2]=cp[k,t,2]*p[(j+(t-1)*8),2,Y[k,(j+(t-1)*8)]];
    cp[k,t,3]=cp[k,t,3]*p[(j+(t-1)*8),3,Y[k,(j+(t-1)*8)]];
    cp[k,t,4]=cp[k,t,4]*p[(j+(t-1)*8),4,Y[k,(j+(t-1)*8)]];}
    }
    }
    
    //Process model : state at t=1 
    for (k in 1:158){
    pz[k,1,1]=unc_pz[1]*cp[k,1,1];//for first year, update probability of not being occupied by either species (kinda like a prior) by the observation history in the first year
    pz[k,1,2]=unc_pz[2]*cp[k,1,2];
    pz[k,1,3]=unc_pz[3]*cp[k,1,3];
    pz[k,1,4]=unc_pz[4]*cp[k,1,4];
    BO[k]=sum(pz[k,1,3:4])/sum(pz[k,1,]);//this is the probability of a BO being present on the site
    NSO[k]=(pz[k,1,2]+pz[k,1,4])/sum(pz[k,1,]);
    }
    
    for (t in 1:21){
    psiNSO_time_mean[t]=sum(NSO)/158;//mean probability of NSO presence
    psiBO_time_mean[t]=sum(BO)/158;//mean probability of NSO presence
    st_psiBO[t]=psiBO_time_mean[t]-0.5;//mean probability of BO presence (subtract 0.5 to center it)
    
    //Process model : extinction/colonization dynamics 
    for (k in 1:158){
    epsNSO_BO=inv_logit(epsNSO_int+epsNSO_b+z_epsNSO_site[k]*sd_epsNSO_site+z_epsNSO_time[t]*sd_epsNSO_time);//probability of extinction of NSO given BO present
    epsNSO_bo=inv_logit(epsNSO_int+z_epsNSO_site[k]*sd_epsNSO_site+z_epsNSO_time[t]*sd_epsNSO_time);//probability of extinction of NSO given BO absent
    gamNSO=inv_logit(gamNSO_int + gamNSO_b*FOR[k,t] +z_gamNSO_site[k]*sd_gamNSO_site+z_gamNSO_time[t]*sd_gamNSO_time);//probability of colonization by NSO is a function of forest cover?
    gamBO_NSO=inv_logit(gamBO_int+gamBO_RF*RF[k,t]+gamBO_autolog*st_psiBO[t]+gamBO_xNSO+ z_gamBO_site[k]*sd_gamBO_site+z_gamBO_time[t]*sd_gamBO_time);//probability of BO colonization given NSO present (think this is a function of riparian forest, and autologistic effects)
    gamBO_nso=inv_logit(gamBO_int+gamBO_RF*RF[k,t]+gamBO_autolog*st_psiBO[t]+ z_gamBO_site[k]*sd_gamBO_site+z_gamBO_time[t]*sd_gamBO_time);
    epsBO_NSO=inv_logit(epsBO_int+epsBO_RF*RF[k,t]+epsBO_autolog*st_psiBO[t]+epsBO_xNSO+z_epsBO_site[k]*sd_epsBO_site+z_epsBO_time[t]*sd_epsBO_time);
    epsBO_nso=inv_logit(epsBO_int+epsBO_RF*RF[k,t]+epsBO_autolog*st_psiBO[t]+z_epsBO_site[k]*sd_epsBO_site+z_epsBO_time[t]*sd_epsBO_time);
    tr[1,1]=(1-gamNSO)*(1-gamBO_nso);
    tr[1,2]=gamNSO*(1-gamBO_nso);
    tr[1,3]=(1-gamNSO)*gamBO_nso;
    tr[1,4]=gamNSO*gamBO_nso;
    tr[2,1]=epsNSO_bo*(1-gamBO_NSO);
    tr[2,2]=(1-epsNSO_bo)*(1-gamBO_NSO);
    tr[2,3]=epsNSO_bo*gamBO_NSO;
    tr[2,4]=(1-epsNSO_bo)*gamBO_NSO;
    tr[3,1]=(1-gamNSO)*epsBO_nso;
    tr[3,2]=gamNSO*epsBO_nso;
    tr[3,3]=(1-gamNSO)*(1-epsBO_nso);
    tr[3,4]=gamNSO*(1-epsBO_nso);
    tr[4,1]=epsNSO_BO*epsBO_NSO;
    tr[4,2]=(1-epsNSO_BO)*epsBO_NSO;
    tr[4,3]=epsNSO_BO*(1-epsBO_NSO);
    tr[4,4]=(1-epsNSO_BO)*(1-epsBO_NSO);
    
    for(j in 1:4){
    temp[j]= pz[k,t,j]*tr[j,1];}
    
    pz[k,(t+1),1]=sum(temp)*cp[k,(t+1),1];
    
    for(j in 1:4){
    temp[j]= pz[k,t,j]*tr[j,2];}
    
    pz[k,(t+1),2]=sum(temp)*cp[k,(t+1),2];
    
    for(j in 1:4){
    temp[j]= pz[k,t,j]*tr[j,3];}
    
    pz[k,(t+1),3]=sum(temp)*cp[k,(t+1),3];
    
    for(j in 1:4){
    temp[j]= pz[k,t,j]*tr[j,4];}
    
    pz[k,(t+1),4]=sum(temp)*cp[k,(t+1),4];
    
    //Determine the percentage of sites occupied by Barred Owls (for use in autologistic extinction/colonization model)
    BO[k]=sum(pz[k,(t+1),3:4])/sum(pz[k,(t+1),]);
    //Determine the percentage of sites occupied by Northern Spotted Owls (for use in multi-level R^2 calculation)
    NSO[k]=(pz[k,(t+1),2]+pz[k,(t+1),4])/sum(pz[k,(t+1),]);
    }
    }
    }
    
    model {
    
    z_gamBO_site~normal(0,1);
    z_epsBO_site~normal(0,1);
    z_gamNSO_site~normal(0,1);
    z_epsNSO_site~normal(0,1);
    z_gamBO_time~normal(0,1);
    z_epsBO_time~normal(0,1);
    z_gamNSO_time~normal(0,1);
    z_epsNSO_time~normal(0,1);
    
    for (k in 1:158) {  
    target +=log(sum(pz[k,22,])); }}
    
    generated quantities{
    
    real psiBO_site[158,21];
    real psiNSO_site[158,21];
    vector<lower=0,upper=1>[158] psiBO_site_mean;
    vector<lower=0,upper=1>[158] psiNSO_site_mean;
    
    vector[158] gamBO_site_pred;
    vector[21] gamBO_time_pred;
    vector[158] epsBO_site_pred;
    vector[21] epsBO_time_pred;
    vector[158] gamNSO_site_pred;
    vector[21] gamNSO_time_pred;
    vector[158] epsNSO_site_pred;
    vector[21] epsNSO_time_pred;
    
    real R2_gamNSO_site;
    real R2_gamNSO_time;
    real R2_epsNSO_site;
    real R2_epsNSO_time;
    real R2_gamBO_site;
    real R2_gamBO_time;
    real R2_epsBO_site;
    real R2_epsBO_time;
    
    for(k in 1:158){
    for (t in 1:21) {
    psiBO_site[k,t] = sum(pz[k,t,3:4])/sum(pz[k,t,]);
    psiNSO_site[k,t] = (pz[k,t,2] + pz[k,t,4])/sum(pz[k,t,]);}
    psiBO_site_mean[k] = mean(psiBO_site[k,1:21]);
    psiNSO_site_mean[k] = mean(psiNSO_site[k,1:21]);}
    
    //excluding autologistic effect on spatial RE because the average autologistic effect for all 21 intervals is the same for each site:
    gamBO_site_pred = gamBO_int+ gamBO_RF*mean_RF_site+gamBO_xNSO*psiNSO_site_mean+ z_gamBO_site*sd_gamBO_site;
    gamBO_time_pred = gamBO_int+ gamBO_RF*mean_RF_time+gamBO_autolog*st_psiBO+gamBO_xNSO*psiNSO_time_mean+ z_gamBO_time*sd_gamBO_time;
    R2_gamBO_site = 1-variance(gamBO_int+z_gamBO_site*sd_gamBO_site)/variance(gamBO_site_pred);
    R2_gamBO_time = 1-variance(gamBO_int+z_gamBO_time*sd_gamBO_time)/variance(gamBO_time_pred);
    
    //excluding autologistic effect on spatial RE because the average autologistic effect for all 21 intervals is the same for each site:
    epsBO_site_pred = epsBO_int+epsBO_RF*mean_RF_site+epsBO_xNSO*psiNSO_site_mean+ z_epsBO_site*sd_epsBO_site;
    epsBO_time_pred = epsBO_int+epsBO_RF*mean_RF_time+epsBO_autolog*st_psiBO+epsBO_xNSO*psiNSO_time_mean+ z_epsBO_time*sd_epsBO_time;
    R2_epsBO_site = 1-variance(epsBO_int+z_epsBO_site*sd_epsBO_site)/variance(epsBO_site_pred);
    R2_epsBO_time = 1-variance(epsBO_int+z_epsBO_time*sd_epsBO_time)/variance(epsBO_time_pred);
    
    gamNSO_site_pred = gamNSO_int+gamNSO_b*mean_FOR_site+ z_gamNSO_site*sd_gamNSO_site;
    gamNSO_time_pred = gamNSO_int+gamNSO_b*mean_FOR_time+ z_gamNSO_time*sd_gamNSO_time;
    R2_gamNSO_site = 1-variance(gamNSO_int+z_gamNSO_site*sd_gamNSO_site)/variance(gamNSO_site_pred);
    R2_gamNSO_time = 1-variance(gamNSO_int+z_gamNSO_time*sd_gamNSO_time)/variance(gamNSO_time_pred);
    
    epsNSO_site_pred =epsNSO_int+epsNSO_b*psiBO_site_mean+z_epsNSO_site*sd_epsNSO_site;
    epsNSO_time_pred =epsNSO_int+epsNSO_b*psiBO_time_mean+z_epsNSO_time*sd_epsNSO_time;
    R2_epsNSO_site = 1-variance(epsNSO_int+z_epsNSO_site*sd_epsNSO_site)/variance(epsNSO_site_pred);
    R2_epsNSO_time = 1-variance(epsNSO_int+z_epsNSO_time*sd_epsNSO_time)/variance(epsNSO_time_pred);
    
    }
    
    
    
    ",fill=TRUE)
sink()  


#-----------------------------------------------------------------------------#
#Run Stan models:

owls_SMRE.data<-list (Y=array(Y,dim=c(dim(Y)[1],dim(Y)[2])),RF=array(as.vector(unlist(RF)),dim=c(dim(RF)[1],dim(RF)[2])),
                      FOR=array(as.vector(unlist(FOR)),dim=c(dim(FOR)[1],dim(FOR)[2])),pX=array(pX,dim=c(dim(pX)[1],dim(pX)[2],dim(pX)[3])),
                      mean_FOR_site= mean_FOR_site, mean_FOR_time=mean_FOR_time, mean_RF_site=mean_RF_site, mean_RF_time=mean_RF_time)

owls_SMRE.par<-c('pBO_b','pNSO_b','psi0_BO','psi0_NSO','gamNSO_b','epsNSO_b','epsNSO_int','gamBO_int','epsBO_int','gamBO_RF','gamBO_autolog', 
                 'gamBO_xNSO', 'epsBO_RF', 'epsBO_autolog','epsBO_xNSO', 'gamNSO_int', 'R2_gamBO_site', 'R2_gamBO_time', 'R2_gamNSO_site','R2_gamNSO_time','R2_epsBO_site', 
                 'R2_epsBO_time', 'R2_epsNSO_site' , 'R2_epsNSO_time')

owls_SMRE.out<- stan(".\\Stan_Marginalized_owls_with_RE.stan",data=owls_SMRE.data,pars=owls_SMRE.par,chains = 3,iter=1000,seed=1)







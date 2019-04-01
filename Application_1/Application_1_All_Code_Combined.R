

# add functions

# need notes about setting time for summary function to work

# add mcmc settings that are described in paper

# add additional code from CY...




###############################################################################
#                                                                     Spring 19
#  Fitting a multi-state version of a CJS model to the RBT data 
#  Discrete JAGS version - fixed time effects
#
#  Notes:
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2jags)

# data.dir = paste0(getwd(), "/Data")  # Need to set directory for data
CH = as.matrix(read.table(file = paste0(data.dir, "/RBT_Capture_History.txt"),
                          header = FALSE, sep = "\t"))
#-----------------------------------------------------------------------------#
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
sink("JAGS_Discrete_Time.jags")
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

# ni <- 1000
# nt <- 1
# nb <- 500

JD.out <- jags.parallel(JD.data, inits = NULL, JD.par, "JAGS_Discrete_Time.jags",
                        n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb,
                        export_obj_names = c("ni", "nt", "nb"))

#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Fitting a multi-state version of a CJS model to the RBT data 
#  Marginalized JAGS version - fixed time effects
#
#  Notes:
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2jags)

# data.dir = paste0(getwd(), "/Data") # Need to set directory for data
CH = as.matrix(read.table(file = paste0(data.dir, "/RBT_Capture_History.txt"),
                          header = FALSE, sep = "\t"))
#-----------------------------------------------------------------------------#
# format data for model fitting
tmpCH = collapse.ch(CH)[[1]]
FR = collapse.ch(CH)[[2]]

# Create vector with occasion of marking
get.first <- function(x) min(which(x != 0))
sumf <- apply(tmpCH, 1, get.first)

S = tmpCH
S[S[,] == 0] = 2

NS = nrow(S)          # number of capture histories 
Nint = ncol(S) - 1    # number of sampling intervals

ones <- FR

zeta = known.state.cjs(tmpCH)

#-----------------------------------------------------------------------------#
sink("JAGS_Marginalized_Time.jags")
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
  
  for(i in 1:NS){
    zeta[i,sumf[i],1] <- 1
    zeta[i,sumf[i],2] <- 0
    
    for(t in sumf[i]:Nint){
      zeta[i,(t+1),1] <- inprod(zeta[i,t,], omega[,1,t]) * rho[1,S[i,(t+1)],t]
      zeta[i,(t+1),2] <- inprod(zeta[i,t,], omega[,2,t]) * rho[2,S[i,(t+1)],t]
    }
    
    lik[i] <- sum(zeta[i,(Nint+1),])
    ones[i] ~ dbin(lik[i], FR[i])
  }
}
    ", fill = TRUE)
sink() 
#-----------------------------------------------------------------------------#
JM.data <- list(NS = NS, Nint = Nint, S = S,
                sumf = sumf, FR = FR, ones = ones)
JM.par <- c('phi', 'p')

# ni <- 1000
# nt <- 1
# nb <- 500

JM.out <- jags.parallel(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_Time.jags",
                        n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb,
                        export_obj_names = c("ni", "nt", "nb"))

#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Fitting a multi-state version of a CJS model to the RBT data 
#  Discrete OpenBUGS version - fixed time effects
#
#  Notes:
#  * OpenBUGS will crash with indexing of omega and rho like other CJS examples,
#    see note below. Same issue as WinBUGS. 
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2OpenBUGS)

# data.dir = paste0(getwd(), "/Data")
CH = as.matrix(read.table(file = paste0(data.dir, "/RBT_Capture_History.txt"),
                          header = FALSE, sep = "\t"))
#-----------------------------------------------------------------------------#
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
sink("OpenBUGS_Discrete_Time.OpenBUGS")
cat("
model{
  
  for(t in 1:Nint){
    phi[t] ~ dunif(0,1)
    p[t] ~ dunif(0,1)
    
    # OpenBUGS crashes if indexed with t as the third dimension
    # omega[1,1,t] <- phi[t]
    # omega[1,2,t] <- (1 - phi[t])
    # omega[2,1,t] <- 0
    # omega[2,2,t] <- 1
    #
    # rho[1,1,t] <- p[t]
    # rho[1,2,t] <- (1 - p[t])
    # rho[2,1,t] <- 0
    # rho[2,2,t] <- 1
    
    omega[1,t,1] <- phi[t]
    omega[1,t,2] <- (1 - phi[t])
    omega[2,t,1] <- 0
    omega[2,t,2] <- 1
    
    rho[1,t,1] <- p[t]
    rho[1,t,2] <- (1 - p[t])
    rho[2,t,1] <- 0
    rho[2,t,2] <- 1
  }
  
  for(i in 1:NY){
    Z[i,indf[i]] <- 1
    
    for(t in indf[i]:Nint){
      # Z[i,(t + 1)] ~ dcat(omega[Z[i, t], , t])
      # Y[i,(t + 1)] ~ dcat(rho[Z[i, (t + 1)], , t])
      
      Z[i,(t+1)] ~ dcat(omega[Z[i, t], t, ])
      Y[i,(t+1)] ~ dcat(rho[Z[i, (t+1)], t, ])
    }
  }
}
    ", fill = TRUE)
sink()    

#-----------------------------------------------------------------------------#
OD.data <- list(NY = NY, Nint = Nint, Y = Y, indf = indf, Z = Z)
OD.par <- c('phi', 'p')

# ni <- 1000
# nt <- 1
# nb <- 500

OD.out <- R2OpenBUGS::bugs(OD.data, inits = NULL, OD.par, "OpenBUGS_Discrete_Time.OpenBUGS",
                           n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb)

#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Fitting a multi-state version of a CJS model to the RBT data 
#  Marginalized OpenBUGS version - fixed time effects
#
#  Notes:
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2OpenBUGS)

# data.dir = paste0(getwd(), "/Data")
CH = as.matrix(read.table(file = paste0(data.dir, "/RBT_Capture_History.txt"),
                          header = FALSE, sep = "\t"))
#-----------------------------------------------------------------------------#
# format data for model fitting
tmpCH = collapse.ch(CH)[[1]]
FR = collapse.ch(CH)[[2]]

# Create vector with occasion of marking
get.first <- function(x) min(which(x != 0))
sumf <- apply(tmpCH, 1, get.first)

S = tmpCH
S[S[,] == 0] = 2

NS = nrow(S)          # number of capture histories 
Nint = ncol(S) - 1    # number of sampling intervals

ones <- FR

zeta = known.state.cjs(tmpCH)

#-----------------------------------------------------------------------------#
sink("OpenBUGS_Marginalized_Time.OpenBUGS")
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
  
  for(i in 1:NS){
    zeta[i,sumf[i],1] <- 1
    zeta[i,sumf[i],2] <- 0
    
    for(t in sumf[i]:Nint){
      zeta[i,(t+1),1] <- inprod(zeta[i,t,], omega[,1,t]) * rho[1,S[i,(t+1)],t]
      zeta[i,(t+1),2] <- inprod(zeta[i,t,], omega[,2,t]) * rho[2,S[i,(t+1)],t]
    }
    
    lik[i] <- sum(zeta[i,(Nint+1),])
    ones[i] ~ dbin(lik[i], FR[i])
  }
}
    ", fill = TRUE)
sink() 
#-----------------------------------------------------------------------------#
OM.data <- list(NS = NS, Nint = Nint, S = S,
                sumf = sumf, FR = FR, ones = ones)
OM.par <- c('phi', 'p')

# ni <- 1000
# nt <- 1
# nb <- 500

# be explicit on the call to 'bugs', as the names are same for both R2WinBUGS and R2OpenBUGS
OM.out <- R2OpenBUGS::bugs(OM.data, inits = NULL, OM.par, "OpenBUGS_Marginalized_Time.OpenBUGS",
               n.chains = 3, n.iter = ni, n.thin = nt)

#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Fitting a multi-state version of a CJS model to the RBT data 
#  Discrete WinBUGS version - fixed time effects
#
#  Notes:
#  * WinBUGS will crash with indexing of omega and rho like other CJS examples,
#    see note below.
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2WinBUGS)

# data.dir = paste0(getwd(), "/Data")
CH = as.matrix(read.table(file = paste0(data.dir, "/RBT_Capture_History.txt"),
                          header = FALSE, sep = "\t"))
#-----------------------------------------------------------------------------#
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
    
    # WinBUGS crashes if indexed with t as the third dimension
    # omega[1,1,t] <- phi[t]
    # omega[1,2,t] <- (1 - phi[t])
    # omega[2,1,t] <- 0
    # omega[2,2,t] <- 1
    # 
    # rho[1,1,t] <- p[t]
    # rho[1,2,t] <- (1 - p[t])
    # rho[2,1,t] <- 0
    # rho[2,2,t] <- 1
    
    omega[1,t,1] <- phi[t]
    omega[1,t,2] <- (1 - phi[t])
    omega[2,t,1] <- 0
    omega[2,t,2] <- 1
    
    rho[1,t,1] <- p[t]
    rho[1,t,2] <- (1 - p[t])
    rho[2,t,1] <- 0
    rho[2,t,2] <- 1
  }
  
  for(i in 1:NY){
    Z[i,indf[i]] <- 1
    
    for(t in indf[i]:Nint){
      # Z[i,(t + 1)] ~ dcat(omega[Z[i, t], , t])
      # Y[i,(t + 1)] ~ dcat(rho[Z[i, (t + 1)], , t])
      
      Z[i,(t+1)] ~ dcat(omega[Z[i, t], t, ])
      Y[i,(t+1)] ~ dcat(rho[Z[i, (t+1)], t, ])
    }
  }
}
    ", fill = TRUE)
sink()    

#-----------------------------------------------------------------------------#
WD.data <- list(NY = NY, Nint = Nint, Y = Y, indf = indf, Z = Z)
WD.par <- c('phi', 'p')

# ni <- 1000
# nt <- 1
# nb <- 500

# be explicit on the call to 'bugs', as the names are same for both R2WinBUGS and R2OpenBUGS
WD.out <- R2WinBUGS::bugs(WD.data, inits = NULL, WD.par, "WinBUGS_Discrete_Time.WinBUGS",
                          n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb, debug = FALSE)

#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Fitting a multi-state version of a CJS model to the RBT data 
#  Marginalized WinBUGS version - fixed time effects
#
#  Notes:
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2WinBUGS)

# data.dir = paste0(getwd(), "/Data")
CH = as.matrix(read.table(file = paste0(data.dir, "/RBT_Capture_History.txt"),
                          header = FALSE, sep = "\t"))
#-----------------------------------------------------------------------------#
# format data for model fitting
tmpCH = collapse.ch(CH)[[1]]
FR = collapse.ch(CH)[[2]]

# Create vector with occasion of marking
get.first <- function(x) min(which(x != 0))
sumf <- apply(tmpCH, 1, get.first)

S = tmpCH
S[S[,] == 0] = 2

NS = nrow(S)          # number of capture histories 
Nint = ncol(S) - 1    # number of sampling intervals

ones <- FR

zeta = known.state.cjs(tmpCH)

#-----------------------------------------------------------------------------#
sink("WinBUGS_Marginalized_Time.WinBUGS")
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
  
  for(i in 1:NS){
    zeta[i,sumf[i],1] <- 1
    zeta[i,sumf[i],2] <- 0
    
    for(t in sumf[i]:Nint){
      zeta[i,(t+1),1] <- inprod(zeta[i,t,], omega[,1,t]) * rho[1,S[i,(t+1)],t]
      zeta[i,(t+1),2] <- inprod(zeta[i,t,], omega[,2,t]) * rho[2,S[i,(t+1)],t]
    }
    
    lik[i] <- sum(zeta[i,(Nint+1),])
    ones[i] ~ dbin(lik[i], FR[i])
  }
}
    ", fill = TRUE)
sink() 
#-----------------------------------------------------------------------------#
WM.data <- list(NS = NS, Nint = Nint, S = S,
                sumf = sumf, FR = FR, ones = ones)
WM.par <- c('phi', 'p')

# ni <- 1000
# nt <- 1
# nb <- 500

# be explicit on the call to 'bugs', as the names are same for both R2WinBUGS and R2OpenBUGS
WM.out <- R2WinBUGS::bugs(WM.data, inits = NULL, WM.par, "WinBUGS_Marginalized_Time.WinBUGS",
               n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb)

#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Fitting a Multi-state version of a CJS model to the RBT data 
#  Marginalized Stan version
#
#  Notes:
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(rstan)
# Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())  

# data.dir = paste0(getwd(), "/Data")
CH = as.matrix(read.table(file = paste0(data.dir, "/RBT_Capture_History.txt"),
                          header = FALSE, sep = "\t"))
#-----------------------------------------------------------------------------#
# format data for model fitting
tmpCH = collapse.ch(CH)[[1]]
sumFR = collapse.ch(CH)[[2]]

# Create vector with occasion of marking
get.first <- function(x) min(which(x != 0))
sumf <- apply(tmpCH, 1, get.first)

sumCH = tmpCH
sumCH[sumCH[,] == 0] = 2

NsumCH = nrow(sumCH)         # number of capture histories 
n.occasions = ncol(sumCH)    # number of sampling occasions

# Catch (for the N versions)
catch = colSums(CH)[2:18]

#-----------------------------------------------------------------------------#
# Some may prefer to work with the Stan model in a seperate tab within Rstudio,
# but the model code is included here for consistency. Note: '//' is for comments
# in Stan

sink("Stan_Marginalized_Time.stan")
cat("
// CJS Model with time varying survival (s) and capture probability (p)

data{
  int<lower = 1> NsumCH;
  int<lower = 1> n_occasions;
  int<lower = 1, upper = 2> sumCH[NsumCH, n_occasions];
  int<lower = 1> sumf[NsumCH];
  int<lower = 1> sumFR[NsumCH];
}

parameters{
  real<lower = 0, upper = 1> s[n_occasions - 1];   // 3 month survivals 
  real<lower = 0, upper = 1> p[n_occasions - 1];   // capture probability
}

transformed parameters{
  simplex[2] tr[2,n_occasions - 1];
  simplex[2] pmat[2,n_occasions - 1];
  
  for(k in 1:n_occasions - 1){
    tr[1,k,1] = s[k];
    tr[1,k,2] = (1 - s[k]);
    tr[2,k,1] = 0;
    tr[2,k,2] = 1;
    
    pmat[1,k,1] = p[k];
    pmat[1,k,2] = (1 - p[k]);
    pmat[2,k,1] = 0;
    pmat[2,k,2] = 1;
  }
}

model{
  vector[2] pz[n_occasions];
  
  p ~ uniform(0,1);
  
  for(i in 1:NsumCH){  
    pz[sumf[i],1] = 1;
    pz[sumf[i],2] = 0;
    
    for(k in (sumf[i] + 1):n_occasions){ 
      pz[k,1] = pz[k-1,1] * tr[1,k-1,1] * pmat[1,k-1,sumCH[i,(k)]];
      pz[k,2] = (pz[k-1, 1] * tr[1,k-1,2] + pz[k-1, 2]) * pmat[2,k-1,sumCH[i,(k)]];
    }  
    
    target += sumFR[i] * log(sum(pz[n_occasions])); 
  }
}

    ", fill = TRUE)
sink()

#-----------------------------------------------------------------------------#
# Time (occasion) specific
sm.params <- c("s", "p")

sm.data <- list(NsumCH = NsumCH, n_occasions = n.occasions, sumCH = sumCH,
                sumf = sumf, sumFR = sumFR)

# MCMC settings
# ni = 1000
# nt = 1
# nb = 500
# nc = 3

# Call Stan from R 
SM.t <- stan("Stan_Marginalized_Time.stan",
              data = sm.data,
              pars = sm.params,
              control = list(adapt_delta = .85),
              chains = nc, iter = ni, thin = nt) 
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Simulation of N-Augmentaion vs. Negative Binomail Abundance Estimates
#
#  Notes:
#  * 
#
###############################################################################
# simulation function
sim <- function(p = 0.1, Nind = 1000, Nvisits = 3, sd = 1){
  set.seed(sd)
  Y <- matrix(NA, nrow = Nind, ncol = Nvisits)
  for(i in 1:Nvisits){
    Y[,i] <- rbinom(Nind, 1, p)
    }
  Y[rowSums(Y) > 0,]
}

Nvisits <- 3
params <- c("N")
#-----------------------------------------------------------------------------#
# jags code - Naug
sink("N_JD.jags")
cat("
model{
  for(i in 1:M){
    Z[i] ~ dbern(omega)
    p_eff[i] <- Z[i] * p
    for(j in 1:Nvisits){
      y_aug[i,j] ~ dbern(p_eff[i])		
    }
  } 
  N <- sum(Z[])
  omega ~ dunif(0, 1)
  p ~ dunif(0, 1)
}
    ",fill = TRUE)
sink()

# JD inits
JD_inits <- function() list(Z = rep(1, nrow(y_aug)))

#-----------------------------------------------------------------------------#
# jags code - marginalized with Negative Binomial used to estimate N
sink("N_JM.jags")
cat("
model{
  for(i in 1:Nobs){
    pz[i] <- pow(p, y[i]) * pow(1 - p, 3 - y[i]) / pstar
    o[i] ~ dbern(pz[i])
  }
  U ~ dnegbin(pstar, Nobs)
  N <- Nobs + U
  p ~ dunif(0, 1)
  pstar <-1 - pow(1 - p, 3)
}
    ",fill = TRUE)
sink()
#-----------------------------------------------------------------------------#
# simulate and fit
OUT_s1 <- array(NA, dim = c(100, 2, 10))
for(j in 1:100){
  Y <- sim(sd = j)
  naug <- 2000
  y_aug <- rbind(Y, array(0, dim = c(naug, Nvisits)))
  win.data <- list(y_aug = y_aug, M = nrow(y_aug), Nvisits = Nvisits)
  win.data2 <-list(y = rowSums(Y), Nobs = nrow(Y), o = rep(1,nrow(Y)))
  t1 <- proc.time()
  sim3_JD <- jags.parallel(win.data, JD_inits, params, "N_JD.jags",
                           n.chains = 3, n.iter = 100000)
  t2<-proc.time()
  sim3_JM <- jags.parallel(win.data2, init=NULL, params, "N_JM.jags",
                           n.chains = 3, n.iter = 1000)
  t3 <- proc.time()
  OUT_s1[j,1,] <- c(sim3_JD$BUGSoutput$summary[1,], (t2 - t1)[3])
  OUT_s1[j,2,] <- c(sim3_JM$BUGSoutput$summary[1,], (t3 - t2)[3])
}

OUT_s2 <- array(NA, dim = c(100, 2, 10))
for(j in 46:100){
  Y <- sim(p = 0.5, sd = j)
  naug <- 1000
  y_aug <- rbind(Y, array(0, dim = c(naug, Nvisits)))
  win.data <- list(y_aug = y_aug, M = nrow(y_aug), Nvisits = Nvisits)
  win.data2 <- list(y = rowSums(Y), Nobs = nrow(Y), o = rep(1, nrow(Y)))
  t1 <- proc.time()
  sim3_JD <- jags.parallel(win.data, JD_inits, params, "N_JD.jags",
                           n.chains = 3, n.iter = 100000)
  t2 <- proc.time()
  sim3_JM <- jags.parallel(win.data2, init=NULL, params, "N_JM.jags",
                           n.chains = 3, n.iter = 1000)
  t3 <- proc.time()
  OUT_s2[j,1,] <- c(sim3_JD$BUGSoutput$summary[1,], (t2 - t1)[3])
  OUT_s2[j,2,] <- c(sim3_JM$BUGSoutput$summary[1,], (t3 - t2)[3])
}

#-----------------------------------------------------------------------------#

# To do:
# add mcmc settings that are described in paper
# add some general notes on the structure of things

###############################################################################
#                                                                     Spring 19
#  Functions to:
#  1). process input data for model fitting 
#  2). process/extract results
#
###############################################################################
#-----------------------------------------------------------------------------#
# Function to create the matrix of latent state z 
known.state.cjs <- function(ch){
  state <- ch
  for(i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,] == 1))
    n2 <- max(which(ch[i,] == 1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state == 0] <- NA
  return(state)
}

#-----------------------------
# Function to create a matrix of initial values for latent state z
cjs.init.z <- function(ch, f){
  for(i in 1:dim(ch)[1]){
    if (sum(ch[i,]) == 1) next
    n2 <- max(which(ch[i,] == 1))
    ch[i,f[i]:n2] <- NA
  }
  for(i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}

#-----------------------------------------------------------------------------#
# Function to collapse capture history matrix into unique histories and the
# frequency of these unique histories. Returns a list: first element is the
# collapsed capture histories, second element is the frequency.

collapse.ch <- function(ch){
  ch.char = apply(ch, 1, function(x) paste(x, collapse = ","))
  
  ch.sum.out = t(sapply(strsplit(names(table(ch.char)), split = ","), as.numeric))
  fr.out = as.numeric(as.vector(table(ch.char)))
  
  return(list(ch.sum.out, fr.out))
}

#-----------------------------------------------------------------------------#
# Function to generate inits for z in JS model with data augmentation. From Kery
# and Schaub 2012

js.multistate.init <- function(ch, nz){
  ch[ch == 2] <- NA
  state <- ch
  for (i in 1:nrow(ch)){
    n1 <- min(which(ch[i,] == 1))
    n2 <- max(which(ch[i,] == 1))
    state[i,n1:n2] <- 2
  }
  state[state == 0] <- NA
  get.first <- function(x) min(which(!is.na(x)))
  get.last <- function(x) max(which(!is.na(x)))   
  f <- apply(state, 1, get.first)
  l <- apply(state, 1, get.last)
  for (i in 1:nrow(ch)){
    state[i,1:f[i]] <- 1
    if(l[i] != ncol(ch)) state[i, (l[i] + 1):ncol(ch)] <- 3
    state[i, f[i]] <- 2
  }   
  state <- rbind(state, matrix(1, ncol = ncol(ch), nrow = nz))
  # change with the book code -- Need first col to be NA !
  state[,1] <- NA  
  return(state)
}

#-----------------------------------------------------------------------------#
# Function to return summary stats for fitted model objects. 
# Added 'ignore' argument, for parms that are set, but your still tracking, but
# want to exclude from the summary (see BNT model). Only for Stan and JAGS. 
# The function takes a list of model fits (WinBUGS, JAGS, Stan) for the first
# argument.  For the timing to work correctly, the time taken to fit a model
# must be appended to the model fit object, as an attribute, for instance:
# t1 <- proc.time()
# model.fit <- jags.parallel(...)
# t2 <- proc.time()
# attr(model.fit, 'time') <- (t2 - t1)[3]

run.times = function(fit.list, ignore = NULL){
  
  out = data.frame(fitter = NA,
                   iterations = NA,
                   min.n.eff = NA,
                   min.n.eff.coda = NA,
                   med.n.eff = NA,
                   med.n.eff.coda = NA,
                   run.time = NA,
                   model = NA,
                   r.hat.count = NA)
  
  for(i in seq_along(fit.list)){
    
    if(any(class(fit.list[[i]]) == "stanfit")){
      
      out[i,]$fitter = "Stan"
      
      # get the n.iter
      out[i,]$iterations = fit.list[[i]]@stan_args[[1]]$iter
      
      if(is.null(ignore) == FALSE){
        # get the n.eff
        tmp.n.eff = rstan::summary(fit.list[[i]])$summary[,"n_eff"]
        n.eff = tmp.n.eff[which(!names(tmp.n.eff) %in% ignore)]
        
        tmp.n.eff.coda = coda::effectiveSize(organize(fit.list[[i]], mcmc.out = TRUE))
        n.eff.coda = tmp.n.eff.coda[which(!names(tmp.n.eff.coda) %in% ignore)]
        
        # count of Rhat over 1.1 (this includes the like/deviance)
        tmp.tmp.r.hat = rstan::summary(fit.list[[i]])$summary[,"Rhat"]
        tmp.r.hat = tmp.tmp.r.hat[which(!names(tmp.tmp.r.hat) %in% ignore)]
      } else {
        # get the n.eff
        n.eff = rstan::summary(fit.list[[i]])$summary[,"n_eff"]
        n.eff.coda = coda::effectiveSize(organize(fit.list[[i]], mcmc.out = TRUE))
        
        # count of Rhat over 1.1 (this includes the like/deviance)
        tmp.r.hat = rstan::summary(fit.list[[i]])$summary[,"Rhat"]
      }
      
      # minus 1 to cut out the likelihood value (only n.eff for parms)
      out[i,]$min.n.eff = min(n.eff[1:(length(n.eff) - 1)])
      
      out[i,]$min.n.eff.coda = min(n.eff.coda[2:length(n.eff.coda)])  # deviance is first here
      
      # median n.eff
      out[i,]$med.n.eff = median(n.eff[1:(length(n.eff) - 1)])
      
      out[i,]$med.n.eff.coda = median(n.eff.coda[2:length(n.eff.coda)])  # deviance is first here
      
      # model name
      out[i,]$model = fit.list[[i]]@model_name
      
      out[i,]$r.hat.count = length(which(tmp.r.hat > 1.1))
      
    }
    
    # had to add 'any' b/c jags in parallel has two classes, one of which will be 'rjags'
    if(any(class(fit.list[[i]]) == "rjags")){  
      out[i,]$fitter = "JAGS"
      
      # get the n.iter
      out[i,]$iterations = fit.list[[i]]$n.iter
      
      if(is.null(ignore) == FALSE){
        # get the n.eff
        tmp.n.eff = fit.list[[i]]$BUGSoutput$summary[,9]
        # n.eff = tmp.n.eff[which(names(tmp.n.eff) != ignore)]
        n.eff = tmp.n.eff[which(!names(tmp.n.eff) %in% ignore)]
        
        tmp.n.eff.coda = coda::effectiveSize(organize(fit.list[[i]], mcmc.out = TRUE))
        # n.eff.coda = tmp.n.eff.coda[which(names(tmp.n.eff.coda) != ignore)]
        n.eff.coda = tmp.n.eff.coda[which(!names(tmp.n.eff.coda) %in% ignore)]
        
        # count of Rhat over 1.1  (this includes the like/deviance)
        tmp.tmp.r.hat = fit.list[[i]]$BUGSoutput$summary[,8]
        # tmp.r.hat = tmp.tmp.r.hat[which(names(tmp.tmp.r.hat) != ignore)]
        tmp.r.hat = tmp.tmp.r.hat[which(!names(tmp.tmp.r.hat) %in% ignore)]
      } else {
        # get the n.eff
        n.eff = fit.list[[i]]$BUGSoutput$summary[,9]
        n.eff.coda = coda::effectiveSize(organize(fit.list[[i]], mcmc.out = TRUE))
        
        # count of Rhat over 1.1  (this includes the like/deviance)
        tmp.r.hat = fit.list[[i]]$BUGSoutput$summary[,8]
      }
      
      # min n.eff - cut out the deviance value (only n.eff for parms)
      out[i,]$min.n.eff = min(n.eff[2:length(n.eff)])
      
      out[i,]$min.n.eff.coda = min(n.eff.coda[2:length(n.eff.coda)])
      
      # median n, eff - cut out the deviance value (only n.eff for parms)
      out[i,]$med.n.eff = median(n.eff[2:length(n.eff)])
      
      out[i,]$med.n.eff.coda = median(n.eff.coda[2:length(n.eff.coda)])
      
      out[i,]$r.hat.count = length(which(tmp.r.hat > 1.1))
      
      # model name
      out[i,]$model = fit.list[[i]]$model.file
    }
    
    if(any(class(fit.list[[i]]) == "bugs")){
      out[i,]$fitter = "bugs"
      
      # get the n.iter
      out[i,]$iterations = fit.list[[i]]$n.iter
      
      # get the n.eff
      n.eff = fit.list[[i]]$summary[,9]
      
      f1 = coda::mcmc.list(lapply(1:fit.list[[i]]$n.chain, function(x) coda::mcmc(fit.list[[i]]$sims.array[,x,])))
      
      n.eff.coda = coda::effectiveSize(f1)
      
      # min n.eff - cut out the deviance value (only n.eff for parms), different than JAGS, deviance is at the end
      out[i,]$min.n.eff = min(n.eff[1:length(n.eff)-1])
      
      out[i,]$min.n.eff.coda = min(n.eff.coda[1:length(n.eff.coda)-1])
      
      # median n.eff
      out[i,]$med.n.eff = median(n.eff[1:length(n.eff)-1])
      
      out[i,]$med.n.eff.coda = median(n.eff.coda[1:length(n.eff.coda)-1])
      
      # model name
      out[i,]$model = fit.list[[i]]$model.file
      
      # count of Rhat over 1.1  (this includes the like/deviance)
      tmp.r.hat = fit.list[[i]]$summary[,8]
      
      out[i,]$r.hat.count = length(which(tmp.r.hat > 1.1))
      
    }
    
    # get time to fit model, stored as an attribute
    out[i,]$run.time = ifelse(is.null(attr(fit.list[[i]], "time")), "NA", attr(fit.list[[i]], "time")) 
  }
  
  # out$efficiency = out$min.n.eff / out$run.time  
  out$efficiency = out$min.n.eff.coda / out$run.time  
  # time required for min n.eff (coda) of 100
  out$time.to.100 = (out$run.time / out$min.n.eff.coda) * 100
  
  return(out)
}

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
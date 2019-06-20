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
#  Application 1
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
#  Application 1
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
#  Application 1
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
#  Application 1
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
#  Application 1
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
#  Application 1
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
#  Application 1
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
#  Application 2
#  Multi-state Capture-Recapture model (Humpback Chub) 
#  JAGS Discrete Version
#
#  Notes:
#  * Need to set directory for data
#  * Note that CH's are formatted so that each row only has the release and
#    subsequent capture information
#
#	   Thus, a CH that is normally : 1 0 1 1 0 0  
#
#	   Would be formatted as such: 1 0 1 0 0 0 (fc = 1, lc =3)
#			       	      0 0 1 1 0 0 (fc = 3, lc = 4)
#			       	      0 0 0 1 0 0 (fc = 4, lc = 6)
#
#		 where fc = release occasion, lc = capture occasion
#
#  * States are defined as such:
#		 1 = small adult captured in LCR
#		 2 = large adult captured in LCR
#		 3 = small adult captured in CR
#		 4 = large adult captured in CR
#		 5 = not observed
#
#	 * Fish are allowed to move between CR and LCR 
#	 * Fish are allowed to grow from small to large adults
#	 * Fish are not allowed to shrink (from large to small adults)
#
###############################################################################
library(R2jags)

# Get data:
dat <- read.csv(".//HBC_data.csv")
#-----------------------------------------------------------------------------#
# Functions:

#	Get first capture occasion:
find.first <- function(x){min(which(x != 5))}

#	Get last capture occasion:
find.last <- function(x){
  ifelse(length(which(x !=5 )) == 1, 27, max(which(x != 5)))
}

#-----------------------------------------------------------------------------#
# format data for model fitting:
sumCH <- dat[,1:27]
sumFR <- dat[,28]
newtag.AO <- dat[,29]

# other data:
season <- c(rep(1:3,v8),v1:2) # season index
LCRs <- c(2, 3, 5, 6, 8, 9, 11:26)
LCRns <- c(1, 4, 7, 10)
CRs <- c(1, 2, 4, 5, 7:26)
CRns <- c(3, 6)
catch <- matrix(c(287, 545, 221, 594, 215, 413, 500, 374, 562,
                  193, 356, 171, 191, 111, 239, 259, 283, 246,
                  102, 129, 154, 54, 63, 48, 36, 41, 35,
                  35, 86, 126, 41, 43, 59, 22, 26, 34), nrow = 4, byrow = TRUE)


# Re-format CH so that lines correspond to individuals (not unique capture histories):
indCH <- numeric()
for(i in 1:length(sumCH[,1])){
  for(j in 1:sumFR[i]){
    indCH <- rbind(indCH, sumCH[i,])
  }
}

indCH <- as.matrix(indCH)
indf <- apply(indCH, 1, find.first)
indl <- apply(indCH, 1, find.last)

# Create CH data matrix where state is known when fish is captured, otherwise NA:
known.state <- indCH
known.state[indCH == 5] <- NA
known.state[,1] <- NA
for(i in 1:dim(indCH)[1]){known.state[i,indf[i]] <- NA}

# Create initial values for z-states (indicates whether fish is a small adult in
# LCR (1), a large adult in LCR (2), a small adult in CR (3), or a large adult
# in CR (4)
# Note code does not allow fish to shrink, so initial values for z-states must
# also not allow for shrinking
z.init <- matrix(NA, nrow = dim(indCH)[1], ncol = 27)
for(j in 1:dim(indCH)[1]){
  z.init[j,(indf[j] + 1)] <- ifelse(indCH[j,(indf[j] + 1)] < 5, NA, indCH[j,indf[j]])
  if(indf[j] < 26 & (indl[j] - indf[j]) > 1){
    for(i in (indf[j]+2):indl[j]){
      size <- c(1, 2, 1, 2, NA)
      sizemax <- max(size[indCH[j,indf[j]:i]], na.rm = TRUE)
      if(sizemax == 1){x <- c(1, 3)}
      if(sizemax == 2){x <- c(2, 4)}
      z.init[j,i] <- ifelse(indCH[j,i] < 5, NA, sample(x, 1))
    }
  }
}

#-----------------------------------------------------------------------------#
sink("JAGS_Discrete_chub.txt")
cat("
    model{
    
    for(j in 1:4){
    s_i[j] ~ dunif(0,1) # 3 month survival rates for small and large adults in LCR and CR
    s[1,j] <- s_i[j]    # 3 month survival for spring to summer interval
    s[2,j] <- s_i[j]    # 3 month survival for summer to fall interval
    s[3,j] <- s_i[j] * s_i[j] # 6 month survival for fall to spring interval
    } 
    
    tau ~ dunif(0,1) # proportion of adults in the observable portion of the LCR
    
    for (i in 1:3){
    g[i,1] ~ dunif(0,1) # LCR growth for three intervals
    g[i,2] ~ dunif(0,1) # CR growth for three intervals
    m[i,1] ~ dunif(0,1) # movement out of LCR by interval for small adults
    m[i,2] ~ dunif(0,1) # movement out of LCR by interval for large adults
    m[i,3] ~ dunif(0,1) # movement into the LCR by interval for small adults
    m[i,4] ~ dunif(0,1) # movement into the LCR by interval for large adults
    
    # tr is an array representing markov transitions for three intervals. some
    # transitions are assumed to be zero (e.g., transitioning from a large to a
    # small adult)
    tr[1,i,1] <- s[i,1] * (1 - g[i,1]) * (1 - m[i,1])
    tr[1,i,2] <- s[i,1] * g[i,1] * (1 - m[i,2])
    tr[1,i,3] <- s[i,1] * (1 - g[i,1]) * m[i,1] * tau
    tr[1,i,4] <- s[i,1] * g[i,1] * m[i,2] * tau
    tr[1,i,5] <- s[i,1] * (1 - g[i,1]) * m[i,1] * (1 - tau)
    tr[1,i,6] <- s[i,1] * g[i,1] * m[i,2] * (1 - tau)
    tr[1,i,7] <- 1 - s[i,1]
    tr[2,i,1] <- 0
    tr[2,i,2] <- s[i,2] * (1 - m[i,2])
    tr[2,i,3] <- 0
    tr[2,i,4] <- s[i,2] * m[i,2] * tau
    tr[2,i,5] <- 0
    tr[2,i,6] <- s[i,2] * m[i,2] * (1 - tau)
    tr[2,i,7] <- 1 - s[i,2]
    tr[3,i,1] <- s[i,3] * (1 - g[i,2]) * m[i,3]
    tr[3,i,2] <- s[i,3] * g[i,2] * m[i,4]
    tr[3,i,3] <- s[i,3] * (1 - g[i,2]) * (1 - m[i,3])
    tr[3,i,4] <- s[i,3] * g[i,2] * (1 - m[i,4])
    tr[3,i,5] <- 0
    tr[3,i,6] <- 0
    tr[3,i,7] <- 1 - s[i,3]
    tr[4,i,1] <- 0
    tr[4,i,2] <- s[i,4] * m[i,4]
    tr[4,i,3] <- 0
    tr[4,i,4] <- s[i,4] * (1 - m[i,4])
    tr[4,i,5] <- 0
    tr[4,i,6] <- 0
    tr[4,i,7] <- 1 - s[i,4]
    tr[5,i,1] <- tr[3,i,1]
    tr[5,i,2] <- tr[3,i,2]
    tr[5,i,3] <- 0
    tr[5,i,4] <- 0
    tr[5,i,5] <- tr[3,i,3]
    tr[5,i,6] <- tr[3,i,4]
    tr[5,i,7] <- tr[3,i,7]
    tr[6,i,1] <- 0
    tr[6,i,2] <- tr[4,i,2]
    tr[6,i,3] <- 0
    tr[6,i,4] <- 0
    tr[6,i,5] <- 0
    tr[6,i,6] <- tr[4,i,4]
    tr[6,i,7] <- tr[4,i,7]
    for(j in 1:6){
    tr[7,i,j] <- 0
    }
    tr[7,i,7] <- 1
    }	
    
    # p is the matrix describing capture probabilities
    # Diagonal values represent capture probabilities for the four states
    # The fifth column represents the probability a fish is not captured
    # Note that for fish in states 5-6 (unobserveable sites in the CR that are not sampled)
    # and for state 7 (dead state) the probability a fish is unobserved is 100%
    for(i in 1:22){
    p[1,LCRs[i],1] ~ dunif(0, 1)
    p[2,LCRs[i],2] ~ dunif(0, 1)
    }
    
    # No sampling in LCR during 1st, 4th, 7th, 10th recap periods (summers of 2009 - 2012)
    for(i in 1:4){
    p[1,LCRns[i],1] <- 0
    p[2,LCRns[i],2] <- 0
    }
    for(i in 1:24){
    p[3,CRs[i],3] ~ dunif(0, 1)
    p[4,CRs[i],4] ~ dunif(0, 1)
    }
    
    # No sampling in CR during 3rd and 6th recap periods (spring of 2010 & 2011)
    for(i in 1:2){
    p[3,CRns[i],3] <- 0
    p[4,CRns[i],4] <- 0
    }
    
    for(i in 1:26){
    p[1,i,2] <- 0
    p[1,i,3] <- 0
    p[1,i,4] <- 0
    p[1,i,5] <- 1 - p[1,i,1]
    p[2,i,1] <- 0
    p[2,i,3] <- 0
    p[2,i,4] <- 0
    p[2,i,5] <- 1 - p[2,i,2]
    p[3,i,1] <- 0
    p[3,i,2] <- 0
    p[3,i,4] <- 0
    p[3,i,5] <- 1 - p[3,i,3]
    p[4,i,1] <- 0
    p[4,i,2] <- 0
    p[4,i,3] <- 0
    p[4,i,5] <- 1 - p[4,i,4]
    p[5,i,1] <- 0
    p[5,i,2] <- 0
    p[5,i,3] <- 0
    p[5,i,4] <- 0
    p[5,i,5] <- 1
    p[6,i,1] <- 0
    p[6,i,2] <- 0
    p[6,i,3] <- 0
    p[6,i,4] <- 0
    p[6,i,5] <- 1
    p[7,i,1] <- 0
    p[7,i,2] <- 0
    p[7,i,3] <- 0
    p[7,i,4] <- 0
    p[7,i,5] <- 1
    }
    
    for(k in 1:NindCH){
    z[k,indf[k]] = indCH[k,indf[k]]
    for(i in (indf[k] + 1):indl[k]){
    z[k,i] ~ dcat(tr[z[k,i - 1], season[i - 1],])
    indCH[k,i] ~ dcat(p[z[k,i], i - 1,])
    }
    }
    }
    
    ",fill=TRUE)
sink() 

#-----------------------------------------------------------------------------#
JD.data <- list(LCRs = LCRs, LCRns = LCRns, CRs = CRs, CRns = CRns,
                season = season, NindCH = dim(indCH)[1], indCH = indCH,
                indf = indf, indl = indl, z = known.state)

JD.par <- c('s_i', 'g', 'm', 'tau')

JD.inits <- function(){
  list(s_i = c(.7, .7, .7, .7), z = array(z.init, dim = c(dim(indCH)[1], 27)))
}

JD.out <- jags.parallel(JD.data, inits = JD.inits, JD.par, ".\\JAGS_Discrete_chub.txt",
                        n.chains = 1, n.iter = 10, export_obj_names = "z.init")

#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Application 2
#  Multi-state Capture-Recapture model (Humpback Chub)
#  JAGS Marginalized Version
#
#  Notes:
#  * Need to set directory for data
#  * See Notes above
###############################################################################
library(R2jags)

# Get data:
dat <- read.csv(".//HBC_data.csv")
#-----------------------------------------------------------------------------#
#  Functions:
#	Get first capture occasion:
find.first <- function(x){min(which(x != 5))}

#	Get last capture occasion:
find.last <- function(x){
  ifelse(length(which(x != 5)) == 1, 27, max(which(x != 5)))
}

#-----------------------------------------------------------------------------#
# format data for model fitting:
sumCH <- as.matrix(dat[,1:27])
sumFR <- dat[,28]
newtag.AO <- dat[,29]

# other data:
season <- c(rep(1:3, 8), 1:2) # season index
LCRs <- c(2, 3, 5, 6, 8, 9, 11:26)
LCRns <- c(1, 4, 7, 10)
CRs <- c(1, 2, 4, 5, 7:26)
CRns <- c(3, 6)
catch <- matrix(c(287, 545, 221, 594, 215, 413, 500, 374, 562,
                  193, 356, 171, 191, 111, 239, 259, 283, 246,
                  102, 129, 154, 54, 63, 48, 36, 41, 35,
                  35, 86, 126, 41, 43, 59, 22, 26, 34), nrow = 4, byrow = TRUE)

# Get first and last capture occasion for line in summarized capture history
fc <- apply(sumCH, 1, find.first)
lc <- apply(sumCH, 1, find.last)

# Use this to calculate abundance in code:
CR_ind <- matrix(0, nrow = 4, ncol = 9)
1 -> CR_ind[3:4,]

#-----------------------------------------------------------------------------#
sink("JAGS_Marginalized_chub.txt")
cat("
    model{
    
    for(j in 1:4){
    s_i[j] ~ dunif(0, 1) # 3 month survival rates for small and large adults in LCR and CR
    s[1,j] <- s_i[j]     # 3 month survival for spring to summer interval
    s[2,j] <- s_i[j]     # 3 month survival for summer to fall interval
    s[3,j] <- s_i[j] * s_i[j] # 6 month survival for fall to spring interval
    } 
    
    tau ~ dunif(0, 1) # proportion of adults in the observable portion of the LCR
    for(i in 1:3){
    g[i,1] ~ dunif(0, 1) # LCR growth for three intervals
    g[i,2] ~ dunif(0, 1) # CR growth for three intervals
    m[i,1] ~ dunif(0, 1) # movement out of LCR by interval for small adults
    m[i,2] ~ dunif(0, 1) # movement out of LCR by interval for large adults
    m[i,3] ~ dunif(0, 1) # movement into the LCR by interval for small adults
    m[i,4] ~ dunif(0, 1) # movement into the LCR by interval for large adults
    
    # tr is an array representing markov transitions for three intervals. some
    # transitions are assumed to be zero (e.g., transitioning from a large to a
    # small adult)
    tr[1,i,1] <- s[i,1] * (1 - g[i,1]) * (1 - m[i,1])
    tr[1,i,2] <- s[i,1] * g[i,1] * (1 - m[i,2])
    tr[1,i,3] <- s[i,1] * (1 - g[i,1]) * m[i,1] * tau
    tr[1,i,4] <- s[i,1] * g[i,1] * m[i,2] * tau
    tr[1,i,5] <- s[i,1] * (1 - g[i,1]) * m[i,1] * (1 - tau)
    tr[1,i,6] <- s[i,1] * g[i,1] * m[i,2] * (1 - tau)
    tr[1,i,7] <- 1 - s[i,1]
    tr[2,i,1] <- 0
    tr[2,i,2] <- s[i,2] * (1 - m[i,2])
    tr[2,i,3] <- 0
    tr[2,i,4] <- s[i,2] * m[i,2] * tau
    tr[2,i,5] <- 0
    tr[2,i,6] <- s[i,2] * m[i,2] * (1 - tau)
    tr[2,i,7] <- 1 - s[i,2]
    tr[3,i,1] <- s[i,3] * (1 - g[i,2]) * m[i,3]
    tr[3,i,2] <- s[i,3] * g[i,2] * m[i,4]
    tr[3,i,3] <- s[i,3] * (1 - g[i,2]) * (1 - m[i,3])
    tr[3,i,4] <- s[i,3] * g[i,2] * (1 - m[i,4])
    tr[3,i,5] <- 0
    tr[3,i,6] <- 0
    tr[3,i,7] <- 1 - s[i,3]
    tr[4,i,1] <- 0
    tr[4,i,2] <- s[i,4] * m[i,4]
    tr[4,i,3] <- 0
    tr[4,i,4] <- s[i,4] * (1 - m[i,4])
    tr[4,i,5] <- 0
    tr[4,i,6] <- 0
    tr[4,i,7] <- 1 - s[i,4]
    tr[5,i,1] <- tr[3,i,1]
    tr[5,i,2] <- tr[3,i,2]
    tr[5,i,3] <- 0
    tr[5,i,4] <- 0
    tr[5,i,5] <- tr[3,i,3]
    tr[5,i,6] <- tr[3,i,4]
    tr[5,i,7] <- tr[3,i,7]
    tr[6,i,1] <- 0
    tr[6,i,2] <- tr[4,i,2]
    tr[6,i,3] <- 0
    tr[6,i,4] <- 0
    tr[6,i,5] <- 0
    tr[6,i,6] <- tr[4,i,4]
    tr[6,i,7] <- tr[4,i,7]
    for (j in 1:6){
    tr[7,i,j] <- 0
    }
    tr[7,i,7] <- 1
    }	
    
    # p is the matrix describing capture probabilities
    # Diagonal values represent capture probabilities for the four states
    # The fifth column represents the probability a fish is not captured
    # Note that for fish in states 5-6 (unobserveable sites in the CR that are not sampled)
    # and for state 7 (dead state) the probability a fish is unobserved is 100%
    for(i in 1:22){
    p[1,LCRs[i],1] ~ dunif(0, 1)
    p[2,LCRs[i],2] ~ dunif(0, 1)
    }
    
    # No sampling in LCR during 1st, 4th, 7th, 10th recap periods (summers of 2009 - 2012)
    for(i in 1:4){
    p[1,LCRns[i],1] <- 0
    p[2,LCRns[i],2] <- 0
    }
    for(i in 1:24){
    p[3,CRs[i],3] ~ dunif(0, 1)
    p[4,CRs[i],4] ~ dunif(0, 1)
    }
    
    # No sampling in CR during 3rd and 6th recap periods (spring of 2010 & 2011)
    for(i in 1:2){
    p[3,CRns[i],3] <- 0
    p[4,CRns[i],4] <- 0
    }
    
    for(i in 1:26){
    p[1,i,2] <- 0
    p[1,i,3] <- 0
    p[1,i,4] <- 0
    p[1,i,5] <- 1 - p[1,i,1]
    p[2,i,1] <- 0
    p[2,i,3] <- 0
    p[2,i,4] <- 0
    p[2,i,5] <- 1 - p[2,i,2]
    p[3,i,1] <- 0
    p[3,i,2] <- 0
    p[3,i,4] <- 0
    p[3,i,5] <- 1 - p[3,i,3]
    p[4,i,1] <- 0
    p[4,i,2] <- 0
    p[4,i,3] <- 0
    p[4,i,5] <- 1 - p[4,i,4]
    p[5,i,1] <- 0
    p[5,i,2] <- 0
    p[5,i,3] <- 0
    p[5,i,4] <- 0
    p[5,i,5] <- 1
    p[6,i,1] <- 0
    p[6,i,2] <- 0
    p[6,i,3] <- 0
    p[6,i,4] <- 0
    p[6,i,5] <- 1
    p[7,i,1] <- 0
    p[7,i,2] <- 0
    p[7,i,3] <- 0
    p[7,i,4] <- 0
    p[7,i,5] <- 1
    }
    
    # The 0.03 is to account for one-time 3% tag loss
    for(k in 1:NsumCH){
    pz[k,sumf[k],1] <- equals(sumCH[k,sumf[k]], 1) * (1 - 0.03 * newtag[k]) 
    pz[k,sumf[k],2] <- equals(sumCH[k,sumf[k]], 2) * (1 - 0.03 * newtag[k])
    pz[k,sumf[k],3] <- equals(sumCH[k,sumf[k]], 3) * (1 - 0.03 * newtag[k])
    pz[k,sumf[k],4] <- equals(sumCH[k,sumf[k]], 4) * (1 - 0.03 * newtag[k])
    pz[k,sumf[k],5] <- 0
    pz[k,sumf[k],6] <- 0
    pz[k,sumf[k],7] <- 0.03 * newtag[k]
    for(i in sumf[k]:(lc[k] - 1)){
    for(j in 1:7){
    pz[k,(i + 1),j] <- inprod(pz[k,i,], tr[,season[i],j]) * p[j,i,sumCH[k,(i + 1)]]
    }
    }
    
    lik[k] <- sum(pz[k,lc[k],])
    one[k] ~ dbin(lik[k], sumFR[k])
    }
    
    # Calculate abundance from catch by simulating from a negative binomial distribution:
    for(i in 1:4){
    for(j in 1:9){
    p_fall[i,j] <- p[i,fall_ind[j],i] * (1 - CR_ind[i,j] * (1 - tau))
    U[i,j] ~ dnegbin(p_fall[i,j], catch[i,j])
    N[i,j] = catch[i,j] + U[i,j] 
    }
    }
    }
    
    ",fill=TRUE)
sink() 

#-----------------------------------------------------------------------------#
JM.data <- list(LCRs = LCRs, LCRns = LCRns, CRs = CRs, CRns = CRns,
                season = season, NsumCH = dim(sumCH)[1],
                sumCH = array(sumCH, dim = c(dim(sumCH)[1], dim(sumCH)[2])),
                newtag = as.vector(newtag.AO), sumf = as.vector(fc),
                sumFR = sumFR, one = sumFR,  
                CR_ind = array(CR_ind, dim = c(dim(CR_ind)[1], dim(CR_ind)[2])),
                lc = as.vector(lc), fall_ind = 1:9 * 3 - 1,
                catch = array(catch, dim = c(dim(catch)[1], dim(catch)[2])))

JM.par <- c('s_i', 'g', 'm', 'tau', 'p_cr', 'p_lcr', 'N')

JM.inits <- function(){list(s_i = c(.7, .7, .7, .7))}

JM.out <- jags.parallel(JM.data, inits = JM.inits, JM.par,
                        ".\\JAGS_Marginalized_chub.txt", n.cluster = 3,
                        n.chains = 1, n.iter = 10)
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Application 2
#  Multi-state Capture-Recapture model (Humpback Chub)
#  Stan Marginalized Version
#
#  Notes:
#  * Need to set directory for data
#  * See Notes above
###############################################################################
library(rstan)

# To run Stan in parallel:
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Get data:
dat <- read.csv(".//HBC_data.csv")
#-----------------------------------------------------------------------------#
# Functions:
#	Get first capture occasion:
find.first <- function(x){min(which(x != 5))}

#	Get last capture occasion:
find.last <- function(x){
  ifelse(length(which(x != 5)) == 1, 27, max(which(x != 5)))
}

#-----------------------------------------------------------------------------#
# format data for model fitting:
sumCH <- as.matrix(dat[,1:27])
sumFR <- dat[,28]
newtag.AO <- dat[,29]

# other data:
season <- c(rep(1:3, 8), 1:2) # season index
LCRs <- c(2, 3, 5, 6, 8, 9, 11:26)
LCRns <- c(1, 4, 7, 10)
CRs <- c(1, 2, 4, 5, 7:26)
CRns <- c(3, 6)
catch <- matrix(c(287, 545, 221, 594, 215, 413, 500, 374, 562,
                  193, 356, 171, 191, 111, 239, 259, 283, 246,
                  102, 129, 154, 54, 63, 48, 36, 41, 35,
                  35, 86, 126, 41, 43, 59, 22, 26, 34), nrow = 4, byrow = TRUE)

# Get first and last capture occasion for line in summarized capture history
fc <- apply(sumCH, 1, find.first)
lc <- apply(sumCH, 1, find.last)

# Use this to calculate abundance in code:
CR_ind <- matrix(0, nrow = 4, ncol = 9)
1 -> CR_ind[3:4,]
#-----------------------------------------------------------------------------#
sink("Stan_Marginalized_chub.stan")
cat("
    data{
    int<lower = 1> NsumCH;
    int<lower = 1, upper = 5> sumCH[NsumCH, 27];
    int<lower = 1, upper = 26> sumf[NsumCH];
    int<lower = 1, upper = 3> season[26];
    int<lower = 1> sumFR[NsumCH];
    int<lower = 1> LCRs [22];
    int<lower = 1> LCRns [4];
    int<lower = 1> CRs [24];
    int<lower = 1> CRns [2];
    int catch_mat[4,9];
    int fall_ind[9];
    int<lower = 0, upper = 1> newtag[NsumCH];
    int<lower = 2, upper = 27> lc[NsumCH];
    }
    
    parameters{
    real<lower = 0, upper = 1> s_i[4]; // 3 month survivals for small and large chub in LCR and CR
    real<lower = 0, upper = 1> tau;    // proportion of CR adults residing in observable location
    vector<lower = 0, upper = 1> [2] g [3]; // seasonal growth for two rivers
    vector<lower = 0, upper = 1> [4] m [3]; // seasonal movement rates for 2 sizes and locations
    vector<lower = 0, upper = 1> [2] p_lcr [22]; // recapture probability in LCR for two size classes
    vector<lower = 0, upper = 1> [2] p_cr [24];  // recapture probability in CR for two size classes
    }
    
    transformed parameters{
    vector<lower = 0, upper = 1> [4] s[3];
    simplex[7] tr[7,3];
    simplex[5] p[7,26];
    
    for(j in 1:4){
    s[1,j] = s_i[j]; // 3 month survival for spring to summer interval
    s[2,j] = s_i[j]; // 3 month survival for summer to fall interval
    s[3,j] = s_i[j] * s_i[j]; // 6 month survival for fall to spring interval
    } 
    
    for(i in 1:3){
    tr[1,i,1] = s[i,1] * (1 - g[i,1]) * (1 - m[i,1]);
    tr[1,i,2] = s[i,1] * g[i,1] * (1 - m[i,2]);
    tr[1,i,3] = s[i,1] * (1 - g[i,1]) * m[i,1] * tau;
    tr[1,i,4] = s[i,1] * g[i,1] * m[i,2] * tau;
    tr[1,i,5] = s[i,1] * (1 - g[i,1]) * m[i,1] * (1 - tau);
    tr[1,i,6] = s[i,1] * g[i,1] * m[i,2] * (1 - tau);
    tr[1,i,7] = 1 - s[i,1];
    tr[2,i,1] = 0;
    tr[2,i,2] = s[i,2] * (1 - m[i,2]);
    tr[2,i,3] = 0;
    tr[2,i,4] = s[i,2] * m[i,2] * tau;
    tr[2,i,5] = 0;
    tr[2,i,6] = s[i,2] * m[i,2] * (1-tau);
    tr[2,i,7] = 1 - s[i,2];
    tr[3,i,1] = s[i,3] * (1 - g[i,2]) * m[i,3];
    tr[3,i,2] = s[i,3] * g[i,2] * m[i,4];
    tr[3,i,3] = s[i,3] * (1 - g[i,2]) * (1 - m[i,3]);
    tr[3,i,4] = s[i,3] * g[i,2] * (1 - m[i,4]);
    tr[3,i,5] = 0;
    tr[3,i,6] = 0;
    tr[3,i,7] = 1 - s[i,3];
    tr[4,i,1] = 0;
    tr[4,i,2] = s[i,4] * m[i,4];
    tr[4,i,3] = 0;
    tr[4,i,4] = s[i,4] * (1 - m[i,4]);
    tr[4,i,5] = 0;
    tr[4,i,6] = 0;
    tr[4,i,7] = 1 - s[i,4];
    tr[5,i,1] = tr[3,i,1];
    tr[5,i,2] = tr[3,i,2];
    tr[5,i,3] = 0;
    tr[5,i,4] = 0;
    tr[5,i,5] = tr[3,i,3];
    tr[5,i,6] = tr[3,i,4];
    tr[5,i,7] = tr[3,i,7];
    tr[6,i,1] = 0;
    tr[6,i,2] = tr[4,i,2];
    tr[6,i,3] = 0;
    tr[6,i,4] = 0;
    tr[6,i,5] = 0;
    tr[6,i,6] = tr[4,i,4];
    tr[6,i,7] = tr[4,i,7];
    for(j in 1:6){
    tr[7,i,j] = 0;
    }
    tr[7,i,7] = 1;
    }
    
    for(i in 1:22){
    p[1,LCRs[i],1] = p_lcr[i,1];
    p[2,LCRs[i],2] = p_lcr[i,2];
    }
    // No sampling in LCR during 1st, 4th, 7th, 10th recap periods (summers of 2009 - 2012)
    for(i in 1:4){
    p[1,LCRns[i],1] = 0;
    p[2,LCRns[i],2] = 0;
    }
    for(i in 1:24){
    p[3,CRs[i],3] = p_cr[i,1];
    p[4,CRs[i],4] = p_cr[i,2];
    }
    // No sampling in CR during 3rd and 6th recap periods (spring of 2010 & 2011)
    for(i in 1:2){
    p[3,CRns[i],3] = 0;
    p[4,CRns[i],4] = 0;
    }
    
    for(i in 1:26){
    p[1,i,2] = 0;
    p[1,i,3] = 0;
    p[1,i,4] = 0;
    p[1,i,5] = 1 - p[1,i,1];
    p[2,i,1] = 0;
    p[2,i,3] = 0;
    p[2,i,4] = 0;
    p[2,i,5] = 1 - p[2,i,2];
    p[3,i,1] = 0;
    p[3,i,2] = 0;
    p[3,i,4] = 0;
    p[3,i,5] = 1 - p[3,i,3];
    p[4,i,1] = 0;
    p[4,i,2] = 0;
    p[4,i,3] = 0;
    p[4,i,5] = 1 - p[4,i,4];
    p[5,i,1] = 0;
    p[5,i,2] = 0;
    p[5,i,3] = 0;
    p[5,i,4] = 0;
    p[5,i,5] = 1;
    p[6,i,1] = 0;
    p[6,i,2] = 0;
    p[6,i,3] = 0;
    p[6,i,4] = 0;
    p[6,i,5] = 1;
    p[7,i,1] = 0;
    p[7,i,2] = 0;
    p[7,i,3] = 0;
    p[7,i,4] = 0;
    p[7,i,5] = 1;
    }
    }
    
    model{
    real temp[7]; 
    vector[7] pz[27]; 
    for(k in 1:NsumCH){ 
    for(j in 1:6){
    pz[sumf[k], j] = (1 - 0.03 * newtag[k]) * (j == sumCH[k, sumf[k]]); 
    } 
    pz[sumf[k], 7] = (0.03 * newtag[k]);  
    for(t in (sumf[k] + 1):lc[k]){ 
    for(i in 1:7){ 
    for(j in 1:7){
    temp[j] = pz[t - 1, j] * tr[j, season[(t - 1)], i] * p[i, t - 1, sumCH[k, t]]; 
    }
    pz[t,i] = sum(temp); 
    } 
    }
    target += sumFR[k] * log(sum(pz[lc[k]])); 
    }
    }
    
    // Code to estimate abundance by simulating from a negative binomial distribution
    generated quantities{
    real ptrans;	
    real<lower = 0> scale_par;
    int U[4,9];
    int N[4,9];
    
    for(i in 1:4){
    for(j in 1:9){
    ptrans = p[i,fall_ind[j],i] * (1 - (i > 2) * (1 - tau));
    scale_par = ptrans / (1 - ptrans);
    U[i,j] = neg_binomial_rng(catch_mat[i,j], scale_par);
    N[i,j] = U[i,j] + catch_mat[i,j];
    }
    }
    }
    
    ",fill=TRUE)
sink() 
#-----------------------------------------------------------------------------#
sm.data<-list(NsumCH = dim(sumCH)[1], newtag = as.vector(newtag.AO),
              sumCH = array(sumCH,  dim = c(dim(sumCH)[1], dim(sumCH)[2])), 
              sumf = fc, season = season, sumFR = sumFR, LCRs = LCRs,
              LCRns = LCRns, CRs = CRs, CRns = CRns,  
              catch_mat = array(catch,  dim = c(dim(catch)[1], dim(catch)[2])),
              lc = lc,  fall_ind = 1:9 * 3 - 1)

sm.inits <- function() list(s_i = runif(4, 0, 1))  
sm.params <- c("s_i", "g", "m", "tau", "p_cr", "p_lcr", "N") 

SM <- stan(".\\Stan_Marginalized_chub.stan", 
           data = sm.data, init = sm.inits, pars = sm.params, 
           chains = 1, iter = 10, thin = 1, 
           seed = 1) 
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Application 2
#  Multi-state Capture-Recapture model (Humpback Chub)
#  Stan Marginalized Version with random effects
#
#  Notes:
#  * Need to set directory for data
#  * See Notes above
###############################################################################
library(rstan)

# To run Stan in parallel:
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Get data:
dat <- read.csv(".//HBC_data.csv")
#-----------------------------------------------------------------------------#
#  Functions:
#	Get first capture occasion:
find.first <- function(x){min(which(x != 5))}

#	Get last capture occasion:
find.last <- function(x){
  ifelse(length(which(x != 5)) == 1, 27, max(which(x != 5)))
}
#-----------------------------------------------------------------------------#
# format data for model fitting:
sumCH <- as.matrix(dat[,1:27])
sumFR <- dat[,28]
newtag.AO <- dat[,29]

# other data:
season <- c(rep(1:3, 8), 1:2) # season index
LCRs <- c(2, 3, 5, 6, 8, 9, 11:26)
LCRns <- c(1, 4, 7, 10)
CRs <- c(1, 2, 4, 5, 7:26)
CRns <- c(3, 6)
catch <- matrix(c(287, 545, 221, 594, 215, 413, 500, 374, 562,
                  193, 356, 171, 191, 111, 239, 259, 283, 246,
                  102, 129, 154, 54, 63, 48, 36, 41, 35,
                  35, 86, 126, 41, 43, 59, 22, 26, 34), nrow = 4, byrow = TRUE)

# Get first and last capture occasion for line in summarized capture history
fc <- apply(sumCH, 1, find.first)
lc <- apply(sumCH, 1, find.last)

# Use this to calculate abundance in code:
CR_ind <- matrix(0, nrow = 4, ncol = 9)
1 -> CR_ind[3:4,]

#-----------------------------------------------------------------------------#
sink("Stan_Marginalized_chub_RE.stan")
cat("
    data{
    int<lower = 1> NsumCH;
    int<lower = 1, upper = 5> sumCH[NsumCH,  27];
    int<lower = 1, upper = 26> sumf[NsumCH];
    int<lower = 1, upper = 3> season[26];
    int<lower = 1> sumFR[NsumCH];
    int<lower = 1> LCRs [22];
    int<lower = 1> LCRns [4];
    int<lower = 1> CRs [24];
    int<lower = 1> CRns [2];
    int catch_mat[4, 9];
    int fall_ind[9];
    int<lower = 2,  upper = 27> lc[NsumCH];
    int<lower = 0,  upper = 1> newtag[NsumCH];
    }
    
    parameters {
    real<lower = -8, upper = 8> mu_ls[4]; // hyperprior on 3 month survivals for small and large chub in LCR and CR
    real<lower = 0> sd_ls; // hyperprior on 3 month survivals for small and large chub in LCR and CR
    real<lower = 0, upper = 1> tau; // proportion of CR adults residing in observable location
    vector<lower = -8, upper = 8> [2] mu_lg [3]; // hyperprior on seasonal growth for two rivers
    vector<lower = -8, upper = 8> [4] mu_lm [3]; // hyperprior on seasonal movement rates for 2 sizes and locations
    real<lower = 0> sd_lg ; // hyperprior on seasonal growth for two rivers
    real<lower = 0> sd_lm ; // hyperprior on seasonal movement rates for 2 sizes and locations
    vector<lower = 0, upper = 1> [2] p_lcr [22]; // recapture probability in LCR for two size classes
    vector<lower = 0, upper = 1> [2] p_cr[24]; // recapture probability in CR for two size classes
    vector [4] z_ls [26]; 
    vector [4] z_lm [26]; 
    vector [2] z_lg [26];
    }
    
    transformed parameters {
    vector<lower = 0, upper = 1> [4] s [26]; // survival for all intervals, 2 locations, and 2 size classes
    vector<lower = 0, upper = 1> [2] g [26]; // growth for all intervals and 2 size classes
    vector<lower = 0, upper = 1> [4] m [26]; // movement for all intervals, 2 locations, and 2 size classes
    simplex[7] tr[7,26];
    simplex[5] p[7,26];
    vector [4] ls [26]; // logit survival for all intervals, 2 locations, and 2 size classes
    vector [2] lg [26]; // logit growth for all intervals and 2 size classes
    vector [4] lm [26]; // logit movement for all intervals, 2 locations, and 2 size classes
    
    for(j in 1:4){
    for(i in 1:26){
    ls[i,j] = mu_ls[j] + z_ls[i,j] * sd_ls;
    }
    }
    for(j in 1:4){
    for(i in 1:26){
    lm[i,j] = mu_lm[season[i],j] + z_lm[i,j] * sd_lm;
    }
    }
    for(j in 1:2){
    for(i in 1:26){
    lg[i,j] = mu_lg[season[i],j] + z_lg[i,j] * sd_lg;
    }
    }   
    
    for(j in 1:2) g[,j] = inv_logit(lg[,j]);
    for(j in 1:4) m[,j] = inv_logit(lm[,j]);
    for(j in 1:4){
    for(i in 1:26){
    s[i,j] = (inv_logit(ls[i,j]))^(season[i]>2 ? 2 : 1);
    }
    }
    
    // state transition matrix:
    for(i in 1:26){
    tr[1,i,1] = s[i,1] * (1 - g[i,1]) * (1 - m[i,1]);
    tr[1,i,2] = s[i,1] * g[i,1] * (1 - m[i,2]);
    tr[1,i,3] = s[i,1] * (1 - g[i,1]) * m[i,1] * tau;
    tr[1,i,4] = s[i,1] * g[i,1] * m[i,2] * tau;
    tr[1,i,5] = s[i,1] * (1 - g[i,1]) * m[i,1] * (1 - tau);
    tr[1,i,6] = s[i,1] * g[i,1] * m[i,2] * (1 - tau);
    tr[1,i,7] = 1 - s[i,1];
    tr[2,i,1] = 0;
    tr[2,i,2] = s[i,2] * (1 - m[i,2]);
    tr[2,i,3] = 0;
    tr[2,i,4] = s[i,2] * m[i,2] * tau;
    tr[2,i,5] = 0;
    tr[2,i,6] = s[i,2] * m[i,2] * (1 - tau);
    tr[2,i,7] = 1 - s[i,2];
    tr[3,i,1] = s[i,3] * (1 - g[i,2]) * m[i,3];
    tr[3,i,2] = s[i,3] * g[i,2] * m[i,4];
    tr[3,i,3] = s[i,3] * (1 - g[i,2]) * (1 - m[i,3]);
    tr[3,i,4] = s[i,3] * g[i,2] * (1 - m[i,4]);
    tr[3,i,5] = 0;
    tr[3,i,6] = 0;
    tr[3,i,7] = 1 - s[i,3];
    tr[4,i,1] = 0;
    tr[4,i,2] = s[i,4] * m[i,4];
    tr[4,i,3] = 0;
    tr[4,i,4] = s[i,4] * (1 - m[i,4]);
    tr[4,i,5] = 0;
    tr[4,i,6] = 0;
    tr[4,i,7] = 1 - s[i,4];
    tr[5,i,1] = tr[3,i,1];
    tr[5,i,2] = tr[3,i,2];
    tr[5,i,3] = 0;
    tr[5,i,4] = 0;
    tr[5,i,5] = tr[3,i,3];
    tr[5,i,6] = tr[3,i,4];
    tr[5,i,7] = tr[3,i,7];
    tr[6,i,1] = 0;
    tr[6,i,2] = tr[4,i,2];
    tr[6,i,3] = 0;
    tr[6,i,4] = 0;
    tr[6,i,5] = 0;
    tr[6,i,6] = tr[4,i,4];
    tr[6,i,7] = tr[4,i,7];
    for(j in 1:6){
    tr[7,i,j] = 0;
    }
    tr[7,i,7] = 1;
    }
    
    // capture probability matrix
    for(i in 1:22){
    p[1,LCRs[i],1] = p_lcr[i,1];
    p[2,LCRs[i],2] = p_lcr[i,2];
    }
    // No sampling in LCR during 1st, 4th, 7th, 10th recap periods (summers of 2009 - 2012)
    for(i in 1:4){
    p[1,LCRns[i],1] = 0;
    p[2,LCRns[i],2] = 0;
    }
    for(i in 1:24){
    p[3,CRs[i],3] = p_cr[i,1];
    p[4,CRs[i],4] = p_cr[i,2];
    }
    // No sampling in CR during 3rd and 6th recap periods (spring of 2010 & 2011)
    for(i in 1:2){
    p[3,CRns[i],3] = 0;
    p[4,CRns[i],4] = 0;
    }
    for(i in 1:26){
    p[1,i,2] = 0;
    p[1,i,3] = 0;
    p[1,i,4] = 0;
    p[1,i,5] = 1 - p[1,i,1]; // prob of fish in state 1 not being captured 
    p[2,i,1] = 0;
    p[2,i,3] = 0;
    p[2,i,4] = 0;
    p[2,i,5] = 1 - p[2,i,2]; // prob of fish in state 2 not being captured 
    p[3,i,1] = 0;
    p[3,i,2] = 0;
    p[3,i,4] = 0;
    p[3,i,5] = 1 - p[3,i,3] ; // prob of fish in state 3 not being captured 
    p[4,i,1] = 0;
    p[4,i,2] = 0;
    p[4,i,3] = 0;
    p[4,i,5] = 1 - p[4,i,4]; // prob of fish in state 4 not being captured 
    p[5,i,1] = 0;
    p[5,i,2] = 0;
    p[5,i,3] = 0;
    p[5,i,4] = 0;
    p[5,i,5] = 1; // prob of fish in state 5 not being captured (set to 100% since these fish are unobserveable) 
    p[6,i,1] = 0;
    p[6,i,2] = 0;
    p[6,i,3] = 0;
    p[6,i,4] = 0;
    p[6,i,5] = 1; // prob of fish in state 6 not being captured (set to 100% since these fish are unobserveable) 
    p[7,i,1] = 0;
    p[7,i,2] = 0;
    p[7,i,3] = 0;
    p[7,i,4] = 0;
    p[7,i,5] = 1; // prob of fish in state 7 not being captured (set to 100% since these fish are dead) 
    }
    }
    
    model{
    real temp[7]; 
    vector[7] pz[27]; 
    
    for(i in 1:26){
    z_ls[i] ~ normal(0, 1);
    z_lg[i] ~ normal(0, 1);
    z_lm[i] ~ normal(0, 1);
    }
    
    mu_ls ~ normal(0, 2);
    
    for(i in 1:3){
    mu_lm[i] ~ normal(0, 2);
    mu_lg[i] ~ normal(0, 2);
    }
    
    for(k in 1:NsumCH){  
    for(j in 1:6){
    pz[sumf[k],j] = (j == sumCH[k, sumf[k]]) * (1 - newtag[k] * 0.03); 
    }
    pz[sumf[k],7] = newtag[k] * 0.03;     
    for(t in (sumf[k] + 1):lc[k]){ 
    for(i in 1:7){ 
    for(j in 1:7){
    temp[j] = pz[t - 1,j] * tr[j,t - 1,i] * p[i,t - 1,sumCH[k,t]];
    }
    pz[t,i] = sum(temp); 
    } 
    } 
    target += (sumFR[k] * log(sum(pz[lc[k]]))); 
    }
    }
    
    //Code to estimate abundance by simulating from a negative binomial distribution
    generated quantities{
    real ptrans;	
    real<lower = 0> scale_par;
    int U[4,9];
    int N[4,9];
    int N_tot[9];
    
    for(i in 1:4){
    for(j in 1:9){
    ptrans = p[i,fall_ind[j],i] * (1 - (i > 2) * (1 - tau));
    scale_par = ptrans / (1 - ptrans);
    U[i,j] = neg_binomial_rng(catch_mat[i,j], scale_par);
    N[i,j] = U[i,j] + catch_mat[i,j];
    }
    }
    
    for(j in 1:9){
    N_tot[j]= sum(N[,j]);
    }
    }
    
    ",fill=TRUE)
sink() 

#-----------------------------------------------------------------------------#
smre.data <- list(NsumCH = dim(sumCH)[1], sumCH = sumCH, sumf = as.vector(fc),
                  season = season, sumFR = sumFR, LCRs = LCRs, LCRns = LCRns,
                  CRs = CRs, CRns = CRns, fall_ind = 1:9 * 3 - 1,
                  catch_mat = array(catch, dim = c(dim(catch)[1], dim(catch)[2])),
                  lc = as.vector(lc), newtag = as.vector(newtag.AO))

smre.inits <- function() {list(mu_ls = rep(0, 4))}

smre.params <- c("s", "g", "m", "mu_ls", "mu_lg", "mu_lm", "sd_ls", "sd_lm", 
                 "sd_lg", "tau", "p_lcr", "p_cr", "N") 

SMRE <- stan(".\\Stan_Marginalized_chub_RE.stan", 
             data = smre.data, init =smre.inits, pars = smre.params, 
             chains = 1, iter = 10, thin = 1, 
             seed = 1) 
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Application 3
#  Dynamic Community Occupancy Models - Sky Islands  
#  Discrete JAGS version 
#
#  Notes:
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2jags)

# data.dir = paste0(getwd(), "/Data/")
dat = read.csv(paste0(data.dir, "ntmb_CM_1991_till_1995.csv"))

# dimensions: species, points, years
nspp = 149
nsites = 92
nyears = 5  # 1991, 1992, 1993, 1994, 1995 
years = seq(1991, 1995, 1)
nsess = 3

# fill the data array
Y = array(NA, dim = c(nspp, nsites, nyears))
for(i in 1:5){
  Y[,,i] = as.matrix(dat[which(dat[,2] == years[i]),3:94])
}

# habitat covariate (dimensions 92 pts) 
hab = c(5L, 5L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 5L, 7L, 7L, 3L, 5L, 5L, 
        5L, 3L, 3L, 5L, 5L, 5L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
        4L, 4L, 3L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 3L, 5L, 3L, 3L, 3L, 3L, 
        3L, 5L, 5L, 5L, 2L, 2L, 2L, 4L, 2L, 3L, 5L, 4L, 5L, 4L, 4L, 3L, 
        4L, 5L, 4L, 5L, 4L, 3L, 3L, 3L, 4L, 4L, 5L, 4L, 4L, 1L, 4L, 3L, 
        5L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 6L, 6L, 7L, 7L, 7L)

#-----------------------------------------------------------------------------#
sink("Cocc_JD.jags")
cat("
    model{
    for(i in 1:(nspp)){
    a1[i] ~ dunif(0, 1)
    alpha1[i] ~ dnorm(alpha1_mu, tau_alpha1)
    
    for(h in 1:7){
    alpha2[i,h] ~ dnorm(0, tau_alpha2)
    }
    
    alpha0[i] ~ dnorm(beta, tau_u)
    mu_eta[i] <- alpha + (rho * sigma_v / sigma_u) * (alpha0[i] - beta)
    beta0[i] ~ dnorm(mu_eta[i], tau_eta)
    logit(p[i]) <- beta0[i] # detection  
    
    for(k in 1:nsites){
    z0[i,k] ~ dbern(a1[i])
    logit(psi[i,k,1]) <- alpha0[i] + alpha1[i] * z0[i,k] + alpha2[i,hab[k]]
    Z[i,k,1] ~ dbern(psi[i,k,1]) # occupancy
    mu_p[i,k,1] <- p[i] * Z[i,k,1] # detection
    Y[i,k,1] ~ dbin(mu_p[i,k,1], nsess)
    
    for(t in 1:(nyears - 1)){
    logit(psi[i,k,(t + 1)]) <- alpha0[i] + alpha1[i] * Z[i,k,t] + alpha2[i,hab[k]]
    Z[i,k,(t + 1)] ~ dbern(psi[i,k,(t + 1)]) # occupancy
    mu_p[i,k,(t + 1)] <- p[i] * Z[i,k,(t + 1)] # detection
    Y[i,k,(t + 1)] ~ dbin(mu_p[i,k,(t + 1)], nsess)
    }
    }
    }
    
    for(k in 1:nsites){  
    for(t in 1:nyears){
    SR[k,t] <- sum(Z[,k,t])
    }
    }
    
    psi_mean ~ dunif(0, 1)
    beta <- log(psi_mean) - log(1 - psi_mean)
    p_mean ~ dunif(0, 1)
    alpha <- log(p_mean) - log(1 - p_mean)
    alpha1_mean ~ dunif(0,1)
    alpha1_mu <- log(alpha1_mean) - log(1 - alpha1_mean)
    sigma_alpha1 ~ dunif(0, 5)
    sigma_alpha2 ~ dunif(0, 5)
    sigma_u ~ dunif(0, 5)
    sigma_v ~ dunif(0, 5)
    tau_alpha1 <- pow(sigma_alpha1, -2)
    tau_alpha2 <- pow(sigma_alpha2, -2)
    tau_u <- pow(sigma_u, -2)
    tau_v <- pow(sigma_v, -2)
    rho ~ dunif(-1, 1)
    tau_eta <- tau_v / (1-pow(rho,2))
    }
    ", fill = TRUE)
sink()   
#-----------------------------------------------------------------------------#
# function to make inits for JAGS 
JD_inits = function() {
  nspp = 149
  nsites = 92
  nyears = 5  # 1991, 1992, 1993, 1994, 1995 
  years = seq(1991, 1995, 1)
  nsess = 3
  
  dat = read.csv(paste0(data.dir, "ntmb_CM_1991_till_1995.csv"))
  
  # fill the data array
  Y = array(NA, dim = c(nspp, nsites, nyears))
  for(i in 1:5){
    Y[,,i] = as.matrix(dat[which(dat[,2] == years[i]),3:94])
  }
  
  a1 = runif(nspp, 0.25, 1) 
  Ztemp = array(rbinom(nspp * 5 * nsites, size = 1, prob = 0.5), dim = c(nspp, nsites, 5))
  z0temp = array(rbinom(nspp * nsites, size = 1, prob = 0.5), dim = c(nspp, nsites))
  
  for(t in 1:5){
    for(i in 1:nspp){
      for(k in 1:nsites){
        if(Ztemp[i,k,t] == 0){
          if(!is.na(Y[i,k,t])){
            if(Y[i,k,t] > 0){
              Ztemp[i,k,t] = 1              
            }
          }
        }
      }
    }
  }
  list(a1 = a1, z0 = z0temp, Z = Ztemp)
}
#-----------------------------------------------------------------------------#
JD_data = list('nspp' = nspp, 'nsites' = nsites, 'nsess' = nsess,
               'nyears' = nyears, 'hab' = hab, 'Y' = Y)

params <- c('psi_mean', 'p_mean', 'rho', 'sigma_v', 'sigma_u', 'sigma_alpha1',
            'sigma_alpha2', 'alpha1_mean', 'SR', 'Z')

JD_Cocc <- jags.parallel(JD_data, inits = JD_inits,
                         params, "Cocc_JD.jags",
                         n.chains = 3, n.iter = 10,
                         export_obj_names = c("data.dir"))
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Application 3
#  Dynamic Community Occupancy Models - Sky Islands  
#  Marginalized JAGS version 
#
#  Notes:
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2jags)

source(paste0(getwd(),"/Functions.R"), chdir = F)

# data.dir = paste0(getwd(), "/Data/")
dat = read.csv(paste0(data.dir, "ntmb_CM_1991_till_1995.csv"))

# dimensions: species, points, years
nspp = 149
nsites = 92
nyears = 5  # 1991, 1992, 1993, 1994, 1995 
years = seq(1991, 1995, 1)
nsess = 3

# fill the data array
Y = array(NA, dim = c(nspp, nsites, nyears))
for(i in 1:5){
  Y[,,i] = as.matrix(dat[which(dat[,2] == years[i]),3:94])
}

visit <- ifelse(is.na(Y[1,,1]) == TRUE, 2, 1)
sY <- ifelse(is.na(Y) == TRUE, 1, Y + 1)

# habitat covariate (dimensions 92 pts) 
hab = c(5L, 5L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 5L, 7L, 7L, 3L, 5L, 5L, 
        5L, 3L, 3L, 5L, 5L, 5L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
        4L, 4L, 3L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 3L, 5L, 3L, 3L, 3L, 3L, 
        3L, 5L, 5L, 5L, 2L, 2L, 2L, 4L, 2L, 3L, 5L, 4L, 5L, 4L, 4L, 3L, 
        4L, 5L, 4L, 5L, 4L, 3L, 3L, 3L, 4L, 4L, 5L, 4L, 4L, 1L, 4L, 3L, 
        5L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 6L, 6L, 7L, 7L, 7L)

#-----------------------------------------------------------------------------#
sink("Cocc_JM.jags")
cat("
    model{
    for(i in 1:(nspp)){
    a1[i] ~ dunif(0, 1)
    z0[i,1] <- a1[i]
    z0[i,2] <- 1 - a1[i]
    alpha1[i] ~ dnorm(alpha1_mu, tau_alpha1)
    for(h in 1:7){
    alpha2[i,h] ~ dnorm(0, tau_alpha2)
    }
    alpha0[i] ~ dnorm(beta, tau_u)
    mu_eta[i] <- alpha + (rho * sigma_v / sigma_u) * (alpha0[i] - beta)
    beta0[i] ~ dnorm(mu_eta[i], tau_eta)
    logit(p[i]) <- beta0[i] # detection
    po[i,1,1,2] = 1;
    po[i,2,1,2] = 1;
    for (j in 1:nsess){
    po[i,1,(j + 1),1] <- 0
    po[i,2,(j + 1),1] <- (p[i] ^ j) * (1 - p[i]) ^ (nsess - j)
    }
    po[i,1,1,1] <- 1
    po[i,2,1,1] <- (1 - p[i]) ^ (nsess)
    for(k in 1:nhab){
    logit(tr[i,k,1,2]) <- alpha0[i] + alpha2[i,k]
    logit(tr[i,k,2,2]) <- alpha0[i] + alpha1[i] + alpha2[i,k]
    tr[i,k,1,1] <- 1 - tr[i,k,1,2]
    tr[i,k,2,1] <- 1 - tr[i,k,2,2]
    }
    }
    
    for(i in 1:NuDH){
    pz[i,1,1] = (z0[spp[i],1] * tr[spp[i],hab2[i],1,1] + z0[spp[i],2] * tr[spp[i],hab2[i],2,1]) *
    po[spp[i],1,sY2[i,1],visit2[i]];
    pz[i,1,2] = (z0[spp[i],1] * tr[spp[i],hab2[i],1,2] + z0[spp[i],2] * tr[spp[i],hab2[i],2,2]) * 
    po[spp[i],2,sY2[i,1],visit2[i]];
    Z[i,1] = pz[i,1,2] / (pz[i,1,1] + pz[i,1,2]);
    for(t in 1:(nyears - 1)){
    pz[i,(t + 1),1] = (pz[i,t,1] * tr[spp[i],hab2[i],1,1] + pz[i,t,2] * tr[spp[i],hab2[i],2,1]) * 
    po[spp[i],1,sY2[i,(t + 1)],1];
    pz[i,(t + 1),2] = (pz[i,t,1] * tr[spp[i],hab2[i],1,2] + pz[i,t,2] * tr[spp[i],hab2[i],2,2]) *
    po[spp[i],2,sY2[i,(t + 1)],1];
    Z[i,(t + 1)] = pz[i,(t + 1),2] / (pz[i,(t + 1),1] + pz[i,(t + 1),2]);
    }
    lik[i] <- sum(pz[i, nyears,])
    fr[i] ~ dbin(lik[i], FR[i])
    }
    
    for (k in 1:nsites){  
    for (t in 1:nyears){
    for (i in 1:nspp){
    Z2[i,k,t] ~ dbern(Z[lookup[i,k],t])
    }
    SR[k,t] <- sum(Z2[1:nspp,k,t])
    }
    }
    
    psi_mean ~ dunif(0, 1)
    beta <- log(psi_mean) - log(1 - psi_mean)
    p_mean ~ dunif(0, 1)
    alpha <- log(p_mean) - log(1 - p_mean)
    alpha1_mean ~ dunif(0, 1)
    alpha1_mu <- log(alpha1_mean) - log(1 - alpha1_mean)
    sigma_alpha1 ~ dunif(0, 5)
    sigma_alpha2 ~ dunif(0, 5)
    sigma_u ~ dunif(0, 5)
    sigma_v ~ dunif(0, 5)
    tau_alpha1 <- pow(sigma_alpha1, -2)
    tau_alpha2 <- pow(sigma_alpha2, -2)
    tau_u <- pow(sigma_u, -2)
    tau_v <- pow(sigma_v, -2)
    rho ~ dunif(-1, 1)
    tau_eta <- tau_v / (1 - pow(rho, 2)) 
    }
    ", fill = TRUE)
sink()   
#-----------------------------------------------------------------------------#
JM_inits = function() {
  nspp = 149
  nsites = 92
  nyears = 5 
  Z2 = array(rbinom(nspp * nsites * nyears, 1, 0.5), dim = c(nspp, nsites, nyears))
  list(Z2 = Z2)
}

pasted <- function(x){
  paste0(x[1], x[2], x[3], x[4], x[5], x[6], x[7], collapse = "")
}

unpaste <- function(x){
  as.numeric(c(substr(x, 1, 1), substr(x, 2, 2), substr(x, 3, 3),
               substr(x, 4, 4), substr(x, 5, 5), substr(x, 6, 6), substr(x, 7, 7)))
}

#-----------------------------------------------------------------------------#
FR <- numeric()
spp <- numeric()
sY2 <- rep(NA, 5)
visit2 <- numeric()
hab2 <- numeric()
lookup <- matrix(NA, nspp, nsites)
temp4 <- 0

for(S in 1:nspp){
  temp <- apply(cbind(sY[S,,], visit,hab), 1, pasted)
  temp2 <- table(temp)
  FR <- c(FR, as.numeric(temp2))
  spp <- c(spp, rep(S, length(temp2)))
  for(j in 1:length(temp2)){
    temp3 <- unpaste(names(temp2[j]))
    sY2 <- rbind(sY2, temp3[1:5])
    visit2 <- c(visit2, temp3[6])
    hab2 <- c(hab2, temp3[7])
  }
  lookup[S,] <- match(temp, names(temp2)) + temp4
  temp4 <- temp4 + length(temp2)
}

sY2 <- sY2[-1,]

#-----------------------------------------------------------------------------#
JM_data =list('nspp' = nspp, 'nsites' = nsites,'nsess' = nsess, 'nyears' = nyears,
              'hab2' = hab2, 'fr' = FR, 'FR' = FR, 'nhab' = 7, 'NuDH' = length(FR),
              'spp' = spp, 'lookup' = lookup, 'sY2' = sY2, 'visit2' = visit2)

params <- c('psi_mean', 'p_mean', 'rho', 'sigma_v', 'sigma_u', 'sigma_alpha1',
            'sigma_alpha2', 'alpha1_mean','SR', 'Z2')

JM_Cocc <- jags.parallel(JM_data, inits = JM_inits, params, 'Cocc_JM.jags',
                         n.chains = 3, n.iter = 10)
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Application 3
#  Dynamic Community Occupancy Models - Sky Islands  
#  Marginalized Stan version 
#
#  Notes:
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())  

# data.dir = paste0(getwd(), "/Data/")
dat = read.csv(paste0(data.dir, "ntmb_CM_1991_till_1995.csv"))

# dimensions: species, points, years
nspp = 149
nsites = 92
nyears = 5  # 1991, 1992, 1993, 1994, 1995 
years = seq(1991, 1995, 1)
nsess = 3

# fill the data array
Y = array(NA, dim = c(nspp, nsites, nyears))
for(i in 1:5){
  Y[,,i] = as.matrix(dat[which(dat[,2] == years[i]),3:94])
}

visit <- ifelse(is.na(Y[1,,1]) == TRUE, 2, 1)
sY <- ifelse(is.na(Y) == TRUE, 1, Y + 1)

# habitat covariate (dimensions 92 pts) 
hab = c(5L, 5L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 5L, 7L, 7L, 3L, 5L, 5L, 
        5L, 3L, 3L, 5L, 5L, 5L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
        4L, 4L, 3L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 3L, 5L, 3L, 3L, 3L, 3L, 
        3L, 5L, 5L, 5L, 2L, 2L, 2L, 4L, 2L, 3L, 5L, 4L, 5L, 4L, 4L, 3L, 
        4L, 5L, 4L, 5L, 4L, 3L, 3L, 3L, 4L, 4L, 5L, 4L, 4L, 1L, 4L, 3L, 
        5L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 6L, 6L, 7L, 7L, 7L)
#-----------------------------------------------------------------------------#
sink("Cocc_SM.stan")
cat("
    data{
    int<lower = 1> nspp;
    int<lower = 1> nsess;
    int<lower = 1> nsites;
    int<lower = 1> nyears;
    int<lower = 1> NuDH; 
    int<lower = 1> nhab; 
    int<lower = 1> sY2 [NuDH,nyears]; 
    int<lower = 1> visit2 [NuDH]; 
    int<lower = 1> hab2 [NuDH]; 
    int<lower = 1> spp [NuDH];
    int<lower = 1> FR [NuDH];
    int<lower = 1> lookup[nspp,nsites];
    }
    
    parameters{
    real<lower = 0, upper = 1> psi_mean;
    real<lower = 0, upper = 1> p_mean;
    real<lower = 0, upper = 1> alpha1_mean;
    real<lower = 0, upper = 5> sigma_alpha1;
    real<lower = 0, upper = 5> sigma_alpha2;
    real<lower = 0, upper = 5> sigma_u;
    real<lower = 0, upper = 5> sigma_v;
    real<lower = -1, upper = 1> rho;
    vector<lower = 0, upper = 1> [nspp] a1;
    vector [nspp] alpha_dev1;
    vector [nspp] alpha_dev0;
    vector [nspp] alpha_dev2 [7];
    real beta_dev0 [nspp];
    }
    
    transformed parameters{
    real beta;
    real alpha;
    real alpha1_mu;
    real sigma_eta;
    real Z [NuDH,nyears];
    real p [nspp];
    simplex [2] z0 [nspp];
    real po [nspp,2,(nsess + 1),2];
    real tr [nspp,nhab,2,2];
    real pz [NuDH,nyears,2];
    
    beta = logit(psi_mean);
    alpha = logit(p_mean);
    alpha1_mu = logit(alpha1_mean);
    sigma_eta = sigma_v * ((1 - rho ^ 2) ^ 0.5);
    
    for(i in 1:(nspp)){
    z0[i,2] = a1[i];
    z0[i,1] = 1 - a1[i];
    po[i,1,1,2] = 1;
    po[i,2,1,2] = 1;
    po[i,1,1,1] = 1;
    p[i] = inv_logit(alpha + rho * sigma_v * alpha_dev0[i] + sigma_eta * beta_dev0[i]);
    for(j in 1:nsess){
    po[i,1,(j + 1),1] = 0;
    po[i,1,(j + 1),2] = 0;
    po[i,2,(j + 1),2] = 0;			
    po[i,2,(j + 1),1] = (p[i] ^ j) * (1 - p[i]) ^ (nsess - j);
    }
    po[i,2,1,1] = (1 - p[i]) ^ (nsess);
    for(k in 1:nhab){
    tr[i,k,2,2] = inv_logit(beta + sigma_u * alpha_dev0[i] + alpha1_mu + sigma_alpha1 * alpha_dev1[i] + sigma_alpha2 * alpha_dev2[k,i]);
    tr[i,k,1,2] = inv_logit(beta + sigma_u * alpha_dev0[i] + sigma_alpha2 * alpha_dev2[k,i]);
    tr[i,k,1,1] = 1 - tr[i,k,1,2];
    tr[i,k,2,1] = 1 - tr[i,k,2,2];
    }
    }
    
    for(i in 1:NuDH){
    pz[i,1,1] = (z0[spp[i],1] * tr[spp[i],hab2[i],1,1] + z0[spp[i],2] * tr[spp[i],hab2[i],2,1]) *
    po[spp[i],1,sY2[i,1],visit2[i]];
    pz[i,1,2] = (z0[spp[i],1] * tr[spp[i],hab2[i],1,2] + z0[spp[i],2] * tr[spp[i],hab2[i],2,2]) *
    po[spp[i],2,sY2[i,1],visit2[i]];
    Z[i,1] = pz[i,1,2] / (pz[i,1,1] + pz[i,1,2]);
    for(t in 1:(nyears - 1)){
    pz[i,(t + 1),1] = (pz[i,t,1] * tr[spp[i],hab2[i],1,1] + pz[i,t,2] * tr[spp[i],hab2[i],2,1]) * 
    po[spp[i],1,sY2[i,(t + 1)],1];
    pz[i,(t + 1),2] = (pz[i,t,1] * tr[spp[i],hab2[i],1,2] + pz[i,t,2] * tr[spp[i],hab2[i],2,2]) * 
    po[spp[i],2,sY2[i,(t + 1)],1];
    Z[i,(t + 1)] = pz[i,(t + 1),2] / (pz[i,(t + 1),1] + pz[i,(t + 1),2]);
    }
    }
    }
    
    model{	
    alpha_dev1 ~ normal(0, 1);
    alpha_dev0 ~ normal(0, 1);
    
    for(h in 1:7){
    alpha_dev2[h] ~ normal(0, 1);
    }
    
    beta_dev0 ~ normal(0, 1);
    
    for(i in 1:NuDH){
    target += FR[i] * log(sum(pz[i,nyears,]));
    }
    }
    
    generated quantities{
    int SR [nsites,nyears];
    int temp [nspp];  
    
    for(k in 1:nsites){  
    for(t in 1:nyears){
    for(i in 1:nspp){
    temp[i] = bernoulli_rng(Z[lookup[i,k],t]);
    }
    SR[k,t] = sum(temp);
    }
    }
    }
    
    ", fill = TRUE)
sink()
#-----------------------------------------------------------------------------#
pasted <- function(x){
  paste0(x[1], x[2], x[3], x[4], x[5], x[6], x[7], collapse = "")
}

unpaste <- function(x){
  as.numeric(c(substr(x, 1, 1), substr(x, 2, 2), substr(x, 3, 3),
               substr(x, 4, 4), substr(x, 5, 5), substr(x, 6, 6), substr(x, 7, 7)))
}

FR <- numeric()
spp <- numeric()
sY2 <- rep(NA, 5)
visit2 <- numeric()
hab2 <- numeric()
lookup <- matrix(NA, nspp, nsites)
temp4 <- 0

for(S in 1:nspp){
  temp <- apply(cbind(sY[S,,], visit,hab), 1, pasted)
  temp2 <- table(temp)
  FR <- c(FR, as.numeric(temp2))
  spp <- c(spp, rep(S, length(temp2)))
  for(j in 1:length(temp2)){
    temp3 <- unpaste(names(temp2[j]))
    sY2 <- rbind(sY2, temp3[1:5])
    visit2 <- c(visit2, temp3[6])
    hab2 <- c(hab2, temp3[7])
  }
  lookup[S,] <- match(temp, names(temp2)) + temp4
  temp4 <- temp4 + length(temp2)
}

sY2 <- sY2[-1,]

#-----------------------------------------------------------------------------#
SM_data <- list(nspp = nspp, nsess = nsess, nsites = nsites, nyears = nyears, 
                hab2 = hab2, sY2 = sY2, visit2 = visit2, nhab = 7, NuDH = length(FR), 
                spp = spp, FR = FR, lookup = lookup)

params <- c('psi_mean', 'p_mean', 'rho', 'sigma_v', 'sigma_u', 'sigma_alpha1',
            'sigma_alpha2', 'alpha1_mean', 'SR')

SM_Cocc <- stan("Cocc_SM.stan",
                data = SM_data,
                pars = params,
                chains = 1, iter = 1000) 
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Application 4
#  Fitting an Integrated Population Model to Brown Trout Data
#  Discrete JAGS version 
#
#  Notes:
#  * The model runs from the fall of 2000 to fall of 2017 on a seasonal basis
#  * We define three size states based on total length in mm 
#    - 0 - 200; 200 - 350; 350 +
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2jags)

# read in data
NO_catch <- read.csv(paste0(getwd(), "/Data/", "NO_catch.csv"), header = TRUE)
AZGF_catch <- read.csv(paste0(getwd(), "/Data/", "AZGF_catch.csv"), header = TRUE)
MR_data <- read.csv(paste0(getwd(), "/Data/", "bCH.csv"), header = FALSE)

# extract/reformat data
bNOc <- NO_catch[,1:3]
NOpasses <- NO_catch$NOpasses
seasNO <- NO_catch$seasNO
spawn <- ifelse(seasNO == 1, 4, 3)
bAZ <- AZGF_catch[,1:3]
ts <- AZGF_catch$ts
AZeff <- AZGF_catch$AZeff
NAZsamps <- length(AZeff)

findlast <- function(x){ifelse(x[24] == 1, 23, max(which(x[1:23] != 4)))}
indlast <- apply(MR_data, 1, findlast)

indCH <- MR_data[,1:23]
NindCH <- length(indlast)
findfirst <- function(x){which(x != 4)[1]}
indf <- apply(indCH, 1, findfirst)

#-----------------------------------------------------------------------------#
sink("JAGS_Discrete.jags")
cat("
    model{
    # lphi is the logit of survival and is given a prior based on the Lorenzen
    # function and the average mass of fish in each size class during each season.
    # variation from the priod mode is determined by an estimated variance
    # parameter (sd_lphi)
    
    lphi[1,1] ~ dnorm(1.08, tau_lphi)
    lphi[1,2] ~ dnorm(1.14, tau_lphi)
    lphi[1,3] ~ dnorm(1.26, tau_lphi)
    lphi[1,4] ~ dnorm(1.38, tau_lphi)
    lphi[2,1] ~ dnorm(1.96, tau_lphi)
    lphi[2,2] ~ dnorm(1.96, tau_lphi)
    lphi[2,3] ~ dnorm(2.02, tau_lphi)
    lphi[2,4] ~ dnorm(2.02, tau_lphi)
    lphi[3,1] ~ dnorm(2.29, tau_lphi)
    lphi[3,2] ~ dnorm(2.29, tau_lphi)
    lphi[3,3] ~ dnorm(2.29, tau_lphi)
    lphi[3,4] ~ dnorm(2.29, tau_lphi)
    tau_lphi <- pow(sd_lphi, -2)
    sd_lphi ~ dunif(0.01,4)
    
    for(j in 1:3){
    for(k in 1:4){
    logit(bphi[j,k]) <- lphi[j,k]
    }
    }
    
    bpsi2 ~ dunif(0,1) # growth of size class 2 fish into size class 3
    
    # define transition matrix that combines survival and growth parameters
    for(i in 1:4){
    btrans[1,3,i] <- 0
    btrans[1,4,i] <- 1 - bphi[1,i]
    btrans[2,1,i] <- 0
    btrans[2,2,i] <- bphi[2,i] * (1 - bpsi2)
    btrans[2,3,i] <- bphi[2,i] * bpsi2
    btrans[2,4,i] <- 1 - bphi[2,i]
    btrans[3,1,i] <- 0
    btrans[3,2,i] <- 0
    btrans[3,3,i] <- bphi[3,i]
    btrans[3,4,i] <- 1 - bphi[3,i]
    btrans[4,1,i] <- 0
    btrans[4,2,i] <- 0
    btrans[4,3,i] <- 0
    btrans[4,4,i] <- 1
    }
    
    # size class one transitions are done separately because growth is allowed to vary between seasons
    for(i in 1:3){
    bpsi1[i] ~ dunif(0,1)
    btrans[1,1,i] <- bphi[1,i] * (1 - bpsi1[i])
    btrans[1,2,i] <- bphi[1,i] * bpsi1[i]
    }
    
    btrans[1,1,4] <- 0
    btrans[1,2,4] <- bphi[1,4]
    
    bpi <- .08 # proportion of brown trout population in NO reach 1
    tau_blp <- pow(sd_blp, -2)
    sd_blp ~ dunif(0.1, 2) # trip to trip deviation in NO pcaps
    
    # mean pcaps per pass on a logit scale for three size classes, plus largest
    # size class during spawning season
    for(i in 1:4){
    mu_blp[i] ~ dnorm(-3, .25) 
    }
    
    # this loop calculates actual per pass pcaps for each trip and modifies based on # of passes  
    for(j in 1:23){
    blp_pass[j,1] ~ dnorm(mu_blp[1], tau_blp)
    blp_pass[j,2] ~ dnorm(mu_blp[2], tau_blp)
    blp_pass[j,3] ~ dnorm(mu_blp[spawn[j]], tau_blp)
    
    for(k in 1:3){
    logit(bp_pass[j,k]) <- blp_pass[j,k]
    bp[j,k,k] <- 1 - pow((1 - bp_pass[j,k]), NOpasses[j])
    bp[j,k,4] <- 1 - bp[j,k,k]
    }
    
    bp[j,1,2] <- 0
    bp[j,1,3] <- 0
    bp[j,2,1] <- 0
    bp[j,2,3] <- 0
    bp[j,3,1] <- 0
    bp[j,3,2] <- 0
    bp[j,4,1] <- 0
    bp[j,4,2] <- 0
    bp[j,4,3] <- 0
    bp[j,4,4] <- 1
    }
    
    for(k in 1:NindCH){
    Z[k,indf[k]] <- indCH[k,indf[k]]    
    for(j in indf[k]:(indlast[k] - 1)){
    Z[k,(j + 1)] ~ dcat(btrans[Z[k,j], ,seasNO[(j + 1)]])
    indCH[k,(j + 1)] ~ dcat(bp[j,Z[k,(j + 1)],])
    }
    }
    
    # calculate offset for each size classes of AZGF effort and calculate expected pcap
    AZadj[1] ~ dnorm(0,1)
    AZadj[2] ~ dnorm(0,1)
    AZadj[3] ~ dnorm(0,1)
    mu_AZ[1] <- mu_blp[1] + AZadj[1]
    mu_AZ[2] <- mu_blp[2] + AZadj[2]
    mu_AZ[3] <- mu_blp[3] + AZadj[3]
    IN[1] <- 0 # initial abundances of size class 1 fish
    IN[2] ~ dunif(0, 1000) # initial abundances of size class 2 fish
    IN[3] ~ dunif(0, 1000) # initial abundances of size class 3 fish
    bN[1,1] <- IN[1]
    bN[1,2] <- IN[2]
    bN[1,3] <- IN[3]
    
    # variance term controlling unexplained variation in reproductive rate (BETA) 
    tau_beta <- pow(sd_beta,-2)
    sd_beta ~ dunif(0.1,4)
    
    # log of the median reproductive rate - i.e., an intercept
    lbeta_0 ~ dunif(-6,0)
    
    # log of the median immigration rate of large brown trout - i.e., the intercept
    mu_I ~ dunif(0,6)
    
    # variance term controlling unexplained variation in immigration
    tau_I <- pow(sd_I,-2)
    sd_I ~ dunif(0.01,3)
    
    # calculate actual immigration in each interval on log scale
    for(j in 1:68){
    I[j] ~ dnorm(mu_I, tau_I)
    }
    
    # calculate latent abundance of brown trout from fall 2000 to end of 2017
    for (j in 1:17){
    for (k in 1:3){
    bN[((j-1) * 4 + k + 1),1] <- btrans[1,1,k] * bN[((j-1) * 4 + k),1]
    bN[((j-1) * 4 + k + 1),2] <- btrans[1,2,k] * bN[((j-1) * 4 + k),1] + btrans[2,2,k] * bN[((j-1) * 4 + k),2]
    bN[((j-1) * 4 + k + 1),3] <- btrans[2,3,k] * bN[((j-1) * 4 + k),2] + btrans[3,3,k] * bN[((j-1) * 4 + k),3] + exp(I[((j-1) * 4 + k)])
    }
    
    # BNT recruits produced in fall as a function weighted sum of adults (wA) and reprodutive rate (Beta) in winter
    wA[j] <- (bN[((j - 1) * 4 + 2),2] + 4 * bN[((j - 1) * 4 + 2),3])
    beta_eps[j] ~ dnorm(0, tau_beta)
    Beta[j] <- exp(lbeta_0 + beta_eps[j]) 
    
    # between summer and fall all bnt graduate to sz 2 & recruits show up
    bN[(1 + j * 4),1] <- wA[j] * Beta[j]
    bN[(1 + j * 4),2] <- btrans[1,2,4] * bN[(j * 4),1] + btrans[2,2,4] * bN[(j * 4),2]
    bN[(1 + j * 4),3] <- btrans[2,3,4] * bN[(j * 4),2] + btrans[3,3,4] * bN[(j * 4),3] + exp(I[(j * 4)])
    }
    
    # 2000 - 2017 AZGF data
    for(j in 1:NAZsamps){
    for(k in 1:3){
    blpAZ[j,k] ~ dnorm(mu_AZ[k], tau_blp)
    logit(bpAZ[j,k]) <- blpAZ[j,k]
    blamAZ[j,k] <- bpAZ[j,k] * bN[ts[j],k] * AZeff[j] / 35 # predicted catch AZ (35 is ~h to do LF, by AZ)
    bAZ[j,k] ~ dpois(blamAZ[j,k])
    }
    }
    
    # 2012 - 2017 NO: starts in april 2012
    for(j in 1:23){
    for(k in 1:3){
    blamNO[j,k] <- bp[j,k,k] * bpi * bN[(j + 46),k]
    bNOc[j,k] ~ dpois(blamNO[j,k])
    }
    }
    
    }
    ", fill = TRUE)
sink()

#-----------------------------------------------------------------------------#
JD_inits <- function(){
  z.init <- matrix(NA, nrow = NindCH, ncol = 23)
  for(i in 1:NindCH){
    z.init[i,indf[i]] <- indCH[i,indf[i]]
    for(j in (indf[i] + 1):indlast[i]){
      z.init[i,j] <- ifelse(indCH[i,j] < 4, indCH[i,j],
                            ifelse(z.init[i,(j - 1)] > 1, z.init[i,(j - 1)],
                                   ifelse(seasNO[j] == 4, 2, 1)))
    }
    z.init[i,indf[i]] <- NA
  }	
  list(Z = array(z.init, dim = c(NindCH, 23)))
}

#-----------------------------------------------------------------------------#
BM_JM.data <- list(NAZsamps = NAZsamps, ts = ts, AZeff = AZeff, bAZ = bAZ,
                   seasNO = seasNO, bNOc = bNOc, NOpasses = NOpasses, 
                   indlast = indlast, indCH = indCH, NindCH = NindCH, indf = indf,
                   spawn = spawn)

BM_JM.par <- c('bphi', 'bpsi1', 'bpsi2', 'mu_blp', 'sd_blp', 'lbeta_0',
               'mu_I', 'I', 'Beta', 'IN', 'AZadj', 'sd_I', 'sd_lphi',
               'sd_beta', 'bN', 'bp_pass')

jags.fit <- jags.parallel(BM_JM.data, inits = JD_inits, BM_JM.par, "JAGS_Discrete.jags",
                          n.chains = 3, n.iter = 10, export_obj_names = c("JD_inits"))
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Application 4
#  Fitting an Integrated Population Model to Brown Trout Data
#  Marginalized JAGS version 
#
#  Notes:
#  * The model runs from the fall of 2000 to fall of 2017 on a seasonal basis
#  * We define three size states based on total length in mm 
#    - 0 - 200; 200 - 350; 350 +
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2jags)

# read in data
NO_catch <- read.csv(paste0(getwd(), "/Data/", "NO_catch.csv"), header = TRUE)
AZGF_catch <- read.csv(paste0(getwd(), "/Data/", "AZGF_catch.csv"), header = TRUE)
MR_data <- read.csv(paste0(getwd(), "/Data/", "bCH.csv"), header = FALSE)

#-----------------------------------------------------------------------------#
# extract/reformat data
bNOc <- NO_catch[,1:3]
NOpasses <- NO_catch$NOpasses
seasNO <- NO_catch$seasNO
spawn <- ifelse(seasNO == 1, 4, 3)
bAZ <- AZGF_catch[,1:3]
ts <- AZGF_catch$ts
AZeff <- AZGF_catch$AZeff
NAZsamps <- length(AZeff)

# capture-recapture data
allCH <- MR_data[,1:23]

bCH = collapse.ch(allCH)[[1]]
FR = collapse.ch(allCH)[[2]]

findlast <- function(x){ifelse(x[23] == 1, 22, max(which(x[1:22] != 4)))}
last <- apply(bCH, 1, findlast)

NCH <- length(last)
findfirst <- function(x){which(x != 4)[1]}
sumf <- apply(bCH, 1, findfirst)

#-----------------------------------------------------------------------------#
sink("JAGS_Marginalized.jags")
cat("
    model{
    # lphi is the logit of survival and is given a prior based on the Lorenzen
    # function and the average mass of fish in each size class during each season.
    # variation from the priod mode is determined by an estimated variance
    # parameter (sd_lphi)
    
    lphi[1,1] ~ dnorm(1.08, tau_lphi)
    lphi[1,2] ~ dnorm(1.14, tau_lphi)
    lphi[1,3] ~ dnorm(1.26, tau_lphi)
    lphi[1,4] ~ dnorm(1.38, tau_lphi)
    lphi[2,1] ~ dnorm(1.96, tau_lphi)
    lphi[2,2] ~ dnorm(1.96, tau_lphi)
    lphi[2,3] ~ dnorm(2.02, tau_lphi)
    lphi[2,4] ~ dnorm(2.02, tau_lphi)
    lphi[3,1] ~ dnorm(2.29, tau_lphi)
    lphi[3,2] ~ dnorm(2.29, tau_lphi)
    lphi[3,3] ~ dnorm(2.29, tau_lphi)
    lphi[3,4] ~ dnorm(2.29, tau_lphi)
    tau_lphi <- pow(sd_lphi, -2)
    sd_lphi ~ dunif(0.01,4)
    
    for(j in 1:3){
    for(k in 1:4){
    logit(bphi[j,k]) <- lphi[j,k]
    }
    }
    
    bpsi2 ~ dunif(0,1) # growth of size class 2 fish into size class 3
    
    # define transition matrix that combines survival and growth parameters
    for(i in 1:4){
    btrans[1,3,i] <- 0
    btrans[1,4,i] <- 1 - bphi[1,i]
    btrans[2,1,i] <- 0
    btrans[2,2,i] <- bphi[2,i] * (1 - bpsi2)
    btrans[2,3,i] <- bphi[2,i] * bpsi2
    btrans[2,4,i] <- 1 - bphi[2,i]
    btrans[3,1,i] <- 0
    btrans[3,2,i] <- 0
    btrans[3,3,i] <- bphi[3,i]
    btrans[3,4,i] <- 1 - bphi[3,i]
    btrans[4,1,i] <- 0
    btrans[4,2,i] <- 0
    btrans[4,3,i] <- 0
    btrans[4,4,i] <- 1
    }
    
    # size class one transitions are done separately because growth is allowed to vary between seasons
    for(i in 1:3){
    bpsi1[i] ~ dunif(0,1)
    btrans[1,1,i] <- bphi[1,i] * (1 - bpsi1[i])
    btrans[1,2,i] <- bphi[1,i] * bpsi1[i]
    }
    
    btrans[1,1,4] <- 0
    btrans[1,2,4] <- bphi[1,4]
    
    bpi <- .08 # proportion of brown trout population in NO reach 1
    tau_blp <- pow(sd_blp, -2)
    sd_blp ~ dunif(0.1, 2) # trip to trip deviation in NO pcaps
    
    for(i in 1:4){
    mu_blp[i] ~ dnorm(-3, .25) # mean pcaps per pass on a logit scale for three size classes, plus largest size class during spawning season
    }
    
    # this loop calculates actual per pass pcaps for each trip and modifies based on # of passes  
    for(j in 1:23){
    # spawn[j] <- 3 + step(-1 * seasNO[j] + 1.1)  # change here to use 'spawn' input
    blp_pass[j,1] ~ dnorm(mu_blp[1], tau_blp)
    blp_pass[j,2] ~ dnorm(mu_blp[2], tau_blp)
    blp_pass[j,3] ~ dnorm(mu_blp[spawn[j]], tau_blp)
    
    for(k in 1:3){
    logit(bp_pass[j,k]) <- blp_pass[j,k]
    bp[j,k,k] <- 1 - pow((1 - bp_pass[j,k]), NOpasses[j])
    bp[j,k,4] <- 1 - bp[j,k,k]
    }
    
    bp[j,1,2] <- 0
    bp[j,1,3] <- 0
    bp[j,2,1] <- 0
    bp[j,2,3] <- 0
    bp[j,3,1] <- 0
    bp[j,3,2] <- 0
    bp[j,4,1] <- 0
    bp[j,4,2] <- 0
    bp[j,4,3] <- 0
    bp[j,4,4] <- 1
    }
    
    for(k in 1:NCH){
    pz[k,sumf[k],1] <- equals(bCH[k,sumf[k]], 1)
    pz[k,sumf[k],2] <- equals(bCH[k,sumf[k]], 2)
    pz[k,sumf[k],3] <- equals(bCH[k,sumf[k]], 3)
    pz[k,sumf[k],4] <- 0
    
    for(j in sumf[k]:(last[k] - 1)){
    for(i in 1:4){
    pz[k,(j + 1),i] <- inprod(pz[k,j,], btrans[,i,seasNO[(j + 1)]]) * bp[j,i,bCH[k,(j + 1)]]
    }
    }
    
    ll[k] <- sum(pz[k, last[k],])
    ones[k] ~ dbin(ll[k], FR[k])
    }
    
    # calculate offset for each size classes of AZGF effort and calculate expected pcap
    AZadj[1] ~ dnorm(0,1)
    AZadj[2] ~ dnorm(0,1)
    AZadj[3] ~ dnorm(0,1)
    mu_AZ[1] <- mu_blp[1] + AZadj[1]
    mu_AZ[2] <- mu_blp[2] + AZadj[2]
    mu_AZ[3] <- mu_blp[3] + AZadj[3]
    IN[1] <- 0 # initial abundances of size class 1 fish
    IN[2] ~ dunif(0, 1000) # initial abundances of size class 2 fish
    IN[3] ~ dunif(0, 1000) # initial abundances of size class 3 fish
    bN[1,1] <- IN[1]
    bN[1,2] <- IN[2]
    bN[1,3] <- IN[3]
    
    # variance term controlling unexplained variation in reproductive rate (BETA) 
    tau_beta <- pow(sd_beta,-2)
    sd_beta ~ dunif(0.1,4)
    
    # log of the median reproductive rate - i.e., an intercept
    lbeta_0 ~ dunif(-6,0)
    
    # log of the median immigration rate of large brown trout - i.e., the intercept
    mu_I ~ dunif(0,6)
    
    # variance term controlling unexplained variation in immigration
    tau_I <- pow(sd_I,-2)
    sd_I ~ dunif(0.01,3)
    
    # calculate actual immigration in each interval on log scale
    for(j in 1:68){
    I[j] ~ dnorm(mu_I, tau_I)
    }
    
    # calculate latent abundance of brown trout from fall 2000 to end of 2017
    for(j in 1:17){
    for(k in 1:3){
    bN[((j-1) * 4 + k + 1),1] <- btrans[1,1,k] * bN[((j-1) * 4 + k),1]
    bN[((j-1) * 4 + k + 1),2] <- btrans[1,2,k] * bN[((j-1) * 4 + k),1] + btrans[2,2,k] * bN[((j-1) * 4 + k),2]
    bN[((j-1) * 4 + k + 1),3] <- btrans[2,3,k] * bN[((j-1) * 4 + k),2] + btrans[3,3,k] * bN[((j-1) * 4 + k),3] + exp(I[((j-1) * 4 + k)])
    }
    
    # BNT recruits produced in fall as a function weighted sum of adults (wA) and reprodutive rate (Beta) in winter
    wA[j] <- (bN[((j - 1) * 4 + 2),2] + 4 * bN[((j - 1) * 4 + 2),3])
    beta_eps[j] ~ dnorm(0, tau_beta)
    Beta[j] <- exp(lbeta_0 + beta_eps[j]) 
    
    # between summer and fall all bnt graduate to sz 2 & recruits show up
    bN[(1 + j * 4),1] <- wA[j] * Beta[j]
    bN[(1 + j * 4),2] <- btrans[1,2,4] * bN[(j * 4),1] + btrans[2,2,4] * bN[(j * 4),2]
    bN[(1 + j * 4),3] <- btrans[2,3,4] * bN[(j * 4),2] + btrans[3,3,4] * bN[(j * 4),3] + exp(I[(j * 4)])
    }
    
    # 2000 - 2017 AZGF data
    for(j in 1:NAZsamps){
    for(k in 1:3){
    blpAZ[j,k] ~ dnorm(mu_AZ[k], tau_blp)
    logit(bpAZ[j,k]) <- blpAZ[j,k]
    blamAZ[j,k] <- bpAZ[j,k] * bN[ts[j],k] * AZeff[j] / 35 # predicted catch AZ (35 is ~h to do LF, by AZ)
    bAZ[j,k] ~ dpois(blamAZ[j,k])
    }
    }
    
    # 2012 - 2017 NO: starts in april 2012
    for(j in 1:23){
    for(k in 1:3){
    blamNO[j,k] <- bp[j,k,k] * bpi * bN[(j + 46),k]
    bNOc[j,k] ~ dpois(blamNO[j,k])
    }
    }
    }
    ",fill=TRUE)
sink()
#-----------------------------------------------------------------------------#
BM_JM.data <- list(NAZsamps = NAZsamps, ts = ts, AZeff = AZeff, bAZ = bAZ,
                   seasNO = seasNO, bNOc = bNOc, NOpasses = NOpasses, ones = FR,
                   FR = FR, last = last, bCH = bCH, NCH = NCH, sumf = sumf,
                   spawn = spawn)

BM_JM.par <- c('bphi', 'bpsi1', 'bpsi2', 'mu_blp', 'sd_blp', 'lbeta_0',
               'mu_I', 'I', 'Beta', 'IN', 'AZadj', 'sd_I', 'sd_lphi',
               'sd_beta', 'bN', 'bp_pass')

jags.fit <- jags.parallel(BM_JM.data, inits = NULL, BM_JM.par, "JAGS_Marginalized.jags",
                          n.chains = 3, n.iter = 10)
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Application 4
#  Fitting an Integrated Population Model to Brown Trout Data
#  Marginalized Stan version 
#
#  Notes:
#  * The model runs from the fall of 2000 to fall of 2017 on a seasonal basis
#  * We define three size states based on total length in mm 
#    - 0 - 200; 200 - 350; 350 +
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

# read in data
NO_catch <- read.csv(paste0(getwd(), "/Data/", "NO_catch.csv"), header = TRUE)
AZGF_catch <- read.csv(paste0(getwd(), "/Data/", "AZGF_catch.csv"), header = TRUE)
MR_data <- read.csv(paste0(getwd(), "/Data/", "bCH.csv"), header = FALSE)

#-----------------------------------------------------------------------------#
# extract/reformat data
bNOc <- as.matrix(NO_catch[,1:3])
NOpasses <- NO_catch$NOpasses
seasNO <- NO_catch$seasNO
spawn <- ifelse(seasNO == 1, 4, 3)
bAZ <- as.matrix(AZGF_catch[,1:3])
ts <- AZGF_catch$ts
AZeff <- AZGF_catch$AZeff
NAZsamps <- length(AZeff)

# capture-recapture data
allCH <- MR_data[,1:23]   

bCH = collapse.ch(allCH)[[1]]
FR = collapse.ch(allCH)[[2]]

findlast <- function(x){ifelse(x[23] == 1, 22, max(which(x[1:22] != 4)))}
last <- apply(bCH, 1, findlast)

NCH <- length(last)
findfirst <- function(x){which(x != 4)[1]}
sumf <- apply(bCH, 1, findfirst)

#-----------------------------------------------------------------------------#
sink("Stan_Marginalized.Stan")
cat("
    // Marginalized Integrated Population Model - Brown Trout Data  
    
    data{
    int NAZsamps;                   // Number of samples AGFD
    int ts[NAZsamps];
    vector [NAZsamps] AZeff;
    int bAZ[NAZsamps, 3];
    int seasNO[23];
    int bNOc[23,3];
    int NOpasses[23];              // Number of passes
    int NCH;                       // Number of capture histories
    int FR[NCH];
    int last[NCH];
    int bCH[NCH, 23];
    int sumf[NCH];
    int spawn[23];               // Indicator for spawning month
    }
    
    parameters{
    matrix[3, 4]lphi;
    real<lower = 0.01, upper = 4> sd_lphi;                  
    real<lower = 0, upper = 1> bpsi2;
    vector<lower = 0, upper = 1>[3] bpsi1;
    real<lower = 0.1, upper = 2> sd_blp;
    vector[4] mu_blp;
    matrix[23, 3]blp_pass;
    vector[3] AZadj;
    vector<lower = 0, upper = 1000>[2] IN;
    real<lower = 0.1, upper = 4> sd_beta;                  
    vector[17] beta_eps;
    real<lower = -6, upper = 0> lbeta_0;                                
    real<lower = 0, upper = 6> mu_I;                                    
    real<lower = 0.01, upper = 3> sd_I;
    vector[68] I;
    matrix[NAZsamps, 3]blpAZ;
    }
    
    transformed parameters{
    matrix[3, 4]bphi;  
    real btrans[4, 4, 4];                                                     
    matrix[23,3]bp_pass;
    real bp[23,4,4];                                                          
    vector[3] mu_AZ;
    matrix[69,3] bN;                                                          
    vector[17] wA;
    vector[17] Beta;
    
    for(j in 1:3){  // remove for no prior version
    for(k in 1:4){
    bphi[j,k] = inv_logit(lphi[j,k]);
    }
    }
    
    // define transition matrix that combines survival and growth parameters
    for(i in 1:4){
    btrans[1,3,i] = 0;
    btrans[1,4,i] = 1 - bphi[1,i];
    btrans[2,1,i] = 0;
    btrans[2,2,i] = bphi[2,i] * (1 - bpsi2);
    btrans[2,3,i] = bphi[2,i] * bpsi2;
    btrans[2,4,i] = 1 - bphi[2,i];
    btrans[3,1,i] = 0;
    btrans[3,2,i] = 0;
    btrans[3,3,i] = bphi[3,i];
    btrans[3,4,i] = 1 - bphi[3,i];
    btrans[4,1,i] = 0;
    btrans[4,2,i] = 0;
    btrans[4,3,i] = 0;
    btrans[4,4,i] = 1;
    }
    
    // size class one transitions are done separately because growth is allowed to vary between seasons
    for(i in 1:3){
    btrans[1,1,i] = bphi[1,i] * (1 - bpsi1[i]);
    btrans[1,2,i] = bphi[1,i] * bpsi1[i];
    }
    
    btrans[1,1,4] = 0;
    btrans[1,2,4] = bphi[1,4];
    
    // this loop calculates actual per pass pcaps for each trip and modifies based on # of passes
    for(j in 1:23){
    for(k in 1:3){
    bp_pass[j,k] = inv_logit(blp_pass[j,k]);
    bp[j,k,k] = 1 - pow((1 - bp_pass[j,k]), NOpasses[j]);                 
    bp[j,k,4] = 1 - bp[j,k,k];
    }
    
    bp[j,1,2] = 0;
    bp[j,1,3] = 0;
    bp[j,2,1] = 0;
    bp[j,2,3] = 0;
    bp[j,3,1] = 0;
    bp[j,3,2] = 0;
    bp[j,4,1] = 0;
    bp[j,4,2] = 0;
    bp[j,4,3] = 0;
    bp[j,4,4] = 1;
    }
    
    // calculate offset for each size classes of AZGF effort and calculate expected pcap
    mu_AZ[1] = mu_blp[1] + AZadj[1];
    mu_AZ[2] = mu_blp[2] + AZadj[2];
    mu_AZ[3] = mu_blp[3] + AZadj[3];
    
    bN[1,1] = 0;                                     
    bN[1,2] = IN[1];
    bN[1,3] = IN[2];
    
    //calculate latent abundance of brown trout from fall 2000 to end of 2017
    for(j in 1:17){
    for(k in 1:3){
    bN[((j-1) * 4 + k + 1),1] = btrans[1,1,k] * bN[((j-1) * 4 + k),1];
    bN[((j-1) * 4 + k + 1),2] = btrans[1,2,k] * bN[((j-1) * 4 + k),1] + btrans[2,2,k] * bN[((j-1) * 4 + k),2];
    bN[((j-1) * 4 + k + 1),3] = btrans[2,3,k] * bN[((j-1) * 4 + k),2] + btrans[3,3,k] * bN[((j-1) * 4 + k),3] + exp(I[((j-1) * 4 + k)]);
    }
    
    // BNT recruits produced in fall as a function weighted sum of adults (wA) and reprodutive rate (Beta) in winter
    wA[j] = (bN[((j - 1) * 4 + 2),2] + 4 * bN[((j - 1) * 4 + 2),3]);
    Beta[j] = exp(lbeta_0 + beta_eps[j]);
    
    // between summer and fall all bnt graduate to sz 2 & recruits show up
    bN[(1 + j * 4),1] = wA[j] * Beta[j];
    bN[(1 + j * 4),2] = btrans[1,2,4] * bN[(j * 4),1] + btrans[2,2,4] * bN[(j * 4),2];
    bN[(1 + j * 4),3] = btrans[2,3,4] * bN[(j * 4),2] + btrans[3,3,4] * bN[(j * 4),3] + exp(I[(j * 4)]);
    }
    }
    
    model{
    real pz[NCH, 23, 4];                                    
    matrix[NAZsamps, 3]bpAZ;
    matrix[NAZsamps, 3]blamAZ;
    matrix[23, 3]blamNO;
    vector[4] temp;
    
    ///////////////////////////
    // lphi is the logit of survival and is given a prior based on the Lorenzen
    // function and the average mass of fish in each size class during each season.
    // variation from the priod mode is determined by an estimated variance
    // parameter (sd_lphi)
    
    lphi[1,1] ~ normal(1.08, sd_lphi);
    lphi[1,2] ~ normal(1.14, sd_lphi);
    lphi[1,3] ~ normal(1.26, sd_lphi);
    lphi[1,4] ~ normal(1.38, sd_lphi);
    lphi[2,1] ~ normal(1.96, sd_lphi);
    lphi[2,2] ~ normal(1.96, sd_lphi);
    lphi[2,3] ~ normal(2.02, sd_lphi);
    lphi[2,4] ~ normal(2.02, sd_lphi);
    lphi[3,1] ~ normal(2.29, sd_lphi);
    lphi[3,2] ~ normal(2.29, sd_lphi);
    lphi[3,3] ~ normal(2.29, sd_lphi);
    lphi[3,4] ~ normal(2.29, sd_lphi);
    ///////////////////////////
    
    // mean pcaps per pass on a logit scale for three size classes, plus largest size class during spawning season
    for(i in 1:4){
    mu_blp[i] ~ normal(-3, 2);                              
    }
    
    // (done above in transformed) this loop calculates actual per pass pcaps for each trip and modifies based on # of passes
    for(j in 1:23){
    blp_pass[j,1] ~ normal(mu_blp[1], sd_blp);
    blp_pass[j,2] ~ normal(mu_blp[2], sd_blp);
    blp_pass[j,3] ~ normal(mu_blp[spawn[j]], sd_blp);
    }
    
    for(k in 1:NCH){
    pz[k,sumf[k],1] = (1 == bCH[k,sumf[k]]);
    pz[k,sumf[k],2] = (2 == bCH[k,sumf[k]]);
    pz[k,sumf[k],3] = (3 == bCH[k,sumf[k]]);
    pz[k,sumf[k],4] = 0;
    for(t in sumf[k]:(last[k] - 1)){
    for(i in 1:4){
    for(j in 1:4){
    temp[j] = pz[k,t,j] * btrans[j,i,seasNO[(t + 1)]] * bp[t,i,bCH[k,(t + 1)]];
    }
    pz[k,(t + 1),i] = sum(temp);
    }
    }
    target += FR[k] * log(sum(pz[k,last[k],]));                             
    }
    
    //////////////////////////
    // (done above in transformed) calculate offset for each size classes of AZGF effort and calculate expected pcap
    AZadj[1] ~ normal(0,1);                                       
    AZadj[2] ~ normal(0,1);
    AZadj[3] ~ normal(0,1);
    
    // calculate actual immigration in each interval on log scale
    for(j in 1:68){                                                 
    I[j] ~ normal(mu_I, sd_I);
    }
    
    for(j in 1:17){                                             
    beta_eps[j] ~ normal(0, sd_beta);
    }
    
    //////////////////////////
    // 2000 - 2017 AZGF data
    for(j in 1:3){
    for(k in 2:3){
    blpAZ[j,k] ~ normal(mu_AZ[k], sd_blp);
    bpAZ[j,k] = inv_logit(blpAZ[j,k]);
    blamAZ[j,k] = bpAZ[j,k] * bN[ts[j],k] * AZeff[j] / 35;
    bAZ[j,k] ~ poisson(blamAZ[j,k]);
    }
    }
    
    for(j in 4:NAZsamps){
    for(k in 1:3){
    blpAZ[j,k] ~ normal(mu_AZ[k], sd_blp);
    bpAZ[j,k] = inv_logit(blpAZ[j,k]);
    blamAZ[j,k] = bpAZ[j,k] * bN[ts[j],k] * AZeff[j] / 35;
    bAZ[j,k] ~ poisson(blamAZ[j,k]);
    }
    }
    
    // 2012 - 2017 NO: starts in april 2012
    for(j in 1:23){
    for(k in 1:3){
    blamNO[j,k] = bp[j,k,k] * 0.08 * bN[(j + 46),k];
    bNOc[j,k] ~ poisson(blamNO[j,k]);
    }
    }
    }
    
    ", fill = TRUE)
sink()
#-----------------------------------------------------------------------------#
sm.data <- list(NAZsamps = NAZsamps, ts = ts, AZeff = AZeff, bAZ = bAZ,
                seasNO = seasNO, bNOc = bNOc, NOpasses = NOpasses, ones = FR,
                FR = FR, last = last, bCH = bCH, NCH = NCH, sumf = sumf,
                spawn = spawn)

sm.params = c('bphi', 'bpsi1', 'bpsi2', 'mu_blp', 'sd_blp', 'lbeta_0',
              'mu_I', 'I', 'Beta', 'IN', 'AZadj', 'sd_I', 'sd_lphi',
              'sd_beta', 'bN', 'bp_pass')

# MCMC settings
ni = 1000
nt = 1
nb = 500
nc = 3

SM.c <- stan("Stan_Marginalized.stan",
             data = sm.data,
             pars = sm.params,
             control = list(max_treedepth = 14, adapt_delta = .85),
             chains = nc, iter = ni, thin = nt, seed = 1) 
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Application 4
#  Fitting an Integrated Population Model to Brown Trout Data
#  Marginalized Stan version - Removes prior on survival & random effects on
#  detection
#
#  Notes:
#  * The model runs from the fall of 2000 to fall of 2017 on a seasonal basis
#  * We define three size states based on total length in mm 
#    - 0 - 200; 200 - 350; 350 +
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

# read in data
NO_catch <- read.csv(paste0(getwd(), "/Data/", "NO_catch.csv"), header = TRUE)
AZGF_catch <- read.csv(paste0(getwd(), "/Data/", "AZGF_catch.csv"), header = TRUE)
MR_data <- read.csv(paste0(getwd(), "/Data/", "bCH.csv"), header = FALSE)

#-----------------------------------------------------------------------------#
# extract/reformat data
bNOc <- as.matrix(NO_catch[,1:3])
NOpasses <- NO_catch$NOpasses
seasNO <- NO_catch$seasNO
spawn <- ifelse(seasNO == 1, 4, 3)
bAZ <- as.matrix(AZGF_catch[,1:3])
ts <- AZGF_catch$ts
AZeff <- AZGF_catch$AZeff
NAZsamps <- length(AZeff)

# capture-recapture data
allCH <- MR_data[,1:23]

bCH = collapse.ch(allCH)[[1]]
FR = collapse.ch(allCH)[[2]]

findlast <- function(x){ifelse(x[23] == 1, 22, max(which(x[1:22] != 4)))}
last <- apply(bCH, 1, findlast)

NCH <- length(last)
findfirst <- function(x){which(x != 4)[1]}
sumf <- apply(bCH, 1, findfirst)

#-----------------------------------------------------------------------------#
sink("nop_nore.stan")
cat("
    data{
    int NAZsamps;                            // Number of samples AGFD
    int ts[NAZsamps];
    vector [NAZsamps] AZeff;
    int bAZ [NAZsamps, 3];
    int seasNO[23];
    int bNOc[23,3];
    int NOpasses[23];
    int NCH;                       // Number of capture histories
    int FR[NCH];
    int last[NCH];
    int bCH[NCH, 23];
    int sumf[NCH];
    int spawn[23];               // Indicator for spawning month
    }
    
    parameters{
    matrix<lower = 0, upper = 1>[3, 4]bphi;
    real<lower = 0, upper = 1> bpsi2;
    vector<lower = 0, upper = 1>[3] bpsi1;
    vector<lower = 0, upper = 1>[4] mu_blp;
    vector<lower = 0, upper = 1>[3] mu_AZ;
    vector<lower = 0, upper = 1000>[2] IN;                              
    vector<lower = -4, upper = 8>[68] I;
    vector<lower = 0, upper = 4>[17] Beta;
    }
    
    transformed parameters{
    real btrans[4, 4, 4];                                                     
    matrix[23,3]bp_pass;
    real bp[23,4,4];                                                          
    matrix[69,3] bN;                                                          
    vector[17] wA;
    
    
    for(j in 1:23){
    bp_pass[j,1]=mu_blp[1];
    bp_pass[j,2]=mu_blp[2];
    bp_pass[j,3]=mu_blp[spawn[j]];
    }
    
    
    // define transition matrix that combines survival and growth parameters
    for(i in 1:4){
    btrans[1,3,i] = 0;
    btrans[1,4,i] = 1 - bphi[1,i];
    btrans[2,1,i] = 0;
    btrans[2,2,i] = bphi[2,i] * (1 - bpsi2);
    btrans[2,3,i] = bphi[2,i] * bpsi2;
    btrans[2,4,i] = 1 - bphi[2,i];
    btrans[3,1,i] = 0;
    btrans[3,2,i] = 0;
    btrans[3,3,i] = bphi[3,i];
    btrans[3,4,i] = 1 - bphi[3,i];
    btrans[4,1,i] = 0;
    btrans[4,2,i] = 0;
    btrans[4,3,i] = 0;
    btrans[4,4,i] = 1;
    }
    
    // size class one transitions are done separately because growth is allowed to vary between seasons
    for(i in 1:3){
    btrans[1,1,i] = bphi[1,i] * (1 - bpsi1[i]);
    btrans[1,2,i] = bphi[1,i] * bpsi1[i];
    }
    
    btrans[1,1,4] = 0;
    btrans[1,2,4] = bphi[1,4];
    
    // this loop calculates actual per pass pcaps for each trip and modifies based on # of passes
    for(j in 1:23){
    for(k in 1:3){
    bp[j,k,k] = 1 - pow((1 - bp_pass[j,k]), NOpasses[j]);                 
    bp[j,k,4] = 1 - bp[j,k,k];
    }
    
    bp[j,1,2] = 0;
    bp[j,1,3] = 0;
    bp[j,2,1] = 0;
    bp[j,2,3] = 0;
    bp[j,3,1] = 0;
    bp[j,3,2] = 0;
    bp[j,4,1] = 0;
    bp[j,4,2] = 0;
    bp[j,4,3] = 0;
    bp[j,4,4] = 1;
    }
    
    // calculate offset for each size classes of AZGF effort and calculate expected pcap
    
    bN[1,1] = 0;                                     
    bN[1,2] = IN[1];
    bN[1,3] = IN[2];
    
    //calculate latent abundance of brown trout from fall 2000 to end of 2017
    for(j in 1:17){
    for(k in 1:3){
    bN[((j-1) * 4 + k + 1),1] = btrans[1,1,k] * bN[((j-1) * 4 + k),1];
    bN[((j-1) * 4 + k + 1),2] = btrans[1,2,k] * bN[((j-1) * 4 + k),1] + btrans[2,2,k] * bN[((j-1) * 4 + k),2];
    bN[((j-1) * 4 + k + 1),3] = btrans[2,3,k] * bN[((j-1) * 4 + k),2] + btrans[3,3,k] * bN[((j-1) * 4 + k),3] + exp(I[((j-1) * 4 + k)]);
    }
    
    // BNT recruits produced in fall as a function weighted sum of adults (wA) and reprodutive rate (Beta) in winter
    wA[j] = (bN[((j - 1) * 4 + 2),2] + 4 * bN[((j - 1) * 4 + 2),3]);
    
    // between summer and fall all bnt graduate to sz 2 & recruits show up
    bN[(1 + j * 4),1] = wA[j] * Beta[j];
    bN[(1 + j * 4),2] = btrans[1,2,4] * bN[(j * 4),1] + btrans[2,2,4] * bN[(j * 4),2];
    bN[(1 + j * 4),3] = btrans[2,3,4] * bN[(j * 4),2] + btrans[3,3,4] * bN[(j * 4),3] + exp(I[(j * 4)]);
    }
    }
    
    model{
    real pz[NCH, 23, 4];                                    
    matrix[NAZsamps, 3]blamAZ;
    matrix[23, 3]blamNO;
    vector[4] temp;
    
    
    // (done above in transformed) this loop calculates actual per pass pcaps for each trip and modifies based on # of passes
    
    for(k in 1:NCH){
    pz[k,sumf[k],1] = (1 == bCH[k,sumf[k]]);
    pz[k,sumf[k],2] = (2 == bCH[k,sumf[k]]);
    pz[k,sumf[k],3] = (3 == bCH[k,sumf[k]]);
    pz[k,sumf[k],4] = 0;
    for(t in sumf[k]:(last[k] - 1)){
    for(i in 1:4){
    for(j in 1:4){
    temp[j] = pz[k,t,j] * btrans[j,i,seasNO[(t + 1)]] * bp[t,i,bCH[k,(t + 1)]];
    }
    pz[k,(t + 1),i] = sum(temp);
    }
    }
    target += FR[k] * log(sum(pz[k,last[k],]));                             
    }
    
    //////////////////////////
    
    // calculate actual immigration in each interval on log scale
    
    //////////////////////////
    // 2000 - 2017 AZGF data
    for(j in 1:3){
    for(k in 2:3){
    blamAZ[j,k] = mu_AZ[k] * bN[ts[j],k] * AZeff[j] / 35;
    bAZ[j,k] ~ poisson(blamAZ[j,k]);
    }
    }
    
    
    for(j in 4:NAZsamps){
    for(k in 1:3){
    blamAZ[j,k] = mu_AZ[k] * bN[ts[j],k] * AZeff[j] / 35;
    bAZ[j,k] ~ poisson(blamAZ[j,k]);
    }
    }
    
    
    // 2012 - 2017 NO: starts in april 2012
    for(j in 1:23){
    for(k in 1:3){
    blamNO[j,k] = bp[j,k,k] * 0.08 * bN[(j + 46),k];
    bNOc[j,k] ~ poisson(blamNO[j,k]);
    }
    }
    }
    
    ", fill = TRUE)
sink()
#-----------------------------------------------------------------------------#
sm.data <- list(NAZsamps = NAZsamps, ts = ts, AZeff = AZeff, bAZ = bAZ,
                seasNO = seasNO, bNOc = bNOc, NOpasses = NOpasses, NCH = NCH,
                FR = FR, last = last, bCH = bCH, sumf = sumf, spawn = spawn)

sm.params <- c('bphi', 'bpsi1', 'bpsi2', "I", "Beta", "IN", "bN", "mu_AZ", "mu_blp")

# MCMC settings
ni = 10
nt = 1
nb = 5
nc = 1

# Call Stan from R 
SM.nop.nore <- stan("nop_nore.stan",
                    data = sm.data,
                    pars = sm.params,
                    control = list(max_treedepth = 14, adapt_delta = .85),
                    chains = nc, iter = ni, thin = nt, seed = 1) 

#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Simulation Code for Results in Table 2 
#
#  Notes:
#  * 
#
###############################################################################
library(R2jags)

# jags code - discrete approach
sink("JD.jags")
cat("
    model{
    Z[1] ~ dbern(psi0)
    Z[2] ~ dbern(Z[1] * (1 - eps) + (1 - Z[1]) * gam)
    Z[3] ~ dbern(Z[2] * (1 - eps) + (1 - Z[2]) * gam)
    for(i in 1:3){
    p_eff[i] <- Z[i] * p
    for(j in 1:2){
    Y[(2 * (i - 1) + j)] ~ dbern(p_eff[i])		
    }
    } 
    }
    ",fill = TRUE)
sink()
#-----------------------------------------------------------------------------#
# run comparisons using equations in ms for example
compare <- function(par, sd = 1){
  set.seed(sd)
  psi0 <- par[1]
  gam <- par[2]
  eps <- par[3]
  p <- par[4]
  out <- numeric()
  out[1] <- psi0 * (1 - eps) + (1 - psi0) * gam
  out[2] <- (out[1] * (1 - p) ^ 2) / (out[1] * (1 - p) ^ 2 + (1 - out[1]))
  out[3] <- ((1 - eps) * ((1 - p) ^ 2)) / ((1 - eps) * ((1 - p) ^ 2) + eps)
  out[4] <- ((1 - eps) * ((1 - p) ^ 2) * (1 - eps) * p ^ 2) / 
    ((1 - eps) * ((1 - p) ^ 2) * (1 - eps) * (p ^ 2) + eps * gam * p ^ 2)
  Y <- c(0, 1, 0, 0, 1, 1)
  data <- list(Y = Y, p = p, psi0 = psi0, gam = gam, eps = eps)
  params <- c("Z")
  JD_inits <- function() list(Z = c(1, rbinom(1, 1, .5), 1))
  fit_JD <- jags(data, JD_inits, params, "JD.jags", n.chains = 100, n.iter = 1000)
  out[5] <- fit_JD$BUGSoutput$summary[2,1]
  return(out)
}

input <- cbind(c(.75, .75, .75, .75), c(.6, .6, .1, .6),
               c(.2, .8, .2, .2), c(0.1, 0.1, 0.1, 0.6))
output <- matrix(NA, ncol = 5, nrow = 4)
for(j in 1:4){
  output[j,] <- compare(input[j,])
}
# write.csv(cbind(input, output), "sim0.csv")
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
###############################################################################
#                                                                     Spring 19
#  Application 5
#  Dynamic two-species occupancy model with autologistic effects
#  (Barred owls and Northern Spotted owls)
#  JAGS Discrete Version
#
#  Notes:
#  * Need to set directory for data
#
###############################################################################
library(R2jags)

#Get data:
owls.dat<-read.csv(".//owls.dat.csv")


#-----------------------------------------------------------------------------#
# format data for model fitting:
nmaxVisit<-8
nYear<-22
pX<-array(0,dim=c(max(owls.dat$Site),nmaxVisit*nYear,7))
Y_BO<-array(0,dim=c(max(owls.dat$Site),nmaxVisit*nYear))
Y_NSO<-array(0,dim=c(max(owls.dat$Site),nmaxVisit*nYear))

#pX : Matrix of detection covariates
#First dimension: Sites
#Second dimension: Number of total visits across years (maximum number of visits per year is 8)
#Third dimension: Covariates
#pX[,,1] = (1/0); whether site was visited 
#pX[,,2] = (1/0); 1 if visited during the daytime
#pX[,,3] = (1/0); 1 if visited during the nightime
#Note if pX[,,2]==0 & pX[,,3]==0 then visited during crepuscular time of day
#pX[,,4] = (1/0); 1 if sampled using Method 1
#pX[,,5] = (1/0); 1 if sampled using Method 2
#Note if pX[,,4]==0 & pX[,,5]==0 then sampled using Method 3
#pX[,,6] = (1/0); Whether visit occurred during the last half of the study
#pX[,,7] = (continuous, standardized); Covariate describing time spent sampling


for(i in 1:length(owls.dat[,1])){
  pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],1]<-owls.dat$Visited[i]
  if(owls.dat$TOD[i]==1){
    pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],2]<-1}
  if(owls.dat$TOD[i]==2){
    pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],3]<-1}
  if(owls.dat$Method[i]==1){
    pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],4]<-1}
  if(owls.dat$Method[i]==2){
    pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],5]<-1}
  pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],7]<-owls.dat$ttot[i]
  Y_BO[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i]]<-owls.dat$BO[i]
  Y_NSO[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i]]<-owls.dat$NSO[i]}

#Dummy variable for last half of the study:
pX[,89:176,6]<-1

#If site is not visited, change all covariates to zero:
for(i in 1:158){for(j in 1:176){ if(pX[i,j,1]==0){pX[i,j,]<-0}}}  

#Extract riparian forest (RF) and forest (FOR) covariates from owl.dat:
RF<-matrix(NA, nrow=dim(pX)[1], ncol=nYear-1)
FOR<-matrix(NA, nrow=dim(pX)[1], ncol=nYear-1)

for(i in 1:dim(RF)[1]){
  for(j in 1:dim(RF)[2]){
    RF[i,j]=unique(owls.dat$RF[which(owls.dat$Year==j & owls.dat$Site==i)])
    FOR[i,j]=unique(owls.dat$FOR[which(owls.dat$Year==j & owls.dat$Site==i)])
  }}


#Combine Y_BO and Y_NSO into a single matrix with four states:
#State 1 = no owls detected
#State 2 = only NSO detected
#State 3 = only BO detected
#State 4 = both BO and NSO detected
Y_BO_temp<-Y_BO
2->Y_BO_temp[Y_BO_temp==1]
Y<-Y_BO_temp+Y_NSO+1

#Matrix 'BO' is whether or not barred owls were detected at each site at least once in each year
BO<-numeric()
for(i in 1:22){BO<-cbind(BO,apply(Y_BO[,(i*8-7):(i*8)],1,sum))}
1->BO[BO>1]
#Matrix 'NSO' is whether or not spotted owls were detected at each site at least once in each year
NSO<-numeric()
for(i in 1:22){NSO<-cbind(NSO,apply(Y_NSO[,(i*8-7):(i*8)],1,sum))}
1->NSO[NSO>1]

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
#  Application 5
#  Dynamic two-species occupancy model with autologistic effects
#  (Barred owls and Northern Spotted owls)
#  JAGS Marginalized Version
#
#  Notes:
#  * Need to set directory for data
#
###############################################################################
library(R2jags)

# Get data
owls.dat<-read.csv(".//owls.dat.csv")

#-----------------------------------------------------------------------------#
# format data for model fitting:
nmaxVisit<-8
nYear<-22
pX<-array(0,dim=c(max(owls.dat$Site),nmaxVisit*nYear,7))
Y_BO<-array(0,dim=c(max(owls.dat$Site),nmaxVisit*nYear))
Y_NSO<-array(0,dim=c(max(owls.dat$Site),nmaxVisit*nYear))

#pX : Matrix of detection covariates
#First dimension: Sites
#Second dimension: Number of total visits across years (maximum number of visits per year is 8)
#Third dimension: Covariates
#pX[,,1] = (1/0); whether site was visited 
#pX[,,2] = (1/0); 1 if visited during the daytime
#pX[,,3] = (1/0); 1 if visited during the nightime
#Note if pX[,,2]==0 & pX[,,3]==0 then visited during crepuscular time of day
#pX[,,4] = (1/0); 1 if sampled using Method 1
#pX[,,5] = (1/0); 1 if sampled using Method 2
#Note if pX[,,4]==0 & pX[,,5]==0 then sampled using Method 3
#pX[,,6] = (1/0); Whether visit occurred during the last half of the study
#pX[,,7] = (continuous, standardized); Covariate describing time spent sampling


for(i in 1:length(owls.dat[,1])){
  pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],1]<-owls.dat$Visited[i]
  if(owls.dat$TOD[i]==1){
    pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],2]<-1}
  if(owls.dat$TOD[i]==2){
    pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],3]<-1}
  if(owls.dat$Method[i]==1){
    pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],4]<-1}
  if(owls.dat$Method[i]==2){
    pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],5]<-1}
  pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],7]<-owls.dat$ttot[i]
  Y_BO[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i]]<-owls.dat$BO[i]
  Y_NSO[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i]]<-owls.dat$NSO[i]}

#Dummy variable for last half of the study:
pX[,89:176,6]<-1

#If site is not visited, change all covariates to zero:
for(i in 1:158){for(j in 1:176){ if(pX[i,j,1]==0){pX[i,j,]<-0}}}  

#Extract riparian forest (RF) and forest (FOR) covariates from owl.dat:
RF<-matrix(NA, nrow=dim(pX)[1], ncol=nYear-1)
FOR<-matrix(NA, nrow=dim(pX)[1], ncol=nYear-1)

for(i in 1:dim(RF)[1]){
  for(j in 1:dim(RF)[2]){
    RF[i,j]=unique(owls.dat$RF[which(owls.dat$Year==j & owls.dat$Site==i)])
    FOR[i,j]=unique(owls.dat$FOR[which(owls.dat$Year==j & owls.dat$Site==i)])
  }}


#Combine Y_BO and Y_NSO into a single matrix with four states:
#State 1 = no owls detected
#State 2 = only NSO detected
#State 3 = only BO detected
#State 4 = both BO and NSO detected
Y_BO_temp<-Y_BO
2->Y_BO_temp[Y_BO_temp==1]
Y<-Y_BO_temp+Y_NSO+1

#Matrix 'BO' is whether or not barred owls were detected at each site at least once in each year
BO<-numeric()
for(i in 1:22){BO<-cbind(BO,apply(Y_BO[,(i*8-7):(i*8)],1,sum))}
1->BO[BO>1]
#Matrix 'NSO' is whether or not spotted owls were detected at each site at least once in each year
NSO<-numeric()
for(i in 1:22){NSO<-cbind(NSO,apply(Y_NSO[,(i*8-7):(i*8)],1,sum))}
1->NSO[NSO>1]

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
    ",fill = TRUE)
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
#                                                                     Spring 19
#  Application 5
#  Dynamic two-species occupancy model with autologistic effects
#  (Barred owls and Northern Spotted owls)
#  Stan Marginalized Version
#
#  Notes:
#  * Need to set directory for data
#
###############################################################################
library(rstan)

# Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) # to run Stan in parallel

#Get data:
owls.dat<-read.csv(".//owls.dat.csv")


#-----------------------------------------------------------------------------#
# format data for model fitting:
nmaxVisit<-8
nYear<-22
pX<-array(0,dim=c(max(owls.dat$Site),nmaxVisit*nYear,7))
Y_BO<-array(0,dim=c(max(owls.dat$Site),nmaxVisit*nYear))
Y_NSO<-array(0,dim=c(max(owls.dat$Site),nmaxVisit*nYear))

#pX : Matrix of detection covariates
#First dimension: Sites
#Second dimension: Number of total visits across years (maximum number of visits per year is 8)
#Third dimension: Covariates
#pX[,,1] = (1/0); whether site was visited 
#pX[,,2] = (1/0); 1 if visited during the daytime
#pX[,,3] = (1/0); 1 if visited during the nightime
#Note if pX[,,2]==0 & pX[,,3]==0 then visited during crepuscular time of day
#pX[,,4] = (1/0); 1 if sampled using Method 1
#pX[,,5] = (1/0); 1 if sampled using Method 2
#Note if pX[,,4]==0 & pX[,,5]==0 then sampled using Method 3
#pX[,,6] = (1/0); Whether visit occurred during the last half of the study
#pX[,,7] = (continuous, standardized); Covariate describing time spent sampling


for(i in 1:length(owls.dat[,1])){
  pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],1]<-owls.dat$Visited[i]
  if(owls.dat$TOD[i]==1){
    pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],2]<-1}
  if(owls.dat$TOD[i]==2){
    pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],3]<-1}
  if(owls.dat$Method[i]==1){
    pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],4]<-1}
  if(owls.dat$Method[i]==2){
    pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],5]<-1}
  pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],7]<-owls.dat$ttot[i]
  Y_BO[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i]]<-owls.dat$BO[i]
  Y_NSO[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i]]<-owls.dat$NSO[i]}

#Dummy variable for last half of the study:
pX[,89:176,6]<-1

#If site is not visited, change all covariates to zero:
for(i in 1:158){for(j in 1:176){ if(pX[i,j,1]==0){pX[i,j,]<-0}}}  

#Extract riparian forest (RF) and forest (FOR) covariates from owl.dat:
RF<-matrix(NA, nrow=dim(pX)[1], ncol=nYear-1)
FOR<-matrix(NA, nrow=dim(pX)[1], ncol=nYear-1)

for(i in 1:dim(RF)[1]){
  for(j in 1:dim(RF)[2]){
    RF[i,j]=unique(owls.dat$RF[which(owls.dat$Year==j & owls.dat$Site==i)])
    FOR[i,j]=unique(owls.dat$FOR[which(owls.dat$Year==j & owls.dat$Site==i)])
  }}


#Combine Y_BO and Y_NSO into a single matrix with four states:
#State 1 = no owls detected
#State 2 = only NSO detected
#State 3 = only BO detected
#State 4 = both BO and NSO detected
Y_BO_temp<-Y_BO
2->Y_BO_temp[Y_BO_temp==1]
Y<-Y_BO_temp+Y_NSO+1

#Matrix 'BO' is whether or not barred owls were detected at each site at least once in each year
BO<-numeric()
for(i in 1:22){BO<-cbind(BO,apply(Y_BO[,(i*8-7):(i*8)],1,sum))}
1->BO[BO>1]
#Matrix 'NSO' is whether or not spotted owls were detected at each site at least once in each year
NSO<-numeric()
for(i in 1:22){NSO<-cbind(NSO,apply(Y_NSO[,(i*8-7):(i*8)],1,sum))}
1->NSO[NSO>1]

#For Stan code:
n_occ<-apply(pX[,1:8,1],1,sum)
start<-1:22*8-7
for(i in 2:length(start)){n_occ<-cbind(n_occ,apply(pX[,start[i]:(start[i]+7),1],1,sum))}

#-----------------------------------------------------------------------------#
# Note: '//' is for comments in Stan
sink("Stan_Marginalized_owls.stan")
cat("
    data{
    int<lower = 1, upper = 4> Y[158,176];
    real RF [158,21];
    real FOR [158,21];
    matrix[176,7] pX[158];
    }
    
    parameters{
    real<lower = 0, upper = 1> psi0_BO;           // initial barred owl (BO) occupancy
    real<lower = 0, upper = 1> psi0_NSO;          // initial northern spotted owl (NSO) occupancy
    real<lower = -5, upper = 5> gamNSO_int [21];  // time varying NSO colonization intercept
    real<lower = -5, upper = 5> gamNSO_b;         // Forest effect on NSO colonization
    real<lower = -5, upper = 5> epsNSO_int;       // constant NSO intercept on extinction
    real<lower = -5, upper = 5> epsNSO_b;         // BO effect on NSO extinction
    real<lower = -5, upper = 5> gamBO_int;        // constant BO intercept on colonization
    real<lower = -5, upper = 5> gamBO_RF;         // Riparian forest effect, 
    real<lower = -15, upper = 15> gamBO_autolog;  // autologistic effect, 
    real<lower = -5, upper = 5> gamBO_xNSO;       // and NSO effect on BO colonization
    real<lower = -5, upper = 5> epsBO_int;        // constant BO intercept on extinction
    real<lower = -5, upper = 5> epsBO_RF;         // Riparian forest effect, 
    real<lower = -15, upper = 15> epsBO_autolog;  //autologistic effect, 
    real<lower = -5, upper = 5> epsBO_xNSO;       // and NSO effect on BO extinction
    vector<lower = -5, upper = 5>[6] pNSO_b;      // intercept, day, night, method1, method3 and BO effects on NSO detection
    vector<lower = -5, upper = 5>[7] pBO_b;       // intercept, day, night, method1, method3, second half, and survey length effects on BO detection
    }
    
    transformed parameters{
    simplex [4] unc_pz; // initial probability of being in any state
    vector<lower = 0, upper = 1>[176] pBO;
    vector<lower = 0, upper = 1>[176] pNSO_BO;
    vector<lower = 0, upper = 1>[176] pNSO_bo;
    real<lower = 0, upper = 1> gamNSO;
    real<lower = 0, upper = 1> epsNSO_BO;
    real<lower = 0, upper = 1> epsNSO_bo;
    real<lower = 0, upper = 1> gamBO_NSO;
    real<lower = 0, upper = 1> gamBO_nso;
    real<lower = 0, upper = 1> epsBO_NSO;
    real<lower = 0, upper = 1> epsBO_nso;
    simplex [4] p [176,4];
    real<lower = 0, upper = 1> cp [158,22,4]; //cumulative probability of observation history (8 occasions a year for 22 years for a total of 176 occasions)
    simplex [4] tr [4];
    real<lower = 0, upper = 1> pz [158,22,4];
    real<lower = 0, upper = 1> BO [158];
    real<lower = -0.5, upper = 0.5> st_psiBO;
    vector[4] temp;
    
    unc_pz[1] = (1 - psi0_BO) * (1 - psi0_NSO); //unc_pz  =  probability of neither species being present in the study are at the start of the study.
    unc_pz[2] = (1 - psi0_BO) * psi0_NSO;
    unc_pz[3] = psi0_BO * (1 - psi0_NSO);
    unc_pz[4] = psi0_BO * psi0_NSO;
    
    for(k in 1:158){
    pBO = pX[k,,1] .* inv_logit(pX[k,,] * pBO_b);
    pNSO_BO = pX[k,,1] .* inv_logit(pX[k,,1:5] * pNSO_b[1:5] + pNSO_b[6]); // NSO detection probability if BO is present
    pNSO_bo = pX[k,,1] .* inv_logit(pX[k,,1:5] * pNSO_b[1:5]); // NSO detection probability if BO is absent
    for(j in 1:176){
    p[j,1,1] = 1; //if two species absent, then of course state of detection will be 1 (i.e., no species will be detected)
    p[j,1,2] = 0;
    p[j,1,3] = 0;
    p[j,1,4] = 0;
    p[j,2,1] = 1 - pNSO_bo[j]; // if NSO present, probability of it not being detected
    p[j,2,2] = pNSO_bo[j];     // if NSO present, probability of it being detected
    p[j,2,3] = 0;
    p[j,2,4] = 0;
    p[j,3,1] = 1 - pBO[j];     // if BO present, probability of it not being detected
    p[j,3,2] = 0;
    p[j,3,3] = pBO[j];         // if BO present, probability of it being detected
    p[j,3,4] = 0;
    p[j,4,1] = 1 - (pBO[j] + pNSO_BO[j] - pBO[j] * pNSO_BO[j]); // I probably would have written this as : (1-pBO)*(1-pNSO_BO)*visited
    p[j,4,2] = pNSO_BO[j] * (1 - pBO[j]); // if both species present, probability of only detecting BO
    p[j,4,3] = (1 - pNSO_BO[j]) * pBO[j]; // if both species present, probability of only detecting NSO
    p[j,4,4] = pNSO_BO[j] * pBO[j]; // probability of detecting both species
    }
    
    // Note: 22 years and 8 occasions per year?  So, this is looking at the first occasion in each year (i.e., t =  1,9,17, etc)
    for (t in 1:22){
    cp[k,t,1] = p[(1+(t-1)*8),1,Y[k,(1+(t-1)*8)]]; //what is the probability that site is in state 1-4 given state is observed in state Y[k,,] 
    cp[k,t,2] = p[(1+(t-1)*8),2,Y[k,(1+(t-1)*8)]];
    cp[k,t,3] = p[(1+(t-1)*8),3,Y[k,(1+(t-1)*8)]];
    cp[k,t,4] = p[(1+(t-1)*8),4,Y[k,(1+(t-1)*8)]];
    for (j in 2:8){
    cp[k,t,1] = cp[k,t,1]*p[(j+(t-1)*8),1,Y[k,(j+(t-1)*8)]]; //update probabilities given subsequent visits within the same year 
    cp[k,t,2] = cp[k,t,2]*p[(j+(t-1)*8),2,Y[k,(j+(t-1)*8)]];
    cp[k,t,3] = cp[k,t,3]*p[(j+(t-1)*8),3,Y[k,(j+(t-1)*8)]];
    cp[k,t,4] = cp[k,t,4]*p[(j+(t-1)*8),4,Y[k,(j+(t-1)*8)]];
    }
    }
    }
    
    for(k in 1:158){
    pz[k,1,1] = unc_pz[1] * cp[k,1,1]; //for first year, update probability of not being occupied by either species (kinda like a prior) by the observation history in the first year
    pz[k,1,2] = unc_pz[2] * cp[k,1,2];
    pz[k,1,3] = unc_pz[3] * cp[k,1,3];
    pz[k,1,4] = unc_pz[4] * cp[k,1,4];
    BO[k] = sum(pz[k,1,3:4]) / sum(pz[k,1,]); //this is the probability of a BO being present on the site
    }
    
    epsNSO_BO = inv_logit(epsNSO_int+epsNSO_b); //probability of extinction of NSO given BO present
    epsNSO_bo = inv_logit(epsNSO_int); //probability of extinction of NSO given BO absent
    
    for (t in 1:21){
    st_psiBO = (sum(BO) / 158 - .5);//mean probability of BO presence (subtract 0.5 to center it)
    for (k in 1:158){
    gamNSO = inv_logit(gamNSO_int[t] + gamNSO_b * FOR[k,t]); //probability of colonization by NSO is a function of forest cover?
    gamBO_NSO = inv_logit(gamBO_int + gamBO_RF * RF[k,t] + gamBO_autolog * st_psiBO + gamBO_xNSO); //probability of BO colonization given NSO present (think this is a function of riparian forest, and autologistic effects)
    gamBO_nso = inv_logit(gamBO_int + gamBO_RF * RF[k,t] + gamBO_autolog * st_psiBO);
    epsBO_NSO = inv_logit(epsBO_int + epsBO_RF * RF[k,t] + epsBO_autolog * st_psiBO + epsBO_xNSO);
    epsBO_nso = inv_logit(epsBO_int + epsBO_RF * RF[k,t] + epsBO_autolog*st_psiBO);
    tr[1,1] = (1 - gamNSO) * (1 - gamBO_nso);
    tr[1,2] = gamNSO * (1 - gamBO_nso);
    tr[1,3] = (1 - gamNSO) * gamBO_nso;
    tr[1,4] = gamNSO * gamBO_nso;
    tr[2,1] = epsNSO_bo * (1 - gamBO_NSO);
    tr[2,2] = (1 - epsNSO_bo) * (1 - gamBO_NSO);
    tr[2,3] = epsNSO_bo * gamBO_NSO;
    tr[2,4] = (1 - epsNSO_bo) * gamBO_NSO;
    tr[3,1] = (1 - gamNSO) * epsBO_nso;
    tr[3,2] = gamNSO * epsBO_nso;
    tr[3,3] = (1 - gamNSO) * (1 - epsBO_nso);
    tr[3,4] = gamNSO * (1 - epsBO_nso);
    tr[4,1] = epsNSO_BO * epsBO_NSO;
    tr[4,2] = (1 - epsNSO_BO) * epsBO_NSO;
    tr[4,3] = epsNSO_BO * (1 - epsBO_NSO);
    tr[4,4] = (1 - epsNSO_BO) * (1 - epsBO_NSO);
    
    for(j in 1:4){
    temp[j] = pz[k,t,j] * tr[j,1];
    }
    pz[k,(t+1),1] = sum(temp) * cp[k,(t+1),1];
    for(j in 1:4){
    temp[j] = pz[k,t,j] * tr[j,2];
    }
    pz[k,(t+1),2] = sum(temp) * cp[k,(t+1),2];
    for(j in 1:4){
    temp[j] = pz[k,t,j] * tr[j,3];
    }
    pz[k,(t+1),3] = sum(temp) * cp[k,(t+1),3];
    for(j in 1:4){
    temp[j] = pz[k,t,j] * tr[j,4];
    }
    pz[k,(t+1),4] = sum(temp) * cp[k,(t+1),4];
    BO[k] = sum(pz[k,(t+1),3:4]) / sum(pz[k,(t+1),]);
    }
    }
    }
    
    model{
    for(k in 1:158) {  
    target += log(sum(pz[k,22,])); 
    }
    }
    ",fill=TRUE)
sink()  

#-----------------------------------------------------------------------------#
owls_SM.data <- list(Y = array(Y, dim = c(dim(Y)[1], dim(Y)[2])),
                     RF = array(as.vector(unlist(RF)), dim = c(dim(RF)[1], dim(RF)[2])),
                     FOR = array(as.vector(unlist(FOR)), dim = c(dim(FOR)[1], dim(FOR)[2])),
                     pX = array(pX, dim = c(dim(pX)[1], dim(pX)[2], dim(pX)[3])))

owls_SM.par <- c('pBO_b', 'pNSO_b', 'psi0_BO', 'psi0_NSO', 'gamNSO_b', 'epsNSO_b',
                 'epsNSO_int', 'gamBO_int', 'epsBO_int', 'gamBO_RF', 'gamBO_autolog', 
                 'gamBO_xNSO', 'epsBO_RF', 'epsBO_autolog', 'epsBO_xNSO', 'gamNSO_int')

owls_SM.out <- stan(".\\Stan_Marginalized_owls.stan", data = owls_SM.data,
                    pars = owls_SM.par, chains = 1, iter = 10, seed = 1)

#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Application 5
#  Dynamic two-species occupancy model with autologistic effects
#  (Barred owls and Northern Spotted owls)
#  Stan Marginalized Version with random effects and multi-level R^2 calculation
#
#  Notes:
#  * Need to set directory for data
#
#############################################################################
library(rstan)

# Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) # to run Stan in parallel

#Get data:
owls.dat<-read.csv(".//owls.dat.csv")


#-----------------------------------------------------------------------------#
# format data for model fitting:
nmaxVisit<-8
nYear<-22
pX<-array(0,dim=c(max(owls.dat$Site),nmaxVisit*nYear,7))
Y_BO<-array(0,dim=c(max(owls.dat$Site),nmaxVisit*nYear))
Y_NSO<-array(0,dim=c(max(owls.dat$Site),nmaxVisit*nYear))

#pX : Matrix of detection covariates
#First dimension: Sites
#Second dimension: Number of total visits across years (maximum number of visits per year is 8)
#Third dimension: Covariates
#pX[,,1] = (1/0); whether site was visited 
#pX[,,2] = (1/0); 1 if visited during the daytime
#pX[,,3] = (1/0); 1 if visited during the nightime
#Note if pX[,,2]==0 & pX[,,3]==0 then visited during crepuscular time of day
#pX[,,4] = (1/0); 1 if sampled using Method 1
#pX[,,5] = (1/0); 1 if sampled using Method 2
#Note if pX[,,4]==0 & pX[,,5]==0 then sampled using Method 3
#pX[,,6] = (1/0); Whether visit occurred during the last half of the study
#pX[,,7] = (continuous, standardized); Covariate describing time spent sampling


for(i in 1:length(owls.dat[,1])){
  pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],1]<-owls.dat$Visited[i]
  if(owls.dat$TOD[i]==1){
    pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],2]<-1}
  if(owls.dat$TOD[i]==2){
    pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],3]<-1}
  if(owls.dat$Method[i]==1){
    pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],4]<-1}
  if(owls.dat$Method[i]==2){
    pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],5]<-1}
  pX[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i],7]<-owls.dat$ttot[i]
  Y_BO[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i]]<-owls.dat$BO[i]
  Y_NSO[owls.dat$Site[i],nmaxVisit*(owls.dat$Year[i]-1)+owls.dat$Visit_num[i]]<-owls.dat$NSO[i]}

#Dummy variable for last half of the study:
pX[,89:176,6]<-1

#If site is not visited, change all covariates to zero:
for(i in 1:158){for(j in 1:176){ if(pX[i,j,1]==0){pX[i,j,]<-0}}}  

#Extract riparian forest (RF) and forest (FOR) covariates from owl.dat:
RF<-matrix(NA, nrow=dim(pX)[1], ncol=nYear-1)
FOR<-matrix(NA, nrow=dim(pX)[1], ncol=nYear-1)

for(i in 1:dim(RF)[1]){
  for(j in 1:dim(RF)[2]){
    RF[i,j]=unique(owls.dat$RF[which(owls.dat$Year==j & owls.dat$Site==i)])
    FOR[i,j]=unique(owls.dat$FOR[which(owls.dat$Year==j & owls.dat$Site==i)])
  }}


#Combine Y_BO and Y_NSO into a single matrix with four states:
#State 1 = no owls detected
#State 2 = only NSO detected
#State 3 = only BO detected
#State 4 = both BO and NSO detected
Y_BO_temp<-Y_BO
2->Y_BO_temp[Y_BO_temp==1]
Y<-Y_BO_temp+Y_NSO+1

#Matrix 'BO' is whether or not barred owls were detected at each site at least once in each year
BO<-numeric()
for(i in 1:22){BO<-cbind(BO,apply(Y_BO[,(i*8-7):(i*8)],1,sum))}
1->BO[BO>1]
#Matrix 'NSO' is whether or not spotted owls were detected at each site at least once in each year
NSO<-numeric()
for(i in 1:22){NSO<-cbind(NSO,apply(Y_NSO[,(i*8-7):(i*8)],1,sum))}
1->NSO[NSO>1]

#For Stan code:
n_occ<-apply(pX[,1:8,1],1,sum)
start<-1:22*8-7
for(i in 2:length(start)){n_occ<-cbind(n_occ,apply(pX[,start[i]:(start[i]+7),1],1,sum))}

#Data for calculating multi-level R^2:
mean_FOR_site<-apply(FOR,1,mean)
mean_FOR_time<-apply(FOR,2,mean)
mean_RF_site<-apply(RF,1,mean)
mean_RF_time<-apply(RF,2,mean)

#-----------------------------------------------------------------------------#
# Note: '//' is for comments in Stan
sink("Stan_Marginalized_owls_with_RE.stan")
cat("
    data{
    int<lower = 1, upper = 4> Y[158,176];
    real RF [158,21];
    real FOR [158,21];
    matrix[176,7] pX[158];
    vector[158] mean_FOR_site; // for multi-level R2 calculation
    vector[158] mean_RF_site;  // for multi-level R2 calculation
    vector[21] mean_FOR_time;  // for multi-level R2 calculation
    vector[21] mean_RF_time;   // for multi-level R2 calculation
    }
    
    parameters{
    real<lower = 0, upper = 1> psi0_BO;
    real<lower = 0, upper = 1> psi0_NSO;
    real<lower = -5, upper = 5> gamNSO_int; // time varying NSO colonization intercept
    real<lower = -5, upper = 5> gamNSO_b;   // Forest effect on NSO colonization
    vector[158] z_gamNSO_site;
    vector[21] z_gamNSO_time;
    real<lower = 0, upper = 5> sd_gamNSO_site;
    real<lower = 0, upper = 5> sd_gamNSO_time;
    
    real<lower = -5, upper = 5> epsNSO_int; // constant NSO intercept on extinction
    real<lower = -5, upper = 5> epsNSO_b;   // BO effect on NSO extinction
    vector[158] z_epsNSO_site;
    vector[21] z_epsNSO_time;
    real<lower = 0, upper = 5> sd_epsNSO_site;
    real<lower = 0, upper = 5> sd_epsNSO_time;
    
    real<lower = -5, upper = 5> gamBO_int;       // constant BO intercept on colonization
    real<lower = -5, upper = 5> gamBO_RF;        // Riparian forest effect, 
    real<lower = -15, upper = 15> gamBO_autolog; // autologistic effect
    real<lower = -5, upper = 5> gamBO_xNSO;      // NSO effect on BO colonization
    vector[158] z_gamBO_site;
    vector[21] z_gamBO_time;
    real<lower = 0, upper = 5> sd_gamBO_site;
    real<lower = 0, upper = 5> sd_gamBO_time;
    
    real<lower = -15,upper = 15> epsBO_int;      // constant BO intercept on extinction
    real<lower = -5, upper = 5> epsBO_RF;        // Riparian forest effect, 
    real<lower = -15, upper = 15> epsBO_autolog; // autologistic effect
    real<lower = -5, upper = 5> epsBO_xNSO;      // NSO effect on BO extinction
    vector[158] z_epsBO_site;
    vector[21] z_epsBO_time;
    real<lower = 0, upper = 5> sd_epsBO_site;
    real<lower = 0, upper = 5> sd_epsBO_time;
    
    vector<lower = -5, upper = 5>[6] pNSO_b; // intercept, day, night, method1, method3 and BO effects on NSO detection
    vector<lower = -5, upper = 5>[7] pBO_b;  // intercept, day, night, method1, method3, second half, and survey length effects on BO detection
    }
    
    transformed parameters{
    simplex [4] unc_pz; // initial probability of being in any state
    vector<lower = 0, upper = 1>[176] pBO;
    vector<lower = 0, upper = 1>[176] pNSO_BO;
    vector<lower = 0, upper = 1>[176] pNSO_bo;
    real<lower = 0, upper = 1> gamNSO;
    real<lower = 0, upper = 1> epsNSO_BO;
    real<lower = 0, upper = 1> epsNSO_bo;
    real<lower = 0, upper = 1> gamBO_NSO;
    real<lower = 0, upper = 1> gamBO_nso;
    real<lower = 0, upper = 1> epsBO_NSO;
    real<lower = 0, upper = 1> epsBO_nso;
    simplex [4] p [176,4];
    
    // cumulative probability of observation history (8 occasions a year 
    // for 22 years for a total of 176 occasions)
    real<lower = 0, upper = 1> cp [158,22,4]; 
    simplex [4] tr [4];
    real<lower = 0, upper = 1> pz [158,22,4];
    real<lower = 0, upper = 1> BO [158];
    real<lower = 0, upper = 1> NSO [158];
    vector<lower = -0.5, upper = 0.5>[21] st_psiBO;
    vector<lower = 0, upper = 1>[21] psiNSO_time_mean;
    vector<lower = 0, upper = 1>[21] psiBO_time_mean;
    vector<lower = 0, upper = 1>[4] temp;
    
    // Unconditional occupancy probability (used to help determine occupancy state at t=1)
    unc_pz[1] = (1 - psi0_BO) * (1 - psi0_NSO); // unc_pz = probability of neither species being present in the study are at the start of the study.
    unc_pz[2] = (1 - psi0_BO) * psi0_NSO;
    unc_pz[3] = psi0_BO * (1 - psi0_NSO);
    unc_pz[4] = psi0_BO * psi0_NSO;
    
    // Observation model
    // 22 years and 8 occasions per year
    for(k in 1:158){
    pBO = pX[k,,1] .* inv_logit(pX[k,,] * pBO_b);
    pNSO_BO = pX[k,,1] .* inv_logit(pX[k,,1:5] * pNSO_b[1:5] + pNSO_b[6]); // NSO detection probability if BO is present
    pNSO_bo = pX[k,,1] .* inv_logit(pX[k,,1:5] * pNSO_b[1:5]);  // NSO detection probability if BO is absent
    
    for(j in 1:176){
    p[j,1,1] = 1; // if two species absent, then of course state of detection will be 1 (i.e., no species will be detected)
    p[j,1,2] = 0;
    p[j,1,3] = 0;
    p[j,1,4] = 0;
    p[j,2,1] = 1 - pNSO_bo[j]; // if NSO present, probability of it not being detected
    p[j,2,2] = pNSO_bo[j];     // if NSO present, probability of it being detected
    p[j,2,3] = 0;
    p[j,2,4] = 0;
    p[j,3,1] = 1 - pBO[j];     // if BO present, probability of it not being detected
    p[j,3,2] = 0;
    p[j,3,3] = pBO[j];         // if BO present, probability of it being detected
    p[j,3,4] = 0;
    p[j,4,1] = 1 - (pBO[j] + pNSO_BO[j] - pBO[j] * pNSO_BO[j]); // I probably would have written this as : (1-pBO)*(1-pNSO_BO)*visited
    p[j,4,2] = pNSO_BO[j] * (1 - pBO[j]); // if both species present, probability of only detecting BO
    p[j,4,3] = (1 - pNSO_BO[j]) * pBO[j]; // if both species present, probability of only detecting NSO
    p[j,4,4] = pNSO_BO[j] * pBO[j]; // probability of detecting both species
    }
    
    for(t in 1:22){
    cp[k,t,1] = p[(1+(t-1)*8),1,Y[k,(1+(t-1)*8)]]; // what is the probability that site is in state 1-4 given state is observed in state Y[k,,] 
    cp[k,t,2] = p[(1+(t-1)*8),2,Y[k,(1+(t-1)*8)]];
    cp[k,t,3] = p[(1+(t-1)*8),3,Y[k,(1+(t-1)*8)]];
    cp[k,t,4] = p[(1+(t-1)*8),4,Y[k,(1+(t-1)*8)]];
    
    for(j in 2:8){
    cp[k,t,1] = cp[k,t,1]*p[(j+(t-1)*8),1,Y[k,(j+(t-1)*8)]]; // update probabilities given subsequent visits within the same year 
    cp[k,t,2] = cp[k,t,2]*p[(j+(t-1)*8),2,Y[k,(j+(t-1)*8)]];
    cp[k,t,3] = cp[k,t,3]*p[(j+(t-1)*8),3,Y[k,(j+(t-1)*8)]];
    cp[k,t,4] = cp[k,t,4]*p[(j+(t-1)*8),4,Y[k,(j+(t-1)*8)]];
    }
    }
    }
    
    // Process model: state at t=1 
    for(k in 1:158){
    pz[k,1,1] = unc_pz[1] * cp[k,1,1]; // for first year, update probability of not being occupied by either species (kinda like a prior) by the observation history in the first year
    pz[k,1,2] = unc_pz[2] * cp[k,1,2];
    pz[k,1,3] = unc_pz[3] * cp[k,1,3];
    pz[k,1,4] = unc_pz[4] * cp[k,1,4];
    BO[k] = sum(pz[k,1,3:4]) / sum(pz[k,1,]); // this is the probability of a BO being present on the site
    NSO[k] = (pz[k,1,2] + pz[k,1,4]) / sum(pz[k,1,]);
    }
    
    for(t in 1:21){
    psiNSO_time_mean[t] = sum(NSO) / 158;   // mean probability of NSO presence
    psiBO_time_mean[t] = sum(BO) / 158;     // mean probability of NSO presence
    st_psiBO[t] = psiBO_time_mean[t] - 0.5; // mean probability of BO presence (subtract 0.5 to center it)
    
    // Process model: extinction / colonization dynamics 
    for(k in 1:158){
    // probability of extinction of NSO given BO present
    epsNSO_BO = inv_logit(epsNSO_int + epsNSO_b + z_epsNSO_site[k] * sd_epsNSO_site + z_epsNSO_time[t] * sd_epsNSO_time);
    // probability of extinction of NSO given BO absent
    epsNSO_bo = inv_logit(epsNSO_int + z_epsNSO_site[k] * sd_epsNSO_site + z_epsNSO_time[t] * sd_epsNSO_time);
    // probability of colonization by NSO is a function of forest cover?
    gamNSO = inv_logit(gamNSO_int + gamNSO_b * FOR[k,t] + z_gamNSO_site[k] * sd_gamNSO_site + z_gamNSO_time[t] * sd_gamNSO_time);
    // probability of BO colonization given NSO present (think this is a function of riparian forest, and autologistic effects)
    gamBO_NSO = inv_logit(gamBO_int + gamBO_RF * RF[k,t] + gamBO_autolog * st_psiBO[t] + gamBO_xNSO + z_gamBO_site[k] * sd_gamBO_site + z_gamBO_time[t] * sd_gamBO_time);
    gamBO_nso = inv_logit(gamBO_int + gamBO_RF * RF[k,t] + gamBO_autolog * st_psiBO[t] + z_gamBO_site[k] * sd_gamBO_site + z_gamBO_time[t] * sd_gamBO_time);
    epsBO_NSO = inv_logit(epsBO_int + epsBO_RF * RF[k,t] + epsBO_autolog * st_psiBO[t] + epsBO_xNSO + z_epsBO_site[k] * sd_epsBO_site + z_epsBO_time[t] * sd_epsBO_time);
    epsBO_nso = inv_logit(epsBO_int + epsBO_RF * RF[k,t] + epsBO_autolog * st_psiBO[t] + z_epsBO_site[k] * sd_epsBO_site + z_epsBO_time[t] * sd_epsBO_time);
    tr[1,1] = (1 - gamNSO) * (1 - gamBO_nso);
    tr[1,2] = gamNSO * (1 - gamBO_nso);
    tr[1,3] = (1 - gamNSO) * gamBO_nso;
    tr[1,4] = gamNSO * gamBO_nso;
    tr[2,1] = epsNSO_bo * (1 - gamBO_NSO);
    tr[2,2] = (1 - epsNSO_bo) * (1 - gamBO_NSO);
    tr[2,3] = epsNSO_bo * gamBO_NSO;
    tr[2,4] = (1 - epsNSO_bo) * gamBO_NSO;
    tr[3,1] = (1 - gamNSO) * epsBO_nso;
    tr[3,2] = gamNSO * epsBO_nso;
    tr[3,3] = (1 - gamNSO) * (1 - epsBO_nso);
    tr[3,4] = gamNSO * (1 - epsBO_nso);
    tr[4,1] = epsNSO_BO * epsBO_NSO;
    tr[4,2] = (1 - epsNSO_BO) * epsBO_NSO;
    tr[4,3] = epsNSO_BO * (1 - epsBO_NSO);
    tr[4,4] = (1 - epsNSO_BO) * (1 - epsBO_NSO);
    
    for(j in 1:4){
    temp[j] = pz[k,t,j] * tr[j,1];
    }
    
    pz[k,(t+1),1] = sum(temp) * cp[k,(t+1),1];
    
    for(j in 1:4){
    temp[j] = pz[k,t,j] * tr[j,2];
    }
    
    pz[k,(t+1),2] = sum(temp) * cp[k,(t+1),2];
    
    for(j in 1:4){
    temp[j] = pz[k,t,j] * tr[j,3];
    }
    
    pz[k,(t+1),3] = sum(temp) * cp[k,(t+1),3];
    
    for(j in 1:4){
    temp[j] = pz[k,t,j] * tr[j,4];
    }
    
    pz[k,(t+1),4] = sum(temp) * cp[k,(t+1),4];
    
    // Determine the percentage of sites occupied by Barred Owls (for use in autologistic extinction/colonization model)
    BO[k] = sum(pz[k,(t+1),3:4]) / sum(pz[k,(t+1),]);
    // Determine the percentage of sites occupied by Northern Spotted Owls (for use in multi-level R^2 calculation)
    NSO[k] = (pz[k,(t+1),2] + pz[k,(t+1),4]) / sum(pz[k,(t+1),]);
    }
    }
    }
    
    model{
    z_gamBO_site ~ normal(0,1);
    z_epsBO_site ~ normal(0,1);
    z_gamNSO_site ~ normal(0,1);
    z_epsNSO_site ~ normal(0,1);
    z_gamBO_time ~ normal(0,1);
    z_epsBO_time ~ normal(0,1);
    z_gamNSO_time ~ normal(0,1);
    z_epsNSO_time ~ normal(0,1);
    
    for(k in 1:158){  
    target += log(sum(pz[k,22,])); 
    }
    }
    
    generated quantities{
    real psiBO_site[158,21];
    real psiNSO_site[158,21];
    vector<lower = 0, upper = 1>[158] psiBO_site_mean;
    vector<lower = 0, upper = 1>[158] psiNSO_site_mean;
    
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
    for(t in 1:21){
    psiBO_site[k,t] = sum(pz[k,t,3:4]) / sum(pz[k,t,]);
    psiNSO_site[k,t] = (pz[k,t,2] + pz[k,t,4]) / sum(pz[k,t,]);
    }
    psiBO_site_mean[k] = mean(psiBO_site[k,1:21]);
    psiNSO_site_mean[k] = mean(psiNSO_site[k,1:21]);
    }
    
    // excluding autologistic effect on spatial RE because the average autologistic 
    // effect for all 21 intervals is the same for each site:
    gamBO_site_pred = gamBO_int + gamBO_RF * mean_RF_site + gamBO_xNSO * psiNSO_site_mean + z_gamBO_site * sd_gamBO_site;
    gamBO_time_pred = gamBO_int + gamBO_RF * mean_RF_time + gamBO_autolog * st_psiBO + gamBO_xNSO * psiNSO_time_mean + z_gamBO_time * sd_gamBO_time;
    R2_gamBO_site = 1 - variance(gamBO_int + z_gamBO_site * sd_gamBO_site) / variance(gamBO_site_pred);
    R2_gamBO_time = 1 - variance(gamBO_int + z_gamBO_time * sd_gamBO_time) / variance(gamBO_time_pred);
    
    // excluding autologistic effect on spatial RE because the average autologistic
    // effect for all 21 intervals is the same for each site:
    epsBO_site_pred = epsBO_int + epsBO_RF * mean_RF_site + epsBO_xNSO * psiNSO_site_mean + z_epsBO_site * sd_epsBO_site;
    epsBO_time_pred = epsBO_int + epsBO_RF * mean_RF_time + epsBO_autolog * st_psiBO + epsBO_xNSO * psiNSO_time_mean + z_epsBO_time * sd_epsBO_time;
    R2_epsBO_site = 1 - variance(epsBO_int + z_epsBO_site * sd_epsBO_site) / variance(epsBO_site_pred);
    R2_epsBO_time = 1 - variance(epsBO_int + z_epsBO_time * sd_epsBO_time) / variance(epsBO_time_pred);
    
    gamNSO_site_pred = gamNSO_int + gamNSO_b * mean_FOR_site + z_gamNSO_site * sd_gamNSO_site;
    gamNSO_time_pred = gamNSO_int + gamNSO_b * mean_FOR_time + z_gamNSO_time * sd_gamNSO_time;
    R2_gamNSO_site = 1 - variance(gamNSO_int + z_gamNSO_site * sd_gamNSO_site) / variance(gamNSO_site_pred);
    R2_gamNSO_time = 1 - variance(gamNSO_int+z_gamNSO_time * sd_gamNSO_time) / variance(gamNSO_time_pred);
    
    epsNSO_site_pred = epsNSO_int + epsNSO_b * psiBO_site_mean + z_epsNSO_site * sd_epsNSO_site;
    epsNSO_time_pred = epsNSO_int + epsNSO_b * psiBO_time_mean + z_epsNSO_time * sd_epsNSO_time;
    R2_epsNSO_site = 1 - variance(epsNSO_int + z_epsNSO_site * sd_epsNSO_site) / variance(epsNSO_site_pred);
    R2_epsNSO_time = 1 - variance(epsNSO_int + z_epsNSO_time * sd_epsNSO_time) / variance(epsNSO_time_pred);
    }
    ", fill = TRUE)
sink()  
#-----------------------------------------------------------------------------#
owls_SMRE.data <- list(Y = array(Y, dim = c(dim(Y)[1], dim(Y)[2])),
                       RF = array(as.vector(unlist(RF)), dim = c(dim(RF)[1],dim(RF)[2])),
                       FOR = array(as.vector(unlist(FOR)), dim = c(dim(FOR)[1], dim(FOR)[2])),
                       pX = array(pX, dim = c(dim(pX)[1], dim(pX)[2], dim(pX)[3])),
                       mean_FOR_site = mean_FOR_site, mean_FOR_time = mean_FOR_time,
                       mean_RF_site = mean_RF_site, mean_RF_time = mean_RF_time)

owls_SMRE.par <- c('pBO_b', 'pNSO_b', 'psi0_BO', 'psi0_NSO', 'gamNSO_b',
                   'epsNSO_b', 'epsNSO_int', 'gamBO_int', 'epsBO_int', 'gamBO_RF',
                   'gamBO_autolog', 'gamBO_xNSO', 'epsBO_RF', 'epsBO_autolog',
                   'epsBO_xNSO', 'gamNSO_int', 'R2_gamBO_site', 'R2_gamBO_time',
                   'R2_gamNSO_site', 'R2_gamNSO_time', 'R2_epsBO_site', 
                   'R2_epsBO_time', 'R2_epsNSO_site' , 'R2_epsNSO_time')

owls_SMRE.out <- stan(".\\Stan_Marginalized_owls_with_RE.stan",
                      data = owls_SMRE.data, pars = owls_SMRE.par,
                      chains = 3, iter = 1000, seed = 1)

#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Application 6
#  N-Occupancy Model
#  Discrete JAGS version 
#
#  Notes:
#   *  Need to set directory for data
#   *  Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2jags)

#-----------------------------------------------------------------------------#
# import, format data for model fitting
nocc.dat<-read.csv(".//nocc.dat.csv")

# extract all BO data

nYears<-26
nVisits<-8
nSites<-158
Ninf<-10
Ns<-c(0:Ninf)

BO<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
day<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
night<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
RF<-array(0,dim=c(max(nocc.dat$Site),nYears-1))


#BO : Matrix of barred owl detections:
# 1 = barred owl detected
# 0 = barred owl not detected
# NA = site not visited

for(i in 1:length(nocc.dat[,1])){
  BO[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-nocc.dat$BO[i]
  if(is.na(nocc.dat$TOD[i])==FALSE & nocc.dat$TOD[i]==1){
    day[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-1}
  if(is.na(nocc.dat$TOD[i])==FALSE & nocc.dat$TOD[i]==3){
    night[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-1}}

day[which(is.na(BO))]<-NA
night[which(is.na(BO))]<-NA
TOD<-2-day+night

#Get riparian forest (RF) covariate:
RF<-matrix(NA, nrow=dim(BO)[1], ncol=nYears-1)

for(i in 1:dim(RF)[1]){
  for(j in 1:dim(RF)[2]){
    RF[i,j]=unique(nocc.dat$RF[which(nocc.dat$Year==j & nocc.dat$Site==i)])
  }}

Y<-numeric()
tod<-numeric()
site<-numeric()
year<-numeric()
for (i in 1:nSites){
  temp<-which(is.na(BO[i,])==FALSE)
  site<-c(site,rep(i,length(temp)))
  Y<-c(Y,as.numeric(BO[i,temp]))
  tod<-c(tod,as.numeric(TOD[i,temp]))
  year<-c(year,floor((temp+7)/nVisits))
}

#-----------------------------------------------------------------------------#
sink("Nocc_JD.jags")
cat("
    model { 
    for (i in 1:nSites) {N[i,1] ~ dpois(lambda)} #latent counts in first year
    for (t in 1:(nYears-1)){
    N_mean[t] <- mean(N[,t])  -1				# mean abundance in prior year
    log(gamma[t]) <- a0 + a1*N_mean[t] 		# expected gains as a function of mean population size
    for (i in 1:nSites){
    logit(omega[i,t]) <- b0 + b1*RF[i,t] # survial probability as a function of riparian forest
    S[i,t] ~ dbin(omega[i,t], N[i,t])    # simulate number that survive
    G[i,t] ~ dpois(gamma[t])             # simulate additions to site population
    N[i,(t+1)] <- S[i,t] + G[i,t]    # add survivors to newcomers
    }}
    for (k in 1:nSamples){
    occ_p[k] <- 1-pow((1-p[tod[k]]),N[site[k],year[k]]) # probability of not detecting any individuals
    Y[k] ~ dbern(occ_p[k])                             
    }
    #Priors
    lambda ~ dunif(0,2) 
    a0 ~ dunif(-4,4) 
    a1 ~ dunif(-4,4) 
    b0 ~ dunif(-4,4) 
    b1 ~ dunif(-4,4) 
    for (k in 1:3) {p[k] ~ dunif(0,1)} # different detection probabilities by tod
    }
    
    ", fill = TRUE)
sink()  

#-----------------------------------------------------------------------------#


JD_data<-list(Y=Y,nYears=nYears,nSites=nSites,nSamples=length(Y),RF=RF,tod=tod,site=site,year=year)

JD_params <- c("a0", "a1", "b0", "b1", "p", "lambda","N_mean")


JD_inits<-function() list(N=matrix(c(rep(5,nSites),rep(NA,nSites*(nYears-1))),nrow=nSites,ncol=nYears),
                          S=matrix(1,nrow=nSites,ncol=(nYears-1)),
                          G=matrix(1,nrow=nSites,ncol=(nYears-1))
)

JD_nocc<-jags.parallel(JD_data, JD_inits, JD_params, model.file="Nocc_JD.txt",n.iter=20000,n.chain=3)
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Application 6
#  N-Occupancy Model
#  Marginalized Stan version 
#
#  Notes:
#   
###############################################################################
library(rstan)

rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

#-----------------------------------------------------------------------------#
# import, format data for model fitting
nocc.dat<-read.csv(".//nocc.dat.csv")

# extract all BO data

nYears<-26
nVisits<-8
nSites<-158
Ninf<-10
Ns<-c(0:Ninf)

BO<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
day<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
night<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
RF<-array(0,dim=c(max(nocc.dat$Site),nYears-1))


#BO : Matrix of barred owl detections:
# 1 = barred owl detected
# 0 = barred owl not detected
# NA = site not visited

for(i in 1:length(nocc.dat[,1])){
  BO[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-nocc.dat$BO[i]
  if(is.na(nocc.dat$TOD[i])==FALSE & nocc.dat$TOD[i]==1){
    day[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-1}
  if(is.na(nocc.dat$TOD[i])==FALSE & nocc.dat$TOD[i]==3){
    night[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-1}}

day[which(is.na(BO))]<-NA
night[which(is.na(BO))]<-NA
TOD<-2-day+night

#Get riparian forest (RF) covariate:
RF<-matrix(NA, nrow=dim(BO)[1], ncol=nYears-1)

for(i in 1:dim(RF)[1]){
  for(j in 1:dim(RF)[2]){
    RF[i,j]=unique(nocc.dat$RF[which(nocc.dat$Year==j & nocc.dat$Site==i)])
  }}



Y<-numeric()
tod<-numeric()
site<-numeric()
year<-numeric()
for (i in 1:nSites){
  temp<-which(is.na(BO[i,])==FALSE)
  site<-c(site,rep(i,length(temp)))
  Y<-c(Y,as.numeric(BO[i,temp]))
  tod<-c(tod,as.numeric(TOD[i,temp]))
  year<-c(year,floor((temp+7)/nVisits))
}





nD0<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was not detected during day samples
nD1<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was detected during day samples
nC0<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was not detected during crepuscular samples
nC1<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was detected during crepuscular samples
nN0<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was not detected during night samples
nN1<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was detected during night samples
for (i in 1:nSites){
  for (t in 1:nYears){
    nD0[i,t]<-length(which(site==i&year==t&tod==1&Y==0))
    nD1[i,t]<-length(which(site==i&year==t&tod==1&Y==1))
    nC0[i,t]<-length(which(site==i&year==t&tod==2&Y==0))
    nC1[i,t]<-length(which(site==i&year==t&tod==2&Y==1))
    nN0[i,t]<-length(which(site==i&year==t&tod==3&Y==0))
    nN1[i,t]<-length(which(site==i&year==t&tod==3&Y==1))
  }}

#-----------------------------------------------------------------------------#

sink("Nocc_SM.stan")
cat("
    data {
    int<lower=1> nYears;
    int<lower=1> nSites;
    int<lower=1> Ninf;
    row_vector [(Ninf+1)] Ns;
    matrix [nSites,(nYears-1)] RF;
    matrix [nSites,nYears] nD0;
    matrix [nSites,nYears] nD1;
    matrix [nSites,nYears] nC0;
    matrix [nSites,nYears] nC1;
    matrix [nSites,nYears] nN0;
    matrix [nSites,nYears] nN1;	
    }
    
    parameters {
    real<lower=0,upper=1> lambda;
    real<lower=0,upper=1> p [3];
    real<lower=-4,upper=4> a0;
    real<lower=-4,upper=4> a1;
    real<lower=-4,upper=4> b0;
    real<lower=-4,upper=4> b1;
    }
    
    transformed parameters {
    real po [3,2,(Ninf+1)];
    row_vector [(Ninf+1)] pN;
    row_vector [(Ninf+1)] cp [nSites,nYears];
    row_vector [(Ninf+1)] pz [nSites,nYears];
    row_vector [(Ninf+1)] spz [nSites,nYears];
    real N_mn [nSites,nYears];
    real N_mean [(nYears-1)];
    real gamma [(nYears-1)];
    matrix [(Ninf+1),(Ninf+1)] gtr [(nYears-1)];
    matrix [(Ninf+1),(Ninf+1)] str [(nYears-1),nSites];
    matrix [nSites,(nYears-1)] omega ;
    
    omega = inv_logit(b0 + b1*RF);
    for (j in 1:(Ninf+1)){
    po[1,1,j]=(1-p[1])^(j-1);
    po[1,2,j]=1-po[1,1,j];
    po[2,1,j]=(1-p[2])^(j-1);
    po[2,2,j]=1-po[2,1,j];
    po[3,1,j]=(1-p[3])^(j-1);
    po[3,2,j]=1-po[3,1,j];
    pN[j]=exp(poisson_lpmf((j-1) | lambda));
    }
    for (i in 1:nSites){
    for (t in 1:nYears){
    for (j in 1:(Ninf+1)){
    cp[i,t,j]=((po[1,1,j])^nD0[i,t])*((po[1,2,j])^nD1[i,t])*((po[2,1,j])^nC0[i,t])*
    ((po[2,2,j])^nC1[i,t])*((po[3,1,j])^nN0[i,t])*((po[3,2,j])^nN1[i,t]);
    }}
    pz[i,1,]=pN .* cp[i,1,];
    spz[i,1,]=pz[i,1,]/sum(pz[i,1,]);
    N_mn[i,1]=sum(Ns .* spz[i,1,]);              
    }
    for(t in 1:(nYears-1)) {
    N_mean[t] = mean(N_mn[,t]) - 1;
    gamma[t] = exp(a0 + a1*N_mean[t]);
    for (j in 2:(Ninf+1)){
    for (k in 1:(j-1)){
    gtr[t,j,k]=0;
    }}
    for (k in 1:(Ninf)) gtr[t,1,k]= exp(poisson_lpmf((k-1) | gamma[t]));
    for (j in 2:(Ninf)){
    for (k in j:(Ninf)) gtr[t,j,k]= exp(poisson_lpmf((k-j) | gamma[t]));
    }
    for (j in 1:(Ninf+1)) gtr[t,j,(Ninf+1)]=1-sum(gtr[t,j,1:Ninf]);
    for(i in 1:nSites) {
    for (j in 1:(Ninf)){
    for (k in (j+1):(Ninf+1)){
    str[t,i,j,k]=0;
    }}
    str[t,i,1,1]=1;
    for (j in 2:(Ninf+1)){
    for (k in 1:j){
    str[t,i,j,k] = exp(binomial_lpmf ((k-1) | (j-1), omega[i,t]));
    }}
    pz[i,(t+1),]=((spz[i,t,]*str[t,i,,])*gtr[t,,]) .* cp[i,t,];
    spz[i,(t+1),]=pz[i,(t+1),]/sum(pz[i,(t+1),]);
    N_mn[i,(t+1)]=sum(Ns .* spz[i,(t+1),]);
    }}}
    
    model {
    for (i in 1:nSites) {
    for (t in 1:nYears){
    target += log(sum(pz[i,t,]));
    }}}
    
    ", fill = TRUE)
sink()
#-----------------------------------------------------------------------------#


SM_data<-list(nYears=nYears,nSites=nSites,RF=RF,nD0=nD0,nD1=nD1,nC0=nC0,nC1=nC1,nN0=nN0,nN1=nN1,Ninf=Ninf,Ns=Ns)

SM_params <- c("a0", "a1", "b0", "b1", "p", "lambda","N_mean")

SM_nocc<-stan("Nocc_SM.stan",data=SM_data,pars=SM_params,iter=10,chains=1) 

###############################################################################
#                                                                     Spring 19
#  Application 6
#  N-Occupancy Model
#  Maximum likelihood version 
#
#  Notes: Includes two different functions
#         function 'Nocc_ML' estimates occupancy by first estimating underlying 
#           site-specific abundance (comparable to Stan and JAGS examples above)
#         function 'occ_ML' estimates occupancy without estimating abundance
#   
###############################################################################


#-----------------------------------------------------------------------------#
# import, format data for model fitting
nocc.dat<-read.csv(".//nocc.dat.csv")

# extract all BO data

nYears<-26
nVisits<-8
nSites<-158
Ninf<-10
Ns<-c(0:Ninf)

BO<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
day<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
night<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
RF<-array(0,dim=c(max(nocc.dat$Site),nYears-1))


#BO : Matrix of barred owl detections:
# 1 = barred owl detected
# 0 = barred owl not detected
# NA = site not visited

for(i in 1:length(nocc.dat[,1])){
  BO[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-nocc.dat$BO[i]
  if(is.na(nocc.dat$TOD[i])==FALSE & nocc.dat$TOD[i]==1){
    day[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-1}
  if(is.na(nocc.dat$TOD[i])==FALSE & nocc.dat$TOD[i]==3){
    night[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-1}}

day[which(is.na(BO))]<-NA
night[which(is.na(BO))]<-NA
TOD<-2-day+night

#Get riparian forest (RF) covariate:
RF<-matrix(NA, nrow=dim(BO)[1], ncol=nYears-1)

for(i in 1:dim(RF)[1]){
  for(j in 1:dim(RF)[2]){
    RF[i,j]=unique(nocc.dat$RF[which(nocc.dat$Year==j & nocc.dat$Site==i)])
  }}



Y<-numeric()
tod<-numeric()
site<-numeric()
year<-numeric()
for (i in 1:nSites){
  temp<-which(is.na(BO[i,])==FALSE)
  site<-c(site,rep(i,length(temp)))
  Y<-c(Y,as.numeric(BO[i,temp]))
  tod<-c(tod,as.numeric(TOD[i,temp]))
  year<-c(year,floor((temp+7)/nVisits))
}




nD0<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was not detected during day samples
nD1<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was detected during day samples
nC0<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was not detected during crepuscular samples
nC1<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was detected during crepuscular samples
nN0<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was not detected during night samples
nN1<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was detected during night samples
for (i in 1:nSites){
  for (t in 1:nYears){
    nD0[i,t]<-length(which(site==i&year==t&tod==1&Y==0))
    nD1[i,t]<-length(which(site==i&year==t&tod==1&Y==1))
    nC0[i,t]<-length(which(site==i&year==t&tod==2&Y==0))
    nC1[i,t]<-length(which(site==i&year==t&tod==2&Y==1))
    nN0[i,t]<-length(which(site==i&year==t&tod==3&Y==0))
    nN1[i,t]<-length(which(site==i&year==t&tod==3&Y==1))
  }}
#-----------------------------------------------------------------------------#

#Maximum likelihood function (N-occupancy):

nocc_ML<-function(par){
  lambda<-par[1]
  a0<-par[2]
  a1<-par[3]
  b0<-par[4]
  b1<-par[5]
  p<-numeric()
  p[1]<-par[6]
  p[2]<-par[7]
  p[3]<-par[8]
  #
  pN<-dpois(0:Ninf,lambda)
  Ns<-0:Ninf
  po<-array(0,dim=c(3,2,(Ninf+1)))
  po[1,1,]<-(1-p[1])^Ns
  po[1,2,]<-1-po[1,1,]
  po[2,1,]<-(1-p[2])^Ns
  po[2,2,]<-1-po[2,1,]
  po[3,1,]<-(1-p[3])^Ns
  po[3,2,]<-1-po[3,1,]
  cp<-array(0,dim=c(nSites,nYears,(Ninf+1)))
  pz<-matrix(NA,nrow=nSites,ncol=(Ninf+1))
  N_mn<-numeric()
  gtr<-matrix(0,nrow=Ninf+1,ncol=Ninf+1)
  str<-matrix(0,nrow=Ninf+1,ncol=Ninf+1)
  llik<-numeric()
  #
  for (t in 1:nYears){
    for (j in 1:(Ninf+1)){
      cp[,t,j]<-((po[1,1,j])^nD0[,t])*((po[1,2,j])^nD1[,t])*((po[2,1,j])^nC0[,t])*
        ((po[2,2,j])^nC1[,t])*((po[3,1,j])^nN0[,t])*((po[3,2,j])^nN1[,t])
    }}
  for (i in 1:nSites){
    pz[i,]=pN*cp[i,1,]
    N_mn[i]<-Ns%*%pz[i,]/sum(pz[i,])
  }
  for(t in 1:(nYears-1)) {
    N_mean <- mean(N_mn)-1
    gamma <- exp(a0 + a1*N_mean)
    for (j in 1:(Ninf+1)){
      gtr[j,] <- dpois((Ns-j+1),gamma)/ppois(Ninf-j+1,gamma)
    }
    for(i in 1:nSites) {
      omega <- plogis(b0 + b1*RF[i,t])
      for (j in 1:(Ninf+1)){
        str[j,] = pz[i,j]*dbinom(Ns,(j-1),omega)
      }
      tr<-colSums(str)%*%gtr
      pz[i,]<-tr*cp[i,t,]
      N_mn[i]<-(Ns%*%pz[i,])/sum(pz[i,])
    }}
  for (i in 1:nSites){
    llik[i]<-log(sum(pz[i,]))
  }		
  -1*sum(llik)
}
#-----------------------------------------------------------------------------#

#Maximum likelihood function (Dynamic occupancy):

occ_ML<-function(par){
  psi<-par[1]
  a0<-par[2]
  a1<-par[3]
  b0<-par[4]
  b1<-par[5]
  p<-numeric()
  p[1]<-par[6]
  p[2]<-par[7]
  p[3]<-par[8]
  #
  po<-array(0,dim=c(3,2,2))
  po[1,1,]<-(1-p[1])^c(0,1)
  po[1,2,]<-1-po[1,1,]
  po[2,1,]<-(1-p[2])^c(0,1)
  po[2,2,]<-1-po[2,1,]
  po[3,1,]<-(1-p[3])^c(0,1)
  po[3,2,]<-1-po[3,1,]
  cp<-array(0,dim=c(nSites,nYears,2))
  pz<-matrix(NA,nrow=nSites,ncol=2)
  psi_mn<-numeric()
  tr<-matrix(0,nrow=2,ncol=2)
  llik<-numeric()
  #
  for (t in 1:nYears){
    for (j in 1:2){
      cp[,t,j]<-((po[1,1,j])^nD0[,t])*((po[1,2,j])^nD1[,t])*((po[2,1,j])^nC0[,t])*
        ((po[2,2,j])^nC1[,t])*((po[3,1,j])^nN0[,t])*((po[3,2,j])^nN1[,t])
    }}
  for (i in 1:nSites){
    pz[i,]=c(1-psi,psi)*cp[i,1,]
    psi_mn[i]<-pz[i,2]/sum(pz[i,])
  }
  for(t in 1:(nYears-1)) {
    psi_mean <- mean(psi_mn)
    gam<-plogis(a0+a1*psi_mean)
    tr[1,]<-c(1-gam,gam)
    for(i in 1:nSites) {
      eps <- plogis(b0 + b1*RF[i,t])
      tr[2,]<-c(eps,1-eps)
      pz[i,]<-(pz[i,]%*%tr)*cp[i,t,]
      psi_mn[i]<-pz[i,2]/sum(pz[i,])
    }}
  for (i in 1:nSites){
    llik[i]<-log(sum(pz[i,]))
  }		
  -1*sum(llik)
}
#-----------------------------------------------------------------------------#

# Run models:
m_Nocc_ML<-optim(rep(0.5,8),nocc_ML,method="BFGS",hessian=TRUE)
m_occ_ML<-optim(rep(0.5,8),occ_ML,method="BFGS",hessian=TRUE)

#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Code for autologistic simulation summarized in Figure 5
#
#  Notes:
#   *  
#
###############################################################################

# creates matrix specifying neighbors of each site 
N <- matrix(0, nrow = 100, ncol = 100)
for(i in 1:100){
  t <- c(i - 2, i - 1, i + 1, i + 2)
  t <- subset(t, t > 0 & t < 101)
  N[i,t] <- 1 / length(t)
}

# returns the xth unique combination of 1's and 0's for N slots
makeZ <- function(x, N){
  i <- 0
  string <- rep(0, N)
  while(x > 0){
    string[N - i] <- x %% 2
    x <- x %/% 2
    i <- i + 1 
  }
  return(string)
}

# create all possible temporal patterns of latent occupancy (potential Z's) over
# six occasions (five intervals)
potZ <- matrix(0, nrow = 2 ^ 6, ncol = 6)
for(j in 1:2 ^ 6){
  potZ[j,] <- makeZ(j,6)
}

pZarray <- array(0, dim = c(64, 6, 2))
for(i in 1:64){
  for(t in 1:6){
    pZarray[i,t,1] <- potZ[i,t]
    pZarray[i,t,2] <- 1-potZ[i,t]
  }
}

# Fully conditional code requires a series of functions to iterate between
# making parameter estimates and conditional estimates, then treating
# conditional estimates as data and rerunning, likelihood for fully conditional
# in first iteration

fc_lik0 <- function(par){
  psi0 <- plogis(par[1])
  eps <- plogis(par[2])
  p <- plogis(par[3])
  g_int <- par[4]
  g <- matrix(0, nrow = 100, ncol = 5)
  for(t in 1:5){
    g[,t] <- plogis(g_int)
  }
  cp <- array(0, dim = c(100, 6, 2))
  zeta <- array(0, dim = c(100, 6, 2))
  po <- matrix(0, nrow = 2, ncol = 2)
  po[1,1] <- p
  po[1,2] <- (1 - p)
  po[2,1] <- 0
  po[2,2] <- 1
  for(t in 1:6){
    cp[,t,1] <- po[1,Y[,(2 * t - 1)]] * po[1,Y[,(2 * t)]]
    cp[,t,2] <- po[2,Y[,(2 * t - 1)]] * po[2,Y[,(2 * t)]]
  }
  zeta[,1,1] <- psi0 * cp[,1,1]
  zeta[,1,2] <- (1 - psi0) * cp[,1,2]
  for(t in 1:5){
    zeta[,(t + 1),1] <- (zeta[,t,1] * (1 - eps) + zeta[,t,2] * g[,t]) * cp[,(t + 1),1]
    zeta[,(t + 1),2] <- (zeta[,t,1] * eps + zeta[,t,2] * (1 - g[,t])) * cp[,(t + 1),2]
  }	
  nll <- 0
  for(k in 1:100){
    nll <- nll - 1 * log(sum(zeta[k,6,]))
  }
  nll
}

# this function creates fully conditional estimates given parameter estimates
fc_pred <- function(par, cX){
  psi0 <- plogis(par[1])
  eps <- plogis(par[2])
  p <- plogis(par[3])
  g_int <- par[4]
  g_auto <- par[5]
  g <- matrix(0, nrow = 100, ncol = 5)
  for(t in 1:5){
    g[,t] <- plogis(g_int + g_auto * (cX[,t] %*% t(N)))
  }
  out <- matrix(NA, nrow = 100, ncol = 5)
  po <- matrix(0, nrow = 2, ncol = 3)
  po[1,1] <- (1 - p) ^ 2
  po[1,2] <- 2 * p * (1 - p)
  po[1,3] <- p ^ 2
  po[2,1] <- 1
  y <- ifelse(Y == 2, 0, 1)
  sy <- 1 + cbind(rowSums(y[,1:2]), rowSums(y[,3:4]), rowSums(y[,5:6]), 
                  rowSums(y[,7:8]), rowSums(y[,9:10]), rowSums(y[,11:12]))
  for(k in 1:100){
    cp <- matrix(0, nrow = 6, ncol = 2)
    for(t in 1:6){
      cp[t,1] <- po[1,sy[k,t]]
      cp[t,2] <- po[2,sy[k,t]]
    }
    unc <- matrix(NA, nrow = 64, ncol = 6)
    unc[,1] <- log(potZ[,1] * psi0 + (1 - potZ[,1]) * (1 - psi0))
    for(t in 1:5){
      unc[,(t + 1)] <- log((1 - potZ[,t]) * (g[k,t] * potZ[,(t + 1)] + (1 - g[k,t]) * (1 - potZ[,(t + 1)])) +
                             potZ[,t] * (eps * (1-potZ[,(t + 1)]) + (1-eps) * potZ[,(t + 1)]))
    }
    UNC <- exp(rowSums(unc))
    tcp <- numeric()
    for(i in 1:64){
      temp <- rowSums(pZarray[i,,] * cp)
      if(sum(ifelse(temp == 0, 1, 0)) > 0){
        (tcp[i] <- 0)
      } else {
        tcp[i] <- UNC[i] * exp(sum(log(temp)))
      }
    }
    for(t in 1:5){
      out[k,t] <- sum(subset(tcp,potZ[,t] == 1)) / sum(tcp)
    }
  }
  return(out)
}

# this function calculates the likelihood after the first iteration and assuming
# cX is present in the workspace
fc_lik <- function(par){
  psi0 <- plogis(par[1])
  eps <- plogis(par[2])
  p <- plogis(par[3])
  g_int <- par[4]
  g_auto <- par[5]
  g <- matrix(0, nrow = 100, ncol = 5)
  for(t in 1:5){
    g[,t] <- plogis(g_int + g_auto * (cX[,t] %*% t(N)))
  }
  cp <- array(0, dim = c(100, 6, 2))
  zeta <- array(0, dim = c(100, 6, 2))
  po <- matrix(0, nrow = 2, ncol = 2)
  po[1,1] <- p
  po[1,2] <- (1 - p)
  po[2,1] <- 0
  po[2,2] <- 1
  for(t in 1:6){
    cp[,t,1] <- po[1,Y[,(2 * t-1)]] * po[1,Y[,(2 * t)]]
    cp[,t,2] <- po[2,Y[,(2 * t-1)]] * po[2,Y[,(2 * t)]]
  }
  zeta[,1,1] <- psi0 * cp[,1,1]
  zeta[,1,2] <- (1 - psi0) * cp[,1,2]
  for(t in 1:5){
    zeta[,(t + 1),1] <- (zeta[,t,1] * (1-eps) + (zeta[,t,2]) * g[,t]) * cp[,(t + 1),1]
    zeta[,(t + 1),2] <- (zeta[,t,1] * eps + (zeta[,t,2]) * (1-g[,t])) * cp[,(t + 1),2]
  }	
  nll <- 0
  for(k in 1:100){
    nll <- nll + -1 * log(sum(zeta[k,6,]))
  }
  nll
}

# function to fit first version of fully conditional and produce first
# conditional estimates
iter0 <- function(){
  m <- optim(c(0, 0, 0, 0), fc_lik0, method = "BFGS")
  tX <- matrix(0, nrow = 100, ncol = 5)
  CX <- fc_pred(c(m$par, 0), cX = tX)
  return(CX)
}
# End of fully conditional code - next comes function for forward conditional
#-----------------------------------------------------------------------------#
# likelihood for forward conditional
fwc_lik <- function(par){
  psi0 <- plogis(par[1])
  eps <- plogis(par[2])
  p <- plogis(par[3])
  g_int <- par[4]
  g_auto <- par[5]
  cp <- array(0, dim = c(100, 6, 2))
  zeta <- array(0, dim = c(100, 6, 2))
  po <- matrix(0, nrow = 2, ncol = 2)
  po[1,1] <- p
  po[1,2] <- (1 - p)
  po[2,1] <- 0
  po[2,2] <- 1
  for(t in 1:6){
    cp[,t,1] <- po[1,Y[,(2 * t - 1)]] * po[1,Y[,(2 * t)]]
    cp[,t,2] <- po[2,Y[,(2 * t - 1)]] * po[2,Y[,(2 * t)]]
  }
  zeta[,1,1] <- psi0 * cp[,1,1]
  zeta[,1,2] <- (1-psi0) * cp[,1,2]
  g <- matrix(0, nrow = 100, ncol = 5)
  psi <- matrix(0, nrow = 100, ncol = 5)
  for(t in 1:5){
    psi[,t] <- zeta[,t,1] / (zeta[,t,1] + zeta[,t,2])
    g[,t] <- plogis(g_int + g_auto * (psi[,t] %*% t(N)))
    zeta[,(t + 1),1] <- (zeta[,t,1] * (1 - eps) + (zeta[,t,2]) * g[,t]) * cp[,(t + 1),1]
    zeta[,(t + 1),2] <- (zeta[,t,1] * eps + (zeta[,t,2]) * (1 - g[,t])) * cp[,(t + 1),2]
  }	
  nll <- 0
  for(k in 1:100){
    nll <- nll + -1 * log(sum(zeta[k,6,]))
  }
  nll
}

#-----------------------------------------------------------------------------#
# code to simulate a dataset
sim <- function(psi0 = c(rep(0, 50), rep(1, 50)), g_int = -1, g_auto = 2,
                eps = 0.5, p = 0.5, Nvisits = 2, Nint = 5, neigh = N){
  Z <- matrix(0, ncol = (1 + Nint), nrow = length(psi0))
  Y <- matrix(0, ncol = (1 + Nint) * Nvisits, nrow = length(psi0))
  Z[,1] <- psi0
  for(k in 1:Nvisits){
    Y[,k] <- rbinom(100, 1, Z[,1] * p)
  }
  for(t in 1:Nint){
    x <- Z[,t] %*% t(neigh)
    gam <- plogis(g_int + g_auto * x)
    Z[,(t + 1)] <- rbinom(100, 1, Z[,t] * (1 - eps) + (1 - Z[,t]) * gam)
    for(k in 1:Nvisits){
      Y[,(t * Nvisits + k)] <- rbinom(100, 1, Z[, (t + 1)] * p)
    }
  }
  return(Y)
}

# create arrays to hold simulations
ss <- array(NA,dim = c(5,2,1000))  # full conditional - random occupancy
ss2 <- array(NA,dim = c(5,2,1000)) # forward conditional - random occupancy
ss3 <- array(NA,dim = c(5,2,1000)) # full conditional - clustered occupancy
ss4 <- array(NA,dim = c(5,2,1000)) # forward conditional - clustered occupancy
#-----------------------------------------------------------------------------#
for(j in 1:1000){
  set.seed(j) # set seed for simulations
  Y <- ifelse(sim(psi0 = rbinom(100, 1, 0.5)) == 0, 2, 1) # data where starting occupancy random
  # fit full conditional through iterative process
  cX <- iter0()
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  # usually converges after a few iterations, but we ran six times to ensure convergence
  ss[,1,j] <- m$par # store results
  ss[,2,j] <- sqrt(diag(solve(m$hessian)))
  # fit forward conditional and store results
  m2 <- optim(rep(0,5), fwc_lik, method = "BFGS", hessian = TRUE)
  ss2[,1,j] <- m2$par
  ss2[,2,j] <- sqrt(diag(solve(m2$hessian)))
  set.seed(j)
  # simulate clustered data
  Y <- ifelse(sim() == 0, 2, 1)
  # fit fully conditional and store
  cX <- iter0()
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  ss3[,1,j] <- m$par
  ss3[,2,j] <- sqrt(diag(solve(m$hessian)))
  # fit forward conditional and store
  m4 <- optim(rep(0,5), fwc_lik, method = "BFGS", hessian = TRUE)
  ss4[,1,j] <- m4$par
  ss4[,2,j] <- sqrt(diag(solve(m4$hessian)))
}

#-----------------------------------------------------------------------------#
sumz <- matrix(NA, ncol = 10, nrow = 4)
truth <- c(0, 0, 0, -1, 2)
for(k in 1:5){
  sumz[1,(2 * k - 1)] <- mean(ss[k,1,] - truth[k])
  sumz[2,(2 * k - 1)] <- mean(ss2[k,1,] - truth[k])
  sumz[3,(2 * k - 1)] <- mean(ss3[k,1,] - truth[k])
  sumz[4,(2 * k - 1)] <- mean(ss4[k,1,] - truth[k])
  sumz[1,(2 * k)] <- mean(ifelse(ss[k,1,] - 1.96 * ss[k,2,] < 
                                   truth[k] & ss[k,1,] + 1.96 * ss[k,2,] > truth[k], 1, 0))
  sumz[2,(2 * k)] <- mean(ifelse(ss2[k,1,] - 1.96 * ss2[k,2,] < 
                                   truth[k] & ss2[k,1,] + 1.96 * ss2[k,2,] > truth[k], 1, 0))
  sumz[3,(2 * k)] <- mean(ifelse(ss3[k,1,] - 1.96 * ss3[k,2,] < 
                                   truth[k] & ss3[k,1,] + 1.96 * ss3[k,2,] > truth[k], 1, 0))
  sumz[4,(2 * k)] <- mean(ifelse(ss4[k,1,] - 1.96 * ss4[k,2,] <
                                   truth[k]&ss4[k,1,] + 1.96 * ss4[k,2,] > truth[k], 1, 0))
}

#-----------------------------------------------------------------------------#
# coverage
cover <- array(NA, dim = c(4, 2, 2))
cover[1,1,1] <- sum(ifelse(ss[4,1,] + 0.67 * ss[4,2,] > (-1) &
                             ss[4,1,] - 0.67 * ss[4,2,] < (-1), 1, 0)) / 1000
cover[1,1,2] <- sum(ifelse(ss[4,1,] + 1.96 * ss[4,2,] > (-1) & 
                             ss[4,1,] - 1.96 * ss[4,2,] < (-1), 1, 0)) / 1000
cover[1,2,1] <- sum(ifelse(ss[5,1,] + 0.67 * ss[5,2,] > (2) & 
                             ss[5,1,] - 0.67 * ss[5,2,] < (2), 1, 0)) / 1000
cover[1,2,2] <- sum(ifelse(ss[5,1,] + 1.96 * ss[5,2,] > (2) & 
                             ss[5,1,] - 1.96 * ss[5,2,] < (2), 1, 0)) / 1000
#
cover[2,1,1] <- sum(ifelse(ss2[4,1,] + 0.67 * ss2[4,2,] > (-1) & 
                             ss2[4,1,] - 0.67 * ss2[4,2,] < (-1), 1, 0)) / 1000
cover[2,1,2] <- sum(ifelse(ss2[4,1,] + 1.96 * ss2[4,2,] > (-1) & 
                             ss2[4,1,] - 1.96 * ss2[4,2,] < (-1), 1, 0)) / 1000
cover[2,2,1] <- sum(ifelse(ss2[5,1,] + 0.67 * ss2[5,2,] > (2) & 
                             ss2[5,1,] - 0.67 * ss2[5,2,] < (2), 1, 0)) / 1000
cover[2,2,2] <- sum(ifelse(ss2[5,1,] + 1.96 * ss2[5,2,] > (2) & 
                             ss2[5,1,] - 1.96 * ss2[5,2,] < (2), 1, 0)) / 1000
#
cover[3,1,1] <- sum(ifelse(ss3[4,1,] + 0.67 * ss3[4,2,] > (-1) & 
                             ss3[4,1,] - 0.67 * ss3[4,2,] < (-1), 1, 0)) / 1000
cover[3,1,2] <- sum(ifelse(ss3[4,1,] + 1.96 * ss3[4,2,] > (-1) & 
                             ss3[4,1,] - 1.96 * ss3[4,2,] < (-1), 1, 0)) / 1000
cover[3,2,1] <- sum(ifelse(ss3[5,1,] + 0.67 * ss3[5,2,] > (2) & 
                             ss3[5,1,] - 0.67 * ss3[5,2,] < (2), 1, 0)) / 1000
cover[3,2,2] <- sum(ifelse(ss3[5,1,] + 1.96 * ss3[5,2,] > (2) & 
                             ss3[5,1,] - 1.96 * ss3[5,2,] < (2), 1, 0)) / 1000
#
cover[4,1,1] <- sum(ifelse(ss4[4,1,] + 0.67 * ss4[4,2,] > (-1) & 
                             ss4[4,1,] - 0.67 * ss4[4,2,] < (-1), 1, 0)) / 1000
cover[4,1,2] <- sum(ifelse(ss4[4,1,] + 1.96 * ss4[4,2,] > (-1) & 
                             ss4[4,1,] - 1.96 * ss4[4,2,] < (-1), 1, 0)) / 1000
cover[4,2,1] <- sum(ifelse(ss4[5,1,] + 0.67 * ss4[5,2,] > (2) & 
                             ss4[5,1,] - 0.67 * ss4[5,2,] < (2), 1, 0)) / 1000
cover[4,2,2] <- sum(ifelse(ss4[5,1,] + 1.96 * ss4[5,2,] > (2) & 
                             ss4[5,1,] - 1.96 * ss4[5,2,] < (2), 1, 0)) / 1000
#-----------------------------------------------------------------------------#
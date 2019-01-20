###############################################################################
#                                                                        Dec 18
#  Functions to:
#  1). process input data for model fitting 
#  2). process/extract results
#
#  Notes:
#  * 
#
#  To do: 
#  *
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

# known.state = known.state.cjs(CH)

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


# known.state = cjs.init.z(iCH, F)

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

# out = collapse.ch(ch)

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
# Notes: this got really slow... speed it up?  --> maybe lapply or mapply

# check out 'gelman.diag' in coda

# require packages for this....

# added 'ignore' argument, for parms that are set, but your still tracking, but
# want to exclude from the summary (see BNT model). Only for Stan and JAGS. 


run.times = function(fit.list, ignore = NULL){
  
  # ignore = c("IN[1]", "bN[1,1]")
  
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
      
      # print(get_elapsed_time(fit))
      # time = get_elapsed_time(fit.list[[i]])
      # out[i,]$run.time = max(time[,1] + time[,2])     # this is different than the system time
      # out[i,]$run.time = attr(fit.list[[i]], "time")  # this doesn't have to be fitter specific, moved below
      
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
  
  return(out)
}
# run.times = function(fit.list){
#   
#   out = data.frame(fitter = NA,
#                    iterations = NA,
#                    min.n.eff = NA,
#                    min.n.eff.coda = NA,
#                    med.n.eff = NA,
#                    med.n.eff.coda = NA,
#                    run.time = NA,
#                    model = NA,
#                    r.hat.count = NA)
#   
#   
#   for(i in seq_along(fit.list)){
#     
#     if(any(class(fit.list[[i]]) == "stanfit")){
#       
#       out[i,]$fitter = "Stan"
#       
#       # get the n.iter
#       out[i,]$iterations = fit.list[[i]]@stan_args[[1]]$iter
#       
#       # get the n.eff
#       n.eff = rstan::summary(fit.list[[i]])$summary[,"n_eff"]
#       
#       n.eff.coda = coda::effectiveSize(organize(fit.list[[i]], mcmc.out = TRUE))
#       
#       # minus 1 to cut out the likelihood value (only n.eff for parms)
#       out[i,]$min.n.eff = min(n.eff[1:(length(n.eff) - 1)])
#       
#       out[i,]$min.n.eff.coda = min(n.eff.coda[2:length(n.eff.coda)])  # deviance is first here
#       
#       # median n.eff
#       out[i,]$med.n.eff = median(n.eff[1:(length(n.eff) - 1)])
#       
#       out[i,]$med.n.eff.coda = median(n.eff.coda[2:length(n.eff.coda)])  # deviance is first here
#       
#       # print(get_elapsed_time(fit))
#       # time = get_elapsed_time(fit.list[[i]])
#       # out[i,]$run.time = max(time[,1] + time[,2])     # this is different than the system time
#       # out[i,]$run.time = attr(fit.list[[i]], "time")  # this doesn't have to be fitter specific, moved below
#       
#       # model name
#       out[i,]$model = fit.list[[i]]@model_name
#       
#       # count of Rhat over 1.1 (this includes the like/deviance)
#       tmp.r.hat = rstan::summary(fit.list[[i]])$summary[,"Rhat"]
#       
#       out[i,]$r.hat.count = length(which(tmp.r.hat > 1.1))
#       
#     }
#     
#     if(any(class(fit.list[[i]]) == "rjags")){  # had to add 'any' b/c jags in parallel has two classes, one of which will be 'rjags'
#       out[i,]$fitter = "JAGS"
#       
#       # get the n.iter
#       out[i,]$iterations = fit.list[[i]]$n.iter
#       
#       # get the n.eff
#       n.eff = fit.list[[i]]$BUGSoutput$summary[,9]
#       
#       n.eff.coda = coda::effectiveSize(organize(fit.list[[i]], mcmc.out = TRUE))
#       
#       # min n.eff - cut out the deviance value (only n.eff for parms)
#       out[i,]$min.n.eff = min(n.eff[2:length(n.eff)])
#       
#       out[i,]$min.n.eff.coda = min(n.eff.coda[2:length(n.eff.coda)])
#       
#       # median n, eff - cut out the deviance value (only n.eff for parms)
#       out[i,]$med.n.eff = median(n.eff[2:length(n.eff)])
#       
#       out[i,]$med.n.eff.coda = median(n.eff.coda[2:length(n.eff.coda)])
#       
#       # model name
#       out[i,]$model = fit.list[[i]]$model.file
#       
#       # count of Rhat over 1.1  (this includes the like/deviance)
#       tmp.r.hat = fit.list[[i]]$BUGSoutput$summary[,8]
#       
#       out[i,]$r.hat.count = length(which(tmp.r.hat > 1.1))
#       
#     }
#     
#     if(any(class(fit.list[[i]]) == "bugs")){
#       out[i,]$fitter = "bugs"
#       
#       # get the n.iter
#       out[i,]$iterations = fit.list[[i]]$n.iter
#       
#       # get the n.eff
#       n.eff = fit.list[[i]]$summary[,9]
#       
#       f1 = coda::mcmc.list(lapply(1:fit.list[[i]]$n.chain, function(x) coda::mcmc(fit.list[[i]]$sims.array[,x,])))
#       
#       n.eff.coda = coda::effectiveSize(f1)
#       
#       # min n.eff - cut out the deviance value (only n.eff for parms), different than JAGS, deviance is at the end
#       out[i,]$min.n.eff = min(n.eff[1:length(n.eff)-1])
#       
#       out[i,]$min.n.eff.coda = min(n.eff.coda[1:length(n.eff.coda)-1])
#       
#       # median n.eff
#       out[i,]$med.n.eff = median(n.eff[1:length(n.eff)-1])
#       
#       out[i,]$med.n.eff.coda = median(n.eff.coda[1:length(n.eff.coda)-1])
#       
#       # model name
#       out[i,]$model = fit.list[[i]]$model.file
#       
#       # count of Rhat over 1.1  (this includes the like/deviance)
#       tmp.r.hat = fit.list[[i]]$summary[,8]
#       
#       out[i,]$r.hat.count = length(which(tmp.r.hat > 1.1))
#       
#     }
#     
#     # get time to fit model, stored as an attribute
#     out[i,]$run.time = ifelse(is.null(attr(fit.list[[i]], "time")), "NA", attr(fit.list[[i]], "time")) 
#   }
#   
#   # out$efficiency = out$min.n.eff / out$run.time  
#   out$efficiency = out$min.n.eff.coda / out$run.time  
#   
#   return(out)
# }

# run.times(fit.list)

#-----------------------------------------------------------------------------#

# this is really slow....


# similar to run.times, run times returns the model run variables, except this
# version will return n.eff for all parms, not just the min

run.times.2 = function(fit.list){
  
  big.out = list()
  
  for(i in seq_along(fit.list)){
    
    if(any(class(fit.list[[i]]) == "stanfit")){
      
      parms = rownames(rstan::summary(fit.list[[i]])$summary)
      
      out = as.data.frame(matrix(ncol = 8, nrow = length(parms)))
      
      colnames(out) = c('fitter',
                        'iterations',
                        'n.eff',
                        'n.eff.2',
                        'run.time',
                        'model',
                        'r.hat',
                        'parm')
      
      out$parm = parms
      
      # get the n.eff
      out$n.eff = rstan::summary(fit.list[[i]])$summary[,"n_eff"]
      
      out$n.eff.2 = coda::effectiveSize(organize(fit.list[[i]], mcmc.out = TRUE))
      
      out$fitter = "Stan"
      
      # get the n.iter
      out$iterations = fit.list[[i]]@stan_args[[1]]$iter
      
      # model name
      out$model = fit.list[[i]]@model_name
      
      # count of Rhat over 1.1 (this includes the like/deviance)
      out$r.hat = rstan::summary(fit.list[[i]])$summary[,"Rhat"]
      
      out$run.time = ifelse(is.null(attr(fit.list[[i]], "time")), "NA", attr(fit.list[[i]], "time")) 
      
    }
    
    if(any(class(fit.list[[i]]) == "rjags")){  # had to add 'any' b/c jags in parallel has two classes, one of which will be 'rjags'
      
      parms = rownames(fit.list[[i]]$BUGSoutput$summary)
      
      out = as.data.frame(matrix(ncol = 8, nrow = length(parms)))
      
      colnames(out) = c('fitter',
                        'iterations',
                        'n.eff',
                        'n.eff.2',
                        'run.time',
                        'model',
                        'r.hat',
                        'parm')
      
      out$parm = parms
      
      # get the n.eff
      out$n.eff = fit.list[[i]]$BUGSoutput$summary[,9]
      
      out$n.eff.2 = coda::effectiveSize(organize(fit.list[[i]], mcmc.out = TRUE))
      
      out$fitter = "JAGS"
      
      # get the n.iter
      out$iterations = fit.list[[i]]$n.iter
      
      # model name
      out$model = fit.list[[i]]$model.file
      
      # count of Rhat over 1.1  (this includes the like/deviance)
      out$r.hat = fit.list[[i]]$BUGSoutput$summary[,8]
      
      out$run.time = ifelse(is.null(attr(fit.list[[i]], "time")), "NA", attr(fit.list[[i]], "time")) 
      
    }
    
    big.out[[i]] = out  
    
  }
  
  
  all.out = do.call('rbind', big.out)
  
  all.out$efficiency = all.out$n.eff / all.out$run.time  
  all.out$efficiency.2 = all.out$n.eff.2 / all.out$run.time  
  
  return(all.out)
}


# run.times.2(fit.list)  

#-----------------------------------------------------------------------------#


organize = function(fit, par.name, mcmc.out = FALSE){
  require(ggmcmc)
  
  # JAGS
  if(class(fit)[1] == "rjags"){
    if(mcmc.out == FALSE){
      tmp = coda::mcmc.list(lapply(1:fit$model$nchain(), function(x) coda::mcmc(fit$BUGSoutput$sims.array[,x,])))
      
      f1 = ggs(tmp, family = par.name)
      
    } else {
      
      f1 = coda::mcmc.list(lapply(1:fit$model$nchain(), function(x) coda::mcmc(fit$BUGSoutput$sims.array[,x,])))
    }
  }
  
  # JAGS Parallel
  if(class(fit)[1] == "rjags.parallel"){
    if(mcmc.out == FALSE){
      tmp = coda::mcmc.list(lapply(1:3, function(x) coda::mcmc(fit$BUGSoutput$sims.array[,x,])))
      
      f1 = ggs(tmp, family = par.name)
      
    } else {
      
      f1 = coda::mcmc.list(lapply(1:3, function(x) coda::mcmc(fit$BUGSoutput$sims.array[,x,])))
    }
  }
  
  # for a stan object
  if(any(class(fit) == "stanfit")){
    if(mcmc.out == FALSE){
      # looks like gss works directly with stan objects now...
      f1 <- ggs(fit, family = par.name)
      # used to do this... but with the new 'family' argument 
      # tmp <- coda::mcmc.list(lapply(1:ncol(fit), function(x) coda::mcmc(as.array(fit)[,x,])))
      # f1 <- ggs(tmp, family = par.name)
      
    } else {
      
      f1 = coda::mcmc.list(lapply(1:ncol(fit), function(x) coda::mcmc(as.array(fit)[,x,])))
      
    }
  }
  
  return(f1)
}


# organize(JD.out, par.name = "p")
# organize(SM.c, par.name = "p")

#-----------------------------------------------------------------------------#

to.cy.format <- function(dat.in) {
  
  dat.out = as.data.frame(matrix(ncol = 9, nrow = nrow(dat.in)))
  names(dat.out) = c("program",
                     "marginalized",
                     "run", 
                     "iterations",
                     "time",
                     "min(n.eff)", 
                     "median(n.eff)",
                     "min(n.eff).coda",
                     "median(n.eff).coda")
  
  dat.out$program = dat.in$fitter
  
  # figure this out...
  # dat.out$marginalized = ifelse(grepl("Marginalized", dat.in$med.n.eff))
  
  dat.out$run = c(1:nrow(dat.in))
  dat.out$iterations = dat.in$iterations
  dat.out$time = dat.in$run.time
  dat.out$`min(n.eff)` = dat.in$min.n.eff
  dat.out$`median(n.eff)` = dat.in$med.n.eff
  dat.out$`min(n.eff).coda` = dat.in$min.n.eff.coda
  dat.out$`median(n.eff).coda` = dat.in$med.n.eff.coda
  
  return(dat.out)
}
#-----------------------------------------------------------------------------#




###############################################################################
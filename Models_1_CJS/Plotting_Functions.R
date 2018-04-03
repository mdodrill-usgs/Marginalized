#-----------------------------------------------------------------------------#
#------------------- trace_plots ---------------------------------------------#

# to do: 
# 1). add in mean bars for right hand side plot 
# 2). for large model objects, figure out a better way to format the data!


# newer version that makes the correct names for parms in a matrix!
make.name = function(par.name, number){
  out = vector()
  if(is.null(number)){
    out = c(par.name)
  } else {
    if(is.list(number)){
      for(i in 1:length(number[[1]])){
        out[i] = c(paste(par.name, "[", number[[1]][i], ",", number[[2]][i], "]", sep = ""))
      }
    } else {
      for(i in 1:length(number)){
        out[i] = c(paste(par.name,"[",number[i],"]", sep = ""))
      }
    }
  }
  return(out)
}


trace_plots = function(fit, par.name, number = NULL, same.scale = FALSE){
  require(ggmcmc)
  require(gridExtra)
  require(ggthemes)
  
  s <- coda::mcmc.list(lapply(1:ncol(fit), function(x) coda::mcmc(as.array(fit)[,x,])))
  # this often thows a warning, don't think this is needed/meaningful, so I turned it off ;)
  s <- suppressWarnings(ggs(s))  
  
  names = make.name(par.name, number)
  
  ltl.s = s[which(s$Parameter %in% names),]
  
  attr(ltl.s, "nParameters") <- length(names)
  attr(ltl.s, "nChains") <- 3
  attr(ltl.s, "nThin") <- 1
  attr(ltl.s, "nBurnin") <- 0
  
  # get the Rhat values 
  r.tmp = rstan::summary(fit)$summary[,"Rhat"]
  r.hat = r.tmp[which(names(r.tmp) %in% names)]
  
  # get the n.eff
  n.tmp = rstan::summary(fit)$summary[,"n_eff"]
  n.eff = n.tmp[which(names(n.tmp) %in% names)]
  
  
  my.y = group_by(ltl.s, Parameter) %>%
    summarize(y.val = quantile(value, .9)) # y.val is the position to place the label (below)
  
  my.y2 = group_by(ltl.s, Parameter) %>%
    summarize(y.val = quantile(value, .2)) # y.val is the position to place the label (below)
  
  r.hat.d = data.frame(Parameter = names(r.hat),
                       Rhat = round(r.hat, 2),
                       Chain = rep(1,length(names)),
                       my.y = my.y$y.val)
  
  n.eff.d = data.frame(Parameter = names(n.eff),
                       Rhat = round(n.eff),
                       Chain = rep(1,length(names)),
                       my.y = my.y2$y.val)
  
  set.1 = ggs_traceplot(ltl.s) +
    geom_label(data = r.hat.d, x = 0, aes(y = my.y, label = Rhat), color = "#cb4b16", fill = "#002b36", size = 3) + # add Rhat values
    geom_label(data = n.eff.d, x = 0, aes(y = my.y, label = Rhat), color = "#cb4b16", fill = "#002b36", size = 3) + # add Rhat values
    
    
    theme_solarized(light = FALSE) +
    theme(strip.background = element_rect(fill = "#002b36", color = "#002b36"),
          strip.text = element_text(color = "#cb4b16"),
          legend.position = "none")
  
  if(same.scale == TRUE){
    set.2 <- ggplot(ltl.s, aes(x = value, colour = as.factor(Chain), fill = as.factor(Chain))) +
      geom_density(alpha = 0.3) + scale_fill_discrete(name = "Chain") + 
      scale_colour_discrete(name = "Chain") +
      facet_wrap(~Parameter, ncol = 1, scales = "free_y") +
      theme_solarized(light = FALSE) + 
      theme(strip.background = element_rect(fill = "#002b36", color = "#002b36"),
            strip.text = element_text(color = "#cb4b16"),
            legend.key = element_rect(fill = "#002b36", color = "#002b36"))
  } else {
    set.2 <- ggplot(ltl.s, aes(x = value, colour = as.factor(Chain), fill = as.factor(Chain))) +
      geom_density(alpha = 0.3) + scale_fill_discrete(name = "Chain") + 
      scale_colour_discrete(name = "Chain") +
      facet_wrap(~Parameter, ncol = 1, scales = "free") +
      theme_solarized(light = FALSE) + 
      theme(strip.background = element_rect(fill = "#002b36", color = "#002b36"),
            strip.text = element_text(color = "#cb4b16"),
            legend.key = element_rect(fill = "#002b36", color = "#002b36"))
  }
  
  grid.arrange(set.1, set.2, ncol = 2)
}


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

# two functions that make quick plots of model parameters for comparing across model runs

# to do: make par.name a list, so that the models can have different names for the same parm.. 


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




q_plot = function(mod.list, par.name){
  require(ggplot2)
  require(ggthemes)
  
  
  out = list()
  
  for(i in seq_along(mod.list)){
    
    out[[i]] = organize(mod.list[[i]], par.name = par.name)
    
    out[[i]]$model = rep(names(mod.list)[i])
    
  }
  
  out.all = do.call('rbind', out)
  
  
  all2 = group_by(out.all, Parameter, model) %>%
    summarize(my.mean = mean(value),
              upper = quantile(value, .95),
              lower = quantile(value, .05))
  
  p = ggplot(all2, aes(y = my.mean, x = Parameter)) +
    geom_point(position = position_dodge(width = 1), aes(color = model), size = 3) +
    geom_errorbar(aes(ymin = lower, ymax = upper, color = model),
                  position = position_dodge(width = 1),
                  width = 0, size = 1) +
    labs(y = "Mean +- 95% CRI")
  
  p + theme_base() + theme(legend.position = c(.15,.9))
  
}



# mod.list = list("name.1" = JM.out, "name.2" = SM.c)
# mod.list = list("name.1" = JM.out, "name.2" = SM.c, "the best" = SM.c)
# q_plot(mod.list, par.name = "p")


#-----------------------------------------------------------------------------#

###############################################################################
#                                                                        Oct 18
#  Fitting a Multi-state version of a CJS model to the RBT data 
#  Marginalized Stan version
#
#  Notes:
#  * 
#
#  To do: 
#  *  Clean up script...
#  *  check to see if this is setting the core
#     options(affinity.list = c(1,2,3,4))     
#
###############################################################################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())  

source(paste0(getwd(),"/Functions.R"), chdir = F)

setwd(paste0(getwd(), '/Application_1'))

data.dir = paste0(getwd(), "/Data")
CH = as.matrix(read.table(file = paste0(data.dir, "/RBT_Capture_History.txt"),
                          header = FALSE, sep = "\t"))
#-----------------------------------------------------------------------------#
# format data for model fitting
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

# Catch (for the N versions)
catch = colSums(CH)

#-----------------------------------------------------------------------------#
# working with the Stan model in a seperate tab... see "Stan_M...stan"

# sink("Stan_M...stan")
# cat("
#     

# Put the model code here...

#     ",fill=TRUE)
# sink()    

#-----------------------------------------------------------------------------#
# Constant 
sm.params <- c("s", "p")

sm.data <- list(NsumCH = NsumCH, n_occasions = n.occasions, sumCH = sumCH,
                sumf = sumf, sumFR = sumFR)

# MCMC settings
ni = 1000
nt = 1
nb = 500
nc = 3


# Call Stan from R 
# t1 <- proc.time()
SM.c <- stan("Stan_Marginalized_Constant.stan",
           # SM <- stan(fit = SM,
           data = sm.data,
           pars = sm.params,
           control = list(adapt_delta = .85),
           # control = list(max_treedepth = 14),
           # control = list(max_treedepth = 15, adapt_delta = .925),
           # chains = nc, iter = ni, warmup = nb, thin = nt, seed = 1) 
           chains = nc, iter = ni, thin = nt, seed = 1) 
# t2 <- proc.time()

SM.c

#-----------------------------------------------------------------------------#
# Constant with abundance 
sm.params <- c("s", "p", "N")

sm.data <- list(NsumCH = NsumCH, n_occasions = n.occasions, sumCH = sumCH,
                sumf = sumf, sumFR = sumFR, Catch = catch)

# MCMC settings
ni = 1000
nt = 1
nb = 500
nc = 3


# Call Stan from R 
# t1 <- proc.time()
SM.c2 <- stan("Stan_Marginalized_Constant_N.stan",
           # SM <- stan(fit = SM,
           data = sm.data,
           pars = sm.params,
           control = list(adapt_delta = .85),
           # control = list(max_treedepth = 14),
           # control = list(max_treedepth = 15, adapt_delta = .925),
           # chains = nc, iter = ni, warmup = nb, thin = nt, seed = 1) 
           chains = nc, iter = ni, thin = nt, seed = 1) 
# t2 <- proc.time()

SM.c2

trace_plots(SM.c2, "N", 1:6, same.scale = TRUE)


#-----------------------------------------------------------------------------#
# Time (occasion) specific parms
sm.params <- c("s", "p")

sm.data <- list(NsumCH = NsumCH, n_occasions = n.occasions, sumCH = sumCH,
                sumf = sumf, sumFR = sumFR)

# MCMC settings
ni = 1000
nt = 1
nb = 500
nc = 3


# Call Stan from R 
# t1 <- proc.time()
# SM.t <- stan("Stan_Marginalized.stan",
# SM.t2 <- stan("Stan_Marginalized.stan",
SM.t3 <- stan("Stan_Marginalized_Time.stan",
           # SM <- stan(fit = SM,
           data = sm.data,
           pars = sm.params,
           control = list(adapt_delta = .85),
           # control = list(max_treedepth = 14),
           # control = list(max_treedepth = 15, adapt_delta = .925),
           # chains = nc, iter = ni, warmup = nb, thin = nt, seed = 1) 
           chains = nc, iter = ni, thin = nt) 
# t2 <- proc.time()

SM.t
SM.t2
SM.t3



q_plot(list("none specified" = SM.t, "logit(p)~normal(0,2)" = SM.t2, "beta(.1,.1)" = SM.t3), par.name = "p")


#-----------------------------------------------------------------------------#
# Time (occasion) specific parms with abundance
sm.params <- c("s", "p", "N")

sm.data <- list(NsumCH = NsumCH, n_occasions = n.occasions, sumCH = sumCH,
                sumf = sumf, sumFR = sumFR, Catch = catch)

# MCMC settings
ni = 1000
nt = 1
nb = 500
nc = 3


# Call Stan from R 
# t1 <- proc.time()
SM.t2 <- stan("Stan_Marginalized_Time_N.stan",
           # SM <- stan(fit = SM,
           data = sm.data,
           pars = sm.params,
           control = list(adapt_delta = .85),
           # control = list(max_treedepth = 14),
           # control = list(max_treedepth = 15, adapt_delta = .925),
           # chains = nc, iter = ni, warmup = nb, thin = nt, seed = 1) 
           chains = nc, iter = ni, thin = nt) 
# t2 <- proc.time()

SM.t2
trace_plots(SM.t2, "N", 1:6, same.scale = TRUE)
trace_plots(SM.t2, "N", 7:13)

#-----------------------------------------------------------------------------#
# variational inference

# m = stan_model("Stan_Marginalized.stan")
# 
# v.SM <- vb(m, data = sm.data, pars = sm.params) 
# 
# mod.list = list("regular" = SM.t3, "V" = v.SM)
# q_plot(mod.list, par.name = "p")

#-----------------------------------------------------------------------------#
# p & s as random effects
# sm.params <- c("s", "p", "mu_s", "sig_s", "mu_p", "sig_p")
sm.params <- c("s", "p", "mean_s", "mean_p")

sm.data <- list(NsumCH = NsumCH, n_occasions = n.occasions, sumCH = sumCH,
                sumf = sumf, sumFR = sumFR)

# MCMC settings
ni = 1000
nt = 1
nb = 500
nc = 3


# Call Stan from R 
# t1 <- proc.time()
SM.re <- stan("Stan_Marginalized_RE.stan",
           # SM <- stan(fit = SM.re,
           data = sm.data,
           pars = sm.params,
           control = list(adapt_delta = .85),
           # control = list(max_treedepth = 14),
           # control = list(max_treedepth = 15, adapt_delta = .925),
           chains = nc, iter = ni, warmup = nb, thin = nt, seed = 1) 
# t2 <- proc.time()

SM.re


q_plot(list("RE" = SM.re, "constant" = SM.c), par.name = "s")


#-----------------------------------------------------------------------------#
# p & s as random effects, with abundance
# sm.params <- c("s", "p", "mu_s", "sig_s", "mu_p", "sig_p")
sm.params <- c("s", "p", "mean_s", "mean_p", "N")

sm.data <- list(NsumCH = NsumCH, n_occasions = n.occasions, sumCH = sumCH,
                sumf = sumf, sumFR = sumFR, Catch = catch)

# MCMC settings
ni = 100
nt = 1
nb = 50
nc = 1


# Call Stan from R 
# t1 <- proc.time()
SM.re2 <- stan("Stan_Marginalized_RE_N.stan",
           # SM <- stan(fit = SM.re,
           data = sm.data,
           pars = sm.params,
           control = list(adapt_delta = .85),
           # control = list(max_treedepth = 14),
           # control = list(max_treedepth = 15, adapt_delta = .925),
           chains = nc, iter = ni, warmup = nb, thin = nt, seed = 1) 
# t2 <- proc.time()

SM.re2

trace_plots(SM.re2, "N", 1:6, same.scale = TRUE)

q_plot(list("dot" = SM.c2, "time" = SM.t2, "RE" = SM.re2), par.name = "N")
#-----------------------------------------------------------------------------#

Rhat = rstan::summary(SM)$summary[,"Rhat"]

bad = Rhat[which(Rhat > 1.2)]
length(bad)
length(Rhat)


windows(xpos = 25, record = T)


trace_plots(SM.re, "s", 1:6, same.scale = TRUE)
trace_plots(SM.re, "s", 7:12, same.scale = TRUE)
trace_plots(SM.re, "s", 13:18, same.scale = TRUE)
trace_plots(SM.re, "p", 1:6, same.scale = TRUE)
trace_plots(SM.re, "p", 7:12, same.scale = TRUE)
trace_plots(SM.re, "p", 13:18, same.scale = TRUE)
trace_plots(SM.re, "mu_p", 1:1)
trace_plots(SM.re, "mu_s", 1:1)

#-----------------------------------------------------------------------------#


q_plot(list("constant" = SM.c, "little t" = SM.t), par.name = "p")
q_plot(list("constant" = SM.c, "little t" = SM.t), par.name = "s")

q_plot(list("constant" = SM.c, "little t" = SM.t, "RE" = SM.re), par.name = "p")
q_plot(list("constant" = SM.c, "little t" = SM.t, "RE" = SM.re), par.name = "s")




#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

n.core = 8 # really n.core * 3

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used 
cllibs <- clusterEvalQ(cl1, c(library(rstan)))

n.runs = 10

# nb = 50
nt = 1
nc = 3

# my.n.iter = c(500, 550)
my.n.iter = seq(0,10000,500)[- 1]
# my.n.iter = seq(5500,10000,500)
# my.n.iter = rep(2000,5)

big.fit.list = list()

seeds <- sample(1:1e5, size = n.runs)   
# seeds <- sample(1:1e5, size = n.runs*length(my.n.iter))   # change how the seed is done, new for every run !

# let all of the clusters have whatever objects are in the workspace
clusterExport(cl = cl1, varlist = ls(), envir = environment())

start.time <- Sys.time()  # start timer

# n = 1
out = list()
#--------------------------------------
out <- foreach(j = seeds) %:% 
  
  foreach(i = my.n.iter) %dopar% {
    
    seed = j
    seed = seed + sample(1:1e5, size = 1)
    
    iter.in = i
    
    t1 <- proc.time()
    
    # my.env = environment()
    
    # JM.re.c = jags.parallel(JM.data, inits = NULL, JM.par, "rbt_JAGS_M_RE.txt",
    #                         n.chains = 3, n.iter = iter.in, export_obj_names = c("iter.in", "seed"),
    #                         jags.seed = seed, envir = my.env)
    
    # SM.re <- stan("Stan_Marginalized_RE.stan",
    # SM.re <- stan(fit = SM.re,     # faster to compile the model above...
    #               data = sm.data,
    #               pars = sm.params,
    #               control = list(adapt_delta = .85),
    #               chains = nc, iter = iter.in, warmup = nb, thin = nt, seed = seed, cores = 3) 
    
    
    # SM.c <- stan(fit = SM.c,     # faster to compile the model above...
    # SM.c <- stan("Stan_Marginalized_Constant.stan",     
    #              data = sm.data,
    #              pars = sm.params,
    #              control = list(adapt_delta = .85),
    #              # chains = nc, iter = iter.in, warmup = nb, thin = nt, seed = seed, cores = 3)
    #              chains = nc, iter = iter.in, thin = nt, seed = seed, cores = 3)
    
    
    SM.t <- stan("Stan_Marginalized_Time.stan",
                 data = sm.data,
                 pars = sm.params,
                 # control = list(adapt_delta = .85),
                 control = list(adapt_enga),
                 # chains = nc, iter = iter.in, warmup = nb, thin = nt, seed = seed, cores = 3)
                 chains = nc, iter = iter.in, thin = nt, seed = seed, cores = 3)
    
    
    
    t2 <- proc.time()
    
    attr(SM.t, 'time') <- (t2 - t1)[3]
    
    SM.t
    
  } 
#--------------------------------------


end.time = Sys.time()
time.taken = end.time - start.time
print(round(time.taken,2))

stopCluster(cl1)  # close the clusters


length(out)

length(out[[1]])


all.out = do.call('c', out)
length(all.out)

tmp = run.times(all.out)

all.stan.m.time = all.out
rm(list=setdiff(ls(), c("all.stan.m.time")))



windows()


p = ggplot(tmp, aes(x = iterations, y = efficiency)) +
    geom_point() +
    geom_boxplot(aes(group = iterations), alpha = 1)
p

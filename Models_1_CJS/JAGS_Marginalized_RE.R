###############################################################################
#                                                                        Feb 18
#        Fitting a state - space version of a CJS model to the RBT data 
#        Marginalized JAGS version - Random effects for p & s
#
#  Notes:
#  * 
#
#  To do: 
#  * Clean up parallel piece at bottom...
#
###############################################################################
setwd(paste0(getwd(), '/Models_1_CJS'))
library(R2jags)

source("RBT_Functions.R", chdir = F)

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

ones <- sumFR

pz = known.state.cjs(tmpCH)

#-----------------------------------------------------------------------------#
sink("JAGS_Marginalized_RE.jags")
cat("
model {
  
  mean.s ~ dunif(0, 1)                     # Prior for mean survival
  mu.s <- log(mean.s / (1 - mean.s))         # Logit transformation
  sigma.s ~ dunif(0, 10)                     # Prior for standard deviation
  tau.s <- pow(sigma.s, -2)
  
  mean.p ~ dunif(0, 1)                     # Prior for mean survival
  mu.p <- log(mean.p / (1 - mean.p))         # Logit transformation
  sigma.p ~ dunif(0, 10)                     # Prior for standard deviation
  tau.p <- pow(sigma.p, -2)
   
  
  for (k in 1:(n.occasions - 1)){
    eps.s[k] ~ dnorm(0, tau.s)
    eps.p[k] ~ dnorm(0, tau.p)
  }
  
  for(k in 1:(n.occasions - 1)){
    
    # s[k] ~ dunif(0,1)
    
    logit(s[k]) <- mu.s + eps.s[k]
    
    tr[1,k,1] <- s[k]
    tr[1,k,2] <- (1 - s[k])
    tr[2,k,1] <- 0
    tr[2,k,2] <- 1
    
    # p[k] ~ dunif(0,1)
    
    logit(p[k]) <- mu.p + eps.p[k]
    
    pmat[1,k,1] <- p[k]
    pmat[1,k,2] <- (1 - p[k])
    pmat[2,k,1] <- 0
    pmat[2,k,2] <- 1
  }
  
  for(i in 1:NsumCH){
    pz[i,sumf[i],1] <- 1
    pz[i,sumf[i],2] <- 0
    
    for(k in (sumf[i] + 1):n.occasions){
      pz[i,k,1] <- pz[i,k-1,1] * tr[1,k-1,1] * pmat[1,k-1,sumCH[i,k]]
      pz[i,k,2] <- inprod(pz[i,k-1,], tr[,k-1,2]) * pmat[2,k-1,sumCH[i,k]]
    }
    
    like[i] <- sum(pz[i,n.occasions,])
    ones[i] ~ dbin(like[i], sumFR[i])
  }
    }
    ", fill = TRUE)
sink() 
#-----------------------------------------------------------------------------#

# JM.inits <- function(){list(s = runif((n.occasions - 1), 0, 1),
#                             p = runif((n.occasions - 1), 0, 1),
#                             z = cjs.init.z(CH, indf))}

JM.data <- list(NsumCH = NsumCH, n.occasions = n.occasions, sumCH = sumCH,
                sumf = sumf, sumFR = sumFR, ones = ones)

JM.par <- c('s', 'p', 'mean.s', 'mean.p')

ni <- 100
nt <- 1
nb <- 50


# t1 <- proc.time()
JM.re <- jags(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_RE.jags",
               n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb)


JM.re.b <- jags(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_RE.jags",
               n.chains = 3, n.iter = ni)


JM.re.c <- jags.parallel(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_RE.jags",
               n.chains = 3, n.iter = ni, export_obj_names = c("ni"))   


# t2 <- proc.time()

print(JM.re, digits = 3)
print(JM.re.b, digits = 3)
print(JM.re.c, digits = 3)

#-----------------------------------------------------------------------------#

q_plot(list("model" = JM.re), par.name = "s")
q_plot(list("model" = JM.re), par.name = "p")

#-----------------------------------------------------------------------------#
 



#-----------------------------------------------------------------------------#

# nb <- 500
all.t1 = proc.time()
n.runs = 5
seeds <- sample(1:1e5, size = n.runs)   

# n.iter = seq(1000,5000,500)
# n.iter = seq(0,5000,500)[- 1]
n.iter = seq(0,2500,500)[- 1]
# n.iter = c(100, 200)

big.fit.list = list()


for(j in 1:n.runs){
  
  seed = seeds[j]
  
  fit.list = list()
  
  for(i in seq_along(n.iter)){
    ni = n.iter[i]
    
    t1 <- proc.time()
    # JM.re <- jags(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_RE.jags",
    #               n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb)
    
    # JM.re.b <- jags(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_RE.jags",
    #                 n.chains = 3, n.iter = ni)
    
    JM.re.c <- jags.parallel(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_RE.jags",
                             n.chains = 3, n.iter = ni, export_obj_names = c("ni"),
                             jags.seed = seed)   
    
    t2 <- proc.time()
    
    attr(JM.re.c, 'time') <- (t2 - t1)[3]
    
    fit.list[[i]] = JM.re.c
  }  
  
  big.fit.list[[j]] = fit.list
}

all.t2 = proc.time()

(all.t2 - all.t1)[3]



# big.list = list()
# 
# big.list[[1]] = fit.list
# big.list[[2]] = fit.list

all.jags = do.call(c, big.fit.list)


# rm(list=setdiff(ls(), "all.jags"))
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

n.core = 8

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used
# cllibs <- clusterEvalQ(cl1, c(library(lme4), library(plyr), library(simr)))
cllibs <- clusterEvalQ(cl1, c(library(R2jags)))
# stopCluster(cl1)


# start.time <- Sys.time()  # start timer


# nb <- 500
all.t1 = proc.time()
n.runs = 1
seeds <- sample(1:1e5, size = n.runs)

# n.iter = seq(1000,5000,500)
n.iter = seq(0,5000,500)[- 1]
# n.iter = c(100, 200)
# n.iter = seq(0,2500,500)[- 1]

big.fit.list = list()


# let all of the clusters have whatever objects are in the workspace
# clusterExport(cl = cl1, varlist = ls())
clusterExport(cl = cl1, varlist = ls(), envir = environment())

start.time <- Sys.time()  # start timer
# the foreach is running the loop in parallel
out = foreach(j = 1:n.runs) %dopar% {
# out = foreach(j = 1:length(n.iter)) %dopar% {
  # for(j in 1:n.runs){

  seed = seeds[j]

  fit.list = list()

  # ni = n.iter[j]

  for(i in seq_along(n.iter)){
  # for(i in 1:n.runs){
    ni = n.iter[i]

    # seed = seeds[i]

    t1 <- proc.time()

    JM.re.c <- jags.parallel(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_RE.jags",
                             n.chains = 3, n.iter = ni, export_obj_names = c("ni"),
                             jags.seed = seed)

    t2 <- proc.time()

    attr(JM.re.c, 'time') <- (t2 - t1)[3]

    fit.list[[i]] = JM.re.c
  }

  big.fit.list[[j]] = fit.list
  # }
}

end.time = Sys.time()
time.taken = end.time - start.time
# message(print(round(time.taken,2)))
print(round(time.taken,2))

all.t2 = proc.time()
stopCluster(cl1)  # close the clusters


# big.list = list()
# 
# big.list[[1]] = fit.list
# big.list[[2]] = fit.list

# all.jags = do.call(c, big.fit.list)


# rm(list=setdiff(ls(), "all.jags"))


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

n.core = 8

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used 
# cllibs <- clusterEvalQ(cl1, c(library(lme4), library(plyr), library(simr)))
cllibs <- clusterEvalQ(cl1, c(library(R2jags)))
# stopCluster(cl1)


# start.time <- Sys.time()  # start timer


# nb <- 500
all.t1 = proc.time()
n.runs = 1

# seeds <- sample(1:1e5, size = n.runs)   

# n.iter = seq(1000,5000,500)
n.iter = seq(0,5000,500)[- 1]
# my.n.iter = c(100)
# n.iter = seq(0,2500,500)[- 1]

big.fit.list = list()

seeds <- sample(1:1e5, size = n.runs*length(my.n.iter))   # change how the seed is done, new for every run !


# let all of the clusters have whatever objects are in the workspace
# clusterExport(cl = cl1, varlist = ls())
# clusterExport(cl = cl1, varlist = ls(), envir = environment())
clusterExport(cl = cl1, varlist = ls())

start.time <- Sys.time()  # start timer
# the foreach is running the loop in parallel
# out = foreach(j = 1:n.runs) %:% {
# out = foreach(j = n.runs) %:% {

n = 1
# out <- foreach(j = 1:n.iter, .combine = 'c') %:% 
out <- foreach(j = n.runs) %:% {
  # foreach(j = n.runs, .combine = 'c') %:% 
  
  # for(j in 1:n.runs){
  
  # seed = seeds[j]
  
  fit.list = list()
  
  # ni = n.iter[j]
  
  # foreach(i = seq_along(n.iter),  .combine = 'c') %dopar% {
  # foreach(i = my.n.iter, .combine = 'c') %dopar% {
  foreach(i = my.n.iter) %dopar% {
  # foreach(i = n.iter) %dopar% {
    # for(i in 1:n.runs){
    # ni = n.iter[i]
    
    seed = seeds[n]
    n = n + 1
    
    # seed = seeds[i]
    
    # t1 <- proc.time()
    
    # print(my.n.iter)
    
    # iter.in = my.n.iter
    # print(iter.in)
    
    # print(i)
    
    iter.in = i
    
    my.env = environment()
    
    # ls(my.env)
    
    # print(iter.in)
    
    # JM.re.c <- jags.parallel(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_RE.jags",
    #                          n.chains = 3, n.iter = ni, export_obj_names = c("ni"),
    #                          jags.seed = seed)
    # jags.parallel(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_RE.jags",
    #               n.chains = 3, n.iter = iter.in, export_obj_names = c("iter.in"),
    #               jags.seed = seed, envir = my.env)
    JM.re.c <- jags.parallel(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_RE.jags",
                             n.chains = 3, n.iter = iter.in, export_obj_names = c("iter.in"),
                             jags.seed = seed, envir = my.env)
    
    # JM.re.b <- jags(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_RE.jags",
    #                 n.chains = 1, n.iter = iter.in)

    
    
    # t2 <- proc.time()
    # 
    # attr(JM.re.c, 'time') <- (t2 - t1)[3]
    # 
    # fit.list[[i]] = JM.re.c
    # 
    JM.re.c
  }  
   # out
  
  # big.fit.list[[j]] = fit.list
  }


end.time = Sys.time()
time.taken = end.time - start.time
# message(print(round(time.taken,2)))
print(round(time.taken,2))

all.t2 = proc.time()
stopCluster(cl1)  # close the clusters


length(out)
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

n.core = 8  # really n.core * 3

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used 
cllibs <- clusterEvalQ(cl1, c(library(R2jags)))

n.runs = 10

# my.n.iter = c(500, 550)
my.n.iter = seq(5500,10000,500)

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
    
    iter.in = i

    t1 <- proc.time()
    
    my.env = environment()

     JM.re.c = jags.parallel(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_RE.jags",
                             n.chains = 3, n.iter = iter.in, export_obj_names = c("iter.in", "seed"),
                             jags.seed = seed, envir = my.env)

     t2 <- proc.time()

     attr(JM.re.c, 'time') <- (t2 - t1)[3]
     
     JM.re.c
     
     } 
#--------------------------------------
   

end.time = Sys.time()
time.taken = end.time - start.time
print(round(time.taken,2))

all.t2 = proc.time()
stopCluster(cl1)  # close the clusters


length(out)

length(out[[1]])


all.out = do.call('c', out)
length(all.out)


# all.jags.m.2 = all.out
# 
# rm(list=setdiff(ls(), "all.jags.m.2"))



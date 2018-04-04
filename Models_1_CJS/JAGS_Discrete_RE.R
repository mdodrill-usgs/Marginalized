###############################################################################
#                                                                        Feb 18
#        Fitting a state - space version of a CJS model to the RBT data 
#        Discrete JAGS version - Random effects for p & s
#
#  Notes:
#  * 
#
#  To do: 
#  * test this on the server...
#
###############################################################################
setwd(paste0(getwd(), '/Models_1_CJS'))
library(R2jags)

source("RBT_Functions.R", chdir = F)

# format data for model fitting
indCH = CH

indCH[indCH[,] == 0] = 2

NindCH = nrow(indCH)         # number of capture histories 
n.occasions = ncol(indCH)    # number of sampling occasions

# Create vector with occasion of marking
get.first <- function(x) min(which(x != 0))
indf <- apply(CH, 1, get.first)

z = known.state.cjs(CH)

#-----------------------------------------------------------------------------#
sink("JAGS_Discrete_RE.jags")
cat("
model {
  
  mean.s ~ dunif(0, 1)                       # Prior for mean survival
  mu.s <- log(mean.s / (1 - mean.s))         # Logit transformation
  sigma.s ~ dunif(0, 10)                     # Prior for standard deviation
  tau.s <- pow(sigma.s, -2)
  
  mean.p ~ dunif(0, 1)                       # Prior for mean survival
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
  
  for(i in 1:NindCH){
    z[i, indf[i]] <- 1
    
    for(k in (indf[i] + 1):n.occasions){
      
      z[i,k] ~ dcat(tr[z[i, k-1], k-1, ])
      
      indCH[i,k] ~ dcat(pmat[z[i, k], k-1, ])
    }
  }
}
    
    ", fill = TRUE)
sink()    
#-----------------------------------------------------------------------------#

# JD.inits <- function(){list(s = runif((n.occasions - 1), 0, 1),
#                             p = runif((n.occasions - 1), 0, 1),
#                             z = cjs.init.z(CH, indf))}

JD.data <- list(NindCH = NindCH, n.occasions = n.occasions, indCH = indCH, indf = indf, z = z)
JD.par <- c('s', 'p', 'mean.s', 'mean.p')

ni <- 10
nt <- 1
nb <- 5

t1 <- proc.time()
JD.re <- jags(JD.data, inits = NULL, JD.par, "JAGS_Discrete_RE.jags",
               n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb)
t2 <- proc.time()

print(JD.re, digits = 3)

# 15 mins to compile on server, wtf?

#-----------------------------------------------------------------------------#
###############################################################################
#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

n.core = 8  # really n.core * 3

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used 
cllibs <- clusterEvalQ(cl1, c(library(R2jags)))

all.t1 = proc.time()
n.runs = 10

# my.n.iter = c(10, 15)
# my.n.iter = seq(0,5000,500)[- c(1,6:10)]
# my.n.iter = c(2500,3000,3500,4000,4500)
my.n.iter = c(7000, 7500, 8000, 8500, 9000, 9500, 10000)


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
    
    # name is wrong here...fix
    JD.re <- jags.parallel(JD.data, inits = NULL, JD.par, "JAGS_Discrete_RE.jags",
                             n.chains = 3, n.iter = iter.in, export_obj_names = c("iter.in", "seed"),
                             jags.seed = seed, envir = my.env)   
    
    t2 <- proc.time()
    
    attr(JD.re, 'time') <- (t2 - t1)[3]
    
    JD.re
    
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



# all.jags.d.re.2 = all.out
# rm(list=setdiff(ls(), "all.jags.d.re.2"))





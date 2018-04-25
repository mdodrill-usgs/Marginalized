###############################################################################
#                                                                      April 18
#        Fitting a state - space version of a CJS model to the RBT data 
#        Discrete JAGS version
#
#  Notes:
#  * Abundance (N) is computed with a Negative Binomial distribution 
#
#  To do: 
#  * 
#
###############################################################################
setwd(paste0(getwd(), '/Models_1_CJS'))
library(R2jags)

source('RBT_Functions.R', chdir = F)

# format data for model fitting
indCH = CH

indCH[indCH[,] == 0] = 2

NindCH = nrow(indCH)         # number of capture histories 
n.occasions = ncol(indCH)    # number of sampling occasions

# Create vector with occasion of marking
get.first <- function(x) min(which(x != 0))
indf <- apply(CH, 1, get.first)

z = known.state.cjs(CH)

catch = colSums(CH)

#-----------------------------------------------------------------------------#
sink("JAGS_Discrete_Constant_N.jags")
cat("
model {
  
  s ~ dunif(0,1)
  tr[1,1] <- s
  tr[1,2] <- (1 - s)
  tr[2,1] <- 0
  tr[2,2] <- 1
  
  p ~ dunif(0,1)
  pmat[1,1] <- p
  pmat[1,2] <- (1 - p)
  pmat[2,1] <- 0
  pmat[2,2] <- 1
  
  for(i in 1:NindCH){
    z[i, indf[i]] <- 1
    
    for(k in (indf[i] + 1):n.occasions){
      
      z[i,k] ~ dcat(tr[z[i, k-1], ])
      
      indCH[i,k] ~ dcat(pmat[z[i, k], ])
    }
  }
  
  #------- Calculate abundance using Negative Binomial 
  # this is the p-cap 
  # ptrans <- p
  
  for(i in 1:n.occasions){
    U[i] ~ dnegbin(p, catch[i])
    N[i] <- U[i] + catch[i]
  }
}
    
    ", fill = TRUE)
sink()    
#-----------------------------------------------------------------------------#


JD.data <- list(NindCH = NindCH, n.occasions = n.occasions, indCH = indCH,
                indf = indf, z = z, catch = catch)
JD.par <- c('s', 'p', 'N')

ni <- 200
nt <- 1
nb <- 100

# t1 <- proc.time()
JD.out <- jags(JD.data, inits = NULL, JD.par, "JAGS_Discrete_Constant_N.jags",
               n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb)
# t2 <- proc.time()

print(JD.out, digits = 3)

#-----------------------------------------------------------------------------#
library(coda)
 
samps <- coda.samples(JD.out)

#-----------------------------------------------------------------------------#
###############################################################################
#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

n.core = 9  # really n.core * 3

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used 
cllibs <- clusterEvalQ(cl1, c(library(R2jags)))

all.t1 = proc.time()
n.runs = 10

# my.n.iter = c(10, 15)
my.n.iter = seq(0,10000,500)[- 1]
# my.n.iter = c(2500,3000,3500,4000,4500)

big.fit.list = list()

seeds <- sample(1:1e5, size = n.runs)   

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
    
    JD.c <- jags.parallel(JD.data, inits = NULL, JD.par, "JAGS_Discrete_Constant_N.jags",
                          n.chains = 3, n.iter = iter.in, export_obj_names = c("iter.in", "seed"),
                          jags.seed = seed, envir = my.env) 
    
    t2 <- proc.time()
    
    attr(JD.c, 'time') <- (t2 - t1)[3]
    
    JD.c
    
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



# all.jags.d.constant.1 = all.out

# rm(list=setdiff(ls(), "all.jags.d.constant.1"))








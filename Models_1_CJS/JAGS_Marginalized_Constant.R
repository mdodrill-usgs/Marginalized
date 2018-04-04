###############################################################################
#                                                                        Feb 18
#        Fitting a state - space version of a CJS model to the RBT data 
#        Marginalized JAGS version
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
sink("JAGS_Marginalized_Constant.jags")
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
  
  for(i in 1:NsumCH){
    pz[i,sumf[i],1] <- 1
    pz[i,sumf[i],2] <- 0
    
    for(k in (sumf[i] + 1):n.occasions){
      pz[i,k,1] <- pz[i,k-1,1] * tr[1,1] * pmat[1,sumCH[i,k]]
      pz[i,k,2] <- inprod(pz[i,k-1,], tr[,2]) * pmat[2,sumCH[i,k]]
    }
    
    like[i] <- sum(pz[i,n.occasions,])
    ones[i] ~ dbin(like[i], sumFR[i])
  }
}
    ", fill = TRUE)
sink() 
#-----------------------------------------------------------------------------#

# 
# JM.inits <- function(){list(s = runif((n.occasions - 1), 0, 1),
#                             p = runif((n.occasions - 1), 0, 1),
#                             z = cjs.init.z(CH, indf))}

JM.data <- list(NsumCH = NsumCH, n.occasions = n.occasions, sumCH = sumCH,
                sumf = sumf, sumFR = sumFR, ones = ones)

JM.par <- c('s', 'p')

ni <- 500
nt <- 1
nb <- 100


# t1 <- proc.time()
JM.out <- jags(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_Constant.jags",
               n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb)

JM.out.b <- jags.parallel(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_Constant.jags",
               n.chains = 3, n.iter = ni, export_obj_names = c("ni"))
# t2 <- proc.time()

print(JM.out, digits = 3)

#-----------------------------------------------------------------------------#

nb <- 500


n.iter = seq(1000,5000,1000)

fit.list = list()

time.list = list()

for(i in 1:5){
  ni = n.iter[i]

  t1 <- proc.time()
  JM.out <- jags(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_Constant.jags",
                 n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb)
  t2 <- proc.time()
  
  time.list[[i]] = t2 - t1
  
  fit.list[[i]] = JM.out
  
}


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

my.n.iter = seq(0,10000,500)[- 1]
# my.n.iter = c(100,200)

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
    
    JM.C <- jags.parallel(JM.data, inits = NULL, JM.par, "JAGS_Marginalized_Constant.jags",
                              n.chains = 3, n.iter = iter.in, export_obj_names = c("iter.in", "seed"),
                          jags.seed = seed, envir = my.env) 
    
    t2 <- proc.time()
    
    attr(JM.C, 'time') <- (t2 - t1)[3]
    
    JM.C
    
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



 # all.jags.m.constant = all.out

# rm(list=setdiff(ls(), "all.jags.m.constant"))









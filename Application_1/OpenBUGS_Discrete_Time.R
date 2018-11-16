###############################################################################
#                                                                        Nov 18
#  Fitting a multi-state version of a CJS model to the RBT data 
#  Discrete OpenBUGS version - fixed time effects
#
#  Notes:
#  * OpenBUGS will crash with indexing of omega and rho like other CJS examples,
#    see note below. Same issue as WinBUGS. 
#
#  To do: 
#  * 
#
###############################################################################
setwd(paste0(getwd(), '/Application_1'))
library(R2OpenBUGS)

source("RBT_Functions.R", chdir = F)

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

ni <- 10
nt <- 1
nb <- 5

t1 <- proc.time()
OD.out <- R2OpenBUGS::bugs(OD.data, inits = NULL, OD.par, "OpenBUGS_Discrete_Time.OpenBUGS",
               n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb)
t2 <- proc.time()

print(OD.out, digits = 3)

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

n.core = 30

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used
cllibs <- clusterEvalQ(cl1, c(library(R2OpenBUGS)))

all.t1 = proc.time()
n.runs = 10

# my.n.iter = c(10, 15)
my.n.iter = seq(0,10000,500)[- 1]
# my.n.iter = rep(500, 1)

# my.n.iter = my.n.iter[1:10]
my.n.iter = my.n.iter[11:20]

big.fit.list = list()

# seeds <- sample(1:1e5, size = n.runs)
seeds <- rep(0, n.runs)

# let all of the clusters have whatever objects are in the workspace
clusterExport(cl = cl1, varlist = ls(), envir = environment())

start.time <- Sys.time()  # start timer

# n = 1
out = list()
#--------------------------------------
out <- foreach(j = seeds, .errorhandling = 'pass') %:%
  
  foreach(i = my.n.iter) %dopar% {
    
    seed = j
    seed = seed + sample(1:14, size = 1) # must be integer between 1 and 14, see ?bugs
    
    iter.in = i
    
    t1 <- proc.time()
    
    result <- tryCatch({
      out = R2OpenBUGS::bugs(OD.data, inits = NULL, OD.par, "OpenBUGS_Discrete_Time.OpenBUGS",
                               n.chains = 3, n.iter = iter.in, n.thin = nt,
                               bugs.seed = seed)
    }, 
      warning = function(war) {
        return('BUGS can return a warning -- What!?')}, 
      error = function(err) {
        return(NULL)} 
    ) # END tryCatch
    
    t2 <- proc.time()
    
    attr(result, 'time') <- (t2 - t1)[3]
    
    return(result)
    
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

tmp = run.times(all.out)

all.Open.d.time.1 = all.out

rm(list=setdiff(ls(), "all.Open.d.time.1"))

save.image("U:/Desktop/Fish_Git/Marginalized/Application_1/working_Runs/Open_D_Time.RData")

#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Fitting a multi-state version of a CJS model to the RBT data 
#  Discrete WinBUGS version - fixed time effects
#
#  Notes:
#  * WinBUGS will crash with indexing of omega and rho like other CJS examples,
#    see note below.
#
#  To do: 
#  * 
#
###############################################################################
library(R2WinBUGS)
setwd("U:/Desktop/Fish_Git/Marginalized")
source(paste0(getwd(),"/Functions.R"), chdir = F)

setwd(paste0(getwd(), '/Application_1'))

data.dir = paste0(getwd(), "/Data")
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

# WD.inits <- function(){list(phi = runif((Nint), .4, .6),
#                             p = runif((Nint), .4, .6),
#                             Z = cjs.init.z(CH, indf))}

ni <- 10
nt <- 1
nb <- 5

t1 <- proc.time()
WD.out <- R2WinBUGS::bugs(WD.data, inits = NULL, WD.par, "WinBUGS_Discrete_Time.WinBUGS",
               n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb, debug = FALSE)
t2 <- proc.time()

print(WD.out, digits = 3)

#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

Sys.time()
n.core = 12

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used
cllibs <- clusterEvalQ(cl1, c(library(R2WinBUGS)))

all.t1 = proc.time()
n.runs = 12

my.n.iter = rep(10000)

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
    # seed = seed + sample(1:14, size = 1) # must be integer between 1 and 14, see ?bugs
    seed = seed + sample(1:1e4, size = 1)
    
    iter.in = i
    
    t1 <- proc.time()
    
    result <- tryCatch({
      out <- R2WinBUGS::bugs(JD.data, inits = NULL, JD.par, "WinBUGS_Discrete_Time.WinBUGS",
                                n.chains = 3, n.iter = iter.in, n.thin = nt,
                             bugs.seed = seed)
      
    }, 
    warning = function(war) {
      return('BUGS can return a warning -- What!?')}, 
    error = function(err) {
      return(NULL)} 
    ) # END tryCatch
    
    t2 <- proc.time()
    
    if(is.null(result) == FALSE){
      attr(result, 'time') <- (t2 - t1)[3]
    }
    
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
tmp

all.bugs.d.time.6 = all.out

rm(list=setdiff(ls(), "all.bugs.d.time.6"))

# save.image("U:/Desktop/Fish_Git/Marginalized/Application_1/working_Runs/BUGS_D_Time_6.RData")

#-----------------------------------------------------------------------------#
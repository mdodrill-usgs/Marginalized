###############################################################################
#                                                                     Spring 19
#  Fitting a multi-state version of a CJS model to the RBT data 
#  Discrete JAGS version - fixed time effects
#
#  Notes:
#  * 
#
###############################################################################
library(R2jags)
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

ni <- 10
nt <- 1
nb <- 5

t1 <- proc.time()
JD.out <- jags(JD.data, inits = NULL, JD.par, "JAGS_Discrete_Time.jags",
               n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb)
t2 <- proc.time()

print(JD.out, digits = 3)
#--------------------------------------

t1 <- proc.time()
JD.out <- jags.parallel(JD.data, inits = NULL, JD.par, "JAGS_Discrete_Time.jags",
                        n.chains = 3, n.iter = ni, n.thin = nt, n.burnin = nb,
                        export_obj_names = c("ni", "nt", "nb"))

t2 <- proc.time()

print(JD.out, digits = 3)


#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

n.core = 5  # really n.core * 3

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used 
cllibs <- clusterEvalQ(cl1, c(library(R2jags)))

all.t1 = proc.time()
n.runs = 10

my.n.iter = c(11000,13000,15000,16000,17000,19000)

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
    seed = seed + sample(1:1e5, size = 1)
    
    iter.in = i
    
    t1 <- proc.time()
    
    my.env = environment()
    
    JD.t <- jags.parallel(JD.data, inits = NULL, JD.par, "JAGS_Discrete_Time.jags",
                          n.chains = 3, n.iter = iter.in, export_obj_names = c("iter.in", "seed"),
                          jags.seed = seed, envir = my.env) 
    
    t2 <- proc.time()
    
    attr(JD.t, 'time') <- (t2 - t1)[3]
    
    JD.t
    
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

all.jags.d.time.5 = all.out

rm(list=setdiff(ls(), "all.jags.d.time.5"))

# save.image("U:/Desktop/Fish_Git/Marginalized/Application_1/working_Runs/JAGS_D_Time_5.RData")
#-----------------------------------------------------------------------------#
# end
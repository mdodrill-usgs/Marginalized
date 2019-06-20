###############################################################################
#                                                                        Dec 18
#  Dynamic Community Occupancy Models - Sky Islands  
#  Marginalized Stan version 
#
#  Notes:
#  * 
#
#  To do: 
#  * format .Stan code (identation!)
#
###############################################################################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())  

setwd("U:\\Desktop\\Fish_Git\\Marginalized")

setwd(paste0(getwd(), '/Application_5'))

data.dir = paste0(getwd(), "/Data/")

# dimensions: species, points, years
nspp = 149
nsites = 92
nyears = 5  # 1991, 1992, 1993, 1994, 1995 
nsess = 3
Y = array(NA, dim = c(nspp, nsites, nyears))

names1 = c('ntmb.CM.1991.csv')  
names2 = c('ntmb.CM.1992.csv')  
names3 = c('ntmb.CM.1993.csv')  
names4 = c('ntmb.CM.1994.csv')  
names5 = c('ntmb.CM.1995.csv')  

names = c(names1, names2, names3, names4, names5)

for(q in 1:5){
  Ymat1 = as.matrix(read.csv(paste0(data.dir, names[q]))[,2:93])
  Y[,,q] = array(c(Ymat1), dim = c(dim(Y)[1], dim(Ymat1)[2]))
}

visit <- ifelse(is.na(Y[1,,1]) == TRUE, 2, 1)
sY <- ifelse(is.na(Y) == TRUE, 1, Y + 1)

# habitat covariate (dimensions 92 pts) 
hab = as.matrix(read.csv(paste0(data.dir, 'habitat.CM.csv'), header = FALSE))[2:93,4] 
hab = as.integer(hab)

#-----------------------------------------------------------------------------#
# see Stan model 'Cocc_SM.stan' 
#-----------------------------------------------------------------------------#
params <- c('psi_mean', 'p_mean', 'rho', 'sigma_v', 'sigma_u', 'sigma_alpha1',
            'sigma_alpha2', 'alpha1_mean','SR')

SM_data = list('nspp' = nspp, 'nsites' = nsites, 'nsess' = nsess,
               'nyears' = nyears, 'hab' = hab,
               'ones' = matrix(1, nrow = nspp, ncol = nsites),
               'sY' = sY, 'visit' = visit)

# ~2.7h for 1000 iter
t1 <- proc.time()
SM_Cocc <- stan("Cocc_SM.stan",
                data = SM_data,
                pars = params,
                chains = 3, iter = 1000) 
t2 <- proc.time()
#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

n.core = 10 # really n.core * 3

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used 
cllibs <- clusterEvalQ(cl1, c(library(rstan)))

n.runs = 10

my.n.iter = c(1000)

big.fit.list = list()

seeds <- sample(1:1e5, size = n.runs)   

# let all of the clusters have whatever objects are in the workspace
clusterExport(cl = cl1, varlist = ls(), envir = environment())

start.time <- Sys.time()  # start timer

out = list()
#--------------------------------------
out <- foreach(j = seeds) %:% 
  
  foreach(i = my.n.iter) %dopar% {
    
    seed = j
    seed = seed + sample(1:1e5, size = 1)
    
    iter.in = i
    
    t1 <- proc.time()
    
    # SM.c <- stan(fit = SM_Cocc,     # faster to compile the model above...
    SM_Cocc <- stan("Cocc_SM.stan",
                    data = SM_data,
                    pars = params,
                    chains = 3, iter = iter.in, seed = seed, cores = 3) 
    
    t2 <- proc.time()
    
    attr(SM_Cocc, 'time') <- (t2 - t1)[3]
    
    SM_Cocc
    
  } 
#--------------------------------------

end.time = Sys.time()
time.taken = end.time - start.time
print(round(time.taken, 2))

stopCluster(cl1)  # close the clusters

length(out)

length(out[[1]])

all.out = do.call('c', out)
length(all.out)

tmp = run.times(all.out)

all.stan.m = all.out
rm(list=setdiff(ls(), c("all.stan.m")))

save.image("U:/Desktop/Fish_Git/Marginalized/Application_5/working_runs/SM_Cocc_1.RData")
#-----------------------------------------------------------------------------#
# end
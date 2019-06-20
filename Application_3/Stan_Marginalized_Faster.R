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
# Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

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
pasted <- function(x){
  paste0(x[1], x[2], x[3], x[4], x[5], x[6], x[7], collapse = "")
}

unpaste <- function(x){
  as.numeric(c(substr(x,1,1), substr(x,2,2), substr(x,3,3),
               substr(x,4,4), substr(x,5,5), substr(x,6,6), substr(x,7,7)))
}

FR <- numeric()
spp <- numeric()
sY2 <- rep(NA,5)
visit2 <- numeric()
hab2 <- numeric()
lookup <- matrix(NA,nspp,nsites)
temp4 <- 0

for(S in 1:nspp){
  temp <- apply(cbind(sY[S,,], visit,hab), 1, pasted)
  temp2 <- table(temp)
  FR <- c(FR, as.numeric(temp2))
  spp <- c(spp, rep(S, length(temp2)))
  for(j in 1:length(temp2)){
    temp3 <- unpaste(names(temp2[j]))
    sY2 <- rbind(sY2, temp3[1:5])
    visit2 <- c(visit2, temp3[6])
    hab2 <- c(hab2, temp3[7])
  }
  lookup[S,] <- match(temp, names(temp2)) + temp4
  temp4 <- temp4 + length(temp2)
}

sY2 <- sY2[-1,]


#-----------------------------------------------------------------------------#
params <- c('psi_mean', 'p_mean', 'rho', 'sigma_v', 'sigma_u', 'sigma_alpha1',
            'sigma_alpha2', 'alpha1_mean','SR')

# SM_data = list('nspp' = nspp, 'nsites' = nsites, 'nsess' = nsess,
#                'nyears' = nyears, 'hab' = hab,
#                'ones' = matrix(1, nrow = nspp, ncol = nsites),
#                'sY' = sY, 'visit' = visit)
SM_data<-list(nspp=nspp,nsess=nsess,nsites=nsites,nyears=nyears,
               hab2=hab2,sY2=sY2,visit2=visit2,nhab=7,NuDH=length(FR),
               spp=spp,FR=FR,lookup=lookup)


# 
t1 <- proc.time()
SM_Cocc <- stan("Cocc_SM2.stan",
                data = SM_data,
                pars = params,
                chains = 1, iter = 1000) 
t2 <- proc.time()
#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

n.core = 4 # really n.core * 3

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used 
cllibs <- clusterEvalQ(cl1, c(library(rstan)))

n.runs = 10

# my.n.iter = c(500, 1000)
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
    SM_Cocc <- stan("Cocc_SM2.stan",
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

all.stan.m.test.7 = all.out
rm(list=setdiff(ls(), c("all.stan.m.test.7")))

save.image("U:/Desktop/Fish_Git/Marginalized/Application_5/working_runs/SM_Cocc2_runs_7_citrix.RData")
#-----------------------------------------------------------------------------#
# end
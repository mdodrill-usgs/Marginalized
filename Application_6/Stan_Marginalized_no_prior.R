###############################################################################
#                                                                        Nov 18
#  Fitting an Integrated Population Model to Brown Trout Data
#  Marginalized Stan version 
#
#  Notes:
#  * The model runs from the fall of 2000 to fall of 2017 on a seasonal basis
#  * We define three size states based on total length in mm 
#    - 0 - 200; 200 - 350; 350 +
#  * No prior on the survival (phi) parameters
#
#  To do: 
#  * 
#
###############################################################################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

source(paste0(getwd(),"/Functions.R"), chdir = F)

setwd(paste0(getwd(), '/Application_6'))

#-----------------------------------------------------------------------------#
# read in data
NO_catch <- read.csv(paste0(getwd(), "/Data/", "NO_catch.csv"), header = TRUE)
AZGF_catch <- read.csv(paste0(getwd(), "/Data/", "AZGF_catch.csv"), header = TRUE)
MR_data <- read.csv(paste0(getwd(), "/Data/", "bCH.csv"), header = FALSE)

# extract/reformat data
bNOc <- as.matrix(NO_catch[,1:3])
NOpasses <- NO_catch$NOpasses
seasNO <- NO_catch$seasNO
spawn <- ifelse(seasNO == 1, 4, 3)
bAZ <- as.matrix(AZGF_catch[,1:3])
ts <- AZGF_catch$ts
AZeff <- AZGF_catch$AZeff
NAZsamps <- length(AZeff)

# capture-recapture data
allCH <- MR_data[,1:23]

bCH = collapse.ch(allCH)[[1]]
FR = collapse.ch(allCH)[[2]]

findlast <- function(x){ifelse(x[23] == 1, 22, max(which(x[1:22] != 4)))}
last <- apply(bCH, 1, findlast)

NCH <- length(last)
findfirst <- function(x){which(x != 4)[1]}
sumf <- apply(bCH, 1, findfirst)

#-----------------------------------------------------------------------------#
# working with the Stan model in a seperate tab... see "Stan_M...stan"

# sink("Stan_M...stan")
# cat("
#     

# Put the model code here...

#     ",fill=TRUE)
# sink()    


#-----------------------------------------------------------------------------#

# sm.params <- c("s", "p")
# 
# sm.data <- list(NsumCH = NsumCH, n_occasions = n.occasions, sumCH = sumCH,
#                 sumf = sumf, sumFR = sumFR)

sm.data <- list(NAZsamps = NAZsamps, ts = ts, AZeff = AZeff, bAZ = bAZ,
                seasNO = seasNO, bNOc = bNOc, NOpasses = NOpasses, ones = FR,
                FR = FR, last = last, bCH = bCH, NCH = NCH, sumf = sumf,
                spawn = spawn)

# BM_JM.par <- c('beta.I', 'bphi', 'bpsi1', 'bpsi2', 'mu.blp', 'sd.blp', "lbeta.0",
#                "mu.I", "I", "Beta", "IN", "AZadj", "sd.I", "sd.lphi", "sd.blp",
#                'sd.beta', "bN")

sm.params = c('bphi', 'bpsi1', 'bpsi2', 'mu_blp', 'sd_blp', 'lbeta_0',
              'mu_I', 'I', 'Beta', 'IN', 'AZadj', 'sd_I',# 'sd_lphi',
              'sd_beta', 'bN', 'bp_pass')

# MCMC settings
ni = 10
nt = 1
nb = 5
nc = 1


# Call Stan from R 
SM.c <- stan("Stan_Marginalized_no_prior.stan",
             data = sm.data,
             pars = sm.params,
             control = list(max_treedepth = 14, adapt_delta = .85),
             chains = nc, iter = ni, thin = nt, seed = 1) 

# stan.fit = SM.c

# rm(list=setdiff(ls(), "stan.fit"))
#-----------------------------------------------------------------------------#
# library(shinystan)
# 
# launch_shinystan(object = SM.c)

#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

n.core = 1 # really n.core * 3

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used
cllibs <- clusterEvalQ(cl1, c(library(rstan)))

n.runs = 1

my.n.iter = c(7500)

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
    
    SM.c <- stan("Stan_Marginalized_no_prior.stan",
                 data = sm.data,
                 pars = sm.params,
                 # control = list(max_treedepth = 14, adapt_delta = .85),
                 chains = 3, iter = iter.in, seed = seed, cores = 3)
    
    
    t2 <- proc.time()
    
    attr(SM.c, 'time') <- (t2 - t1)[3]
    
    SM.c
    
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

all.stan.m.no.prior = all.out
rm(list=setdiff(ls(), c("all.stan.m.no.prior")))

save.image("U:/Desktop/Fish_Git/Marginalized/Application_6/working_runs/Stan_run_no_prior.RData")
#-----------------------------------------------------------------------------#
# end

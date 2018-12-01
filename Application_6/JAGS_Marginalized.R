###############################################################################
#                                                                        Nov 18
#  Fitting an Integrated Population Model to Brown Trout Data
#  Marginalized JAGS version 
#
#  Notes:
#  * The model runs from the fall of 2000 to fall of 2017 on a seasonal basis
#  * We define three size states based on total length in mm 
#    - 0 - 200; 200 - 350; 350 +
#
#  To do: 
#  * 
#
###############################################################################
setwd(paste0(getwd(), '/Application_6'))
library(R2jags)

#-----------------------------------------------------------------------------#
# read in data
NO_catch <- read.csv(paste0(getwd(), "/Data/", "NO_catch.csv"), header = TRUE)
AZGF_catch <- read.csv(paste0(getwd(), "/Data/", "AZGF_catch.csv"), header = TRUE)
MR_data <- read.csv(paste0(getwd(), "/Data/", "bCH.csv"), header = FALSE)

# extract/reformat data
bNOc <- NO_catch[,1:3]
NOpasses <- NO_catch$NOpasses
seasNO <- NO_catch$seasNO
bAZ <- AZGF_catch[,1:3]
ts <- AZGF_catch$ts
AZeff <- AZGF_catch$AZeff
NAZsamps <- length(AZeff)

findlast <- function(x){ifelse(x[24] == 1, 23, max(which(x[1:23] != 4)))}
last <- apply(MR_data, 1, findlast)

bCH <- MR_data[,1:23]
NCH <- length(last)
FR <- rep(1, NCH)
findfirst <- function(x){which(x != 4)[1]}
sumf <- apply(bCH, 1, findfirst)

#-----------------------------------------------------------------------------#
sink("JAGS_Marginalized.jags")
cat("
model {
  # lphi is the logit of survival and is given a prior based on the Lorenzen
  # function and the average mass of fish in each size class during each season.
  # variation from the priod mode is determined by an estimated variance
  # parameter (sd_lphi)
  
  lphi[1,1] ~ dnorm(1.08, tau_lphi)
  lphi[1,2] ~ dnorm(1.14, tau_lphi)
  lphi[1,3] ~ dnorm(1.26, tau_lphi)
  lphi[1,4] ~ dnorm(1.38, tau_lphi)
  lphi[2,1] ~ dnorm(1.96, tau_lphi)
  lphi[2,2] ~ dnorm(1.96, tau_lphi)
  lphi[2,3] ~ dnorm(2.02, tau_lphi)
  lphi[2,4] ~ dnorm(2.02, tau_lphi)
  lphi[3,1] ~ dnorm(2.29, tau_lphi)
  lphi[3,2] ~ dnorm(2.29, tau_lphi)
  lphi[3,3] ~ dnorm(2.29, tau_lphi)
  lphi[3,4] ~ dnorm(2.29, tau_lphi)
  tau_lphi <- pow(sd_lphi, -2)
  sd_lphi ~ dunif(0.01,4)
  
  for(j in 1:3){
    for(k in 1:4){
      logit(bphi[j,k]) <- lphi[j,k]
    }
  }
  
  bpsi2 ~ dunif(0,1) # growth os size class 2 fish into size class 3
  
  # define transition matrix that combines survival and growth parameters
  for(i in 1:4){
    btrans[1,3,i] <- 0
    btrans[1,4,i] <- 1 - bphi[1,i]
    btrans[2,1,i] <- 0
    btrans[2,2,i] <- bphi[2,i] * (1 - bpsi2)
    btrans[2,3,i] <- bphi[2,i] * bpsi2
    btrans[2,4,i] <- 1 - bphi[2,i]
    btrans[3,1,i] <- 0
    btrans[3,2,i] <- 0
    btrans[3,3,i] <- bphi[3,i]
    btrans[3,4,i] <- 1 - bphi[3,i]
    btrans[4,1,i] <- 0
    btrans[4,2,i] <- 0
    btrans[4,3,i] <- 0
    btrans[4,4,i] <- 1
  }
  
  # size class one transitions are done separately because growth is allowed to vary between seasons
  for(i in 1:3){
    bpsi1[i] ~ dunif(0,1)
    btrans[1,1,i] <- bphi[1,i] * (1 - bpsi1[i])
    btrans[1,2,i] <- bphi[1,i] * bpsi1[i]
  }
  
  btrans[1,1,4] <- 0
  btrans[1,2,4] <- bphi[1,4]
  
  bpi <- .08 # proportion of brown trout population in NO reach 1
  tau_blp <- pow(sd_blp, -2)
  sd_blp ~ dunif(0.1, 2) # trip to trip deviation in NO pcaps
  
  for(i in 1:4){
    mu_blp[i] ~ dnorm(-3, .25) # mean pcaps per pass on a logit scale for three size classes, plus largest size class during spawning season
  }
  
  # this loop calculates actual per pass pcaps for each trip and modifies based on # of passes  
  for(j in 1:23){
    spawn[j] <- 3 + step(-1 * seasNO[j] + 1.1)  # change here to use 'spawn' input
    blp_pass[j,1] ~ dnorm(mu_blp[1], tau_blp)
    blp_pass[j,2] ~ dnorm(mu_blp[2], tau_blp)
    blp_pass[j,3] ~ dnorm(mu_blp[spawn[j]], tau_blp)
    
    for(k in 1:3){
      logit(bp_pass[j,k]) <- blp_pass[j,k]
      bp[j,k,k] <- 1 - pow((1 - bp_pass[j,k]), NOpasses[j])
      bp[j,k,4] <- 1 - bp[j,k,k]
    }
    
    bp[j,1,2] <- 0
    bp[j,1,3] <- 0
    bp[j,2,1] <- 0
    bp[j,2,3] <- 0
    bp[j,3,1] <- 0
    bp[j,3,2] <- 0
    bp[j,4,1] <- 0
    bp[j,4,2] <- 0
    bp[j,4,3] <- 0
    bp[j,4,4] <- 1
  }
  
  for(k in 1:NCH){
    pz[k,sumf[k],1] <- equals(bCH[k,sumf[k]], 1)
    pz[k,sumf[k],2] <- equals(bCH[k,sumf[k]], 2)
    pz[k,sumf[k],3] <- equals(bCH[k,sumf[k]], 3)
    pz[k,sumf[k],4] <- 0
    
    for(j in sumf[k]:(last[k] - 1)){
      for(i in 1:4){
        pz[k,(j + 1),i] <- inprod(pz[k,j,], btrans[,i,seasNO[(j + 1)]]) * bp[j,i,bCH[k,(j + 1)]]
      }
    }
    
    ll[k] <- sum(pz[k, last[k],])
    one[k] ~ dbin(ll[k], FR[k])
  }
  
  # calculate offset for each size classes of AZGF effort and calculate expected pcap
  AZadj[1] ~ dnorm(0,1)
  AZadj[2] ~ dnorm(0,1)
  AZadj[3] ~ dnorm(0,1)
  mu_AZ[1] <- mu_blp[1] + AZadj[1]
  mu_AZ[2] <- mu_blp[2] + AZadj[2]
  mu_AZ[3] <- mu_blp[3] + AZadj[3]
  IN[1] <- 0 # initial abundances of size class 1 fish
  IN[2] ~ dunif(0, 1000) # initial abundances of size class 2 fish
  IN[3] ~ dunif(0, 1000) # initial abundances of size class 3 fish
  bN[1,1] <- IN[1]
  bN[1,2] <- IN[2]
  bN[1,3] <- IN[3]
  
  # variance term controlling unexplained variation in reproductive rate (BETA) 
  tau_beta <- pow(sd_beta,-2)
  sd_beta ~ dunif(0.1,4)
  
  # log of the median reproductive rate - i.e., an intercept
  lbeta_0 ~ dunif(-6,0)
  
  # log of the median immigration rate of large brown trout - i.e., the intercept
  mu_I ~ dunif(0,6)
  
  # variance term controlling unexplained variation in immigration
  tau_I <- pow(sd_I,-2)
  sd_I ~ dunif(0.01,3)
  
  # calculate actual immigration in each interval on log scale
  for(j in 1:68){
    I[j] ~ dnorm(mu_I, tau_I)
  }
  
  # calculate latent abundance of brown trout from fall 2000 to end of 2017
  for (j in 1:17){
    for (k in 1:3){
      bN[((j-1) * 4 + k + 1),1] <- btrans[1,1,k] * bN[((j-1) * 4 + k),1]
      bN[((j-1) * 4 + k + 1),2] <- btrans[1,2,k] * bN[((j-1) * 4 + k),1] + btrans[2,2,k] * bN[((j-1) * 4 + k),2]
      bN[((j-1) * 4 + k + 1),3] <- btrans[2,3,k] * bN[((j-1) * 4 + k),2] + btrans[3,3,k] * bN[((j-1) * 4 + k),3] + exp(I[((j-1) * 4 + k)])
    }
    
    # BNT recruits produced in fall as a function weighted sum of adults (wA) and reprodutive rate (Beta) in winter
    wA[j] <- (bN[((j - 1) * 4 + 2),2] + 4 * bN[((j - 1) * 4 + 2),3])
    beta_eps[j] ~ dnorm(0, tau_beta)
    Beta[j] <- exp(lbeta_0 + beta_eps[j]) 
    
    # between summer and fall all bnt graduate to sz 2 & recruits show up
    bN[(1 + j * 4),1] <- wA[j] * Beta[j]
    bN[(1 + j * 4),2] <- btrans[1,2,4] * bN[(j * 4),1] + btrans[2,2,4] * bN[(j * 4),2]
    bN[(1 + j * 4),3] <- btrans[2,3,4] * bN[(j * 4),2] + btrans[3,3,4] * bN[(j * 4),3] + exp(I[(j * 4)])
  }
  
  # 2000 - 2017 AZGF data
  for(j in 1:NAZsamps){
    for(k in 1:3){
      blpAZ[j,k] ~ dnorm(mu_AZ[k], tau_blp)
      logit(bpAZ[j,k]) <- blpAZ[j,k]
      blamAZ[j,k] <- bpAZ[j,k] * bN[ts[j],k] * AZeff[j] / 35 # predicted catch AZ (35 is ~h to do LF, by AZ)
      bAZ[j,k] ~ dpois(blamAZ[j,k])
    }
  }
  
  # 2012 - 2017 NO: starts in april 2012
  for(j in 1:23){
    for(k in 1:3){
      blamNO[j,k] <- bp[j,k,k] * bpi * bN[(j + 46),k]
      bNOc[j,k] ~ dpois(blamNO[j,k])
    }
  }
  
}
    ",fill=TRUE)
sink()

#-----------------------------------------------------------------------------#
BM_JM.data <- list(NAZsamps = NAZsamps, ts = ts, AZeff = AZeff, bAZ = bAZ,
                   seasNO = seasNO, bNOc = bNOc, NOpasses = NOpasses, ones = FR,
                   FR = FR, last = last, bCH = bCH, NCH = NCH, sumf = sumf)

BM_JM.par <- c('bphi', 'bpsi1', 'bpsi2', 'mu_blp', 'sd_blp', "lbeta_0",
             "mu_I", "I", "Beta", "IN", "AZadj", "sd_I", "sd_lphi",
             'sd_beta', "bN", "bp_pass")

# BM_JM.par = c('blp_pass')

ni <- 20000
# ni <- 1000

t1 <- proc.time()
jags.fit <- jags.parallel(BM_JM.data, inits = NULL, BM_JM.par, "JAGS_Marginalized.jags",
                       n.chains = 3, n.iter = ni, export_obj_names = c("ni"))
t2 <- proc.time()
#-----------------------------------------------------------------------------#

# rm(list=setdiff(ls(), "jags.fit"))
#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

n.core = 3  # really n.core * 3

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used 
cllibs <- clusterEvalQ(cl1, c(library(R2jags)))

all.t1 = proc.time()
n.runs = 10

# my.n.iter = c(10, 15)
# my.n.iter = seq(0,10000,500)[- 1]
# my.n.iter = c(12000,14000,18000,20000)
my.n.iter = c(11000,13000,15000,16000,17000,19000)
# my.n.iter = c(20000)

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

# tmp = run.times(all.out)

all.jags.d.time.3 = all.out

rm(list=setdiff(ls(), "all.jags.d.time.3"))

save.image("U:/Desktop/Fish_Git/Marginalized/Application_1/working_Runs/JAGS_D_Time_3.RData")
#-----------------------------------------------------------------------------#
# end








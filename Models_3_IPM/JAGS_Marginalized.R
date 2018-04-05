###############################################################################
#                                                                      March 18
#        Fitting an Integrated Population Model to Brown Trout Data
#        Marginalized JAGS version 
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
setwd('C:\\Users\\mdodrill\\Desktop\\Fish_Git\\marginalized_2\\Models_3_IPM')
library(R2jags)

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
  # lphi is the logit of survival and is given a prior based on the lorenzen
  # function and the average mass of fish in each size class during each season.
  # variation from the priod mode is determined by an estimated variance
  # parameter (sd.lphi)
  
  lphi[1,1] ~ dnorm(1.08, tau.lphi)
  lphi[1,2] ~ dnorm(1.14, tau.lphi)
  lphi[1,3] ~ dnorm(1.26, tau.lphi)
  lphi[1,4] ~ dnorm(1.38, tau.lphi)
  lphi[2,1] ~ dnorm(1.96, tau.lphi)
  lphi[2,2] ~ dnorm(1.96, tau.lphi)
  lphi[2,3] ~ dnorm(2.02, tau.lphi)
  lphi[2,4] ~ dnorm(2.02, tau.lphi)
  lphi[3,1] ~ dnorm(2.29, tau.lphi)
  lphi[3,2] ~ dnorm(2.29, tau.lphi)
  lphi[3,3] ~ dnorm(2.29, tau.lphi)
  lphi[3,4] ~ dnorm(2.29, tau.lphi)
  tau.lphi <- pow(sd.lphi, -2)
  sd.lphi ~ dunif(0.01,4)
  
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
  tau.blp <- pow(sd.blp, -2)
  sd.blp ~ dunif(0.1, 2) # trip to trip deviation in NO pcaps
  
  for(i in 1:4){
    mu.blp[i] ~ dnorm(-3, .25) # mean pcaps per pass on a logit scale for three size classes, plus largest size class during spawning season
  }
  
  for(j in 1:23){
    #this loop calculates actual per pass pcaps for each trip and modifies based on # of passes
    spawn[j] <- 3 + step(-1 * seasNO[j] + 1.1)
    blp_pass[j,1] ~ dnorm(mu.blp[1], tau.blp)
    blp_pass[j,2] ~ dnorm(mu.blp[2], tau.blp)
    blp_pass[j,3] ~ dnorm(mu.blp[spawn[j]], tau.blp)
    
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
  mu.AZ[1] <- mu.blp[1] + AZadj[1]
  mu.AZ[2] <- mu.blp[2] + AZadj[2]
  mu.AZ[3] <- mu.blp[3] + AZadj[3]
  IN[1] <- 0 # initial abundances of size class 1 fish
  IN[2] ~ dunif(0, 1000)# initial abundances of size class 2 fish
  IN[3] ~ dunif(0, 1000)# initial abundances of size class 3 fish
  bN[1,1] <- IN[1]
  bN[1,2] <- IN[2]
  bN[1,3] <- IN[3]
  
  # variance term controlling unexplained variation in reproductive rate (BETA) 
  tau.beta <- pow(sd.beta,-2)
  sd.beta ~ dunif(0.1,4)
  
  # log of the median reproductive rate - i.e., an intercept
  lbeta.0 ~ dunif(-6,0)
  
  # log of the median immigration rate of large brown trout - i.e., the intercept
  mu.I ~ dunif(0,6)
  
  # variance term controlling unexplained variation in immigration
  tau.I <- pow(sd.I,-2)
  sd.I ~ dunif(0.01,3)
  
  # calculate actual immigration in each interval on log scale
  for(j in 1:68){
    I[j] ~ dnorm(mu.I, tau.I)
  }
  
  # calculate latent abundance of brown trout from fall 2000 to end of 2017
  for (j in 1:17){
    for (k in 1:3){
      bN[((j-1) * 4 + k + 1),1] <- btrans[1,1,k] * bN[((j-1) * 4 + k),1]
      bN[((j-1) * 4 + k + 1),2] <- btrans[1,2,k] * bN[((j-1) * 4 + k),1] + btrans[2,2,k] * bN[((j-1) * 4 + k),2]
      bN[((j-1) * 4 + k + 1),3] <- btrans[2,3,k] * bN[((j-1) * 4 + k),2] + btrans[3,3,k] * bN[((j-1) * 4 + k),3] + exp(I[((j-1) * 4 + k)])
    }
    
    # BNT eggs produced in winter as a function weighted sum of adults (wA) and reprodutive rate (Beta)
    wA[j] <- (bN[((j - 1) * 4 + 2),2] + 4 * bN[((j - 1) * 4 + 2),3])
    beta.eps[j] ~ dnorm(0, tau.beta)
    Beta[j] <- exp(lbeta.0 + beta.eps[j])
    
    # between summer and fall all bnt graduate to sz 2 & recruits show up
    bN[(1 + j * 4),1] <- wA[j] * Beta[j]
    bN[(1 + j * 4),2] <- btrans[1,2,4] * bN[(j * 4),1] + btrans[2,2,4] * bN[(j * 4),2]
    bN[(1 + j * 4),3] <- btrans[2,3,4] * bN[(j * 4),2] + btrans[3,3,4] * bN[(j * 4),3] + exp(I[(j * 4)])
  }
  
  # 2000 - 2017 AZGF data
  for(j in 1:NAZsamps){
    for(k in 1:3){
      blpAZ[j,k] ~ dnorm(mu.AZ[k], tau.blp)
      logit(bpAZ[j,k]) <- blpAZ[j,k]
      blamAZ[j,k] <- bpAZ[j,k] * bN[ts[j],k] * AZeff[j] / 35
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

BM_JM.par <- c('beta.I', 'bphi', 'bpsi1', 'bpsi2', 'mu.blp', 'sd.blp', "lbeta.0",
             "mu.I", "I", "Beta", "IN", "AZadj", "sd.I", "sd.lphi", "sd.blp",
             'sd.beta', "bN")

# ni <- 10000
ni <- 1000

BM_JM <- jags.parallel(BM_JM.data, inits = NULL, BM_JM.par, "JAGS_Marginalized.jags",
                       n.chains = 3, n.iter = ni, export_obj_names = c("ni"))

print(BM_JM, digits = 3)

#-----------------------------------------------------------------------------#
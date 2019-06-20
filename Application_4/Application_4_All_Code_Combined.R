
# do I need any of the main functions for this....?


###############################################################################
#                                                                     Spring 19
#  Fitting an Integrated Population Model to Brown Trout Data
#  Discrete JAGS version 
#
#  Notes:
#  * The model runs from the fall of 2000 to fall of 2017 on a seasonal basis
#  * We define three size states based on total length in mm 
#    - 0 - 200; 200 - 350; 350 +
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2jags)

# read in data
NO_catch <- read.csv(paste0(getwd(), "/Data/", "NO_catch.csv"), header = TRUE)
AZGF_catch <- read.csv(paste0(getwd(), "/Data/", "AZGF_catch.csv"), header = TRUE)
MR_data <- read.csv(paste0(getwd(), "/Data/", "bCH.csv"), header = FALSE)

# extract/reformat data
bNOc <- NO_catch[,1:3]
NOpasses <- NO_catch$NOpasses
seasNO <- NO_catch$seasNO
spawn <- ifelse(seasNO == 1, 4, 3)
bAZ <- AZGF_catch[,1:3]
ts <- AZGF_catch$ts
AZeff <- AZGF_catch$AZeff
NAZsamps <- length(AZeff)

findlast <- function(x){ifelse(x[24] == 1, 23, max(which(x[1:23] != 4)))}
indlast <- apply(MR_data, 1, findlast)

indCH <- MR_data[,1:23]
NindCH <- length(indlast)
findfirst <- function(x){which(x != 4)[1]}
indf <- apply(indCH, 1, findfirst)

#-----------------------------------------------------------------------------#
sink("JAGS_Discrete.jags")
cat("
model{
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
  
  bpsi2 ~ dunif(0,1) # growth of size class 2 fish into size class 3
  
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
  
  # mean pcaps per pass on a logit scale for three size classes, plus largest
  # size class during spawning season
  for(i in 1:4){
    mu_blp[i] ~ dnorm(-3, .25) 
  }
  
  # this loop calculates actual per pass pcaps for each trip and modifies based on # of passes  
  for(j in 1:23){
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
  
  for(k in 1:NindCH){
    Z[k,indf[k]] <- indCH[k,indf[k]]    
    for(j in indf[k]:(indlast[k] - 1)){
      Z[k,(j + 1)] ~ dcat(btrans[Z[k,j], ,seasNO[(j + 1)]])
      indCH[k,(j + 1)] ~ dcat(bp[j,Z[k,(j + 1)],])
    }
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
    ", fill = TRUE)
sink()

#-----------------------------------------------------------------------------#
JD_inits <- function(){
  z.init <- matrix(NA, nrow = NindCH, ncol = 23)
  for(i in 1:NindCH){
    z.init[i,indf[i]] <- indCH[i,indf[i]]
    for(j in (indf[i] + 1):indlast[i]){
      z.init[i,j] <- ifelse(indCH[i,j] < 4, indCH[i,j],
                            ifelse(z.init[i,(j - 1)] > 1, z.init[i,(j - 1)],
                                   ifelse(seasNO[j] == 4, 2, 1)))
    }
    z.init[i,indf[i]] <- NA
  }	
  list(Z = array(z.init, dim = c(NindCH, 23)))
}

#-----------------------------------------------------------------------------#
BM_JM.data <- list(NAZsamps = NAZsamps, ts = ts, AZeff = AZeff, bAZ = bAZ,
                   seasNO = seasNO, bNOc = bNOc, NOpasses = NOpasses, 
                   indlast = indlast, indCH = indCH, NindCH = NindCH, indf = indf,
                   spawn = spawn)

BM_JM.par <- c('bphi', 'bpsi1', 'bpsi2', 'mu_blp', 'sd_blp', 'lbeta_0',
               'mu_I', 'I', 'Beta', 'IN', 'AZadj', 'sd_I', 'sd_lphi',
               'sd_beta', 'bN', 'bp_pass')

jags.fit <- jags.parallel(BM_JM.data, inits = JD_inits, BM_JM.par, "JAGS_Discrete.jags",
                          n.chains = 3, n.iter = 10, export_obj_names = c("JD_inits"))
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Fitting an Integrated Population Model to Brown Trout Data
#  Marginalized JAGS version 
#
#  Notes:
#  * The model runs from the fall of 2000 to fall of 2017 on a seasonal basis
#  * We define three size states based on total length in mm 
#    - 0 - 200; 200 - 350; 350 +
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2jags)

# read in data
NO_catch <- read.csv(paste0(getwd(), "/Data/", "NO_catch.csv"), header = TRUE)
AZGF_catch <- read.csv(paste0(getwd(), "/Data/", "AZGF_catch.csv"), header = TRUE)
MR_data <- read.csv(paste0(getwd(), "/Data/", "bCH.csv"), header = FALSE)

#-----------------------------------------------------------------------------#
# extract/reformat data
bNOc <- NO_catch[,1:3]
NOpasses <- NO_catch$NOpasses
seasNO <- NO_catch$seasNO
spawn <- ifelse(seasNO == 1, 4, 3)
bAZ <- AZGF_catch[,1:3]
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
sink("JAGS_Marginalized.jags")
cat("
model{
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
  
  bpsi2 ~ dunif(0,1) # growth of size class 2 fish into size class 3
  
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
    # spawn[j] <- 3 + step(-1 * seasNO[j] + 1.1)  # change here to use 'spawn' input
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
    ones[k] ~ dbin(ll[k], FR[k])
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
  for(j in 1:17){
    for(k in 1:3){
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
                   FR = FR, last = last, bCH = bCH, NCH = NCH, sumf = sumf,
                   spawn = spawn)

BM_JM.par <- c('bphi', 'bpsi1', 'bpsi2', 'mu_blp', 'sd_blp', 'lbeta_0',
               'mu_I', 'I', 'Beta', 'IN', 'AZadj', 'sd_I', 'sd_lphi',
               'sd_beta', 'bN', 'bp_pass')

jags.fit <- jags.parallel(BM_JM.data, inits = NULL, BM_JM.par, "JAGS_Marginalized.jags",
                          n.chains = 3, n.iter = 10)
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Fitting an Integrated Population Model to Brown Trout Data
#  Marginalized Stan version 
#
#  Notes:
#  * The model runs from the fall of 2000 to fall of 2017 on a seasonal basis
#  * We define three size states based on total length in mm 
#    - 0 - 200; 200 - 350; 350 +
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

# read in data
NO_catch <- read.csv(paste0(getwd(), "/Data/", "NO_catch.csv"), header = TRUE)
AZGF_catch <- read.csv(paste0(getwd(), "/Data/", "AZGF_catch.csv"), header = TRUE)
MR_data <- read.csv(paste0(getwd(), "/Data/", "bCH.csv"), header = FALSE)

#-----------------------------------------------------------------------------#
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
sink("Stan_Marginalized.Stan")
cat("
// Marginalized Integrated Population Model - Brown Trout Data  
    
data{
  int NAZsamps;                   // Number of samples AGFD
  int ts[NAZsamps];
  vector [NAZsamps] AZeff;
  int bAZ[NAZsamps, 3];
  int seasNO[23];
  int bNOc[23,3];
  int NOpasses[23];              // Number of passes
  int NCH;                       // Number of capture histories
  int FR[NCH];
  int last[NCH];
  int bCH[NCH, 23];
  int sumf[NCH];
  int spawn[23];               // Indicator for spawning month
}

parameters{
  matrix[3, 4]lphi;
  real<lower = 0.01, upper = 4> sd_lphi;                  
  real<lower = 0, upper = 1> bpsi2;
  vector<lower = 0, upper = 1>[3] bpsi1;
  real<lower = 0.1, upper = 2> sd_blp;
  vector[4] mu_blp;
  matrix[23, 3]blp_pass;
  vector[3] AZadj;
  vector<lower = 0, upper = 1000>[2] IN;
  real<lower = 0.1, upper = 4> sd_beta;                  
  vector[17] beta_eps;
  real<lower = -6, upper = 0> lbeta_0;                                
  real<lower = 0, upper = 6> mu_I;                                    
  real<lower = 0.01, upper = 3> sd_I;
  vector[68] I;
  matrix[NAZsamps, 3]blpAZ;
}

transformed parameters{
  matrix[3, 4]bphi;  
  real btrans[4, 4, 4];                                                     
  matrix[23,3]bp_pass;
  real bp[23,4,4];                                                          
  vector[3] mu_AZ;
  matrix[69,3] bN;                                                          
  vector[17] wA;
  vector[17] Beta;
  
  for(j in 1:3){  // remove for no prior version
    for(k in 1:4){
      bphi[j,k] = inv_logit(lphi[j,k]);
    }
  }
  
  // define transition matrix that combines survival and growth parameters
  for(i in 1:4){
    btrans[1,3,i] = 0;
    btrans[1,4,i] = 1 - bphi[1,i];
    btrans[2,1,i] = 0;
    btrans[2,2,i] = bphi[2,i] * (1 - bpsi2);
    btrans[2,3,i] = bphi[2,i] * bpsi2;
    btrans[2,4,i] = 1 - bphi[2,i];
    btrans[3,1,i] = 0;
    btrans[3,2,i] = 0;
    btrans[3,3,i] = bphi[3,i];
    btrans[3,4,i] = 1 - bphi[3,i];
    btrans[4,1,i] = 0;
    btrans[4,2,i] = 0;
    btrans[4,3,i] = 0;
    btrans[4,4,i] = 1;
  }
  
  // size class one transitions are done separately because growth is allowed to vary between seasons
  for(i in 1:3){
    btrans[1,1,i] = bphi[1,i] * (1 - bpsi1[i]);
    btrans[1,2,i] = bphi[1,i] * bpsi1[i];
  }
  
  btrans[1,1,4] = 0;
  btrans[1,2,4] = bphi[1,4];
  
  // this loop calculates actual per pass pcaps for each trip and modifies based on # of passes
  for(j in 1:23){
    for(k in 1:3){
      bp_pass[j,k] = inv_logit(blp_pass[j,k]);
      bp[j,k,k] = 1 - pow((1 - bp_pass[j,k]), NOpasses[j]);                 
      bp[j,k,4] = 1 - bp[j,k,k];
    }
    
    bp[j,1,2] = 0;
    bp[j,1,3] = 0;
    bp[j,2,1] = 0;
    bp[j,2,3] = 0;
    bp[j,3,1] = 0;
    bp[j,3,2] = 0;
    bp[j,4,1] = 0;
    bp[j,4,2] = 0;
    bp[j,4,3] = 0;
    bp[j,4,4] = 1;
  }
  
  // calculate offset for each size classes of AZGF effort and calculate expected pcap
  mu_AZ[1] = mu_blp[1] + AZadj[1];
  mu_AZ[2] = mu_blp[2] + AZadj[2];
  mu_AZ[3] = mu_blp[3] + AZadj[3];
  
  bN[1,1] = 0;                                     
  bN[1,2] = IN[1];
  bN[1,3] = IN[2];
  
  //calculate latent abundance of brown trout from fall 2000 to end of 2017
  for(j in 1:17){
    for(k in 1:3){
      bN[((j-1) * 4 + k + 1),1] = btrans[1,1,k] * bN[((j-1) * 4 + k),1];
      bN[((j-1) * 4 + k + 1),2] = btrans[1,2,k] * bN[((j-1) * 4 + k),1] + btrans[2,2,k] * bN[((j-1) * 4 + k),2];
      bN[((j-1) * 4 + k + 1),3] = btrans[2,3,k] * bN[((j-1) * 4 + k),2] + btrans[3,3,k] * bN[((j-1) * 4 + k),3] + exp(I[((j-1) * 4 + k)]);
    }
    
    // BNT recruits produced in fall as a function weighted sum of adults (wA) and reprodutive rate (Beta) in winter
    wA[j] = (bN[((j - 1) * 4 + 2),2] + 4 * bN[((j - 1) * 4 + 2),3]);
    Beta[j] = exp(lbeta_0 + beta_eps[j]);
    
    // between summer and fall all bnt graduate to sz 2 & recruits show up
    bN[(1 + j * 4),1] = wA[j] * Beta[j];
    bN[(1 + j * 4),2] = btrans[1,2,4] * bN[(j * 4),1] + btrans[2,2,4] * bN[(j * 4),2];
    bN[(1 + j * 4),3] = btrans[2,3,4] * bN[(j * 4),2] + btrans[3,3,4] * bN[(j * 4),3] + exp(I[(j * 4)]);
  }
}

model{
  real pz[NCH, 23, 4];                                    
  matrix[NAZsamps, 3]bpAZ;
  matrix[NAZsamps, 3]blamAZ;
  matrix[23, 3]blamNO;
  vector[4] temp;
  
  ///////////////////////////
  // lphi is the logit of survival and is given a prior based on the Lorenzen
  // function and the average mass of fish in each size class during each season.
  // variation from the priod mode is determined by an estimated variance
  // parameter (sd_lphi)
  
  lphi[1,1] ~ normal(1.08, sd_lphi);
  lphi[1,2] ~ normal(1.14, sd_lphi);
  lphi[1,3] ~ normal(1.26, sd_lphi);
  lphi[1,4] ~ normal(1.38, sd_lphi);
  lphi[2,1] ~ normal(1.96, sd_lphi);
  lphi[2,2] ~ normal(1.96, sd_lphi);
  lphi[2,3] ~ normal(2.02, sd_lphi);
  lphi[2,4] ~ normal(2.02, sd_lphi);
  lphi[3,1] ~ normal(2.29, sd_lphi);
  lphi[3,2] ~ normal(2.29, sd_lphi);
  lphi[3,3] ~ normal(2.29, sd_lphi);
  lphi[3,4] ~ normal(2.29, sd_lphi);
  ///////////////////////////
    
  // mean pcaps per pass on a logit scale for three size classes, plus largest size class during spawning season
  for(i in 1:4){
    mu_blp[i] ~ normal(-3, 2);                              
  }
  
  // (done above in transformed) this loop calculates actual per pass pcaps for each trip and modifies based on # of passes
  for(j in 1:23){
    blp_pass[j,1] ~ normal(mu_blp[1], sd_blp);
    blp_pass[j,2] ~ normal(mu_blp[2], sd_blp);
    blp_pass[j,3] ~ normal(mu_blp[spawn[j]], sd_blp);
  }
  
  for(k in 1:NCH){
    pz[k,sumf[k],1] = (1 == bCH[k,sumf[k]]);
    pz[k,sumf[k],2] = (2 == bCH[k,sumf[k]]);
    pz[k,sumf[k],3] = (3 == bCH[k,sumf[k]]);
    pz[k,sumf[k],4] = 0;
    for(t in sumf[k]:(last[k] - 1)){
      for(i in 1:4){
        for(j in 1:4){
          temp[j] = pz[k,t,j] * btrans[j,i,seasNO[(t + 1)]] * bp[t,i,bCH[k,(t + 1)]];
        }
        pz[k,(t + 1),i] = sum(temp);
      }
    }
    target += FR[k] * log(sum(pz[k,last[k],]));                             
  }
  
  //////////////////////////
  // (done above in transformed) calculate offset for each size classes of AZGF effort and calculate expected pcap
  AZadj[1] ~ normal(0,1);                                       
  AZadj[2] ~ normal(0,1);
  AZadj[3] ~ normal(0,1);
  
  // calculate actual immigration in each interval on log scale
  for(j in 1:68){                                                 
    I[j] ~ normal(mu_I, sd_I);
  }
  
  for(j in 1:17){                                             
    beta_eps[j] ~ normal(0, sd_beta);
  }
  
  //////////////////////////
  // 2000 - 2017 AZGF data
  for(j in 1:3){
    for(k in 2:3){
      blpAZ[j,k] ~ normal(mu_AZ[k], sd_blp);
      bpAZ[j,k] = inv_logit(blpAZ[j,k]);
      blamAZ[j,k] = bpAZ[j,k] * bN[ts[j],k] * AZeff[j] / 35;
      bAZ[j,k] ~ poisson(blamAZ[j,k]);
    }
  }
  
  for(j in 4:NAZsamps){
    for(k in 1:3){
      blpAZ[j,k] ~ normal(mu_AZ[k], sd_blp);
      bpAZ[j,k] = inv_logit(blpAZ[j,k]);
      blamAZ[j,k] = bpAZ[j,k] * bN[ts[j],k] * AZeff[j] / 35;
      bAZ[j,k] ~ poisson(blamAZ[j,k]);
    }
  }
  
  // 2012 - 2017 NO: starts in april 2012
  for(j in 1:23){
    for(k in 1:3){
      blamNO[j,k] = bp[j,k,k] * 0.08 * bN[(j + 46),k];
      bNOc[j,k] ~ poisson(blamNO[j,k]);
    }
  }
}
    
    ", fill = TRUE)
sink()
#-----------------------------------------------------------------------------#
sm.data <- list(NAZsamps = NAZsamps, ts = ts, AZeff = AZeff, bAZ = bAZ,
                seasNO = seasNO, bNOc = bNOc, NOpasses = NOpasses, ones = FR,
                FR = FR, last = last, bCH = bCH, NCH = NCH, sumf = sumf,
                spawn = spawn)

sm.params = c('bphi', 'bpsi1', 'bpsi2', 'mu_blp', 'sd_blp', 'lbeta_0',
              'mu_I', 'I', 'Beta', 'IN', 'AZadj', 'sd_I', 'sd_lphi',
              'sd_beta', 'bN', 'bp_pass')

# MCMC settings
ni = 1000
nt = 1
nb = 500
nc = 3

SM.c <- stan("Stan_Marginalized.stan",
             data = sm.data,
             pars = sm.params,
             control = list(max_treedepth = 14, adapt_delta = .85),
             chains = nc, iter = ni, thin = nt, seed = 1) 
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Fitting an Integrated Population Model to Brown Trout Data
#  Marginalized Stan version - Removes prior on survival & random effects on
#  detection
#
#  Notes:
#  * The model runs from the fall of 2000 to fall of 2017 on a seasonal basis
#  * We define three size states based on total length in mm 
#    - 0 - 200; 200 - 350; 350 +
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

# read in data
NO_catch <- read.csv(paste0(getwd(), "/Data/", "NO_catch.csv"), header = TRUE)
AZGF_catch <- read.csv(paste0(getwd(), "/Data/", "AZGF_catch.csv"), header = TRUE)
MR_data <- read.csv(paste0(getwd(), "/Data/", "bCH.csv"), header = FALSE)

#-----------------------------------------------------------------------------#
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
sink("nop_nore.stan")
cat("
data{
  int NAZsamps;                            // Number of samples AGFD
  int ts[NAZsamps];
  vector [NAZsamps] AZeff;
  int bAZ [NAZsamps, 3];
  int seasNO[23];
  int bNOc[23,3];
  int NOpasses[23];
  int NCH;                       // Number of capture histories
  int FR[NCH];
  int last[NCH];
  int bCH[NCH, 23];
  int sumf[NCH];
  int spawn[23];               // Indicator for spawning month
}

parameters{
  matrix<lower = 0, upper = 1>[3, 4]bphi;
  real<lower = 0, upper = 1> bpsi2;
  vector<lower = 0, upper = 1>[3] bpsi1;
  vector<lower = 0, upper = 1>[4] mu_blp;
  vector<lower = 0, upper = 1>[3] mu_AZ;
  vector<lower = 0, upper = 1000>[2] IN;                              
  vector<lower = -4, upper = 8>[68] I;
  vector<lower = 0, upper = 4>[17] Beta;
}

transformed parameters{
  real btrans[4, 4, 4];                                                     
  matrix[23,3]bp_pass;
  real bp[23,4,4];                                                          
  matrix[69,3] bN;                                                          
  vector[17] wA;
  
  
  for(j in 1:23){
    bp_pass[j,1]=mu_blp[1];
    bp_pass[j,2]=mu_blp[2];
    bp_pass[j,3]=mu_blp[spawn[j]];
  }
  
  
  // define transition matrix that combines survival and growth parameters
  for(i in 1:4){
    btrans[1,3,i] = 0;
    btrans[1,4,i] = 1 - bphi[1,i];
    btrans[2,1,i] = 0;
    btrans[2,2,i] = bphi[2,i] * (1 - bpsi2);
    btrans[2,3,i] = bphi[2,i] * bpsi2;
    btrans[2,4,i] = 1 - bphi[2,i];
    btrans[3,1,i] = 0;
    btrans[3,2,i] = 0;
    btrans[3,3,i] = bphi[3,i];
    btrans[3,4,i] = 1 - bphi[3,i];
    btrans[4,1,i] = 0;
    btrans[4,2,i] = 0;
    btrans[4,3,i] = 0;
    btrans[4,4,i] = 1;
  }
  
  // size class one transitions are done separately because growth is allowed to vary between seasons
  for(i in 1:3){
    btrans[1,1,i] = bphi[1,i] * (1 - bpsi1[i]);
    btrans[1,2,i] = bphi[1,i] * bpsi1[i];
  }
  
  btrans[1,1,4] = 0;
  btrans[1,2,4] = bphi[1,4];
  
  // this loop calculates actual per pass pcaps for each trip and modifies based on # of passes
  for(j in 1:23){
    for(k in 1:3){
      bp[j,k,k] = 1 - pow((1 - bp_pass[j,k]), NOpasses[j]);                 
      bp[j,k,4] = 1 - bp[j,k,k];
    }
    
    bp[j,1,2] = 0;
    bp[j,1,3] = 0;
    bp[j,2,1] = 0;
    bp[j,2,3] = 0;
    bp[j,3,1] = 0;
    bp[j,3,2] = 0;
    bp[j,4,1] = 0;
    bp[j,4,2] = 0;
    bp[j,4,3] = 0;
    bp[j,4,4] = 1;
  }
  
  // calculate offset for each size classes of AZGF effort and calculate expected pcap
  
  bN[1,1] = 0;                                     
  bN[1,2] = IN[1];
  bN[1,3] = IN[2];
  
  //calculate latent abundance of brown trout from fall 2000 to end of 2017
  for(j in 1:17){
    for(k in 1:3){
      bN[((j-1) * 4 + k + 1),1] = btrans[1,1,k] * bN[((j-1) * 4 + k),1];
      bN[((j-1) * 4 + k + 1),2] = btrans[1,2,k] * bN[((j-1) * 4 + k),1] + btrans[2,2,k] * bN[((j-1) * 4 + k),2];
      bN[((j-1) * 4 + k + 1),3] = btrans[2,3,k] * bN[((j-1) * 4 + k),2] + btrans[3,3,k] * bN[((j-1) * 4 + k),3] + exp(I[((j-1) * 4 + k)]);
    }
    
    // BNT recruits produced in fall as a function weighted sum of adults (wA) and reprodutive rate (Beta) in winter
    wA[j] = (bN[((j - 1) * 4 + 2),2] + 4 * bN[((j - 1) * 4 + 2),3]);
    
    // between summer and fall all bnt graduate to sz 2 & recruits show up
    bN[(1 + j * 4),1] = wA[j] * Beta[j];
    bN[(1 + j * 4),2] = btrans[1,2,4] * bN[(j * 4),1] + btrans[2,2,4] * bN[(j * 4),2];
    bN[(1 + j * 4),3] = btrans[2,3,4] * bN[(j * 4),2] + btrans[3,3,4] * bN[(j * 4),3] + exp(I[(j * 4)]);
  }
}

model{
  real pz[NCH, 23, 4];                                    
  matrix[NAZsamps, 3]blamAZ;
  matrix[23, 3]blamNO;
  vector[4] temp;
  
  
  // (done above in transformed) this loop calculates actual per pass pcaps for each trip and modifies based on # of passes
  
  for(k in 1:NCH){
    pz[k,sumf[k],1] = (1 == bCH[k,sumf[k]]);
    pz[k,sumf[k],2] = (2 == bCH[k,sumf[k]]);
    pz[k,sumf[k],3] = (3 == bCH[k,sumf[k]]);
    pz[k,sumf[k],4] = 0;
    for(t in sumf[k]:(last[k] - 1)){
      for(i in 1:4){
        for(j in 1:4){
          temp[j] = pz[k,t,j] * btrans[j,i,seasNO[(t + 1)]] * bp[t,i,bCH[k,(t + 1)]];
        }
        pz[k,(t + 1),i] = sum(temp);
      }
    }
    target += FR[k] * log(sum(pz[k,last[k],]));                             
  }
  
  //////////////////////////
    
    // calculate actual immigration in each interval on log scale
  
  //////////////////////////
    // 2000 - 2017 AZGF data
  for(j in 1:3){
    for(k in 2:3){
      blamAZ[j,k] = mu_AZ[k] * bN[ts[j],k] * AZeff[j] / 35;
      bAZ[j,k] ~ poisson(blamAZ[j,k]);
    }
  }
  
  
  for(j in 4:NAZsamps){
    for(k in 1:3){
      blamAZ[j,k] = mu_AZ[k] * bN[ts[j],k] * AZeff[j] / 35;
      bAZ[j,k] ~ poisson(blamAZ[j,k]);
    }
  }
  
  
  // 2012 - 2017 NO: starts in april 2012
  for(j in 1:23){
    for(k in 1:3){
      blamNO[j,k] = bp[j,k,k] * 0.08 * bN[(j + 46),k];
      bNOc[j,k] ~ poisson(blamNO[j,k]);
    }
  }
}

    ", fill = TRUE)
sink()
#-----------------------------------------------------------------------------#
sm.data <- list(NAZsamps = NAZsamps, ts = ts, AZeff = AZeff, bAZ = bAZ,
                seasNO = seasNO, bNOc = bNOc, NOpasses = NOpasses, NCH = NCH,
                FR = FR, last = last, bCH = bCH, sumf = sumf, spawn = spawn)

sm.params <- c('bphi', 'bpsi1', 'bpsi2', "I", "Beta", "IN", "bN", "mu_AZ", "mu_blp")

# MCMC settings
ni = 10
nt = 1
nb = 5
nc = 1

# Call Stan from R 
SM.nop.nore <- stan("nop_nore.stan",
                    data = sm.data,
                    pars = sm.params,
                    control = list(max_treedepth = 14, adapt_delta = .85),
                    chains = nc, iter = ni, thin = nt, seed = 1) 

#-----------------------------------------------------------------------------#
# End
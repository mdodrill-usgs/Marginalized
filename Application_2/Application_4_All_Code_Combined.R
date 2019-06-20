###############################################################################
#                                                                     Spring 19
#  Multi-state Capture-Recapture model (Humpback Chub) 
#  JAGS Discrete Version
#
#  Notes:
#  * Need to set directory for data
#  * Note that CH's are formatted so that each row only has the release and
#    subsequent capture information
#
#	   Thus, a CH that is normally : 1 0 1 1 0 0  
#
#	   Would be formatted as such: 1 0 1 0 0 0 (fc = 1, lc =3)
#			       	      0 0 1 1 0 0 (fc = 3, lc = 4)
#			       	      0 0 0 1 0 0 (fc = 4, lc = 6)
#
#		 where fc = release occasion, lc = capture occasion
#
#  * States are defined as such:
#		 1 = small adult captured in LCR
#		 2 = large adult captured in LCR
#		 3 = small adult captured in CR
#		 4 = large adult captured in CR
#		 5 = not observed
#
#	 * Fish are allowed to move between CR and LCR 
#	 * Fish are allowed to grow from small to large adults
#	 * Fish are not allowed to shrink (from large to small adults)
#
###############################################################################
library(R2jags)

# Get data:
dat <- read.csv(".//HBC_data.csv")
#-----------------------------------------------------------------------------#
# Functions:

#	Get first capture occasion:
find.first <- function(x){min(which(x != 5))}

#	Get last capture occasion:
find.last <- function(x){
  ifelse(length(which(x !=5 )) == 1, 27, max(which(x != 5)))
  }

#-----------------------------------------------------------------------------#
# format data for model fitting:
sumCH <- dat[,1:27]
sumFR <- dat[,28]
newtag.AO <- dat[,29]

# other data:
season <- c(rep(1:3,v8),v1:2) # season index
LCRs <- c(2, 3, 5, 6, 8, 9, 11:26)
LCRns <- c(1, 4, 7, 10)
CRs <- c(1, 2, 4, 5, 7:26)
CRns <- c(3, 6)
catch <- matrix(c(287, 545, 221, 594, 215, 413, 500, 374, 562,
                  193, 356, 171, 191, 111, 239, 259, 283, 246,
                  102, 129, 154, 54, 63, 48, 36, 41, 35,
                  35, 86, 126, 41, 43, 59, 22, 26, 34), nrow = 4, byrow = TRUE)


# Re-format CH so that lines correspond to individuals (not unique capture histories):
indCH <- numeric()
for(i in 1:length(sumCH[,1])){
  for(j in 1:sumFR[i]){
    indCH <- rbind(indCH, sumCH[i,])
  }
}

indCH <- as.matrix(indCH)
indf <- apply(indCH, 1, find.first)
indl <- apply(indCH, 1, find.last)

# Create CH data matrix where state is known when fish is captured, otherwise NA:
known.state <- indCH
known.state[indCH == 5] <- NA
known.state[,1] <- NA
for(i in 1:dim(indCH)[1]){known.state[i,indf[i]] <- NA}

# Create initial values for z-states (indicates whether fish is a small adult in
# LCR (1), a large adult in LCR (2), a small adult in CR (3), or a large adult
# in CR (4)
# Note code does not allow fish to shrink, so initial values for z-states must
# also not allow for shrinking
z.init <- matrix(NA, nrow = dim(indCH)[1], ncol = 27)
for(j in 1:dim(indCH)[1]){
  z.init[j,(indf[j] + 1)] <- ifelse(indCH[j,(indf[j] + 1)] < 5, NA, indCH[j,indf[j]])
  if(indf[j] < 26 & (indl[j] - indf[j]) > 1){
    for(i in (indf[j]+2):indl[j]){
      size <- c(1, 2, 1, 2, NA)
      sizemax <- max(size[indCH[j,indf[j]:i]], na.rm = TRUE)
      if(sizemax == 1){x <- c(1, 3)}
      if(sizemax == 2){x <- c(2, 4)}
      z.init[j,i] <- ifelse(indCH[j,i] < 5, NA, sample(x, 1))
    }
  }
}

#-----------------------------------------------------------------------------#
sink("JAGS_Discrete_chub.txt")
cat("
model{
  
  for(j in 1:4){
    s_i[j] ~ dunif(0,1) # 3 month survival rates for small and large adults in LCR and CR
    s[1,j] <- s_i[j]    # 3 month survival for spring to summer interval
    s[2,j] <- s_i[j]    # 3 month survival for summer to fall interval
    s[3,j] <- s_i[j] * s_i[j] # 6 month survival for fall to spring interval
  } 
  
  tau ~ dunif(0,1) # proportion of adults in the observable portion of the LCR
  
  for (i in 1:3){
    g[i,1] ~ dunif(0,1) # LCR growth for three intervals
    g[i,2] ~ dunif(0,1) # CR growth for three intervals
    m[i,1] ~ dunif(0,1) # movement out of LCR by interval for small adults
    m[i,2] ~ dunif(0,1) # movement out of LCR by interval for large adults
    m[i,3] ~ dunif(0,1) # movement into the LCR by interval for small adults
    m[i,4] ~ dunif(0,1) # movement into the LCR by interval for large adults
    
    # tr is an array representing markov transitions for three intervals. some
    # transitions are assumed to be zero (e.g., transitioning from a large to a
    # small adult)
    tr[1,i,1] <- s[i,1] * (1 - g[i,1]) * (1 - m[i,1])
    tr[1,i,2] <- s[i,1] * g[i,1] * (1 - m[i,2])
    tr[1,i,3] <- s[i,1] * (1 - g[i,1]) * m[i,1] * tau
    tr[1,i,4] <- s[i,1] * g[i,1] * m[i,2] * tau
    tr[1,i,5] <- s[i,1] * (1 - g[i,1]) * m[i,1] * (1 - tau)
    tr[1,i,6] <- s[i,1] * g[i,1] * m[i,2] * (1 - tau)
    tr[1,i,7] <- 1 - s[i,1]
    tr[2,i,1] <- 0
    tr[2,i,2] <- s[i,2] * (1 - m[i,2])
    tr[2,i,3] <- 0
    tr[2,i,4] <- s[i,2] * m[i,2] * tau
    tr[2,i,5] <- 0
    tr[2,i,6] <- s[i,2] * m[i,2] * (1 - tau)
    tr[2,i,7] <- 1 - s[i,2]
    tr[3,i,1] <- s[i,3] * (1 - g[i,2]) * m[i,3]
    tr[3,i,2] <- s[i,3] * g[i,2] * m[i,4]
    tr[3,i,3] <- s[i,3] * (1 - g[i,2]) * (1 - m[i,3])
    tr[3,i,4] <- s[i,3] * g[i,2] * (1 - m[i,4])
    tr[3,i,5] <- 0
    tr[3,i,6] <- 0
    tr[3,i,7] <- 1 - s[i,3]
    tr[4,i,1] <- 0
    tr[4,i,2] <- s[i,4] * m[i,4]
    tr[4,i,3] <- 0
    tr[4,i,4] <- s[i,4] * (1 - m[i,4])
    tr[4,i,5] <- 0
    tr[4,i,6] <- 0
    tr[4,i,7] <- 1 - s[i,4]
    tr[5,i,1] <- tr[3,i,1]
    tr[5,i,2] <- tr[3,i,2]
    tr[5,i,3] <- 0
    tr[5,i,4] <- 0
    tr[5,i,5] <- tr[3,i,3]
    tr[5,i,6] <- tr[3,i,4]
    tr[5,i,7] <- tr[3,i,7]
    tr[6,i,1] <- 0
    tr[6,i,2] <- tr[4,i,2]
    tr[6,i,3] <- 0
    tr[6,i,4] <- 0
    tr[6,i,5] <- 0
    tr[6,i,6] <- tr[4,i,4]
    tr[6,i,7] <- tr[4,i,7]
    for(j in 1:6){
      tr[7,i,j] <- 0
    }
    tr[7,i,7] <- 1
  }	
  
  # p is the matrix describing capture probabilities
  # Diagonal values represent capture probabilities for the four states
  # The fifth column represents the probability a fish is not captured
  # Note that for fish in states 5-6 (unobserveable sites in the CR that are not sampled)
  # and for state 7 (dead state) the probability a fish is unobserved is 100%
  for(i in 1:22){
    p[1,LCRs[i],1] ~ dunif(0, 1)
    p[2,LCRs[i],2] ~ dunif(0, 1)
  }
  
  # No sampling in LCR during 1st, 4th, 7th, 10th recap periods (summers of 2009 - 2012)
  for(i in 1:4){
    p[1,LCRns[i],1] <- 0
    p[2,LCRns[i],2] <- 0
  }
  for(i in 1:24){
    p[3,CRs[i],3] ~ dunif(0, 1)
    p[4,CRs[i],4] ~ dunif(0, 1)
  }
  
  # No sampling in CR during 3rd and 6th recap periods (spring of 2010 & 2011)
  for(i in 1:2){
    p[3,CRns[i],3] <- 0
    p[4,CRns[i],4] <- 0
  }
  
  for(i in 1:26){
    p[1,i,2] <- 0
    p[1,i,3] <- 0
    p[1,i,4] <- 0
    p[1,i,5] <- 1 - p[1,i,1]
    p[2,i,1] <- 0
    p[2,i,3] <- 0
    p[2,i,4] <- 0
    p[2,i,5] <- 1 - p[2,i,2]
    p[3,i,1] <- 0
    p[3,i,2] <- 0
    p[3,i,4] <- 0
    p[3,i,5] <- 1 - p[3,i,3]
    p[4,i,1] <- 0
    p[4,i,2] <- 0
    p[4,i,3] <- 0
    p[4,i,5] <- 1 - p[4,i,4]
    p[5,i,1] <- 0
    p[5,i,2] <- 0
    p[5,i,3] <- 0
    p[5,i,4] <- 0
    p[5,i,5] <- 1
    p[6,i,1] <- 0
    p[6,i,2] <- 0
    p[6,i,3] <- 0
    p[6,i,4] <- 0
    p[6,i,5] <- 1
    p[7,i,1] <- 0
    p[7,i,2] <- 0
    p[7,i,3] <- 0
    p[7,i,4] <- 0
    p[7,i,5] <- 1
  }
  
  for(k in 1:NindCH){
    z[k,indf[k]] = indCH[k,indf[k]]
    for(i in (indf[k] + 1):indl[k]){
      z[k,i] ~ dcat(tr[z[k,i - 1], season[i - 1],])
      indCH[k,i] ~ dcat(p[z[k,i], i - 1,])
    }
  }
}

    ",fill=TRUE)
sink() 

#-----------------------------------------------------------------------------#
JD.data <- list(LCRs = LCRs, LCRns = LCRns, CRs = CRs, CRns = CRns,
                season = season, NindCH = dim(indCH)[1], indCH = indCH,
                indf = indf, indl = indl, z = known.state)

JD.par <- c('s_i', 'g', 'm', 'tau')

JD.inits <- function(){
  list(s_i = c(.7, .7, .7, .7), z = array(z.init, dim = c(dim(indCH)[1], 27)))
  }

JD.out <- jags.parallel(JD.data, inits = JD.inits, JD.par, ".\\JAGS_Discrete_chub.txt",
                        n.chains = 1, n.iter = 10, export_obj_names = "z.init")

#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
# Multi-state Capture-Recapture model (Humpback Chub)
# JAGS Marginalized Version
#
#  Notes:
#  * Need to set directory for data
#  * See Notes above
###############################################################################
library(R2jags)

# Get data:
dat <- read.csv(".//HBC_data.csv")
#-----------------------------------------------------------------------------#
#  Functions:
#	Get first capture occasion:
find.first <- function(x){min(which(x != 5))}

#	Get last capture occasion:
find.last <- function(x){
  ifelse(length(which(x != 5)) == 1, 27, max(which(x != 5)))
  }

#-----------------------------------------------------------------------------#
# format data for model fitting:
sumCH <- as.matrix(dat[,1:27])
sumFR <- dat[,28]
newtag.AO <- dat[,29]

# other data:
season <- c(rep(1:3, 8), 1:2) # season index
LCRs <- c(2, 3, 5, 6, 8, 9, 11:26)
LCRns <- c(1, 4, 7, 10)
CRs <- c(1, 2, 4, 5, 7:26)
CRns <- c(3, 6)
catch <- matrix(c(287, 545, 221, 594, 215, 413, 500, 374, 562,
                 193, 356, 171, 191, 111, 239, 259, 283, 246,
                 102, 129, 154, 54, 63, 48, 36, 41, 35,
                 35, 86, 126, 41, 43, 59, 22, 26, 34), nrow = 4, byrow = TRUE)

# Get first and last capture occasion for line in summarized capture history
fc <- apply(sumCH, 1, find.first)
lc <- apply(sumCH, 1, find.last)

# Use this to calculate abundance in code:
CR_ind <- matrix(0, nrow = 4, ncol = 9)
1 -> CR_ind[3:4,]

#-----------------------------------------------------------------------------#
sink("JAGS_Marginalized_chub.txt")
cat("
model{
  
  for(j in 1:4){
    s_i[j] ~ dunif(0, 1) # 3 month survival rates for small and large adults in LCR and CR
    s[1,j] <- s_i[j]     # 3 month survival for spring to summer interval
    s[2,j] <- s_i[j]     # 3 month survival for summer to fall interval
    s[3,j] <- s_i[j] * s_i[j] # 6 month survival for fall to spring interval
  } 
  
  tau ~ dunif(0, 1) # proportion of adults in the observable portion of the LCR
  for(i in 1:3){
    g[i,1] ~ dunif(0, 1) # LCR growth for three intervals
    g[i,2] ~ dunif(0, 1) # CR growth for three intervals
    m[i,1] ~ dunif(0, 1) # movement out of LCR by interval for small adults
    m[i,2] ~ dunif(0, 1) # movement out of LCR by interval for large adults
    m[i,3] ~ dunif(0, 1) # movement into the LCR by interval for small adults
    m[i,4] ~ dunif(0, 1) # movement into the LCR by interval for large adults
    
    # tr is an array representing markov transitions for three intervals. some
    # transitions are assumed to be zero (e.g., transitioning from a large to a
    # small adult)
    tr[1,i,1] <- s[i,1] * (1 - g[i,1]) * (1 - m[i,1])
    tr[1,i,2] <- s[i,1] * g[i,1] * (1 - m[i,2])
    tr[1,i,3] <- s[i,1] * (1 - g[i,1]) * m[i,1] * tau
    tr[1,i,4] <- s[i,1] * g[i,1] * m[i,2] * tau
    tr[1,i,5] <- s[i,1] * (1 - g[i,1]) * m[i,1] * (1 - tau)
    tr[1,i,6] <- s[i,1] * g[i,1] * m[i,2] * (1 - tau)
    tr[1,i,7] <- 1 - s[i,1]
    tr[2,i,1] <- 0
    tr[2,i,2] <- s[i,2] * (1 - m[i,2])
    tr[2,i,3] <- 0
    tr[2,i,4] <- s[i,2] * m[i,2] * tau
    tr[2,i,5] <- 0
    tr[2,i,6] <- s[i,2] * m[i,2] * (1 - tau)
    tr[2,i,7] <- 1 - s[i,2]
    tr[3,i,1] <- s[i,3] * (1 - g[i,2]) * m[i,3]
    tr[3,i,2] <- s[i,3] * g[i,2] * m[i,4]
    tr[3,i,3] <- s[i,3] * (1 - g[i,2]) * (1 - m[i,3])
    tr[3,i,4] <- s[i,3] * g[i,2] * (1 - m[i,4])
    tr[3,i,5] <- 0
    tr[3,i,6] <- 0
    tr[3,i,7] <- 1 - s[i,3]
    tr[4,i,1] <- 0
    tr[4,i,2] <- s[i,4] * m[i,4]
    tr[4,i,3] <- 0
    tr[4,i,4] <- s[i,4] * (1 - m[i,4])
    tr[4,i,5] <- 0
    tr[4,i,6] <- 0
    tr[4,i,7] <- 1 - s[i,4]
    tr[5,i,1] <- tr[3,i,1]
    tr[5,i,2] <- tr[3,i,2]
    tr[5,i,3] <- 0
    tr[5,i,4] <- 0
    tr[5,i,5] <- tr[3,i,3]
    tr[5,i,6] <- tr[3,i,4]
    tr[5,i,7] <- tr[3,i,7]
    tr[6,i,1] <- 0
    tr[6,i,2] <- tr[4,i,2]
    tr[6,i,3] <- 0
    tr[6,i,4] <- 0
    tr[6,i,5] <- 0
    tr[6,i,6] <- tr[4,i,4]
    tr[6,i,7] <- tr[4,i,7]
    for (j in 1:6){
      tr[7,i,j] <- 0
    }
    tr[7,i,7] <- 1
  }	
  
  # p is the matrix describing capture probabilities
  # Diagonal values represent capture probabilities for the four states
  # The fifth column represents the probability a fish is not captured
  # Note that for fish in states 5-6 (unobserveable sites in the CR that are not sampled)
  # and for state 7 (dead state) the probability a fish is unobserved is 100%
  for(i in 1:22){
    p[1,LCRs[i],1] ~ dunif(0, 1)
    p[2,LCRs[i],2] ~ dunif(0, 1)
  }
  
  # No sampling in LCR during 1st, 4th, 7th, 10th recap periods (summers of 2009 - 2012)
  for(i in 1:4){
    p[1,LCRns[i],1] <- 0
    p[2,LCRns[i],2] <- 0
  }
  for(i in 1:24){
    p[3,CRs[i],3] ~ dunif(0, 1)
    p[4,CRs[i],4] ~ dunif(0, 1)
  }
  
  # No sampling in CR during 3rd and 6th recap periods (spring of 2010 & 2011)
  for(i in 1:2){
    p[3,CRns[i],3] <- 0
    p[4,CRns[i],4] <- 0
  }
  
  for(i in 1:26){
    p[1,i,2] <- 0
    p[1,i,3] <- 0
    p[1,i,4] <- 0
    p[1,i,5] <- 1 - p[1,i,1]
    p[2,i,1] <- 0
    p[2,i,3] <- 0
    p[2,i,4] <- 0
    p[2,i,5] <- 1 - p[2,i,2]
    p[3,i,1] <- 0
    p[3,i,2] <- 0
    p[3,i,4] <- 0
    p[3,i,5] <- 1 - p[3,i,3]
    p[4,i,1] <- 0
    p[4,i,2] <- 0
    p[4,i,3] <- 0
    p[4,i,5] <- 1 - p[4,i,4]
    p[5,i,1] <- 0
    p[5,i,2] <- 0
    p[5,i,3] <- 0
    p[5,i,4] <- 0
    p[5,i,5] <- 1
    p[6,i,1] <- 0
    p[6,i,2] <- 0
    p[6,i,3] <- 0
    p[6,i,4] <- 0
    p[6,i,5] <- 1
    p[7,i,1] <- 0
    p[7,i,2] <- 0
    p[7,i,3] <- 0
    p[7,i,4] <- 0
    p[7,i,5] <- 1
  }
  
  # The 0.03 is to account for one-time 3% tag loss
  for(k in 1:NsumCH){
    pz[k,sumf[k],1] <- equals(sumCH[k,sumf[k]], 1) * (1 - 0.03 * newtag[k]) 
    pz[k,sumf[k],2] <- equals(sumCH[k,sumf[k]], 2) * (1 - 0.03 * newtag[k])
    pz[k,sumf[k],3] <- equals(sumCH[k,sumf[k]], 3) * (1 - 0.03 * newtag[k])
    pz[k,sumf[k],4] <- equals(sumCH[k,sumf[k]], 4) * (1 - 0.03 * newtag[k])
    pz[k,sumf[k],5] <- 0
    pz[k,sumf[k],6] <- 0
    pz[k,sumf[k],7] <- 0.03 * newtag[k]
    for(i in sumf[k]:(lc[k] - 1)){
      for(j in 1:7){
        pz[k,(i + 1),j] <- inprod(pz[k,i,], tr[,season[i],j]) * p[j,i,sumCH[k,(i + 1)]]
      }
    }
    
    lik[k] <- sum(pz[k,lc[k],])
    one[k] ~ dbin(lik[k], sumFR[k])
  }
  
  # Calculate abundance from catch by simulating from a negative binomial distribution:
  for(i in 1:4){
    for(j in 1:9){
      p_fall[i,j] <- p[i,fall_ind[j],i] * (1 - CR_ind[i,j] * (1 - tau))
      U[i,j] ~ dnegbin(p_fall[i,j], catch[i,j])
      N[i,j] = catch[i,j] + U[i,j] 
    }
  }
}

    ",fill=TRUE)
sink() 

#-----------------------------------------------------------------------------#
JM.data <- list(LCRs = LCRs, LCRns = LCRns, CRs = CRs, CRns = CRns,
                season = season, NsumCH = dim(sumCH)[1],
                sumCH = array(sumCH, dim = c(dim(sumCH)[1], dim(sumCH)[2])),
                newtag = as.vector(newtag.AO), sumf = as.vector(fc),
                sumFR = sumFR, one = sumFR,  
                CR_ind = array(CR_ind, dim = c(dim(CR_ind)[1], dim(CR_ind)[2])),
                lc = as.vector(lc), fall_ind = 1:9 * 3 - 1,
                catch = array(catch, dim = c(dim(catch)[1], dim(catch)[2])))

JM.par <- c('s_i', 'g', 'm', 'tau', 'p_cr', 'p_lcr', 'N')

JM.inits <- function(){list(s_i = c(.7, .7, .7, .7))}

JM.out <- jags.parallel(JM.data, inits = JM.inits, JM.par,
                        ".\\JAGS_Marginalized_chub.txt", n.cluster = 3,
                        n.chains = 1, n.iter = 10)
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
# Multi-state Capture-Recapture model (Humpback Chub)
# Stan Marginalized Version
#
#  Notes:
#  * Need to set directory for data
#  * See Notes above
###############################################################################
library(rstan)

# To run Stan in parallel:
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Get data:
dat <- read.csv(".//HBC_data.csv")
#-----------------------------------------------------------------------------#
# Functions:
#	Get first capture occasion:
find.first <- function(x){min(which(x != 5))}

#	Get last capture occasion:
find.last <- function(x){
  ifelse(length(which(x != 5)) == 1, 27, max(which(x != 5)))
  }

#-----------------------------------------------------------------------------#
# format data for model fitting:
sumCH <- as.matrix(dat[,1:27])
sumFR <- dat[,28]
newtag.AO <- dat[,29]

# other data:
season <- c(rep(1:3, 8), 1:2) # season index
LCRs <- c(2, 3, 5, 6, 8, 9, 11:26)
LCRns <- c(1, 4, 7, 10)
CRs <- c(1, 2, 4, 5, 7:26)
CRns <- c(3, 6)
catch <- matrix(c(287, 545, 221, 594, 215, 413, 500, 374, 562,
                 193, 356, 171, 191, 111, 239, 259, 283, 246,
                 102, 129, 154, 54, 63, 48, 36, 41, 35,
                 35, 86, 126, 41, 43, 59, 22, 26, 34), nrow = 4, byrow = TRUE)

# Get first and last capture occasion for line in summarized capture history
fc <- apply(sumCH, 1, find.first)
lc <- apply(sumCH, 1, find.last)

# Use this to calculate abundance in code:
CR_ind <- matrix(0, nrow = 4, ncol = 9)
1 -> CR_ind[3:4,]
#-----------------------------------------------------------------------------#
sink("Stan_Marginalized_chub.stan")
cat("
data{
  int<lower = 1> NsumCH;
  int<lower = 1, upper = 5> sumCH[NsumCH, 27];
  int<lower = 1, upper = 26> sumf[NsumCH];
  int<lower = 1, upper = 3> season[26];
  int<lower = 1> sumFR[NsumCH];
  int<lower = 1> LCRs [22];
  int<lower = 1> LCRns [4];
  int<lower = 1> CRs [24];
  int<lower = 1> CRns [2];
  int catch_mat[4,9];
  int fall_ind[9];
  int<lower = 0, upper = 1> newtag[NsumCH];
  int<lower = 2, upper = 27> lc[NsumCH];
}

parameters{
  real<lower = 0, upper = 1> s_i[4]; // 3 month survivals for small and large chub in LCR and CR
  real<lower = 0, upper = 1> tau;    // proportion of CR adults residing in observable location
  vector<lower = 0, upper = 1> [2] g [3]; // seasonal growth for two rivers
  vector<lower = 0, upper = 1> [4] m [3]; // seasonal movement rates for 2 sizes and locations
  vector<lower = 0, upper = 1> [2] p_lcr [22]; // recapture probability in LCR for two size classes
  vector<lower = 0, upper = 1> [2] p_cr [24];  // recapture probability in CR for two size classes
}

transformed parameters{
  vector<lower = 0, upper = 1> [4] s[3];
  simplex[7] tr[7,3];
  simplex[5] p[7,26];
  
  for(j in 1:4){
    s[1,j] = s_i[j]; // 3 month survival for spring to summer interval
    s[2,j] = s_i[j]; // 3 month survival for summer to fall interval
    s[3,j] = s_i[j] * s_i[j]; // 6 month survival for fall to spring interval
  } 
  
  for(i in 1:3){
    tr[1,i,1] = s[i,1] * (1 - g[i,1]) * (1 - m[i,1]);
    tr[1,i,2] = s[i,1] * g[i,1] * (1 - m[i,2]);
    tr[1,i,3] = s[i,1] * (1 - g[i,1]) * m[i,1] * tau;
    tr[1,i,4] = s[i,1] * g[i,1] * m[i,2] * tau;
    tr[1,i,5] = s[i,1] * (1 - g[i,1]) * m[i,1] * (1 - tau);
    tr[1,i,6] = s[i,1] * g[i,1] * m[i,2] * (1 - tau);
    tr[1,i,7] = 1 - s[i,1];
    tr[2,i,1] = 0;
    tr[2,i,2] = s[i,2] * (1 - m[i,2]);
    tr[2,i,3] = 0;
    tr[2,i,4] = s[i,2] * m[i,2] * tau;
    tr[2,i,5] = 0;
    tr[2,i,6] = s[i,2] * m[i,2] * (1-tau);
    tr[2,i,7] = 1 - s[i,2];
    tr[3,i,1] = s[i,3] * (1 - g[i,2]) * m[i,3];
    tr[3,i,2] = s[i,3] * g[i,2] * m[i,4];
    tr[3,i,3] = s[i,3] * (1 - g[i,2]) * (1 - m[i,3]);
    tr[3,i,4] = s[i,3] * g[i,2] * (1 - m[i,4]);
    tr[3,i,5] = 0;
    tr[3,i,6] = 0;
    tr[3,i,7] = 1 - s[i,3];
    tr[4,i,1] = 0;
    tr[4,i,2] = s[i,4] * m[i,4];
    tr[4,i,3] = 0;
    tr[4,i,4] = s[i,4] * (1 - m[i,4]);
    tr[4,i,5] = 0;
    tr[4,i,6] = 0;
    tr[4,i,7] = 1 - s[i,4];
    tr[5,i,1] = tr[3,i,1];
    tr[5,i,2] = tr[3,i,2];
    tr[5,i,3] = 0;
    tr[5,i,4] = 0;
    tr[5,i,5] = tr[3,i,3];
    tr[5,i,6] = tr[3,i,4];
    tr[5,i,7] = tr[3,i,7];
    tr[6,i,1] = 0;
    tr[6,i,2] = tr[4,i,2];
    tr[6,i,3] = 0;
    tr[6,i,4] = 0;
    tr[6,i,5] = 0;
    tr[6,i,6] = tr[4,i,4];
    tr[6,i,7] = tr[4,i,7];
    for(j in 1:6){
      tr[7,i,j] = 0;
    }
    tr[7,i,7] = 1;
  }
  
  for(i in 1:22){
    p[1,LCRs[i],1] = p_lcr[i,1];
    p[2,LCRs[i],2] = p_lcr[i,2];
  }
  // No sampling in LCR during 1st, 4th, 7th, 10th recap periods (summers of 2009 - 2012)
  for(i in 1:4){
    p[1,LCRns[i],1] = 0;
    p[2,LCRns[i],2] = 0;
  }
  for(i in 1:24){
    p[3,CRs[i],3] = p_cr[i,1];
    p[4,CRs[i],4] = p_cr[i,2];
  }
  // No sampling in CR during 3rd and 6th recap periods (spring of 2010 & 2011)
  for(i in 1:2){
    p[3,CRns[i],3] = 0;
    p[4,CRns[i],4] = 0;
  }

  for(i in 1:26){
    p[1,i,2] = 0;
    p[1,i,3] = 0;
    p[1,i,4] = 0;
    p[1,i,5] = 1 - p[1,i,1];
    p[2,i,1] = 0;
    p[2,i,3] = 0;
    p[2,i,4] = 0;
    p[2,i,5] = 1 - p[2,i,2];
    p[3,i,1] = 0;
    p[3,i,2] = 0;
    p[3,i,4] = 0;
    p[3,i,5] = 1 - p[3,i,3];
    p[4,i,1] = 0;
    p[4,i,2] = 0;
    p[4,i,3] = 0;
    p[4,i,5] = 1 - p[4,i,4];
    p[5,i,1] = 0;
    p[5,i,2] = 0;
    p[5,i,3] = 0;
    p[5,i,4] = 0;
    p[5,i,5] = 1;
    p[6,i,1] = 0;
    p[6,i,2] = 0;
    p[6,i,3] = 0;
    p[6,i,4] = 0;
    p[6,i,5] = 1;
    p[7,i,1] = 0;
    p[7,i,2] = 0;
    p[7,i,3] = 0;
    p[7,i,4] = 0;
    p[7,i,5] = 1;
  }
}

model{
  real temp[7]; 
  vector[7] pz[27]; 
  for(k in 1:NsumCH){ 
    for(j in 1:6){
      pz[sumf[k], j] = (1 - 0.03 * newtag[k]) * (j == sumCH[k, sumf[k]]); 
    } 
    pz[sumf[k], 7] = (0.03 * newtag[k]);  
    for(t in (sumf[k] + 1):lc[k]){ 
      for(i in 1:7){ 
        for(j in 1:7){
          temp[j] = pz[t - 1, j] * tr[j, season[(t - 1)], i] * p[i, t - 1, sumCH[k, t]]; 
        }
        pz[t,i] = sum(temp); 
      } 
    }
    target += sumFR[k] * log(sum(pz[lc[k]])); 
  }
}

// Code to estimate abundance by simulating from a negative binomial distribution
generated quantities{
  real ptrans;	
  real<lower = 0> scale_par;
  int U[4,9];
  int N[4,9];
  
  for(i in 1:4){
    for(j in 1:9){
      ptrans = p[i,fall_ind[j],i] * (1 - (i > 2) * (1 - tau));
      scale_par = ptrans / (1 - ptrans);
      U[i,j] = neg_binomial_rng(catch_mat[i,j], scale_par);
      N[i,j] = U[i,j] + catch_mat[i,j];
    }
  }
}
    
    ",fill=TRUE)
sink() 
#-----------------------------------------------------------------------------#
sm.data<-list(NsumCH = dim(sumCH)[1], newtag = as.vector(newtag.AO),
              sumCH = array(sumCH,  dim = c(dim(sumCH)[1], dim(sumCH)[2])), 
              sumf = fc, season = season, sumFR = sumFR, LCRs = LCRs,
              LCRns = LCRns, CRs = CRs, CRns = CRns,  
              catch_mat = array(catch,  dim = c(dim(catch)[1], dim(catch)[2])),
              lc = lc,  fall_ind = 1:9 * 3 - 1)

sm.inits <- function() list(s_i = runif(4, 0, 1))  
sm.params <- c("s_i", "g", "m", "tau", "p_cr", "p_lcr", "N") 

SM <- stan(".\\Stan_Marginalized_chub.stan", 
           data = sm.data, init = sm.inits, pars = sm.params, 
           chains = 1, iter = 10, thin = 1, 
           seed = 1) 
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
# Multi-state Capture-Recapture model (Humpback Chub)
# Stan Marginalized Version with random effects
#
#  Notes:
#  * Need to set directory for data
#  * See Notes above
###############################################################################
library(rstan)

# To run Stan in parallel:
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Get data:
dat <- read.csv(".//HBC_data.csv")
#-----------------------------------------------------------------------------#
#  Functions:
#	Get first capture occasion:
find.first <- function(x){min(which(x != 5))}

#	Get last capture occasion:
find.last <- function(x){
  ifelse(length(which(x != 5)) == 1, 27, max(which(x != 5)))
  }
#-----------------------------------------------------------------------------#
# format data for model fitting:
sumCH <- as.matrix(dat[,1:27])
sumFR <- dat[,28]
newtag.AO <- dat[,29]

# other data:
season <- c(rep(1:3, 8), 1:2) # season index
LCRs <- c(2, 3, 5, 6, 8, 9, 11:26)
LCRns <- c(1, 4, 7, 10)
CRs <- c(1, 2, 4, 5, 7:26)
CRns <- c(3, 6)
catch <- matrix(c(287, 545, 221, 594, 215, 413, 500, 374, 562,
                 193, 356, 171, 191, 111, 239, 259, 283, 246,
                 102, 129, 154, 54, 63, 48, 36, 41, 35,
                 35, 86, 126, 41, 43, 59, 22, 26, 34), nrow = 4, byrow = TRUE)

# Get first and last capture occasion for line in summarized capture history
fc <- apply(sumCH, 1, find.first)
lc <- apply(sumCH, 1, find.last)

# Use this to calculate abundance in code:
CR_ind <- matrix(0, nrow = 4, ncol = 9)
1 -> CR_ind[3:4,]

#-----------------------------------------------------------------------------#
sink("Stan_Marginalized_chub_RE.stan")
cat("
data{
  int<lower = 1> NsumCH;
  int<lower = 1, upper = 5> sumCH[NsumCH,  27];
  int<lower = 1, upper = 26> sumf[NsumCH];
  int<lower = 1, upper = 3> season[26];
  int<lower = 1> sumFR[NsumCH];
  int<lower = 1> LCRs [22];
  int<lower = 1> LCRns [4];
  int<lower = 1> CRs [24];
  int<lower = 1> CRns [2];
  int catch_mat[4, 9];
  int fall_ind[9];
  int<lower = 2,  upper = 27> lc[NsumCH];
  int<lower = 0,  upper = 1> newtag[NsumCH];
}

parameters {
  real<lower = -8, upper = 8> mu_ls[4]; // hyperprior on 3 month survivals for small and large chub in LCR and CR
  real<lower = 0> sd_ls; // hyperprior on 3 month survivals for small and large chub in LCR and CR
  real<lower = 0, upper = 1> tau; // proportion of CR adults residing in observable location
  vector<lower = -8, upper = 8> [2] mu_lg [3]; // hyperprior on seasonal growth for two rivers
  vector<lower = -8, upper = 8> [4] mu_lm [3]; // hyperprior on seasonal movement rates for 2 sizes and locations
  real<lower = 0> sd_lg ; // hyperprior on seasonal growth for two rivers
  real<lower = 0> sd_lm ; // hyperprior on seasonal movement rates for 2 sizes and locations
  vector<lower = 0, upper = 1> [2] p_lcr [22]; // recapture probability in LCR for two size classes
  vector<lower = 0, upper = 1> [2] p_cr[24]; // recapture probability in CR for two size classes
  vector [4] z_ls [26]; 
  vector [4] z_lm [26]; 
  vector [2] z_lg [26];
}

transformed parameters {
  vector<lower = 0, upper = 1> [4] s [26]; // survival for all intervals, 2 locations, and 2 size classes
  vector<lower = 0, upper = 1> [2] g [26]; // growth for all intervals and 2 size classes
  vector<lower = 0, upper = 1> [4] m [26]; // movement for all intervals, 2 locations, and 2 size classes
  simplex[7] tr[7,26];
  simplex[5] p[7,26];
  vector [4] ls [26]; // logit survival for all intervals, 2 locations, and 2 size classes
  vector [2] lg [26]; // logit growth for all intervals and 2 size classes
  vector [4] lm [26]; // logit movement for all intervals, 2 locations, and 2 size classes
  
  for(j in 1:4){
    for(i in 1:26){
      ls[i,j] = mu_ls[j] + z_ls[i,j] * sd_ls;
    }
  }
  for(j in 1:4){
    for(i in 1:26){
      lm[i,j] = mu_lm[season[i],j] + z_lm[i,j] * sd_lm;
    }
  }
  for(j in 1:2){
    for(i in 1:26){
      lg[i,j] = mu_lg[season[i],j] + z_lg[i,j] * sd_lg;
    }
  }   
  
  for(j in 1:2) g[,j] = inv_logit(lg[,j]);
  for(j in 1:4) m[,j] = inv_logit(lm[,j]);
  for(j in 1:4){
    for(i in 1:26){
      s[i,j] = (inv_logit(ls[i,j]))^(season[i]>2 ? 2 : 1);
    }
  }
  
  // state transition matrix:
    for(i in 1:26){
      tr[1,i,1] = s[i,1] * (1 - g[i,1]) * (1 - m[i,1]);
      tr[1,i,2] = s[i,1] * g[i,1] * (1 - m[i,2]);
      tr[1,i,3] = s[i,1] * (1 - g[i,1]) * m[i,1] * tau;
      tr[1,i,4] = s[i,1] * g[i,1] * m[i,2] * tau;
      tr[1,i,5] = s[i,1] * (1 - g[i,1]) * m[i,1] * (1 - tau);
      tr[1,i,6] = s[i,1] * g[i,1] * m[i,2] * (1 - tau);
      tr[1,i,7] = 1 - s[i,1];
      tr[2,i,1] = 0;
      tr[2,i,2] = s[i,2] * (1 - m[i,2]);
      tr[2,i,3] = 0;
      tr[2,i,4] = s[i,2] * m[i,2] * tau;
      tr[2,i,5] = 0;
      tr[2,i,6] = s[i,2] * m[i,2] * (1 - tau);
      tr[2,i,7] = 1 - s[i,2];
      tr[3,i,1] = s[i,3] * (1 - g[i,2]) * m[i,3];
      tr[3,i,2] = s[i,3] * g[i,2] * m[i,4];
      tr[3,i,3] = s[i,3] * (1 - g[i,2]) * (1 - m[i,3]);
      tr[3,i,4] = s[i,3] * g[i,2] * (1 - m[i,4]);
      tr[3,i,5] = 0;
      tr[3,i,6] = 0;
      tr[3,i,7] = 1 - s[i,3];
      tr[4,i,1] = 0;
      tr[4,i,2] = s[i,4] * m[i,4];
      tr[4,i,3] = 0;
      tr[4,i,4] = s[i,4] * (1 - m[i,4]);
      tr[4,i,5] = 0;
      tr[4,i,6] = 0;
      tr[4,i,7] = 1 - s[i,4];
      tr[5,i,1] = tr[3,i,1];
      tr[5,i,2] = tr[3,i,2];
      tr[5,i,3] = 0;
      tr[5,i,4] = 0;
      tr[5,i,5] = tr[3,i,3];
      tr[5,i,6] = tr[3,i,4];
      tr[5,i,7] = tr[3,i,7];
      tr[6,i,1] = 0;
      tr[6,i,2] = tr[4,i,2];
      tr[6,i,3] = 0;
      tr[6,i,4] = 0;
      tr[6,i,5] = 0;
      tr[6,i,6] = tr[4,i,4];
      tr[6,i,7] = tr[4,i,7];
      for(j in 1:6){
        tr[7,i,j] = 0;
      }
      tr[7,i,7] = 1;
    }
  
  // capture probability matrix
  for(i in 1:22){
    p[1,LCRs[i],1] = p_lcr[i,1];
    p[2,LCRs[i],2] = p_lcr[i,2];
  }
  // No sampling in LCR during 1st, 4th, 7th, 10th recap periods (summers of 2009 - 2012)
  for(i in 1:4){
    p[1,LCRns[i],1] = 0;
    p[2,LCRns[i],2] = 0;
  }
  for(i in 1:24){
    p[3,CRs[i],3] = p_cr[i,1];
    p[4,CRs[i],4] = p_cr[i,2];
  }
  // No sampling in CR during 3rd and 6th recap periods (spring of 2010 & 2011)
  for(i in 1:2){
    p[3,CRns[i],3] = 0;
    p[4,CRns[i],4] = 0;
  }
  for(i in 1:26){
    p[1,i,2] = 0;
    p[1,i,3] = 0;
    p[1,i,4] = 0;
    p[1,i,5] = 1 - p[1,i,1]; // prob of fish in state 1 not being captured 
    p[2,i,1] = 0;
    p[2,i,3] = 0;
    p[2,i,4] = 0;
    p[2,i,5] = 1 - p[2,i,2]; // prob of fish in state 2 not being captured 
    p[3,i,1] = 0;
    p[3,i,2] = 0;
    p[3,i,4] = 0;
    p[3,i,5] = 1 - p[3,i,3] ; // prob of fish in state 3 not being captured 
    p[4,i,1] = 0;
    p[4,i,2] = 0;
    p[4,i,3] = 0;
    p[4,i,5] = 1 - p[4,i,4]; // prob of fish in state 4 not being captured 
    p[5,i,1] = 0;
    p[5,i,2] = 0;
    p[5,i,3] = 0;
    p[5,i,4] = 0;
    p[5,i,5] = 1; // prob of fish in state 5 not being captured (set to 100% since these fish are unobserveable) 
    p[6,i,1] = 0;
    p[6,i,2] = 0;
    p[6,i,3] = 0;
    p[6,i,4] = 0;
    p[6,i,5] = 1; // prob of fish in state 6 not being captured (set to 100% since these fish are unobserveable) 
    p[7,i,1] = 0;
    p[7,i,2] = 0;
    p[7,i,3] = 0;
    p[7,i,4] = 0;
    p[7,i,5] = 1; // prob of fish in state 7 not being captured (set to 100% since these fish are dead) 
  }
}

model{
  real temp[7]; 
  vector[7] pz[27]; 
  
  for(i in 1:26){
    z_ls[i] ~ normal(0, 1);
    z_lg[i] ~ normal(0, 1);
    z_lm[i] ~ normal(0, 1);
  }
  
  mu_ls ~ normal(0, 2);
  
  for(i in 1:3){
    mu_lm[i] ~ normal(0, 2);
    mu_lg[i] ~ normal(0, 2);
  }
  
  for(k in 1:NsumCH){  
    for(j in 1:6){
      pz[sumf[k],j] = (j == sumCH[k, sumf[k]]) * (1 - newtag[k] * 0.03); 
    }
    pz[sumf[k],7] = newtag[k] * 0.03;     
    for(t in (sumf[k] + 1):lc[k]){ 
      for(i in 1:7){ 
        for(j in 1:7){
          temp[j] = pz[t - 1,j] * tr[j,t - 1,i] * p[i,t - 1,sumCH[k,t]];
        }
        pz[t,i] = sum(temp); 
      } 
    } 
    target += (sumFR[k] * log(sum(pz[lc[k]]))); 
  }
}

//Code to estimate abundance by simulating from a negative binomial distribution
generated quantities{
  real ptrans;	
  real<lower = 0> scale_par;
  int U[4,9];
  int N[4,9];
  int N_tot[9];
  
  for(i in 1:4){
    for(j in 1:9){
      ptrans = p[i,fall_ind[j],i] * (1 - (i > 2) * (1 - tau));
      scale_par = ptrans / (1 - ptrans);
      U[i,j] = neg_binomial_rng(catch_mat[i,j], scale_par);
      N[i,j] = U[i,j] + catch_mat[i,j];
    }
  }
  
  for(j in 1:9){
    N_tot[j]= sum(N[,j]);
  }
}

    ",fill=TRUE)
sink() 

#-----------------------------------------------------------------------------#
smre.data <- list(NsumCH = dim(sumCH)[1], sumCH = sumCH, sumf = as.vector(fc),
                  season = season, sumFR = sumFR, LCRs = LCRs, LCRns = LCRns,
                  CRs = CRs, CRns = CRns, fall_ind = 1:9 * 3 - 1,
                  catch_mat = array(catch, dim = c(dim(catch)[1], dim(catch)[2])),
                  lc = as.vector(lc), newtag = as.vector(newtag.AO))

smre.inits <- function() {list(mu_ls = rep(0, 4))}

smre.params <- c("s", "g", "m", "mu_ls", "mu_lg", "mu_lm", "sd_ls", "sd_lm", 
                 "sd_lg", "tau", "p_lcr", "p_cr", "N") 

SMRE <- stan(".\\Stan_Marginalized_chub_RE.stan", 
             data = smre.data, init =smre.inits, pars = smre.params, 
             chains = 1, iter = 10, thin = 1, 
             seed = 1) 
#-----------------------------------------------------------------------------#
# End
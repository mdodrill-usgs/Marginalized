###############################################################################
#                                                                     Spring 19
#  Combination Occupancy and N-Mixture Model
#  Discrete JAGS version 
#
#  Notes:
#   *  Need to set directory for data
#   *  Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2jags)

# data.dir = paste0(getwd(), "/Data/")

#-----------------------------------------------------------------------------#
# import, format data for model fitting
raw.data <- read.csv(paste0(data.dir, "barred_owl_survey_data.csv"), header = TRUE) 

# extract all BO data
nYears <- 22
tBO <- as.matrix(raw.data[,grep("BO", colnames(raw.data))])
nSites <- dim(tBO)[1]
Y <- tBO[,241:246]
visited <- ifelse(is.na(Y) == TRUE, 0, 1)
Y <- ifelse(is.na(Y) == TRUE, 0, Y)
maxY <- matrix(0, nrow = nSites, ncol = 2)

for(t in 1:2){
  maxY[,t] <- apply(Y[,(3  *  (t - 1) + 1):( 3  *  t)], 1, max)
}

num_det <- matrix(0, nrow = nSites, ncol = 20)
num_trial <- matrix(0, nrow = nSites, ncol = 20)
num_nodet <- matrix(0, nrow = nSites, ncol = 20)

for(t in 1:20){
  num_trial[,t] <- rowSums(ifelse(is.na(tBO[,((t - 1)  *  12 + 1):(t  *  12)]) == TRUE, 0, 1))
  num_det[,t] <- rowSums(ifelse(is.na(tBO[,((t - 1)  *  12 + 1):(t  *  12)]) == TRUE, 0,
                                tBO[,((t - 1)  *  12 + 1):(t  *  12)]))
  num_nodet[,t] <- num_trial[,t] - num_det[,t]
}

maxdet <- max(num_det)
maxnodet <- max(num_nodet)

#-----------------------------------------------------------------------------#
# Pull out and standardize covariates
# PHAB = percent older riparian growth forest)
PHAB <- as.matrix(raw.data[,grep("PH", colnames(raw.data))]) 
X <-  (PHAB - mean(PHAB)) / sd(PHAB)
tOV <- as.matrix(raw.data[,grep("OV", colnames(raw.data))])
visited <- ifelse(tOV == 0, 0, 1)
OV <- ifelse(tOV == 0, 0, log(tOV))
Ninf = 15

C <- array(0, dim = c(nSites, 2, (Ninf + 1)))
for(i in 1:nSites){
  for(j in 1:(Ninf + 1)){
    C[i,1,j] <- choose((j-1), Y[i,1])  *  choose((j-1), Y[i,2])  *  choose((j-1), Y[i,3])
    C[i,2,j] <- choose((j-1), Y[i,4])  *  choose((j-1), Y[i,5])  *  choose((j-1), Y[i,6])
  }
}

iC <- ifelse(C == 0, 0, 1)
lC <- ifelse(C == 0, 0, log(C))

#-----------------------------------------------------------------------------#
sink("Combo_JD.jags")
cat("
model {
  for(i in 1:nSites){
    N[i,1] ~ dpois(lambda)
  }
  
  for(t in 1:(nYears-1)){
    N_mean[t] <- mean(N[,t]) - 1
    log(gamma[t]) <- g0 + g1  *  N_mean[t] + g2  *  N_mean[t]  *  N_mean[t]
    for(i in 1:nSites){
      logit(omega[i,t]) <- b0 + b1  *  X[i,t]
      S[i,t] ~ dbin(omega[i,t], N[i,t])
      G[i,t] ~ dpois(gamma[t])
      N[i,(t+1)] <- S[i,t] + G[i,t] 
    }
  }
  
  for(i in 1:nSites){
    for(t in 1:20){
      p_site[i,t] <- 1 - pow((1 - p_occ), N[i,t]) 
      num_det[i,t] ~ dbin(p_site[i,t], num_trial[i,t])
    }
    
    for(t in 1:2){
      for(k in 1:3){
        cloglog(tpc[i,(3  *  (t - 1) + k)]) <- a0 + OV[i,(3  *  (t - 1) + k)]
        pc[i,(3  *  (t - 1) + k)] <- tpc[i,(3  *  (t - 1) + k)]  *  visited[i,(3  *  (t - 1) + k)]
        Y[i,(3  *  (t - 1) + k)] ~ dbin(pc[i,(3  *  (t - 1) + k)], N[i,(t + 20)])
      }
    }
  }
  
  
  lambda ~ dunif(0, 1) # initial abundance
  p_occ ~ dunif(0, 1)  # detection
  b0 ~ dunif(-3, 3)    # intercept on survival
  b1 ~ dunif(-1, 1)    # slope on survival
  a0 ~ dunif(-3, 3)    # intercept on effort (p_count)
  g0 ~ dunif(-3, 3)    # intercept on gamma
  g1 ~ dunif(-1, 1)    # slope on mean(N)
  g2 ~ dunif(-1, 1)    # squared term on mean(N)
}
    ", fill = TRUE)
sink()   
#-----------------------------------------------------------------------------#
JD_data <- list(Y = Y, visited = visited, X = X, OV = OV, nYears = nYears,
                nSites = nSites, num_det = num_det, num_trial = num_trial)

params <- c("a0", "b0", "b1", "g0", "g1", "g2", "p_occ", "lambda", "N_mean")

JD_inits <- function() list(N = matrix(c(rep(20, nSites), rep(NA, nSites  *  (nYears - 1))),
                                       nrow = nSites, ncol = nYears),
                          S = matrix(10, nrow = nSites, ncol = nYears - 1),
                          G = matrix(10, nrow = nSites, ncol = nYears - 1))

JD_combo <- jags(JD_data, JD_inits, params, model.file = "combo_JD.txt", n.iter = 10)
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Combination Occupancy and N-Mixture Model
#  Marginalized JAGS version 
#
#  Notes:
#   *  Need to set directory for data
#
###############################################################################
library(R2jags)

# data.dir = paste0(getwd(), "/Data/")
#-----------------------------------------------------------------------------#
# import, format data for model fitting
raw.data <- read.csv(paste0(data.dir, "barred_owl_survey_data.csv"), header = TRUE) 

# extract all BO data
nYears <- 22
tBO <- as.matrix(raw.data[,grep("BO", colnames(raw.data))])
nSites <- dim(tBO)[1]
Y <- tBO[,241:246]
visited <- ifelse(is.na(Y) == TRUE, 0, 1)
Y <- ifelse(is.na(Y) == TRUE, 0, Y)

maxY <- matrix(0, nrow = nSites, ncol = 2)
for(t in 1:2){
  maxY[,t] <- apply(Y[,(3 * (t - 1) + 1):(3 * t)], 1, max)
}

num_det <- matrix(0, nrow = nSites, ncol = 20)
num_trial <- matrix(0, nrow = nSites, ncol = 20)
num_nodet <- matrix(0, nrow = nSites, ncol = 20)

for(t in 1:20){
  num_trial[,t] <- rowSums(ifelse(is.na(tBO[,((t - 1) * 12 + 1):(t * 12)]) == TRUE, 0, 1))
  num_det[,t] <- rowSums(ifelse(is.na(tBO[,((t - 1) * 12 + 1):(t * 12)]) == TRUE, 0,
                                tBO[,((t - 1) * 12 + 1):(t * 12)]))
  num_nodet[,t] <- num_trial[,t] - num_det[,t]
}

maxdet <- max(num_det)
maxnodet <- max(num_nodet)
#-----------------------------------------------------------------------------#
# Pull out and standardize covariates
# PHAB = percent older riparian growth forest)
PHAB <- as.matrix(raw.data[,grep("PH", colnames(raw.data))]) 
X <- (PHAB - mean(PHAB)) / sd(PHAB)
tOV <- as.matrix(raw.data[,grep("OV", colnames(raw.data))])
visited <- ifelse(tOV == 0, 0, 1)
OV <- ifelse(tOV == 0, 0, log(tOV))
Ninf = 15

C <- array(0, dim = c(nSites, 2, (Ninf + 1)))
for(i in 1:nSites){
  for(j in 1:(Ninf + 1)){
    C[i,1,j] <- choose((j - 1), Y[i,1]) * choose((j - 1), Y[i,2]) * choose((j - 1), Y[i,3])
    C[i,2,j] <- choose((j - 1), Y[i,4]) * choose((j - 1), Y[i,5]) * choose((j - 1), Y[i,6])
  }
}

iC <- ifelse(C == 0, 0, 1)
lC <- ifelse(C == 0, 0, log(C))

#-----------------------------------------------------------------------------#
sink("Combo_JM.jags")
cat("
model {
  for(j in 1:(Ninf + 1)){
    po[j,2] = 1 - exp((j - 1) * log(1 - p_occ))
    po[j,1] = 1 - po[j,2]
    pN[j] = dpois((j - 1), lambda)
    Ns[j] = (j - 1)
  }
  
  for(i in 1:nSites){
    for(t in 1:20){
      for(j in 1:(Ninf + 1)){
        cp[i,t,j] <- (po[j,2] ^ num_det[i,t]) * (po[j,1] ^ num_nodet[i,t])
      }
    }
    
    for(t in 1:2){ 
      for(k in 1:3){
        cloglog(tpc[i,(3 * (t - 1) + k)]) <- a0 + OV[i,(3 * (t - 1) + k)]
        pc[i,(3 * (t - 1) + k)] <- tpc[i,(3 * (t - 1) + k)] * visited[i,(3 * (t - 1) + k)]
      }
      
      for(j in 1:(Ninf + 1)){
        cp[i,(t + 20),j] <- dbin(Y[i,(3 * t-2)], pc[i,(3 * t - 2)], (j - 1)) * 
          dbin(Y[i,(3 * t - 1)], pc[i,(3 * t - 1)], (j - 1)) *
          dbin(Y[i,(3 * t)], pc[i,(3 * t)], (j - 1))
      }
    }
  }
  
  for(i in 1:nSites){
    for(j in 1:(Ninf + 1)){
      pz[i,1,j] = pN[j] * cp[i,1,j]
    }
    N_mn[i,1] <- inprod(Ns,pz[i,1,]) / sum(pz[i,1,])
  }
  
  for(t in 1:(nYears - 1)){
    N_mean[t] <- mean(N_mn[1:nSites,t]) - 1
    log(gamma[t]) <- g0 + g1 * N_mean[t] + g2 * N_mean[t] * N_mean[t]
    for(j in 1:(Ninf + 1)){
      for(k in 1:(Ninf + 1)){
        t_gtr[t,j,k] = dpois((k - j), gamma[t])
      }
      for(k in 1:(Ninf + 1)){
        gtr[t,j,k] = t_gtr[t,j,k] / sum(t_gtr[t,j,])
      }
    }
    for(i in 1:nSites){
      logit(omega[i,t]) = b0 + b1 * X[i,t]
      for(j in 1:(Ninf + 1)){
        for(k in 1:(Ninf + 1)){
          str[i,t,j,k] = pz[i,t,j] * dbin((k - 1), omega[i,t], (j - 1))
        }
      }
      for(j in 1:(Ninf + 1)){
        for(k in 1:(Ninf + 1)){
          tr[i,t,j,k] <- sum(str[i,t,j:(Ninf + 1),j]) * gtr[t,j,k]
        }
      }
      for(j in 1:(Ninf + 1)){
        pz[i,(t + 1),j] <- sum(tr[i,t,1:j,j]) * cp[i,t,j]
      }
      N_mn[i,(t + 1)] <- inprod(Ns, pz[i,(t + 1),]) / sum(pz[i,(t + 1),])
    }
  }
  
  for(i in 1:nSites){
    lik[i] <- sum(pz[i,22,])
    ones[i] ~ dbern(lik[i])
  }		
  
  lambda ~ dunif(0, 1) # initial abundance
  p_occ ~ dunif(0, 1)  # detection
  b0 ~ dunif(-3, 3)    # intercept on survival
  b1 ~ dunif(-1, 1)    # slope on survival
  a0 ~ dunif(-3, 3)    # intercept on effort (p_count)
  g0 ~ dunif(-3, 3)    # intercept on gamma
  g1 ~ dunif(-1, 1)    # slope on mean(N)
  g2 ~ dunif(-1, 1)    # squared term on mean(N)
}
    ", fill = TRUE)
sink()   
#-----------------------------------------------------------------------------#
JM_data <- list(Y = Y, visited = visited, X = X, OV = OV, nYears = nYears,
                nSites = nSites, num_det = num_det, num_nodet = num_nodet, 
                Ninf = Ninf, ones = rep(1, nSites), maxdet = maxdet,
                maxnodet = maxnodet)

params <- c("a0", "b0", "b1", "g0", "g1", "g2", "p_occ", "lambda", "N_mean")


JM_combo  <-  jags.parallel(JM_data, inits = NULL, params, 'Combo_JM.jags',
                         n.chains = 3, n.iter = 10)
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                        Jan 19
#  Combination Occupancy and N-Mixture Model
#  Marginalized Stan version 
#
#  Notes:
#  * 
#
#  To do: 
#  * format .Stan code (identation!)
#  * clean up and format this code
#
###############################################################################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())  

# data.dir = paste0(getwd(), "/Data/")
#-----------------------------------------------------------------------------#
# import, format data for model fitting
raw.data <- read.csv(paste0(data.dir, "barred_owl_survey_data.csv"), header = TRUE) 

# extract all BO data
nYears <- 22
tBO <- as.matrix(raw.data[,grep("BO", colnames(raw.data))])
nSites <- dim(tBO)[1]
Y <- tBO[,241:246]
visited <- ifelse(is.na(Y) == TRUE, 0, 1)
Y <- ifelse(is.na(Y) == TRUE, 0, Y)

maxY <- matrix(0, nrow = nSites, ncol = 2)
for(t in 1:2){
  maxY[,t] <- apply(Y[,(3 * (t - 1) + 1):(3 * t)], 1, max)
}

num_det <- matrix(0, nrow = nSites, ncol = 20)
num_trial <- matrix(0, nrow = nSites, ncol = 20)
num_nodet <- matrix(0, nrow = nSites, ncol = 20)

for(t in 1:20){
  num_trial[,t] <- rowSums(ifelse(is.na(tBO[,((t - 1) * 12 + 1):(t * 12)]) == TRUE, 0, 1))
  num_det[,t] <- rowSums(ifelse(is.na(tBO[,((t - 1) * 12 + 1):(t * 12)]) == TRUE, 0, 
                                tBO[,((t - 1) * 12 + 1):(t * 12)]))
  num_nodet[,t] <- num_trial[,t] - num_det[,t]
}

maxdet <- max(num_det)
maxnodet <- max(num_nodet)

#-----------------------------------------------------------------------------#
# Pull out and standardize covariates
# PHAB = percent older riparian growth forest)
PHAB <- as.matrix(raw.data[,grep("PH", colnames(raw.data))]) 
X <- (PHAB - mean(PHAB)) / sd(PHAB)
tOV <- as.matrix(raw.data[,grep("OV", colnames(raw.data))])
visited <- ifelse(tOV == 0, 0, 1)
OV <- ifelse(tOV == 0, 0, log(tOV))
Ninf = 15

C <- array(0, dim = c(nSites, 2, (Ninf + 1)))
for(i in 1:nSites){
  for(j in 1:(Ninf + 1)){
    C[i,1,j] <- choose((j - 1), Y[i,1]) * choose((j - 1), Y[i,2]) * choose((j - 1), Y[i,3])
    C[i,2,j] <- choose((j - 1), Y[i,4]) * choose((j - 1), Y[i,5]) * choose((j - 1), Y[i,6])
  }
}

iC <- ifelse(C == 0, 0, 1)
lC <- ifelse(C == 0, 0, log(C))

#-----------------------------------------------------------------------------#
sink("Combo_SM.stan")
cat("

data{
  int<lower = 1> nYears;
  int<lower = 1> nSites;
  int<lower = 1> Ninf;
  int<lower = 0> Y [nSites,6];
  int<lower = 0> visited [nSites,6];
  matrix [nSites,6] OV;
  matrix [nSites,(nYears - 1)] X;
  int<lower = 0> num_det [nSites,(nYears - 2)];
  int<lower = 0> num_nodet [nSites,(nYears - 2)];
  int<lower = 0> maxdet;
  int<lower = 0> maxnodet;
  int<lower = 0,upper = 1> iC [nSites,2,(Ninf + 1)];
  real<lower = 0> lC [nSites,2,(Ninf + 1)];
  row_vector [(Ninf + 1)] Ns; 
}

parameters{
  real<lower = 0,upper = 1> lambda;
  real<lower = 0,upper = 1> p_occ;
  real<lower = -3,upper = 3> b0;
  real<lower = -1,upper = 1> b1;
  real<lower = -3,upper = 3> a0;
  real<lower = -3,upper = 3> g0;
  real<lower = -1,upper = 1> g1;
  real<lower = -1,upper = 1> g2;
}

transformed parameters {
  simplex [2] po [(Ninf + 1)];
  row_vector<lower = 0, upper = 1> [(Ninf + 1)] sumdet [(maxdet + 1),(maxnodet + 1)];
  row_vector<lower = 0, upper = 1> [(Ninf + 1)] pN;
  row_vector<lower = 0, upper = 1> [(Ninf + 1)] cp [nSites,nYears];
  matrix<lower = 0, upper = 1> [nSites,6] pc;
  row_vector<lower = 0, upper = 1> [(Ninf + 1)] pz [nSites,nYears];
  real<lower = 0> N_mn [nSites,nYears];
  real N_mean [(nYears - 1)];
  real<lower = 0> gamma [(nYears - 1)];
  matrix<lower = 0, upper = 1> [(Ninf + 1),(Ninf + 1)] gtr;
  matrix<lower = 0, upper = 1> [(Ninf + 1),(Ninf + 1)] str;
  matrix<lower = 0, upper = 1> [nSites,(nYears - 1)] omega ;
  
  pc = (1 - exp( - 1 * exp(a0 + OV)));
  omega = inv_logit(b0 + b1 * X);
  
  for(j in 1:(Ninf + 1)){
    po[j,2] = 1 - exp((j - 1) * log(1 - p_occ));
    po[j,1] = 1 - po[j,2];
    pN[j] = exp(poisson_lpmf((j - 1) | lambda));
  }
  
  for(d in 1:(maxdet + 1)){
    for(n in 1:(maxnodet + 1)){
      for(j in 1:(Ninf + 1)){
        sumdet[d,n,j] = (po[j,2] ^ (d - 1)) * (po[j,1] ^ (n - 1));
      }
    }
  }
  
  for(i in 1:nSites){
    for(t in 1:20){
      cp[i,t,] = sumdet[(num_det[i,t] + 1),(num_nodet[i,t] + 1),];
    }
    pz[i,1,] = pN .* cp[i,1,];
    N_mn[i,1] = sum(Ns .* pz[i,1,]) / sum(pz[i,1,]);              
    for(t in 1:2){
      for(j in 1:(Ninf + 1)){
        cp[i,(t + 20),j] = iC[i,t,j] * 
          (exp(lC[i,t,j] + visited[i,(3 * t - 2)] *
                 (Y[i,(3 * t - 2)] * log(pc[i,(3 * t - 2)]) +
                    (j - 1 - Y[i,(3 * t - 2)]) * log(1 - pc[i,(3 * t - 2)])) +
                 visited[i,(3 * t - 1)] * (Y[i,(3 * t - 1)] * log(pc[i,(3 * t - 1)]) +
                                             (j - 1 - Y[i,(3 * t - 1)]) * log(1 - pc[i,(3 * t - 1)])) +
                 visited[i,(3 * t)] * (Y[i,(3 * t)] * log(pc[i,(3 * t)]) +
                                         (j - 1 - Y[i,(3 * t)]) * log(1 - pc[i,(3 * t)]))));
      }
    }
  }
  
  for(j in 2:(Ninf + 1)){
    for(k in 1:(j - 1)){
      gtr[j,k] = 0;
    }
  }
  
  for(j in 1:(Ninf)){
    for(k in (j + 1):(Ninf + 1)){
      str[j,k] = 0;
    }
  }
  
  str[1,1] = 1;
  for(t in 1:(nYears - 1)){
    N_mean[t] = mean(N_mn[,t]) - 1;
    gamma[t] = exp(g0 + g1 * N_mean[t] + g2 * N_mean[t] * N_mean[t]);
    for(k in 1:(Ninf + 1)){
      gtr[1,k] = exp(poisson_lpmf((k - 1) | gamma[t]) - poisson_lcdf(Ninf | gamma[t]));
    }
    for(j in 2:(Ninf + 1)){
      for(k in j:(Ninf + 1)){
        gtr[j,k] = exp(poisson_lpmf((k - j) | gamma[t]) - poisson_lcdf((Ninf + 1 - j) | gamma[t]));
      }
    }
    for(i in 1:nSites){
      for(j in 2:(Ninf + 1)){
        for(k in 1:j){
          str[j,k] = exp(binomial_lpmf((k - 1) | (j - 1), omega[i,t]));
        }
      }
      pz[i,(t + 1),] = ((pz[i,t,] * str) * gtr) .* cp[i,t,];
      N_mn[i,(t + 1)] = sum(Ns .* pz[i,(t + 1),]) / sum(pz[i,(t + 1),]);
    }
  }
}

model{
  for(i in 1:nSites){
    target += log(sum(pz[i,22,]));
  }
}

    ", fill = TRUE)
sink()
#-----------------------------------------------------------------------------#
SM_data <- list(Y = Y, visited = visited, X = X, OV = OV, nYears = nYears,
              nSites = nSites, Ns = c(0:Ninf), num_det = num_det,
              num_nodet = num_nodet, Ninf = Ninf, maxdet = maxdet,
              maxnodet = maxnodet, iC = iC, lC = lC)

params <- c("a0", "b0", "b1", "g0", "g1", "g2", "p_occ","lambda", "N_mean")

SM_Cocc <- stan("Combo_SM.stan",
                data = SM_data,
                pars = params,
                chains = 1, iter = 10) 
#-----------------------------------------------------------------------------#
# End
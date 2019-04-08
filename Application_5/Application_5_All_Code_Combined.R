
# Need to add functions...?

###############################################################################
#                                                                     Spring 19
#  Dynamic Community Occupancy Models - Sky Islands  
#  Discrete JAGS version 
#
#  Notes:
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2jags)

# data.dir = paste0(getwd(), "/Data/")
dat = read.csv(paste0(data.dir, "ntmb_CM_1991_till_1995.csv"))

# dimensions: species, points, years
nspp = 149
nsites = 92
nyears = 5  # 1991, 1992, 1993, 1994, 1995 
years = seq(1991, 1995, 1)
nsess = 3

# fill the data array
Y = array(NA, dim = c(nspp, nsites, nyears))
for(i in 1:5){
  Y[,,i] = as.matrix(dat[which(dat[,2] == years[i]),3:94])
}

# habitat covariate (dimensions 92 pts) 
hab = c(5L, 5L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 5L, 7L, 7L, 3L, 5L, 5L, 
        5L, 3L, 3L, 5L, 5L, 5L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
        4L, 4L, 3L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 3L, 5L, 3L, 3L, 3L, 3L, 
        3L, 5L, 5L, 5L, 2L, 2L, 2L, 4L, 2L, 3L, 5L, 4L, 5L, 4L, 4L, 3L, 
        4L, 5L, 4L, 5L, 4L, 3L, 3L, 3L, 4L, 4L, 5L, 4L, 4L, 1L, 4L, 3L, 
        5L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 6L, 6L, 7L, 7L, 7L)

#-----------------------------------------------------------------------------#
sink("Cocc_JD.jags")
cat("
model{
  for(i in 1:(nspp)){
    a1[i] ~ dunif(0, 1)
    alpha1[i] ~ dnorm(alpha1_mu, tau_alpha1)
    
    for(h in 1:7){
      alpha2[i,h] ~ dnorm(0, tau_alpha2)
    }
    
    alpha0[i] ~ dnorm(beta, tau_u)
    mu_eta[i] <- alpha + (rho * sigma_v / sigma_u) * (alpha0[i] - beta)
    beta0[i] ~ dnorm(mu_eta[i], tau_eta)
    logit(p[i]) <- beta0[i] # detection  
    
    for(k in 1:nsites){
      z0[i,k] ~ dbern(a1[i])
      logit(psi[i,k,1]) <- alpha0[i] + alpha1[i] * z0[i,k] + alpha2[i,hab[k]]
      Z[i,k,1] ~ dbern(psi[i,k,1]) # occupancy
      mu_p[i,k,1] <- p[i] * Z[i,k,1] # detection
      Y[i,k,1] ~ dbin(mu_p[i,k,1], nsess)
      
      for(t in 1:(nyears - 1)){
        logit(psi[i,k,(t + 1)]) <- alpha0[i] + alpha1[i] * Z[i,k,t] + alpha2[i,hab[k]]
        Z[i,k,(t + 1)] ~ dbern(psi[i,k,(t + 1)]) # occupancy
        mu_p[i,k,(t + 1)] <- p[i] * Z[i,k,(t + 1)] # detection
        Y[i,k,(t + 1)] ~ dbin(mu_p[i,k,(t + 1)], nsess)
      }
    }
  }
  
  for(k in 1:nsites){  
    for(t in 1:nyears){
      SR[k,t] <- sum(Z[,k,t])
    }
  }
  
  psi_mean ~ dunif(0, 1)
  beta <- log(psi_mean) - log(1 - psi_mean)
  p_mean ~ dunif(0, 1)
  alpha <- log(p_mean) - log(1 - p_mean)
  alpha1_mean ~ dunif(0,1)
  alpha1_mu <- log(alpha1_mean) - log(1 - alpha1_mean)
  sigma_alpha1 ~ dunif(0, 5)
  sigma_alpha2 ~ dunif(0, 5)
  sigma_u ~ dunif(0, 5)
  sigma_v ~ dunif(0, 5)
  tau_alpha1 <- pow(sigma_alpha1, -2)
  tau_alpha2 <- pow(sigma_alpha2, -2)
  tau_u <- pow(sigma_u, -2)
  tau_v <- pow(sigma_v, -2)
  rho ~ dunif(-1, 1)
  tau_eta <- tau_v / (1-pow(rho,2))
}
    ", fill = TRUE)
sink()   
#-----------------------------------------------------------------------------#
# function to make inits for JAGS 
JD_inits = function() {
  nspp = 149
  nsites = 92
  nyears = 5  # 1991, 1992, 1993, 1994, 1995 
  years = seq(1991, 1995, 1)
  nsess = 3
  
  dat = read.csv(paste0(data.dir, "ntmb_CM_1991_till_1995.csv"))
  
  # fill the data array
  Y = array(NA, dim = c(nspp, nsites, nyears))
  for(i in 1:5){
    Y[,,i] = as.matrix(dat[which(dat[,2] == years[i]),3:94])
  }
  
  a1 = runif(nspp, 0.25, 1) 
  Ztemp = array(rbinom(nspp * 5 * nsites, size = 1, prob = 0.5), dim = c(nspp, nsites, 5))
  z0temp = array(rbinom(nspp * nsites, size = 1, prob = 0.5), dim = c(nspp, nsites))
  
  for(t in 1:5){
    for(i in 1:nspp){
      for(k in 1:nsites){
        if(Ztemp[i,k,t] == 0){
          if(!is.na(Y[i,k,t])){
            if(Y[i,k,t] > 0){
              Ztemp[i,k,t] = 1              
            }
          }
        }
      }
    }
  }
  list(a1 = a1, z0 = z0temp, Z = Ztemp)
}
#-----------------------------------------------------------------------------#
JD_data = list('nspp' = nspp, 'nsites' = nsites, 'nsess' = nsess,
               'nyears' = nyears, 'hab' = hab, 'Y' = Y)

params <- c('psi_mean', 'p_mean', 'rho', 'sigma_v', 'sigma_u', 'sigma_alpha1',
            'sigma_alpha2', 'alpha1_mean', 'SR', 'Z')

JD_Cocc <- jags.parallel(JD_data, inits = JD_inits,
                         params, "Cocc_JD.jags",
                         n.chains = 3, n.iter = 10,
                         export_obj_names = c("data.dir"))
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Dynamic Community Occupancy Models - Sky Islands  
#  Marginalized JAGS version 
#
#  Notes:
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2jags)

source(paste0(getwd(),"/Functions.R"), chdir = F)

# data.dir = paste0(getwd(), "/Data/")
dat = read.csv(paste0(data.dir, "ntmb_CM_1991_till_1995.csv"))

# dimensions: species, points, years
nspp = 149
nsites = 92
nyears = 5  # 1991, 1992, 1993, 1994, 1995 
years = seq(1991, 1995, 1)
nsess = 3

# fill the data array
Y = array(NA, dim = c(nspp, nsites, nyears))
for(i in 1:5){
  Y[,,i] = as.matrix(dat[which(dat[,2] == years[i]),3:94])
}

visit <- ifelse(is.na(Y[1,,1]) == TRUE, 2, 1)
sY <- ifelse(is.na(Y) == TRUE, 1, Y + 1)

# habitat covariate (dimensions 92 pts) 
hab = c(5L, 5L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 5L, 7L, 7L, 3L, 5L, 5L, 
        5L, 3L, 3L, 5L, 5L, 5L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
        4L, 4L, 3L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 3L, 5L, 3L, 3L, 3L, 3L, 
        3L, 5L, 5L, 5L, 2L, 2L, 2L, 4L, 2L, 3L, 5L, 4L, 5L, 4L, 4L, 3L, 
        4L, 5L, 4L, 5L, 4L, 3L, 3L, 3L, 4L, 4L, 5L, 4L, 4L, 1L, 4L, 3L, 
        5L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 6L, 6L, 7L, 7L, 7L)

#-----------------------------------------------------------------------------#
sink("Cocc_JM.jags")
cat("
model{
  for(i in 1:(nspp)){
    a1[i] ~ dunif(0, 1)
    z0[i,1] <- a1[i]
    z0[i,2] <- 1 - a1[i]
    alpha1[i] ~ dnorm(alpha1_mu, tau_alpha1)
    for(h in 1:7){
      alpha2[i,h] ~ dnorm(0, tau_alpha2)
    }
    alpha0[i] ~ dnorm(beta, tau_u)
    mu_eta[i] <- alpha + (rho * sigma_v / sigma_u) * (alpha0[i] - beta)
    beta0[i] ~ dnorm(mu_eta[i], tau_eta)
    logit(p[i]) <- beta0[i] # detection
    po[i,1,1,2] = 1;
    po[i,2,1,2] = 1;
    for (j in 1:nsess){
      po[i,1,(j + 1),1] <- 0
      po[i,2,(j + 1),1] <- (p[i] ^ j) * (1 - p[i]) ^ (nsess - j)
    }
    po[i,1,1,1] <- 1
    po[i,2,1,1] <- (1 - p[i]) ^ (nsess)
    for(k in 1:nhab){
      logit(tr[i,k,1,2]) <- alpha0[i] + alpha2[i,k]
      logit(tr[i,k,2,2]) <- alpha0[i] + alpha1[i] + alpha2[i,k]
      tr[i,k,1,1] <- 1 - tr[i,k,1,2]
      tr[i,k,2,1] <- 1 - tr[i,k,2,2]
    }
  }
  
  for(i in 1:NuDH){
    pz[i,1,1] = (z0[spp[i],1] * tr[spp[i],hab2[i],1,1] + z0[spp[i],2] * tr[spp[i],hab2[i],2,1]) *
      po[spp[i],1,sY2[i,1],visit2[i]];
    pz[i,1,2] = (z0[spp[i],1] * tr[spp[i],hab2[i],1,2] + z0[spp[i],2] * tr[spp[i],hab2[i],2,2]) * 
      po[spp[i],2,sY2[i,1],visit2[i]];
    Z[i,1] = pz[i,1,2] / (pz[i,1,1] + pz[i,1,2]);
    for(t in 1:(nyears - 1)){
      pz[i,(t + 1),1] = (pz[i,t,1] * tr[spp[i],hab2[i],1,1] + pz[i,t,2] * tr[spp[i],hab2[i],2,1]) * 
        po[spp[i],1,sY2[i,(t + 1)],1];
      pz[i,(t + 1),2] = (pz[i,t,1] * tr[spp[i],hab2[i],1,2] + pz[i,t,2] * tr[spp[i],hab2[i],2,2]) *
        po[spp[i],2,sY2[i,(t + 1)],1];
      Z[i,(t + 1)] = pz[i,(t + 1),2] / (pz[i,(t + 1),1] + pz[i,(t + 1),2]);
    }
    lik[i] <- sum(pz[i, nyears,])
    fr[i] ~ dbin(lik[i], FR[i])
  }
  
  for (k in 1:nsites){  
    for (t in 1:nyears){
      for (i in 1:nspp){
        Z2[i,k,t] ~ dbern(Z[lookup[i,k],t])
      }
      SR[k,t] <- sum(Z2[1:nspp,k,t])
    }
  }
  
  psi_mean ~ dunif(0, 1)
  beta <- log(psi_mean) - log(1 - psi_mean)
  p_mean ~ dunif(0, 1)
  alpha <- log(p_mean) - log(1 - p_mean)
  alpha1_mean ~ dunif(0, 1)
  alpha1_mu <- log(alpha1_mean) - log(1 - alpha1_mean)
  sigma_alpha1 ~ dunif(0, 5)
  sigma_alpha2 ~ dunif(0, 5)
  sigma_u ~ dunif(0, 5)
  sigma_v ~ dunif(0, 5)
  tau_alpha1 <- pow(sigma_alpha1, -2)
  tau_alpha2 <- pow(sigma_alpha2, -2)
  tau_u <- pow(sigma_u, -2)
  tau_v <- pow(sigma_v, -2)
  rho ~ dunif(-1, 1)
  tau_eta <- tau_v / (1 - pow(rho, 2)) 
}
    ", fill = TRUE)
sink()   
#-----------------------------------------------------------------------------#
JM_inits = function() {
  nspp = 149
  nsites = 92
  nyears = 5 
  Z2 = array(rbinom(nspp * nsites * nyears, 1, 0.5), dim = c(nspp, nsites, nyears))
  list(Z2 = Z2)
}

pasted <- function(x){
  paste0(x[1], x[2], x[3], x[4], x[5], x[6], x[7], collapse = "")
}

unpaste <- function(x){
  as.numeric(c(substr(x, 1, 1), substr(x, 2, 2), substr(x, 3, 3),
               substr(x, 4, 4), substr(x, 5, 5), substr(x, 6, 6), substr(x, 7, 7)))
}

#-----------------------------------------------------------------------------#
FR <- numeric()
spp <- numeric()
sY2 <- rep(NA, 5)
visit2 <- numeric()
hab2 <- numeric()
lookup <- matrix(NA, nspp, nsites)
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
JM_data =list('nspp' = nspp, 'nsites' = nsites,'nsess' = nsess, 'nyears' = nyears,
              'hab2' = hab2, 'fr' = FR, 'FR' = FR, 'nhab' = 7, 'NuDH' = length(FR),
              'spp' = spp, 'lookup' = lookup, 'sY2' = sY2, 'visit2' = visit2)

params <- c('psi_mean', 'p_mean', 'rho', 'sigma_v', 'sigma_u', 'sigma_alpha1',
            'sigma_alpha2', 'alpha1_mean','SR', 'Z2')

JM_Cocc <- jags.parallel(JM_data, inits = JM_inits, params, 'Cocc_JM.jags',
                         n.chains = 3, n.iter = 10)
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Dynamic Community Occupancy Models - Sky Islands  
#  Marginalized Stan version 
#
#  Notes:
#  * Need to set directory for data
#  * Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())  

# data.dir = paste0(getwd(), "/Data/")
dat = read.csv(paste0(data.dir, "ntmb_CM_1991_till_1995.csv"))

# dimensions: species, points, years
nspp = 149
nsites = 92
nyears = 5  # 1991, 1992, 1993, 1994, 1995 
years = seq(1991, 1995, 1)
nsess = 3

# fill the data array
Y = array(NA, dim = c(nspp, nsites, nyears))
for(i in 1:5){
  Y[,,i] = as.matrix(dat[which(dat[,2] == years[i]),3:94])
}

visit <- ifelse(is.na(Y[1,,1]) == TRUE, 2, 1)
sY <- ifelse(is.na(Y) == TRUE, 1, Y + 1)

# habitat covariate (dimensions 92 pts) 
hab = c(5L, 5L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 5L, 7L, 7L, 3L, 5L, 5L, 
        5L, 3L, 3L, 5L, 5L, 5L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
        4L, 4L, 3L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 3L, 5L, 3L, 3L, 3L, 3L, 
        3L, 5L, 5L, 5L, 2L, 2L, 2L, 4L, 2L, 3L, 5L, 4L, 5L, 4L, 4L, 3L, 
        4L, 5L, 4L, 5L, 4L, 3L, 3L, 3L, 4L, 4L, 5L, 4L, 4L, 1L, 4L, 3L, 
        5L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 6L, 6L, 7L, 7L, 7L)
#-----------------------------------------------------------------------------#
sink("Cocc_SM.stan")
cat("
data{
  int<lower = 1> nspp;
  int<lower = 1> nsess;
  int<lower = 1> nsites;
  int<lower = 1> nyears;
  int<lower = 1> NuDH; 
  int<lower = 1> nhab; 
  int<lower = 1> sY2 [NuDH,nyears]; 
  int<lower = 1> visit2 [NuDH]; 
  int<lower = 1> hab2 [NuDH]; 
  int<lower = 1> spp [NuDH];
  int<lower = 1> FR [NuDH];
  int<lower = 1> lookup[nspp,nsites];
}

parameters{
  real<lower = 0, upper = 1> psi_mean;
  real<lower = 0, upper = 1> p_mean;
  real<lower = 0, upper = 1> alpha1_mean;
  real<lower = 0, upper = 5> sigma_alpha1;
  real<lower = 0, upper = 5> sigma_alpha2;
  real<lower = 0, upper = 5> sigma_u;
  real<lower = 0, upper = 5> sigma_v;
  real<lower = -1, upper = 1> rho;
  vector<lower = 0, upper = 1> [nspp] a1;
  vector [nspp] alpha_dev1;
  vector [nspp] alpha_dev0;
  vector [nspp] alpha_dev2 [7];
  real beta_dev0 [nspp];
}

transformed parameters{
  real beta;
  real alpha;
  real alpha1_mu;
  real sigma_eta;
  real Z [NuDH,nyears];
  real p [nspp];
  simplex [2] z0 [nspp];
  real po [nspp,2,(nsess + 1),2];
  real tr [nspp,nhab,2,2];
  real pz [NuDH,nyears,2];
  
  beta = logit(psi_mean);
  alpha = logit(p_mean);
  alpha1_mu = logit(alpha1_mean);
  sigma_eta = sigma_v * ((1 - rho ^ 2) ^ 0.5);
  
  for(i in 1:(nspp)){
    z0[i,2] = a1[i];
    z0[i,1] = 1 - a1[i];
    po[i,1,1,2] = 1;
    po[i,2,1,2] = 1;
    po[i,1,1,1] = 1;
    p[i] = inv_logit(alpha + rho * sigma_v * alpha_dev0[i] + sigma_eta * beta_dev0[i]);
    for(j in 1:nsess){
      po[i,1,(j + 1),1] = 0;
      po[i,1,(j + 1),2] = 0;
      po[i,2,(j + 1),2] = 0;			
      po[i,2,(j + 1),1] = (p[i] ^ j) * (1 - p[i]) ^ (nsess - j);
    }
    po[i,2,1,1] = (1 - p[i]) ^ (nsess);
    for(k in 1:nhab){
      tr[i,k,2,2] = inv_logit(beta + sigma_u * alpha_dev0[i] + alpha1_mu + sigma_alpha1 * alpha_dev1[i] + sigma_alpha2 * alpha_dev2[k,i]);
      tr[i,k,1,2] = inv_logit(beta + sigma_u * alpha_dev0[i] + sigma_alpha2 * alpha_dev2[k,i]);
      tr[i,k,1,1] = 1 - tr[i,k,1,2];
      tr[i,k,2,1] = 1 - tr[i,k,2,2];
    }
  }
  
  for(i in 1:NuDH){
    pz[i,1,1] = (z0[spp[i],1] * tr[spp[i],hab2[i],1,1] + z0[spp[i],2] * tr[spp[i],hab2[i],2,1]) *
      po[spp[i],1,sY2[i,1],visit2[i]];
    pz[i,1,2] = (z0[spp[i],1] * tr[spp[i],hab2[i],1,2] + z0[spp[i],2] * tr[spp[i],hab2[i],2,2]) *
      po[spp[i],2,sY2[i,1],visit2[i]];
    Z[i,1] = pz[i,1,2] / (pz[i,1,1] + pz[i,1,2]);
    for(t in 1:(nyears - 1)){
      pz[i,(t + 1),1] = (pz[i,t,1] * tr[spp[i],hab2[i],1,1] + pz[i,t,2] * tr[spp[i],hab2[i],2,1]) * 
        po[spp[i],1,sY2[i,(t + 1)],1];
      pz[i,(t + 1),2] = (pz[i,t,1] * tr[spp[i],hab2[i],1,2] + pz[i,t,2] * tr[spp[i],hab2[i],2,2]) * 
        po[spp[i],2,sY2[i,(t + 1)],1];
      Z[i,(t + 1)] = pz[i,(t + 1),2] / (pz[i,(t + 1),1] + pz[i,(t + 1),2]);
    }
  }
}

model{	
  alpha_dev1 ~ normal(0, 1);
  alpha_dev0 ~ normal(0, 1);
  
  for(h in 1:7){
    alpha_dev2[h] ~ normal(0, 1);
  }
  
  beta_dev0 ~ normal(0, 1);
  
  for(i in 1:NuDH){
    target += FR[i] * log(sum(pz[i,nyears,]));
  }
}

generated quantities{
  int SR [nsites,nyears];
  int temp [nspp];  
  
  for(k in 1:nsites){  
    for(t in 1:nyears){
      for(i in 1:nspp){
        temp[i] = bernoulli_rng(Z[lookup[i,k],t]);
      }
      SR[k,t] = sum(temp);
    }
  }
}

    ", fill = TRUE)
sink()
#-----------------------------------------------------------------------------#
pasted <- function(x){
  paste0(x[1], x[2], x[3], x[4], x[5], x[6], x[7], collapse = "")
}

unpaste <- function(x){
  as.numeric(c(substr(x, 1, 1), substr(x, 2, 2), substr(x, 3, 3),
               substr(x, 4, 4), substr(x, 5, 5), substr(x, 6, 6), substr(x, 7, 7)))
}

FR <- numeric()
spp <- numeric()
sY2 <- rep(NA, 5)
visit2 <- numeric()
hab2 <- numeric()
lookup <- matrix(NA, nspp, nsites)
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
SM_data <- list(nspp = nspp, nsess = nsess, nsites = nsites, nyears = nyears, 
              hab2 = hab2, sY2 = sY2, visit2 = visit2, nhab = 7, NuDH = length(FR), 
              spp = spp, FR = FR, lookup = lookup)

params <- c('psi_mean', 'p_mean', 'rho', 'sigma_v', 'sigma_u', 'sigma_alpha1',
            'sigma_alpha2', 'alpha1_mean', 'SR')

SM_Cocc <- stan("Cocc_SM.stan",
                data = SM_data,
                pars = params,
                chains = 1, iter = 1000) 
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Simulation Code for Figure 4
#
#  Notes:
#  * 
#
###############################################################################
library(R2jags)

# jags code - discrete approach
sink("JD.jags")
cat("
    model{
    Z[1] ~ dbern(psi0)
    Z[2] ~ dbern(Z[1] * (1 - eps) + (1 - Z[1]) * gam)
    Z[3] ~ dbern(Z[2] * (1 - eps) + (1 - Z[2]) * gam)
    for(i in 1:3){
    p_eff[i] <- Z[i] * p
    for(j in 1:2){
    Y[(2 * (i - 1) + j)] ~ dbern(p_eff[i])		
    }
    } 
    }
    ",fill = TRUE)
sink()
#-----------------------------------------------------------------------------#
# run comparisons using equations in ms for example
compare <- function(par, sd = 1){
  set.seed(sd)
  psi0 <- par[1]
  gam <- par[2]
  eps <- par[3]
  p <- par[4]
  out <- numeric()
  out[1] <- psi0 * (1 - eps) + (1 - psi0) * gam
  out[2] <- (out[1] * (1 - p) ^ 2) / (out[1] * (1 - p) ^ 2 + (1 - out[1]))
  out[3] <- ((1 - eps) * ((1 - p) ^ 2)) / ((1 - eps) * ((1 - p) ^ 2) + eps)
  out[4] <- ((1 - eps) * ((1 - p) ^ 2) * (1 - eps) * p ^ 2) / 
    ((1 - eps) * ((1 - p) ^ 2) * (1 - eps) * (p ^ 2) + eps * gam * p ^ 2)
  Y <- c(0, 1, 0, 0, 1, 1)
  data <- list(Y = Y, p = p, psi0 = psi0, gam = gam, eps = eps)
  params <- c("Z")
  JD_inits <- function() list(Z = c(1, rbinom(1, 1, .5), 1))
  fit_JD <- jags(data, JD_inits, params, "JD.jags", n.chains = 100, n.iter = 1000)
  out[5] <- fit_JD$BUGSoutput$summary[2,1]
  return(out)
}

input <- cbind(c(.75, .75, .75, .75), c(.6, .6, .1, .6),
               c(.2, .8, .2, .2), c(0.1, 0.1, 0.1, 0.6))
output <- matrix(NA, ncol = 5, nrow = 4)
for(j in 1:4){
  output[j,] <- compare(input[j,])
}
# write.csv(cbind(input, output), "sim0.csv")
#-----------------------------------------------------------------------------#
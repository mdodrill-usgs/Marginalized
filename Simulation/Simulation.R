###############################################################################
#                                                                     Spring 19
#  Simulation of N-Augmentaion vs. Negative Binomail Abundance Estimates
#
#  Notes:
#  * 
#
###############################################################################
# simulation function
sim <- function(p = 0.1, Nind = 1000, Nvisits = 3, sd = 1){
  set.seed(sd)
  Y <- matrix(NA, nrow = Nind, ncol = Nvisits)
  for(i in 1:Nvisits){
    Y[,i] <- rbinom(Nind, 1, p)
  }
  Y[rowSums(Y) > 0,]
}

Nvisits <- 3
params <- c("N")
#-----------------------------------------------------------------------------#
# jags code - Naug
sink("N_JD.jags")
cat("
    model{
    for(i in 1:M){
    Z[i] ~ dbern(omega)
    p_eff[i] <- Z[i] * p
    for(j in 1:Nvisits){
    y_aug[i,j] ~ dbern(p_eff[i])		
    }
    } 
    N <- sum(Z[])
    omega ~ dunif(0, 1)
    p ~ dunif(0, 1)
    }
    ",fill = TRUE)
sink()

# JD inits
JD_inits <- function() list(Z = rep(1, nrow(y_aug)))

#-----------------------------------------------------------------------------#
# jags code - marginalized with Negative Binomial used to estimate N
sink("N_JM.jags")
cat("
    model{
    for(i in 1:Nobs){
    pz[i] <- pow(p, y[i]) * pow(1 - p, 3 - y[i]) / pstar
    o[i] ~ dbern(pz[i])
    }
    U ~ dnegbin(pstar, Nobs)
    N <- Nobs + U
    p ~ dunif(0, 1)
    pstar <-1 - pow(1 - p, 3)
    }
    ",fill = TRUE)
sink()
#-----------------------------------------------------------------------------#
# simulate and fit
OUT_s1 <- array(NA, dim = c(100, 2, 10))
for(j in 1:100){
  Y <- sim(sd = j)
  naug <- 2000
  y_aug <- rbind(Y, array(0, dim = c(naug, Nvisits)))
  win.data <- list(y_aug = y_aug, M = nrow(y_aug), Nvisits = Nvisits)
  win.data2 <-list(y = rowSums(Y), Nobs = nrow(Y), o = rep(1,nrow(Y)))
  t1 <- proc.time()
  sim3_JD <- jags.parallel(win.data, JD_inits, params, "N_JD.jags",
                           n.chains = 3, n.iter = 100000)
  t2<-proc.time()
  sim3_JM <- jags.parallel(win.data2, init=NULL, params, "N_JM.jags",
                           n.chains = 3, n.iter = 1000)
  t3 <- proc.time()
  OUT_s1[j,1,] <- c(sim3_JD$BUGSoutput$summary[1,], (t2 - t1)[3])
  OUT_s1[j,2,] <- c(sim3_JM$BUGSoutput$summary[1,], (t3 - t2)[3])
}

OUT_s2 <- array(NA, dim = c(100, 2, 10))
for(j in 46:100){
  Y <- sim(p = 0.5, sd = j)
  naug <- 1000
  y_aug <- rbind(Y, array(0, dim = c(naug, Nvisits)))
  win.data <- list(y_aug = y_aug, M = nrow(y_aug), Nvisits = Nvisits)
  win.data2 <- list(y = rowSums(Y), Nobs = nrow(Y), o = rep(1, nrow(Y)))
  t1 <- proc.time()
  sim3_JD <- jags.parallel(win.data, JD_inits, params, "N_JD.jags",
                           n.chains = 3, n.iter = 100000)
  t2 <- proc.time()
  sim3_JM <- jags.parallel(win.data2, init=NULL, params, "N_JM.jags",
                           n.chains = 3, n.iter = 1000)
  t3 <- proc.time()
  OUT_s2[j,1,] <- c(sim3_JD$BUGSoutput$summary[1,], (t2 - t1)[3])
  OUT_s2[j,2,] <- c(sim3_JM$BUGSoutput$summary[1,], (t3 - t2)[3])
}

#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Code for autologistic simulation summarized in Figure 5
#
#  Notes:
#   *  
#
###############################################################################

# creates matrix specifying neighbors of each site 
N <- matrix(0, nrow = 100, ncol = 100)
for(i in 1:100){
  t <- c(i - 2, i - 1, i + 1, i + 2)
  t <- subset(t, t > 0 & t < 101)
  N[i,t] <- 1 / length(t)
}

# returns the xth unique combination of 1's and 0's for N slots
makeZ <- function(x, N){
  i <- 0
  string <- rep(0, N)
  while(x > 0){
    string[N - i] <- x %% 2
    x <- x %/% 2
    i <- i + 1 
  }
  return(string)
}

# create all possible temporal patterns of latent occupancy (potential Z's) over
# six occasions (five intervals)
potZ <- matrix(0, nrow = 2 ^ 6, ncol = 6)
for(j in 1:2 ^ 6){
  potZ[j,] <- makeZ(j,6)
}

pZarray <- array(0, dim = c(64, 6, 2))
for(i in 1:64){
  for(t in 1:6){
    pZarray[i,t,1] <- potZ[i,t]
    pZarray[i,t,2] <- 1-potZ[i,t]
  }
}

# Fully conditional code requires a series of functions to iterate between
# making parameter estimates and conditional estimates, then treating
# conditional estimates as data and rerunning, likelihood for fully conditional
# in first iteration

fc_lik0 <- function(par){
  psi0 <- plogis(par[1])
  eps <- plogis(par[2])
  p <- plogis(par[3])
  g_int <- par[4]
  g <- matrix(0, nrow = 100, ncol = 5)
  for(t in 1:5){
    g[,t] <- plogis(g_int)
  }
  cp <- array(0, dim = c(100, 6, 2))
  zeta <- array(0, dim = c(100, 6, 2))
  po <- matrix(0, nrow = 2, ncol = 2)
  po[1,1] <- p
  po[1,2] <- (1 - p)
  po[2,1] <- 0
  po[2,2] <- 1
  for(t in 1:6){
    cp[,t,1] <- po[1,Y[,(2 * t - 1)]] * po[1,Y[,(2 * t)]]
    cp[,t,2] <- po[2,Y[,(2 * t - 1)]] * po[2,Y[,(2 * t)]]
  }
  zeta[,1,1] <- psi0 * cp[,1,1]
  zeta[,1,2] <- (1 - psi0) * cp[,1,2]
  for(t in 1:5){
    zeta[,(t + 1),1] <- (zeta[,t,1] * (1 - eps) + zeta[,t,2] * g[,t]) * cp[,(t + 1),1]
    zeta[,(t + 1),2] <- (zeta[,t,1] * eps + zeta[,t,2] * (1 - g[,t])) * cp[,(t + 1),2]
  }	
  nll <- 0
  for(k in 1:100){
    nll <- nll - 1 * log(sum(zeta[k,6,]))
  }
  nll
}

# this function creates fully conditional estimates given parameter estimates
fc_pred <- function(par, cX){
  psi0 <- plogis(par[1])
  eps <- plogis(par[2])
  p <- plogis(par[3])
  g_int <- par[4]
  g_auto <- par[5]
  g <- matrix(0, nrow = 100, ncol = 5)
  for(t in 1:5){
    g[,t] <- plogis(g_int + g_auto * (cX[,t] %*% t(N)))
  }
  out <- matrix(NA, nrow = 100, ncol = 5)
  po <- matrix(0, nrow = 2, ncol = 3)
  po[1,1] <- (1 - p) ^ 2
  po[1,2] <- 2 * p * (1 - p)
  po[1,3] <- p ^ 2
  po[2,1] <- 1
  y <- ifelse(Y == 2, 0, 1)
  sy <- 1 + cbind(rowSums(y[,1:2]), rowSums(y[,3:4]), rowSums(y[,5:6]), 
                  rowSums(y[,7:8]), rowSums(y[,9:10]), rowSums(y[,11:12]))
  for(k in 1:100){
    cp <- matrix(0, nrow = 6, ncol = 2)
    for(t in 1:6){
      cp[t,1] <- po[1,sy[k,t]]
      cp[t,2] <- po[2,sy[k,t]]
    }
    unc <- matrix(NA, nrow = 64, ncol = 6)
    unc[,1] <- log(potZ[,1] * psi0 + (1 - potZ[,1]) * (1 - psi0))
    for(t in 1:5){
      unc[,(t + 1)] <- log((1 - potZ[,t]) * (g[k,t] * potZ[,(t + 1)] + (1 - g[k,t]) * (1 - potZ[,(t + 1)])) +
                             potZ[,t] * (eps * (1-potZ[,(t + 1)]) + (1-eps) * potZ[,(t + 1)]))
    }
    UNC <- exp(rowSums(unc))
    tcp <- numeric()
    for(i in 1:64){
      temp <- rowSums(pZarray[i,,] * cp)
      if(sum(ifelse(temp == 0, 1, 0)) > 0){
        (tcp[i] <- 0)
      } else {
        tcp[i] <- UNC[i] * exp(sum(log(temp)))
      }
    }
    for(t in 1:5){
      out[k,t] <- sum(subset(tcp,potZ[,t] == 1)) / sum(tcp)
    }
  }
  return(out)
}

# this function calculates the likelihood after the first iteration and assuming
# cX is present in the workspace
fc_lik <- function(par){
  psi0 <- plogis(par[1])
  eps <- plogis(par[2])
  p <- plogis(par[3])
  g_int <- par[4]
  g_auto <- par[5]
  g <- matrix(0, nrow = 100, ncol = 5)
  for(t in 1:5){
    g[,t] <- plogis(g_int + g_auto * (cX[,t] %*% t(N)))
  }
  cp <- array(0, dim = c(100, 6, 2))
  zeta <- array(0, dim = c(100, 6, 2))
  po <- matrix(0, nrow = 2, ncol = 2)
  po[1,1] <- p
  po[1,2] <- (1 - p)
  po[2,1] <- 0
  po[2,2] <- 1
  for(t in 1:6){
    cp[,t,1] <- po[1,Y[,(2 * t-1)]] * po[1,Y[,(2 * t)]]
    cp[,t,2] <- po[2,Y[,(2 * t-1)]] * po[2,Y[,(2 * t)]]
  }
  zeta[,1,1] <- psi0 * cp[,1,1]
  zeta[,1,2] <- (1 - psi0) * cp[,1,2]
  for(t in 1:5){
    zeta[,(t + 1),1] <- (zeta[,t,1] * (1-eps) + (zeta[,t,2]) * g[,t]) * cp[,(t + 1),1]
    zeta[,(t + 1),2] <- (zeta[,t,1] * eps + (zeta[,t,2]) * (1-g[,t])) * cp[,(t + 1),2]
  }	
  nll <- 0
  for(k in 1:100){
    nll <- nll + -1 * log(sum(zeta[k,6,]))
  }
  nll
}

# function to fit first version of fully conditional and produce first
# conditional estimates
iter0 <- function(){
  m <- optim(c(0, 0, 0, 0), fc_lik0, method = "BFGS")
  tX <- matrix(0, nrow = 100, ncol = 5)
  CX <- fc_pred(c(m$par, 0), cX = tX)
  return(CX)
}
# End of fully conditional code - next comes function for forward conditional
#-----------------------------------------------------------------------------#
# likelihood for forward conditional
fwc_lik <- function(par){
  psi0 <- plogis(par[1])
  eps <- plogis(par[2])
  p <- plogis(par[3])
  g_int <- par[4]
  g_auto <- par[5]
  cp <- array(0, dim = c(100, 6, 2))
  zeta <- array(0, dim = c(100, 6, 2))
  po <- matrix(0, nrow = 2, ncol = 2)
  po[1,1] <- p
  po[1,2] <- (1 - p)
  po[2,1] <- 0
  po[2,2] <- 1
  for(t in 1:6){
    cp[,t,1] <- po[1,Y[,(2 * t - 1)]] * po[1,Y[,(2 * t)]]
    cp[,t,2] <- po[2,Y[,(2 * t - 1)]] * po[2,Y[,(2 * t)]]
  }
  zeta[,1,1] <- psi0 * cp[,1,1]
  zeta[,1,2] <- (1-psi0) * cp[,1,2]
  g <- matrix(0, nrow = 100, ncol = 5)
  psi <- matrix(0, nrow = 100, ncol = 5)
  for(t in 1:5){
    psi[,t] <- zeta[,t,1] / (zeta[,t,1] + zeta[,t,2])
    g[,t] <- plogis(g_int + g_auto * (psi[,t] %*% t(N)))
    zeta[,(t + 1),1] <- (zeta[,t,1] * (1 - eps) + (zeta[,t,2]) * g[,t]) * cp[,(t + 1),1]
    zeta[,(t + 1),2] <- (zeta[,t,1] * eps + (zeta[,t,2]) * (1 - g[,t])) * cp[,(t + 1),2]
  }	
  nll <- 0
  for(k in 1:100){
    nll <- nll + -1 * log(sum(zeta[k,6,]))
  }
  nll
}

#-----------------------------------------------------------------------------#
# code to simulate a dataset
sim <- function(psi0 = c(rep(0, 50), rep(1, 50)), g_int = -1, g_auto = 2,
                eps = 0.5, p = 0.5, Nvisits = 2, Nint = 5, neigh = N){
  Z <- matrix(0, ncol = (1 + Nint), nrow = length(psi0))
  Y <- matrix(0, ncol = (1 + Nint) * Nvisits, nrow = length(psi0))
  Z[,1] <- psi0
  for(k in 1:Nvisits){
    Y[,k] <- rbinom(100, 1, Z[,1] * p)
  }
  for(t in 1:Nint){
    x <- Z[,t] %*% t(neigh)
    gam <- plogis(g_int + g_auto * x)
    Z[,(t + 1)] <- rbinom(100, 1, Z[,t] * (1 - eps) + (1 - Z[,t]) * gam)
    for(k in 1:Nvisits){
      Y[,(t * Nvisits + k)] <- rbinom(100, 1, Z[, (t + 1)] * p)
    }
  }
  return(Y)
}

# create arrays to hold simulations
ss <- array(NA,dim = c(5,2,1000))  # full conditional - random occupancy
ss2 <- array(NA,dim = c(5,2,1000)) # forward conditional - random occupancy
ss3 <- array(NA,dim = c(5,2,1000)) # full conditional - clustered occupancy
ss4 <- array(NA,dim = c(5,2,1000)) # forward conditional - clustered occupancy
#-----------------------------------------------------------------------------#
for(j in 1:1000){
  set.seed(j) # set seed for simulations
  Y <- ifelse(sim(psi0 = rbinom(100, 1, 0.5)) == 0, 2, 1) # data where starting occupancy random
  # fit full conditional through iterative process
  cX <- iter0()
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  # usually converges after a few iterations, but we ran six times to ensure convergence
  ss[,1,j] <- m$par # store results
  ss[,2,j] <- sqrt(diag(solve(m$hessian)))
  # fit forward conditional and store results
  m2 <- optim(rep(0,5), fwc_lik, method = "BFGS", hessian = TRUE)
  ss2[,1,j] <- m2$par
  ss2[,2,j] <- sqrt(diag(solve(m2$hessian)))
  set.seed(j)
  # simulate clustered data
  Y <- ifelse(sim() == 0, 2, 1)
  # fit fully conditional and store
  cX <- iter0()
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  cX <- fc_pred(m$par, cX = cX)
  m <- optim(rep(0,5), fc_lik, method = "BFGS", hessian = TRUE)
  ss3[,1,j] <- m$par
  ss3[,2,j] <- sqrt(diag(solve(m$hessian)))
  # fit forward conditional and store
  m4 <- optim(rep(0,5), fwc_lik, method = "BFGS", hessian = TRUE)
  ss4[,1,j] <- m4$par
  ss4[,2,j] <- sqrt(diag(solve(m4$hessian)))
}

#-----------------------------------------------------------------------------#
sumz <- matrix(NA, ncol = 10, nrow = 4)
truth <- c(0, 0, 0, -1, 2)
for(k in 1:5){
  sumz[1,(2 * k - 1)] <- mean(ss[k,1,] - truth[k])
  sumz[2,(2 * k - 1)] <- mean(ss2[k,1,] - truth[k])
  sumz[3,(2 * k - 1)] <- mean(ss3[k,1,] - truth[k])
  sumz[4,(2 * k - 1)] <- mean(ss4[k,1,] - truth[k])
  sumz[1,(2 * k)] <- mean(ifelse(ss[k,1,] - 1.96 * ss[k,2,] < 
                                   truth[k] & ss[k,1,] + 1.96 * ss[k,2,] > truth[k], 1, 0))
  sumz[2,(2 * k)] <- mean(ifelse(ss2[k,1,] - 1.96 * ss2[k,2,] < 
                                   truth[k] & ss2[k,1,] + 1.96 * ss2[k,2,] > truth[k], 1, 0))
  sumz[3,(2 * k)] <- mean(ifelse(ss3[k,1,] - 1.96 * ss3[k,2,] < 
                                   truth[k] & ss3[k,1,] + 1.96 * ss3[k,2,] > truth[k], 1, 0))
  sumz[4,(2 * k)] <- mean(ifelse(ss4[k,1,] - 1.96 * ss4[k,2,] <
                                   truth[k]&ss4[k,1,] + 1.96 * ss4[k,2,] > truth[k], 1, 0))
}

#-----------------------------------------------------------------------------#
# coverage
cover <- array(NA, dim = c(4, 2, 2))
cover[1,1,1] <- sum(ifelse(ss[4,1,] + 0.67 * ss[4,2,] > (-1) &
                             ss[4,1,] - 0.67 * ss[4,2,] < (-1), 1, 0)) / 1000
cover[1,1,2] <- sum(ifelse(ss[4,1,] + 1.96 * ss[4,2,] > (-1) & 
                             ss[4,1,] - 1.96 * ss[4,2,] < (-1), 1, 0)) / 1000
cover[1,2,1] <- sum(ifelse(ss[5,1,] + 0.67 * ss[5,2,] > (2) & 
                             ss[5,1,] - 0.67 * ss[5,2,] < (2), 1, 0)) / 1000
cover[1,2,2] <- sum(ifelse(ss[5,1,] + 1.96 * ss[5,2,] > (2) & 
                             ss[5,1,] - 1.96 * ss[5,2,] < (2), 1, 0)) / 1000
#
cover[2,1,1] <- sum(ifelse(ss2[4,1,] + 0.67 * ss2[4,2,] > (-1) & 
                             ss2[4,1,] - 0.67 * ss2[4,2,] < (-1), 1, 0)) / 1000
cover[2,1,2] <- sum(ifelse(ss2[4,1,] + 1.96 * ss2[4,2,] > (-1) & 
                             ss2[4,1,] - 1.96 * ss2[4,2,] < (-1), 1, 0)) / 1000
cover[2,2,1] <- sum(ifelse(ss2[5,1,] + 0.67 * ss2[5,2,] > (2) & 
                             ss2[5,1,] - 0.67 * ss2[5,2,] < (2), 1, 0)) / 1000
cover[2,2,2] <- sum(ifelse(ss2[5,1,] + 1.96 * ss2[5,2,] > (2) & 
                             ss2[5,1,] - 1.96 * ss2[5,2,] < (2), 1, 0)) / 1000
#
cover[3,1,1] <- sum(ifelse(ss3[4,1,] + 0.67 * ss3[4,2,] > (-1) & 
                             ss3[4,1,] - 0.67 * ss3[4,2,] < (-1), 1, 0)) / 1000
cover[3,1,2] <- sum(ifelse(ss3[4,1,] + 1.96 * ss3[4,2,] > (-1) & 
                             ss3[4,1,] - 1.96 * ss3[4,2,] < (-1), 1, 0)) / 1000
cover[3,2,1] <- sum(ifelse(ss3[5,1,] + 0.67 * ss3[5,2,] > (2) & 
                             ss3[5,1,] - 0.67 * ss3[5,2,] < (2), 1, 0)) / 1000
cover[3,2,2] <- sum(ifelse(ss3[5,1,] + 1.96 * ss3[5,2,] > (2) & 
                             ss3[5,1,] - 1.96 * ss3[5,2,] < (2), 1, 0)) / 1000
#
cover[4,1,1] <- sum(ifelse(ss4[4,1,] + 0.67 * ss4[4,2,] > (-1) & 
                             ss4[4,1,] - 0.67 * ss4[4,2,] < (-1), 1, 0)) / 1000
cover[4,1,2] <- sum(ifelse(ss4[4,1,] + 1.96 * ss4[4,2,] > (-1) & 
                             ss4[4,1,] - 1.96 * ss4[4,2,] < (-1), 1, 0)) / 1000
cover[4,2,1] <- sum(ifelse(ss4[5,1,] + 0.67 * ss4[5,2,] > (2) & 
                             ss4[5,1,] - 0.67 * ss4[5,2,] < (2), 1, 0)) / 1000
cover[4,2,2] <- sum(ifelse(ss4[5,1,] + 1.96 * ss4[5,2,] > (2) & 
                             ss4[5,1,] - 1.96 * ss4[5,2,] < (2), 1, 0)) / 1000
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  Simulation Code for Results in Table 2 
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
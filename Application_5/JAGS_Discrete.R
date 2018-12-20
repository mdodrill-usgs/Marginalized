###############################################################################
#                                                                        Dec 18
#  Dynamic Community Occupancy Models - Sky Islands  
#  Discrete JAGS version 
#
#  Notes:
#  * 
#
#  To do: 
#  * 
#
###############################################################################
library(R2jags)
# source(paste0(getwd(),"/Functions.R"), chdir = F)

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
  Ymat1 = as.matrix(read.csv(paste0(data.dir,names[q]))[,2:93])
  Y[,,q]=array(c(Ymat1),dim=c(dim(Y)[1],dim(Ymat1)[2]))
}

# habitat covariate (dimensions 92 pts) 
hab = as.matrix(read.csv(paste0(data.dir, 'habitat.CM.csv'), header = FALSE))[2:93,4] 
hab = as.integer(hab)
#-----------------------------------------------------------------------------#
sink("Cocc_JD.jags")
cat("
model {
  for(i in 1:(nspp)){
    a1[i] ~ dunif(0,1)
    alpha1[i] ~ dnorm(alpha1_mean, tau_alpha1)
    
    for(h in 1:7){
      alpha2[i,h] ~ dnorm(0, tau_alpha2)
    }
    
    alpha0[i] ~ dnorm(beta, tau_u)
    mu_eta[i] <- alpha + (rho*sigma_v / sigma_u) * (alpha0[i] - beta)
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
  
  psi_mean ~ dunif(0,1)
  beta <- log(psi_mean) - log(1 - psi_mean)
  p_mean ~ dunif(0,1)
  alpha <- log(p_mean) - log(1-p_mean)
  alpha1_mean ~ dunif(0,1)
  alpha1_mu <- log(alpha1_mean) - log(1-alpha1_mean)
  sigma_alpha1 ~ dunif(0,5)
  sigma_alpha2 ~ dunif(0,5)
  sigma_u ~ dunif(0,5)
  sigma_v ~ dunif(0,5)
  tau_alpha1 <- pow(sigma_alpha1,-2)
  tau_alpha2 <- pow(sigma_alpha2,-2)
  tau_u <- pow(sigma_u,-2)
  tau_v <- pow(sigma_v,-2)
  rho ~ dunif(-1,1)
  tau_eta <- tau_v / (1-pow(rho,2))
}
    ", fill = TRUE)
sink()   
#-----------------------------------------------------------------------------#
# tried to cut this down, but JAGS doesn't like it
# JD_inits = function(nspp, nsites, nyears, nsess, Y){
#   a1 = runif(nspp, 0.25, 1)
#   Ztemp = array(rbinom(nspp * 5 * nsites, size = 1, prob = 0.5), dim = c(nspp, nsites, 5))
#   z0temp = array(rbinom(nspp * nsites, size = 1, prob = 0.5), dim = c(nspp, nsites))
# 
#   for(t in 1:5){
#     for(i in 1:nspp){
#       for(k in 1:nsites){
#         if(Ztemp[i,k,t] == 0){
#           if(!is.na(Y[i,k,t])){
#             if(Y[i,k,t] > 0){
#               Ztemp[i,k,t] = 1
#             }
#           }
#         }
#       }
#     }
#   }
# 
#   list(a1 = a1, z0 = z0temp, Z = Ztemp)
# }
JD_inits = function() {
  nspp = 149
  nsites = 92
  nyears = 5  
  nsess = 3
  Y = array(NA, dim = c(nspp, nsites, nyears))
  
  names1 = c('ntmb.CM.1991.csv')  
  names2 = c('ntmb.CM.1992.csv')  
  names3 = c('ntmb.CM.1993.csv')  
  names4 = c('ntmb.CM.1994.csv')  
  names5 = c('ntmb.CM.1995.csv')  
  
  names = c(names1, names2, names3, names4, names5)
  
  for (q in 1:5){
    Ymat1 = as.matrix(read.csv(paste0(data.dir, names[q]))[,2:93])
    Y[,,q] = array(c(Ymat1), dim = c(dim(Y)[1], dim(Ymat1)[2]))
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
params <- c('psi_mean', 'p_mean', 'rho', 'sigma_v', 'sigma_u', 'sigma_alpha1',
            'sigma_alpha2', 'alpha1_mean','SR')

JD_data = list('nspp' = nspp, 'nsites' = nsites, 'nsess' = nsess,
               'nyears' = nyears, 'hab' = hab, 'Y' = Y)

# ~8h for 75000
t1 <- proc.time()
JD_Cocc <- jags.parallel(JD_data, inits = JD_inits,
                         params, "Cocc_JD.jags",
                         n.chains = 3, n.iter = 10,
                         export_obj_names = c("data.dir"))
t2 <- proc.time()
#-----------------------------------------------------------------------------#

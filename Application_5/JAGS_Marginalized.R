###############################################################################
#                                                                        Dec 18
#  Dynamic Community Occupancy Models - Sky Islands  
#  Marginalized JAGS version 
#
#  Notes:
#  * 
#
#  To do: 
#  * Clean up and format better!
#
###############################################################################
library(R2jags)
# source(paste0(getwd(),"/Functions.R"), chdir = F)
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
  Ymat1 = as.matrix(read.csv(paste0(data.dir,names[q]))[,2:93])
  Y[,,q]=array(c(Ymat1),dim=c(dim(Y)[1],dim(Ymat1)[2]))
}

visit <- ifelse(is.na(Y[1,,1]) == TRUE, 2, 1)
sY <- ifelse(is.na(Y) == TRUE, 1, Y + 1)

# habitat covariate (dimensions 92 pts) 
hab = as.matrix(read.csv(paste0(data.dir, 'habitat.CM.csv'), header = FALSE))[2:93,4] 
hab = as.integer(hab)

#-----------------------------------------------------------------------------#
sink("Cocc_JM.jags")
cat("
model {
  for (i in 1:(nspp)) {
    a1[i] ~ dunif(0,1)
    z0[i,1]<-a1[i]
    z0[i,2]<-1-a1[i]
    alpha1[i] ~ dnorm(alpha1_mean,tau_alpha1)
    for (h in 1:7){alpha2[i,h]~dnorm(0,tau_alpha2)}    
    alpha0[i] ~ dnorm(beta, tau_u)
    mu_eta[i] <- alpha + (rho*sigma_v/sigma_u)*(alpha0[i] - beta)
    beta0[i] ~ dnorm(mu_eta[i], tau_eta)
    logit(p[i]) <- beta0[i] #detection
    po[i,1,1,2]=1;
    po[i,2,1,2]=1;
    po[i,1,2,2]=0;
    po[i,2,2,2]=0;
    for (j in 1:nsess){
      po[i,1,(j+1),1]<-0
      po[i,2,(j+1),1]<-(p[i]^j)*(1-p[i])^(nsess-j)
    }
    po[i,1,1,1]<-1
    po[i,2,1,1]<-(1-p[i])^(nsess)
    for (k in 1:nsites) {
      logit(tr[i,k,1,2])<-alpha0[i]+alpha2[i,hab[k]]
      logit(tr[i,k,2,2])<-alpha0[i]+alpha1[i]+alpha2[i,hab[k]]
      tr[i,k,1,1]<-1-tr[i,k,1,2]
      tr[i,k,2,1]<-1-tr[i,k,2,2]
      pz[i,k,1,1]<-inprod(z0[i,],tr[i,k,,1])*po[i,1,sY[i,k,1],visit[k]]
      pz[i,k,1,2]<-inprod(z0[i,],tr[i,k,,2])*po[i,2,sY[i,k,1],visit[k]]
      Z[i,k,1]<-pz[i,k,1,2]/(pz[i,k,1,1]+pz[i,k,1,2])
      Z2[i,k,1]~dbern(Z[i,k,1])
      for (t in 1:(nyears-1)){
        pz[i,k,(t+1),1]<-inprod(pz[i,k,t,],tr[i,k,,1])*po[i,1,sY[i,k,(t+1)],1]
        pz[i,k,(t+1),2]<-inprod(pz[i,k,t,],tr[i,k,,2])*po[i,2,sY[i,k,(t+1)],1]
        Z[i,k,(t+1)]<-pz[i,k,(t+1),2]/(pz[i,k,(t+1),1]+pz[i,k,(t+1),2])
        Z2[i,k,(t+1)]~dbern(Z[i,k,(t+1)])
      }
      lik[i,k]<-sum(pz[i,k,nyears,])
      ones[i,k]~dbern(lik[i,k])
    }}
  for (k in 1:nsites){  
    for (t in 1:nyears){
      SR[k,t] <- sum(Z2[1:nspp,k,t])
    }}
  psi_mean ~ dunif(0,1)
  beta <- log(psi_mean) - log(1-psi_mean)
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
  tau_eta <- tau_v/(1-pow(rho,2)) 
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
#-----------------------------------------------------------------------------#
params <- c('psi_mean', 'p_mean', 'rho', 'sigma_v', 'sigma_u', 'sigma_alpha1',
            'sigma_alpha2', 'alpha1_mean','SR')

JM_data = list('nspp' = nspp, 'nsites' = nsites, 'nsess' = nsess,
               'nyears' = nyears, 'hab' = hab, 'ones' = matrix(1, nrow = nspp, ncol = nsites),
               'sY' = sY, 'visit' = visit)

t1 <- proc.time()
JM_Cocc<- jags.parallel(JM_data, inits = JM_inits, params, 'Cocc_JM.jags',
                        n.chains = 3, n.iter=10)
t2 <- proc.time()
#-----------------------------------------------------------------------------#

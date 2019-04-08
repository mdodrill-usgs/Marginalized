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
setwd("U:\\Desktop\\Fish_Git\\Marginalized")

source(paste0(getwd(),"/Functions.R"), chdir = F)

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
sink("Cocc_JM_2.jags")
cat("
  model {
  for (i in 1:(nspp)) {
    a1[i] ~ dunif(0,1)
    z0[i,1]<-a1[i]
    z0[i,2]<-1-a1[i]
    alpha1[i] ~ dnorm(alpha1_mu,tau_alpha1)
    for (h in 1:7){alpha2[i,h]~dnorm(0,tau_alpha2)}    
    alpha0[i] ~ dnorm(beta, tau_u)
    mu_eta[i] <- alpha + (rho*sigma_v/sigma_u)*(alpha0[i] - beta)
    beta0[i] ~ dnorm(mu_eta[i], tau_eta)
    logit(p[i]) <- beta0[i] #detection
    po[i,1,1,2]=1;
    po[i,2,1,2]=1;
    for (j in 1:nsess){
    po[i,1,(j+1),1]<-0
    po[i,2,(j+1),1]<-(p[i]^j)*(1-p[i])^(nsess-j)
    }
    po[i,1,1,1]<-1
    po[i,2,1,1]<-(1-p[i])^(nsess)
    for (k in 1:nhab) {
    logit(tr[i,k,1,2])<-alpha0[i]+alpha2[i,k]
    logit(tr[i,k,2,2])<-alpha0[i]+alpha1[i]+alpha2[i,k]
    tr[i,k,1,1]<-1-tr[i,k,1,2]
    tr[i,k,2,1]<-1-tr[i,k,2,2]
    }}
    for (i in 1:NuDH){
    pz[i,1,1]=(z0[spp[i],1]*tr[spp[i],hab2[i],1,1]+z0[spp[i],2]*tr[spp[i],hab2[i],2,1])*po[spp[i],1,sY2[i,1],visit2[i]];
    pz[i,1,2]=(z0[spp[i],1]*tr[spp[i],hab2[i],1,2]+z0[spp[i],2]*tr[spp[i],hab2[i],2,2])*po[spp[i],2,sY2[i,1],visit2[i]];
    Z[i,1]=pz[i,1,2]/(pz[i,1,1]+pz[i,1,2]);
    for (t in 1:(nyears-1)){
    pz[i,(t+1),1]=(pz[i,t,1]*tr[spp[i],hab2[i],1,1]+pz[i,t,2]*tr[spp[i],hab2[i],2,1])*po[spp[i],1,sY2[i,(t+1)],1];
    pz[i,(t+1),2]=(pz[i,t,1]*tr[spp[i],hab2[i],1,2]+pz[i,t,2]*tr[spp[i],hab2[i],2,2])*po[spp[i],2,sY2[i,(t+1)],1];
    Z[i,(t+1)]=pz[i,(t+1),2]/(pz[i,(t+1),1]+pz[i,(t+1),2]);
    }
    lik[i]<-sum(pz[i,nyears,])
    fr[i]~dbin(lik[i],FR[i])
    }
    for (k in 1:nsites){  
    for (t in 1:nyears){
    for (i in 1:nspp){Z2[i,k,t]~dbern(Z[lookup[i,k],t])}
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

pasted <- function(x){
  paste0(x[1], x[2], x[3], x[4], x[5], x[6], x[7], collapse = "")
}

unpaste <- function(x){
  as.numeric(c(substr(x,1,1), substr(x,2,2), substr(x,3,3),
               substr(x,4,4), substr(x,5,5), substr(x,6,6), substr(x,7,7)))
}

#-----------------------------------------------------------------------------#
FR <- numeric()
spp <- numeric()
sY2 <- rep(NA,5)
visit2 <- numeric()
hab2 <- numeric()
lookup <- matrix(NA,nspp,nsites)
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
params <- c('psi_mean', 'p_mean', 'rho', 'sigma_v', 'sigma_u', 'sigma_alpha1',
            'sigma_alpha2', 'alpha1_mean','SR', 'Z2')

# JM_data = list('nspp' = nspp, 'nsites' = nsites, 'nsess' = nsess,
#                'nyears' = nyears, 'hab' = hab, 'ones' = matrix(1, nrow = nspp, ncol = nsites),
#                'sY' = sY, 'visit' = visit)

JM_data =list('nspp' = nspp, 'nsites' = nsites,'nsess' = nsess, 'nyears' = nyears,
               'hab2' = hab2, 'fr' = FR, 'FR' = FR, 'nhab' = 7, 'NuDH' = length(FR),
               'spp' = spp, 'lookup' = lookup, 'sY2' = sY2, 'visit2' = visit2)


# 
t1 <- proc.time()
JM_Cocc <- jags.parallel(JM_data, inits = JM_inits, params, 'Cocc_JM_2.jags',
                        n.chains = 3, n.iter = 10)
t2 <- proc.time()
#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

n.core = 3  # really n.core * 3

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used 
cllibs <- clusterEvalQ(cl1, c(library(R2jags)))

all.t1 = proc.time()
n.runs = 3

my.n.iter = c(10000)

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
    
    my.env = environment()
    
    JM_Cocc<- jags.parallel(JM_data, inits = JM_inits,
                            params, 'Cocc_JM_2.jags',
                            n.chains = 3, n.iter = iter.in,
                            export_obj_names = c("iter.in", "seed"),
                            jags.seed = seed, envir = my.env)
    
    t2 <- proc.time()
    
    attr(JM_Cocc, 'time') <- (t2 - t1)[3]
    
    JM_Cocc
    
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
# 
all.jags.m.5 = all.out
# 
rm(list=setdiff(ls(), "all.jags.m.5"))
# 
save.image("U:/Desktop/Fish_Git/Marginalized/Application_5/working_runs/Cocc_JM_2_runs_5_test.RData")
#-----------------------------------------------------------------------------#
# end
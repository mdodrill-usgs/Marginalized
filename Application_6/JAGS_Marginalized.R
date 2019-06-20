###############################################################################
#                                                                        Dec 18
#  Combination Occupancy and N-Mixture Model
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

setwd(paste0(getwd(), '/Application_3'))

data.dir = paste0(getwd(), "/Data/")

# import data
raw.data <- read.csv(paste0(data.dir, "barred_owl_survey_data.csv"), header = TRUE) 
# extract all BO data
nYears<-22
tBO<-as.matrix(raw.data[,grep("BO",colnames(raw.data))])
nSites<-dim(tBO)[1]
Y<-tBO[,241:246]
visited<-ifelse(is.na(Y)==TRUE,0,1)
Y<-ifelse(is.na(Y)==TRUE,0,Y)
maxY<-matrix(0,nrow=nSites,ncol=2)
for (t in 1:2){
  maxY[,t]<-apply(Y[,(3*(t-1)+1):(3*t)],1,max)}
num_det<-matrix(0,nrow=nSites,ncol=20)
num_trial<-matrix(0,nrow=nSites,ncol=20)
num_nodet<-matrix(0,nrow=nSites,ncol=20)
for (t in 1:20){
  num_trial[,t]<-rowSums(ifelse(is.na(tBO[,((t-1)*12+1):(t*12)])==TRUE,0,1))
  num_det[,t]<-rowSums(ifelse(is.na(tBO[,((t-1)*12+1):(t*12)])==TRUE,0,tBO[,((t-1)*12+1):(t*12)]))
  num_nodet[,t]<-num_trial[,t]-num_det[,t]
}
maxdet<-max(num_det)
maxnodet<-max(num_nodet)
# Pull out and standardize covariates
PHAB<-as.matrix(raw.data[,grep("PH",colnames(raw.data))]) # PHAB = percent older riparian growth forest)
X<- (PHAB - mean(PHAB)) / sd(PHAB)
tOV<-as.matrix(raw.data[,grep("OV",colnames(raw.data))])
visited<-ifelse(tOV==0,0,1)
OV<-ifelse(tOV==0,0,log(tOV))
params<-c("a0","b0","b1","g0","g1","g2","p_occ","lambda","N_mean")
Ninf=15
C<-array(0,dim=c(nSites,2,(Ninf+1)))
for (i in 1:nSites){
  for (j in 1:(Ninf+1)){
    C[i,1,j]<-choose((j-1),Y[i,1])*choose((j-1),Y[i,2])*choose((j-1),Y[i,3])
    C[i,2,j]<-choose((j-1),Y[i,4])*choose((j-1),Y[i,5])*choose((j-1),Y[i,6])
  }}
iC<-ifelse(C==0,0,1)
lC<-ifelse(C==0,0,log(C))

#-----------------------------------------------------------------------------#
sink("Combo_JM.jags")
cat("
model {
  for (j in 1:(Ninf+1)){
    po[j,2]=1-exp((j-1)*log(1-p_occ))
    po[j,1]=1-po[j,2]
    pN[j]=dpois((j-1),lambda)
    Ns[j]=(j-1)
  }
  for (i in 1:nSites){
    for (t in 1:20){
      for (j in 1:(Ninf+1)){
        cp[i,t,j]<-(po[j,2]^num_det[i,t])*(po[j,1]^num_nodet[i,t])
      }}
    for (t in 1:2){ 
      for (k in 1:3){
        cloglog(tpc[i,(3*(t-1)+k)]) <- a0 + OV[i,(3*(t-1)+k)]
        pc[i,(3*(t-1)+k)]<-tpc[i,(3*(t-1)+k)]*visited[i,(3*(t-1)+k)]
      }
      for (j in 1:(Ninf+1)){
        cp[i,(t+20),j]<-dbin(Y[i,(3*t-2)],pc[i,(3*t-2)],(j-1))*
          dbin(Y[i,(3*t-1)],pc[i,(3*t-1)],(j-1))*dbin(Y[i,(3*t)],pc[i,(3*t)],(j-1))
      }}}
  for (i in 1:nSites){
    for (j in 1:(Ninf+1)){
      pz[i,1,j]=pN[j]*cp[i,1,j]
    }
    N_mn[i,1]<-inprod(Ns,pz[i,1,])/sum(pz[i,1,])
  }
  for(t in 1:(nYears-1)) {
    N_mean[t] <- mean(N_mn[1:nSites,t]) - 1
    log(gamma[t]) <- g0 + g1*N_mean[t] + g2*N_mean[t]*N_mean[t]
    for (j in 1:(Ninf+1)){
      for (k in 1:(Ninf+1)){t_gtr[t,j,k]=dpois((k-j),gamma[t])}
      for (k in 1:(Ninf+1)){gtr[t,j,k]=t_gtr[t,j,k]/sum(t_gtr[t,j,])}
    }
    for(i in 1:nSites) {
      logit(omega[i,t]) = b0 + b1*X[i,t]
      for (j in 1:(Ninf+1)){
        for (k in 1:(Ninf+1)){
          str[i,t,j,k] = pz[i,t,j]*dbin((k-1),omega[i,t],(j-1))
        }}
      for (j in 1:(Ninf+1)){
        for (k in 1:(Ninf+1)){
          tr[i,t,j,k]<-sum(str[i,t,j:(Ninf+1),j])*gtr[t,j,k]
        }}
      for (j in 1:(Ninf+1)){
        pz[i,(t+1),j]<-sum(tr[i,t,1:j,j])*cp[i,t,j]
      }
      N_mn[i,(t+1)]<-inprod(Ns,pz[i,(t+1),])/sum(pz[i,(t+1),])
    }}
  
  for (i in 1:nSites){
    lik[i]<-sum(pz[i,22,])
    ones[i]~dbern(lik[i])
  }		
  lambda ~ dunif(0,1) # initial abundance
  p_occ ~ dunif(0,1) # detection
  b0 ~ dunif(-3,3) # intercept on survival
  b1 ~ dunif(-1,1) # slope on survival
  a0 ~ dunif(-3,3) # intercept on effort (p_count)
  g0 ~ dunif(-3,3) # intercept on gamma
  g1 ~ dunif(-1,1) # slope on mean(N)
  g2 ~ dunif(-1,1) # squared term on mean(N)
}
", fill = TRUE)
sink()   
#-----------------------------------------------------------------------------#
params<-c("a0","b0","b1","g0","g1","g2","p_occ","lambda","N_mean")

JM_data<-list(Y=Y,visited=visited,X=X,OV=OV,nYears=nYears,nSites=nSites,
              num_det=num_det,num_nodet=num_nodet,Ninf=Ninf,ones=rep(1,nSites),maxdet=maxdet,maxnodet=maxnodet)


# JM_combo<-jags(JM_data, inits=NULL, params, model.file="combo_JM.txt",n.iter=10000)

t1 <- proc.time()
JM_combo <- jags.parallel(JM_data, inits = NULL, params, 'Combo_JM.jags',
                         n.chains = 3, n.iter = 10)
t2 <- proc.time()

#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

n.core = 5  # really n.core * 3

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used 
cllibs <- clusterEvalQ(cl1, c(library(R2jags)))

all.t1 = proc.time()
n.runs = 5

my.n.iter = c(20000)

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
    
    JM_combo<- jags.parallel(JM_data, inits = NULL,
                            params, 'Combo_JM.jags',
                            n.chains = 3, n.iter = iter.in,
                            export_obj_names = c("iter.in", "seed"),
                            jags.seed = seed, envir = my.env)
    
    t2 <- proc.time()
    
    attr(JM_combo, 'time') <- (t2 - t1)[3]
    
    JM_combo
    
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

tmp = run.times(all.out)
tmp

all.jags.m.3 = all.out
  
rm(list=setdiff(ls(), "all.jags.m.3"))
  
save.image("U:/Desktop/Fish_Git/Marginalized/Application_3/working_runs/Combo2_JM_runs_3.RData")
#-----------------------------------------------------------------------------#
# end



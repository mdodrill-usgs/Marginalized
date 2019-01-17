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
# see Stan model 'Combo2_SM.stan' 
#-----------------------------------------------------------------------------#
params<-c("a0","b0","b1","g0","g1","g2","p_occ","lambda","N_mean")

SM_data<-list(Y=Y,visited=visited,X=X,OV=OV,nYears=nYears,nSites=nSites,Ns=c(0:Ninf),
              num_det=num_det,num_nodet=num_nodet,Ninf=Ninf,maxdet=maxdet,maxnodet=maxnodet,iC=iC,lC=lC)


t1 <- proc.time()
SM_Cocc <- stan("Combo2_SM.stan",
                data = SM_data,
                pars = params,
                chains = 1, iter = 10) 
t2 <- proc.time()
#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

n.core = 10 # really n.core * 3

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used 
cllibs <- clusterEvalQ(cl1, c(library(rstan)))

n.runs = 10

my.n.iter = c(1000)

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
    
    # SM.c <- stan(fit = SM_Cocc,     # faster to compile the model above...
    SM_Cocc <- stan("Combo2_SM.stan",
                    data = SM_data,
                    pars = params,
                    chains = 3, iter = iter.in, seed = seed, cores = 3) 
    
    t2 <- proc.time()
    
    attr(SM_Cocc, 'time') <- (t2 - t1)[3]
    
    SM_Cocc
    
  } 
#--------------------------------------

end.time = Sys.time()
time.taken = end.time - start.time
print(round(time.taken, 2))

stopCluster(cl1)  # close the clusters

length(out)

length(out[[1]])

all.out = do.call('c', out)
length(all.out)

tmp = run.times(all.out)

all.stan.m = all.out
rm(list=setdiff(ls(), c("all.stan.m")))

save.image("U:/Desktop/Fish_Git/Marginalized/Application_3/working_runs/Combo2_SM_runs_1.RData")
#-----------------------------------------------------------------------------#
# end
###############################################################################
#                                                                        Dec 18
#  Combination Occupancy and N-Mixture Model
#  Discrete JAGS version 
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
sink("Combo_JD.jags")
cat("
    model {
	for(i in 1:nSites) {N[i,1] ~ dpois(lambda)}	
	for(t in 1:(nYears-1)) {
		N_mean[t] <- mean(N[,t]) - 1
		log(gamma[t]) <- g0 + g1*N_mean[t] + g2*N_mean[t]*N_mean[t]
		for(i in 1:nSites) {
			logit(omega[i,t]) <- b0 + b1*X[i,t]
			S[i,t] ~ dbin(omega[i,t], N[i,t])
			G[i,t] ~ dpois(gamma[t])
			N[i,(t+1)] <- S[i,t] + G[i,t] 
			}}
	for(i in 1:nSites) {
		for (t in 1:20){
			p_site[i,t] <- 1-pow((1-p_occ),N[i,t]) 
			num_det[i,t] ~ dbin(p_site[i,t],num_trial[i,t])
			}
		for (t in 1:2){
			for (k in 1:3){
				cloglog(tpc[i,(3*(t-1)+k)]) <- a0 + OV[i,(3*(t-1)+k)]
				pc[i,(3*(t-1)+k)]<-tpc[i,(3*(t-1)+k)]*visited[i,(3*(t-1)+k)]
				Y[i,(3*(t-1)+k)]~dbin(pc[i,(3*(t-1)+k)],N[i,(t+20)])
				}}}
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
JD_data<-list(Y=Y,visited=visited,X=X,OV=OV,nYears=nYears,nSites=nSites,num_det=num_det,num_trial=num_trial)
JD_inits<-function() list(N=matrix(c(rep(20,nSites),rep(NA,nSites*(nYears-1))),nrow=nSites,ncol=nYears),
                          S=matrix(10,nrow=nSites,ncol=nYears-1),
                          G=matrix(10,nrow=nSites,ncol=nYears-1))
# JD_combo<-jags(JD_data, JD_inits, params, model.file="combo_JD.txt",n.iter=100000)
JD_combo<-jags(JD_data, JD_inits, params, model.file="combo_JD.txt",n.iter=10)


#-----------------------------------------------------------------------------#
library(foreach)
library(doParallel)

n.core = 4  # really n.core * 3

cl1 = makeCluster(n.core) # number of cores you want to use
registerDoParallel(cl1)

# make sure each cluster has the packages used 
cllibs <- clusterEvalQ(cl1, c(library(R2jags)))

all.t1 = proc.time()
n.runs = 10

# my.n.iter = c(10000, 20000, 30000)
my.n.iter = c(20000) # ~2h

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
    
    JD_combo <- jags.parallel(JD_data, inits = JD_inits,
                            params, 'Combo_JD.jags',
                            n.chains = 3, n.iter = iter.in,
                            export_obj_names = c("iter.in", "seed", "JD_inits"),
                            jags.seed = seed, envir = my.env)
    
    t2 <- proc.time()
    
    attr(JD_combo, 'time') <- (t2 - t1)[3]
    
    JD_combo
    
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
# 
all.jags.d.2 = all.out
# 
rm(list=setdiff(ls(), "all.jags.d.2"))
# 
save.image("U:/Desktop/Fish_Git/Marginalized/Application_3/working_runs/Combo2_JD_runs_2.RData")
#-----------------------------------------------------------------------------#
# end
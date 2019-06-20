###############################################################################
#                                                                     Spring 19
#  N-Occupancy Model
#  Discrete JAGS version 
#
#  Notes:
#   *  Need to set directory for data
#   *  Need to supply JAGS/WinBUGS/Stan settings
#
###############################################################################
library(R2jags)

#-----------------------------------------------------------------------------#
# import, format data for model fitting
nocc.dat<-read.csv(".//nocc.dat.csv")

# extract all BO data

nYears<-26
nVisits<-8
nSites<-158
Ninf<-10
Ns<-c(0:Ninf)

BO<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
day<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
night<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
RF<-array(0,dim=c(max(nocc.dat$Site),nYears-1))


#BO : Matrix of barred owl detections:
	# 1 = barred owl detected
	# 0 = barred owl not detected
	# NA = site not visited
	
for(i in 1:length(nocc.dat[,1])){
	BO[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-nocc.dat$BO[i]
   if(is.na(nocc.dat$TOD[i])==FALSE & nocc.dat$TOD[i]==1){
	day[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-1}
   if(is.na(nocc.dat$TOD[i])==FALSE & nocc.dat$TOD[i]==3){
	night[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-1}}

day[which(is.na(BO))]<-NA
night[which(is.na(BO))]<-NA
TOD<-2-day+night

#Get riparian forest (RF) covariate:
RF<-matrix(NA, nrow=dim(BO)[1], ncol=nYears-1)

for(i in 1:dim(RF)[1]){
	for(j in 1:dim(RF)[2]){
			RF[i,j]=unique(nocc.dat$RF[which(nocc.dat$Year==j & nocc.dat$Site==i)])
					}}

Y<-numeric()
tod<-numeric()
site<-numeric()
year<-numeric()
for (i in 1:nSites){
		temp<-which(is.na(BO[i,])==FALSE)
		site<-c(site,rep(i,length(temp)))
		Y<-c(Y,as.numeric(BO[i,temp]))
		tod<-c(tod,as.numeric(TOD[i,temp]))
		year<-c(year,floor((temp+7)/nVisits))
		}

#-----------------------------------------------------------------------------#
sink("Nocc_JD.jags")
cat("
model { 
	for (i in 1:nSites) {N[i,1] ~ dpois(lambda)} #latent counts in first year
    for (t in 1:(nYears-1)){
		N_mean[t] <- mean(N[,t])  -1				# mean abundance in prior year
		log(gamma[t]) <- a0 + a1*N_mean[t] 		# expected gains as a function of mean population size
		for (i in 1:nSites){
			logit(omega[i,t]) <- b0 + b1*RF[i,t] # survial probability as a function of riparian forest
			S[i,t] ~ dbin(omega[i,t], N[i,t])    # simulate number that survive
			G[i,t] ~ dpois(gamma[t])             # simulate additions to site population
			N[i,(t+1)] <- S[i,t] + G[i,t]    # add survivors to newcomers
			}}
	for (k in 1:nSamples){
		occ_p[k] <- 1-pow((1-p[tod[k]]),N[site[k],year[k]]) # probability of not detecting any individuals
		Y[k] ~ dbern(occ_p[k])                             
		}
	#Priors
	lambda ~ dunif(0,2) 
	a0 ~ dunif(-4,4) 
	a1 ~ dunif(-4,4) 
	b0 ~ dunif(-4,4) 
	b1 ~ dunif(-4,4) 
	for (k in 1:3) {p[k] ~ dunif(0,1)} # different detection probabilities by tod
	}

    ", fill = TRUE)
sink()  

#-----------------------------------------------------------------------------#


JD_data<-list(Y=Y,nYears=nYears,nSites=nSites,nSamples=length(Y),RF=RF,tod=tod,site=site,year=year)

JD_params <- c("a0", "a1", "b0", "b1", "p", "lambda","N_mean")


JD_inits<-function() list(N=matrix(c(rep(5,nSites),rep(NA,nSites*(nYears-1))),nrow=nSites,ncol=nYears),
                       S=matrix(1,nrow=nSites,ncol=(nYears-1)),
                       G=matrix(1,nrow=nSites,ncol=(nYears-1))
                       )

JD_nocc<-jags.parallel(JD_data, JD_inits, JD_params, model.file="Nocc_JD.txt",n.iter=20000,n.chain=3)
#-----------------------------------------------------------------------------#
###############################################################################
#                                                                     Spring 19
#  N-Occupancy Model
#  Marginalized Stan version 
#
#  Notes:
#   
###############################################################################
library(rstan)

rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

#-----------------------------------------------------------------------------#
# import, format data for model fitting
nocc.dat<-read.csv(".//nocc.dat.csv")

# extract all BO data

nYears<-26
nVisits<-8
nSites<-158
Ninf<-10
Ns<-c(0:Ninf)

BO<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
day<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
night<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
RF<-array(0,dim=c(max(nocc.dat$Site),nYears-1))


#BO : Matrix of barred owl detections:
	# 1 = barred owl detected
	# 0 = barred owl not detected
	# NA = site not visited
	
for(i in 1:length(nocc.dat[,1])){
	BO[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-nocc.dat$BO[i]
   if(is.na(nocc.dat$TOD[i])==FALSE & nocc.dat$TOD[i]==1){
	day[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-1}
   if(is.na(nocc.dat$TOD[i])==FALSE & nocc.dat$TOD[i]==3){
	night[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-1}}

day[which(is.na(BO))]<-NA
night[which(is.na(BO))]<-NA
TOD<-2-day+night

#Get riparian forest (RF) covariate:
RF<-matrix(NA, nrow=dim(BO)[1], ncol=nYears-1)

for(i in 1:dim(RF)[1]){
	for(j in 1:dim(RF)[2]){
			RF[i,j]=unique(nocc.dat$RF[which(nocc.dat$Year==j & nocc.dat$Site==i)])
					}}



Y<-numeric()
tod<-numeric()
site<-numeric()
year<-numeric()
for (i in 1:nSites){
		temp<-which(is.na(BO[i,])==FALSE)
		site<-c(site,rep(i,length(temp)))
		Y<-c(Y,as.numeric(BO[i,temp]))
		tod<-c(tod,as.numeric(TOD[i,temp]))
		year<-c(year,floor((temp+7)/nVisits))
		}





nD0<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was not detected during day samples
nD1<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was detected during day samples
nC0<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was not detected during crepuscular samples
nC1<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was detected during crepuscular samples
nN0<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was not detected during night samples
nN1<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was detected during night samples
for (i in 1:nSites){
  for (t in 1:nYears){
    nD0[i,t]<-length(which(site==i&year==t&tod==1&Y==0))
    nD1[i,t]<-length(which(site==i&year==t&tod==1&Y==1))
    nC0[i,t]<-length(which(site==i&year==t&tod==2&Y==0))
    nC1[i,t]<-length(which(site==i&year==t&tod==2&Y==1))
    nN0[i,t]<-length(which(site==i&year==t&tod==3&Y==0))
    nN1[i,t]<-length(which(site==i&year==t&tod==3&Y==1))
  }}

#-----------------------------------------------------------------------------#

sink("Nocc_SM.stan")
cat("
data {
	int<lower=1> nYears;
	int<lower=1> nSites;
	int<lower=1> Ninf;
	row_vector [(Ninf+1)] Ns;
	matrix [nSites,(nYears-1)] RF;
	matrix [nSites,nYears] nD0;
	matrix [nSites,nYears] nD1;
	matrix [nSites,nYears] nC0;
	matrix [nSites,nYears] nC1;
	matrix [nSites,nYears] nN0;
	matrix [nSites,nYears] nN1;	
	}

parameters {
	real<lower=0,upper=1> lambda;
	real<lower=0,upper=1> p [3];
	real<lower=-4,upper=4> a0;
	real<lower=-4,upper=4> a1;
	real<lower=-4,upper=4> b0;
	real<lower=-4,upper=4> b1;
	}

transformed parameters {
	real po [3,2,(Ninf+1)];
	row_vector [(Ninf+1)] pN;
	row_vector [(Ninf+1)] cp [nSites,nYears];
	row_vector [(Ninf+1)] pz [nSites,nYears];
	row_vector [(Ninf+1)] spz [nSites,nYears];
	real N_mn [nSites,nYears];
	real N_mean [(nYears-1)];
	real gamma [(nYears-1)];
	matrix [(Ninf+1),(Ninf+1)] gtr [(nYears-1)];
	matrix [(Ninf+1),(Ninf+1)] str [(nYears-1),nSites];
	matrix [nSites,(nYears-1)] omega ;
	
	omega = inv_logit(b0 + b1*RF);
	for (j in 1:(Ninf+1)){
		po[1,1,j]=(1-p[1])^(j-1);
		po[1,2,j]=1-po[1,1,j];
		po[2,1,j]=(1-p[2])^(j-1);
		po[2,2,j]=1-po[2,1,j];
		po[3,1,j]=(1-p[3])^(j-1);
		po[3,2,j]=1-po[3,1,j];
		pN[j]=exp(poisson_lpmf((j-1) | lambda));
		}
	for (i in 1:nSites){
		for (t in 1:nYears){
			for (j in 1:(Ninf+1)){
				cp[i,t,j]=((po[1,1,j])^nD0[i,t])*((po[1,2,j])^nD1[i,t])*((po[2,1,j])^nC0[i,t])*
					((po[2,2,j])^nC1[i,t])*((po[3,1,j])^nN0[i,t])*((po[3,2,j])^nN1[i,t]);
				}}
		pz[i,1,]=pN .* cp[i,1,];
		spz[i,1,]=pz[i,1,]/sum(pz[i,1,]);
		N_mn[i,1]=sum(Ns .* spz[i,1,]);              
		}
	for(t in 1:(nYears-1)) {
		N_mean[t] = mean(N_mn[,t]) - 1;
		gamma[t] = exp(a0 + a1*N_mean[t]);
		for (j in 2:(Ninf+1)){
			for (k in 1:(j-1)){
				gtr[t,j,k]=0;
				}}
		for (k in 1:(Ninf)) gtr[t,1,k]= exp(poisson_lpmf((k-1) | gamma[t]));
		for (j in 2:(Ninf)){
			for (k in j:(Ninf)) gtr[t,j,k]= exp(poisson_lpmf((k-j) | gamma[t]));
			}
		for (j in 1:(Ninf+1)) gtr[t,j,(Ninf+1)]=1-sum(gtr[t,j,1:Ninf]);
		for(i in 1:nSites) {
			for (j in 1:(Ninf)){
				for (k in (j+1):(Ninf+1)){
					str[t,i,j,k]=0;
					}}
			str[t,i,1,1]=1;
			for (j in 2:(Ninf+1)){
				for (k in 1:j){
					str[t,i,j,k] = exp(binomial_lpmf ((k-1) | (j-1), omega[i,t]));
						}}
			pz[i,(t+1),]=((spz[i,t,]*str[t,i,,])*gtr[t,,]) .* cp[i,t,];
			spz[i,(t+1),]=pz[i,(t+1),]/sum(pz[i,(t+1),]);
			N_mn[i,(t+1)]=sum(Ns .* spz[i,(t+1),]);
			}}}

model {
	for (i in 1:nSites) {
		for (t in 1:nYears){
			target += log(sum(pz[i,t,]));
		}}}

    ", fill = TRUE)
sink()
#-----------------------------------------------------------------------------#


SM_data<-list(nYears=nYears,nSites=nSites,RF=RF,nD0=nD0,nD1=nD1,nC0=nC0,nC1=nC1,nN0=nN0,nN1=nN1,Ninf=Ninf,Ns=Ns)

SM_params <- c("a0", "a1", "b0", "b1", "p", "lambda","N_mean")

SM_nocc<-stan("Nocc_SM.stan",data=SM_data,pars=SM_params,iter=10,chains=1) 

###############################################################################
#                                                                     Spring 19
#  N-Occupancy Model
#  Maximum likelihood version 
#
#  Notes: Includes two different functions
#         function 'Nocc_ML' estimates occupancy by first estimating underlying site-specific abundance (comparable to Stan and JAGS examples above)
#         function 'occ_ML' estimates occupancy without estimating abundance
#   
###############################################################################


#-----------------------------------------------------------------------------#
# import, format data for model fitting
nocc.dat<-read.csv(".//nocc.dat.csv")

# extract all BO data

nYears<-26
nVisits<-8
nSites<-158
Ninf<-10
Ns<-c(0:Ninf)

BO<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
day<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
night<-array(0,dim=c(max(nocc.dat$Site),nVisits*nYears))
RF<-array(0,dim=c(max(nocc.dat$Site),nYears-1))


#BO : Matrix of barred owl detections:
# 1 = barred owl detected
# 0 = barred owl not detected
# NA = site not visited

for(i in 1:length(nocc.dat[,1])){
  BO[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-nocc.dat$BO[i]
  if(is.na(nocc.dat$TOD[i])==FALSE & nocc.dat$TOD[i]==1){
    day[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-1}
  if(is.na(nocc.dat$TOD[i])==FALSE & nocc.dat$TOD[i]==3){
    night[nocc.dat$Site[i],nVisits*(nocc.dat$Year[i]-1)+nocc.dat$Visit_num[i]]<-1}}

day[which(is.na(BO))]<-NA
night[which(is.na(BO))]<-NA
TOD<-2-day+night

#Get riparian forest (RF) covariate:
RF<-matrix(NA, nrow=dim(BO)[1], ncol=nYears-1)

for(i in 1:dim(RF)[1]){
  for(j in 1:dim(RF)[2]){
    RF[i,j]=unique(nocc.dat$RF[which(nocc.dat$Year==j & nocc.dat$Site==i)])
  }}



Y<-numeric()
tod<-numeric()
site<-numeric()
year<-numeric()
for (i in 1:nSites){
  temp<-which(is.na(BO[i,])==FALSE)
  site<-c(site,rep(i,length(temp)))
  Y<-c(Y,as.numeric(BO[i,temp]))
  tod<-c(tod,as.numeric(TOD[i,temp]))
  year<-c(year,floor((temp+7)/nVisits))
}




nD0<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was not detected during day samples
nD1<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was detected during day samples
nC0<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was not detected during crepuscular samples
nC1<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was detected during crepuscular samples
nN0<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was not detected during night samples
nN1<-matrix(0,nrow=nSites,ncol=nYears) #matrix describing number of times BO was detected during night samples
for (i in 1:nSites){
  for (t in 1:nYears){
    nD0[i,t]<-length(which(site==i&year==t&tod==1&Y==0))
    nD1[i,t]<-length(which(site==i&year==t&tod==1&Y==1))
    nC0[i,t]<-length(which(site==i&year==t&tod==2&Y==0))
    nC1[i,t]<-length(which(site==i&year==t&tod==2&Y==1))
    nN0[i,t]<-length(which(site==i&year==t&tod==3&Y==0))
    nN1[i,t]<-length(which(site==i&year==t&tod==3&Y==1))
  }}
#-----------------------------------------------------------------------------#

#Maximum likelihood function (N-occupancy):

nocc_ML<-function(par){
  lambda<-par[1]
  a0<-par[2]
  a1<-par[3]
  b0<-par[4]
  b1<-par[5]
  p<-numeric()
  p[1]<-par[6]
  p[2]<-par[7]
  p[3]<-par[8]
  #
  pN<-dpois(0:Ninf,lambda)
  Ns<-0:Ninf
  po<-array(0,dim=c(3,2,(Ninf+1)))
  po[1,1,]<-(1-p[1])^Ns
  po[1,2,]<-1-po[1,1,]
  po[2,1,]<-(1-p[2])^Ns
  po[2,2,]<-1-po[2,1,]
  po[3,1,]<-(1-p[3])^Ns
  po[3,2,]<-1-po[3,1,]
  cp<-array(0,dim=c(nSites,nYears,(Ninf+1)))
  pz<-matrix(NA,nrow=nSites,ncol=(Ninf+1))
  N_mn<-numeric()
  gtr<-matrix(0,nrow=Ninf+1,ncol=Ninf+1)
  str<-matrix(0,nrow=Ninf+1,ncol=Ninf+1)
  llik<-numeric()
  #
  for (t in 1:nYears){
    for (j in 1:(Ninf+1)){
      cp[,t,j]<-((po[1,1,j])^nD0[,t])*((po[1,2,j])^nD1[,t])*((po[2,1,j])^nC0[,t])*
        ((po[2,2,j])^nC1[,t])*((po[3,1,j])^nN0[,t])*((po[3,2,j])^nN1[,t])
    }}
  for (i in 1:nSites){
    pz[i,]=pN*cp[i,1,]
    N_mn[i]<-Ns%*%pz[i,]/sum(pz[i,])
  }
  for(t in 1:(nYears-1)) {
    N_mean <- mean(N_mn)-1
    gamma <- exp(a0 + a1*N_mean)
    for (j in 1:(Ninf+1)){
      gtr[j,] <- dpois((Ns-j+1),gamma)/ppois(Ninf-j+1,gamma)
    }
    for(i in 1:nSites) {
      omega <- plogis(b0 + b1*RF[i,t])
      for (j in 1:(Ninf+1)){
        str[j,] = pz[i,j]*dbinom(Ns,(j-1),omega)
      }
      tr<-colSums(str)%*%gtr
      pz[i,]<-tr*cp[i,t,]
      N_mn[i]<-(Ns%*%pz[i,])/sum(pz[i,])
    }}
  for (i in 1:nSites){
    llik[i]<-log(sum(pz[i,]))
  }		
  -1*sum(llik)
}
#-----------------------------------------------------------------------------#

#Maximum likelihood function (Dynamic occupancy):

occ_ML<-function(par){
  psi<-par[1]
  a0<-par[2]
  a1<-par[3]
  b0<-par[4]
  b1<-par[5]
  p<-numeric()
  p[1]<-par[6]
  p[2]<-par[7]
  p[3]<-par[8]
  #
  po<-array(0,dim=c(3,2,2))
  po[1,1,]<-(1-p[1])^c(0,1)
  po[1,2,]<-1-po[1,1,]
  po[2,1,]<-(1-p[2])^c(0,1)
  po[2,2,]<-1-po[2,1,]
  po[3,1,]<-(1-p[3])^c(0,1)
  po[3,2,]<-1-po[3,1,]
  cp<-array(0,dim=c(nSites,nYears,2))
  pz<-matrix(NA,nrow=nSites,ncol=2)
  psi_mn<-numeric()
  tr<-matrix(0,nrow=2,ncol=2)
  llik<-numeric()
  #
  for (t in 1:nYears){
    for (j in 1:2){
      cp[,t,j]<-((po[1,1,j])^nD0[,t])*((po[1,2,j])^nD1[,t])*((po[2,1,j])^nC0[,t])*
        ((po[2,2,j])^nC1[,t])*((po[3,1,j])^nN0[,t])*((po[3,2,j])^nN1[,t])
    }}
  for (i in 1:nSites){
    pz[i,]=c(1-psi,psi)*cp[i,1,]
    psi_mn[i]<-pz[i,2]/sum(pz[i,])
  }
  for(t in 1:(nYears-1)) {
    psi_mean <- mean(psi_mn)
    gam<-plogis(a0+a1*psi_mean)
    tr[1,]<-c(1-gam,gam)
    for(i in 1:nSites) {
      eps <- plogis(b0 + b1*RF[i,t])
      tr[2,]<-c(eps,1-eps)
      pz[i,]<-(pz[i,]%*%tr)*cp[i,t,]
      psi_mn[i]<-pz[i,2]/sum(pz[i,])
    }}
  for (i in 1:nSites){
    llik[i]<-log(sum(pz[i,]))
  }		
  -1*sum(llik)
}
#-----------------------------------------------------------------------------#

#Run models:

m_Nocc_ML<-optim(rep(0.5,8),nocc_ML,method="BFGS",hessian=TRUE)

m_occ_ML<-optim(rep(0.5,8),occ_ML,method="BFGS",hessian=TRUE)




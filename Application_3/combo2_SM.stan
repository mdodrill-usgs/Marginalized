data {
	int<lower=1> nYears;
	int<lower=1> nSites;
	int<lower=1> Ninf;
	int<lower=0> Y [nSites,6];
	int<lower=0> visited [nSites,6];
	matrix [nSites,6] OV;
	matrix [nSites,(nYears-1)] X;
	int<lower=0> num_det [nSites,(nYears-2)];
	int<lower=0> num_nodet [nSites,(nYears-2)];
	int<lower=0> maxdet;
	int<lower=0> maxnodet;
	int<lower=0,upper=1> iC [nSites,2,(Ninf+1)];
	real<lower=0> lC [nSites,2,(Ninf+1)];
	row_vector [(Ninf+1)] Ns; 
	}

parameters {
	real<lower=0,upper=1> lambda;
	real<lower=0,upper=1> p_occ;
	real<lower=-3,upper=3> b0;
	real<lower=-1,upper=1> b1;
	real<lower=-3,upper=3> a0;
	real<lower=-3,upper=3> g0;
	real<lower=-1,upper=1> g1;
	real<lower=-1,upper=1> g2;
	}

transformed parameters {
	simplex [2] po [(Ninf+1)];
	row_vector<lower=0,upper=1> [(Ninf+1)] sumdet [(maxdet+1),(maxnodet+1)];
	row_vector<lower=0,upper=1> [(Ninf+1)] pN;
	row_vector<lower=0,upper=1> [(Ninf+1)] cp [nSites,nYears];
	matrix<lower=0,upper=1> [nSites,6] pc;
	row_vector<lower=0,upper=1> [(Ninf+1)] pz [nSites,nYears];
	real<lower=0> N_mn [nSites,nYears];
	real N_mean [(nYears-1)];
	real<lower=0> gamma [(nYears-1)];
	matrix<lower=0,upper=1> [(Ninf+1),(Ninf+1)] gtr;
	matrix<lower=0,upper=1> [(Ninf+1),(Ninf+1)] str;
	matrix<lower=0,upper=1> [nSites,(nYears-1)] omega ;
	
	pc=(1-exp(-1*exp(a0 + OV)));
	omega = inv_logit(b0 + b1*X);
	for (j in 1:(Ninf+1)){
		po[j,2]=1-exp((j-1)*log(1-p_occ));
		po[j,1]=1-po[j,2];
		pN[j]=exp(poisson_lpmf((j-1) | lambda));
		}
	for (d in 1:(maxdet+1)){
		for (n in 1:(maxnodet+1)){
			for (j in 1:(Ninf+1)){
				sumdet[d,n,j]=(po[j,2]^(d-1))*(po[j,1]^(n-1));
					}}}				
	for (i in 1:nSites){
		for (t in 1:20){
				cp[i,t,]=sumdet[(num_det[i,t]+1),(num_nodet[i,t]+1),];
						}
		pz[i,1,]=pN .* cp[i,1,];
		N_mn[i,1]=sum(Ns .* pz[i,1,])/sum(pz[i,1,]);              
		for (t in 1:2){
			for (j in 1:(Ninf+1)){
				cp[i,(t+20),j] = iC[i,t,j]*(exp(lC[i,t,j]+
				visited[i,(3*t-2)]*(Y[i,(3*t-2)]*log(pc[i,(3*t-2)])+
				(j-1-Y[i,(3*t-2)])*log(1-pc[i,(3*t-2)]))+
				visited[i,(3*t-1)]*(Y[i,(3*t-1)]*log(pc[i,(3*t-1)])+
				(j-1-Y[i,(3*t-1)])*log(1-pc[i,(3*t-1)]))+
				visited[i,(3*t)]*(Y[i,(3*t)]*log(pc[i,(3*t)])+
				(j-1-Y[i,(3*t)])*log(1-pc[i,(3*t)]))));
				}}}
	for (j in 2:(Ninf+1)){
		for (k in 1:(j-1)){
			gtr[j,k]=0;
			}}
	for (j in 1:(Ninf)){
		for (k in (j+1):(Ninf+1)){
				str[j,k]=0;
				}}
	str[1,1]=1;
	for(t in 1:(nYears-1)) {
		N_mean[t] = mean(N_mn[,t]) - 1;
		gamma[t] = exp(g0 + g1*N_mean[t] + g2*N_mean[t]*N_mean[t]);
		for (k in 1:(Ninf+1)) gtr[1,k]= exp(poisson_lpmf((k-1) | gamma[t])- poisson_lcdf(Ninf | gamma[t]));
		for (j in 2:(Ninf+1)){
			for (k in j:(Ninf+1)) gtr[j,k]= exp(poisson_lpmf((k-j) | gamma[t])- poisson_lcdf((Ninf+1-j) | gamma[t]));
			}
		for(i in 1:nSites) {
			for (j in 2:(Ninf+1)){
				for (k in 1:j){
					str[j,k] = exp(binomial_lpmf ((k-1) | (j-1), omega[i,t]));
						}}
			pz[i,(t+1),]=((pz[i,t,]*str)*gtr) .* cp[i,t,];
			N_mn[i,(t+1)]=sum(Ns .* pz[i,(t+1),])/sum(pz[i,(t+1),]);
			}}}

model {
	for (i in 1:nSites) {
		target += log(sum(pz[i,22,]));
		}}

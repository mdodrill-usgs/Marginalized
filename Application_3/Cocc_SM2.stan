
data {
	int<lower=1> nspp;
	int<lower=1> nsess;
	int<lower=1> nsites;
	int<lower=1> nyears;
	int<lower=1> NuDH; 
	int<lower=1> nhab; 
	int<lower=1> sY2 [NuDH,nyears]; 
	int<lower=1> visit2 [NuDH]; 
	int<lower=1> hab2 [NuDH]; 
	int<lower=1> spp [NuDH];
	int<lower=1> FR [NuDH];
	int<lower=1> lookup[nspp,nsites];
	}

parameters {
	real<lower=0,upper=1> psi_mean;
	real<lower=0,upper=1> p_mean;
	real<lower=0,upper=1> alpha1_mean;
	real<lower=0,upper=5> sigma_alpha1;
	real<lower=0,upper=5> sigma_alpha2;
	real<lower=0,upper=5> sigma_u;
	real<lower=0,upper=5> sigma_v;
	real<lower=-1,upper=1> rho;
	vector<lower=0,upper=1> [nspp] a1;
	vector [nspp] alpha_dev1;
	vector [nspp] alpha_dev0;
	vector [nspp] alpha_dev2 [7];
	real beta_dev0 [nspp];
	}

transformed parameters {
	real beta;
	real alpha;
	real alpha1_mu;
	real sigma_eta;
	real Z [NuDH,nyears];
	real p [nspp];
	simplex [2] z0 [nspp];
	real po [nspp,2,(nsess+1),2];
	real tr [nspp,nhab,2,2];
	real pz [NuDH,nyears,2];

	beta = logit(psi_mean);
	alpha = logit(p_mean);
	alpha1_mu = logit(alpha1_mean);
	sigma_eta=sigma_v*((1-rho^2)^.5);
    for (i in 1:(nspp)) {
		z0[i,2]=a1[i];
		z0[i,1]=1-a1[i];
		po[i,1,1,2]=1;
		po[i,2,1,2]=1;
		po[i,1,1,1]=1;
		p[i] = inv_logit(alpha + rho*sigma_v*alpha_dev0[i]+sigma_eta*beta_dev0[i]);
		for (j in 1:nsess){
			po[i,1,(j+1),1]=0;
			po[i,1,(j+1),2]=0;
			po[i,2,(j+1),2]=0;			
			po[i,2,(j+1),1]=(p[i]^j)*(1-p[i])^(nsess-j);
			}
		po[i,2,1,1]=(1-p[i])^(nsess);
		for (k in 1:nhab) {
			tr[i,k,2,2]=inv_logit(beta+sigma_u*alpha_dev0[i]+alpha1_mu+sigma_alpha1*alpha_dev1[i]+sigma_alpha2*alpha_dev2[k,i]);
			tr[i,k,1,2]=inv_logit(beta+sigma_u*alpha_dev0[i]+sigma_alpha2*alpha_dev2[k,i]);
			tr[i,k,1,1]=1-tr[i,k,1,2];
			tr[i,k,2,1]=1-tr[i,k,2,2];
			}}
	for (i in 1:NuDH){
		pz[i,1,1]=(z0[spp[i],1]*tr[spp[i],hab2[i],1,1]+z0[spp[i],2]*tr[spp[i],hab2[i],2,1])*po[spp[i],1,sY2[i,1],visit2[i]];
			pz[i,1,2]=(z0[spp[i],1]*tr[spp[i],hab2[i],1,2]+z0[spp[i],2]*tr[spp[i],hab2[i],2,2])*po[spp[i],2,sY2[i,1],visit2[i]];
			Z[i,1]=pz[i,1,2]/(pz[i,1,1]+pz[i,1,2]);
			for (t in 1:(nyears-1)){
				pz[i,(t+1),1]=(pz[i,t,1]*tr[spp[i],hab2[i],1,1]+pz[i,t,2]*tr[spp[i],hab2[i],2,1])*po[spp[i],1,sY2[i,(t+1)],1];
				pz[i,(t+1),2]=(pz[i,t,1]*tr[spp[i],hab2[i],1,2]+pz[i,t,2]*tr[spp[i],hab2[i],2,2])*po[spp[i],2,sY2[i,(t+1)],1];
				Z[i,(t+1)]=pz[i,(t+1),2]/(pz[i,(t+1),1]+pz[i,(t+1),2]);
				}}}

model {	
	alpha_dev1 ~ normal(0,1);
	alpha_dev0 ~ normal(0,1);
	for (h in 1:7) alpha_dev2[h]~normal(0,1);
	beta_dev0 ~ normal(0,1);
	
	for (i in 1:NuDH) {
		target += FR[i]*log(sum(pz[i,nyears,]));
			}}

generated quantities{
	int SR [nsites,nyears];
	int temp [nspp];  

	for (k in 1:nsites){  
		for (t in 1:nyears){
			for (i in 1:nspp) temp[i]=bernoulli_rng(Z[lookup[i,k],t]);
				SR[k,t] = sum(temp);
    }}}

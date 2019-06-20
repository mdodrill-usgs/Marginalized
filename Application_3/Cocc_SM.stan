// Dynamic Community Occupancy Models - Sky Islands 
data {
	int<lower=1> nspp;
	int<lower=1> nsess;
	int<lower=1> nsites;
	int<lower=1> nyears;
	int<lower=1> sY [nspp,nsites,nyears];
	int<lower=1> visit [nsites];
	int<lower=1> hab [nsites];
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
	real Z [nspp,nsites,nyears];
	real p;
	real mu_eta;
	real z0 [2];
	real po [2,(nsess+1),2];
	real tr [2,2];
	vector [2] pz [nspp,nsites,nyears];

	beta = logit(psi_mean);
	alpha = logit(p_mean);
	alpha1_mu = logit(alpha1_mean);
	sigma_eta = sigma_v * ((1 - rho ^ 2) ^ .5);
	po[1,1,2] = 1;
	po[2,1,2] = 1;
	po[1,2,2] = 0;
	po[2,2,2] = 0;
	po[1,1,1] = 1;
	
    for(i in 1:(nspp)){
		mu_eta = alpha + rho * sigma_v * alpha_dev0[i];
		p = inv_logit(mu_eta + sigma_eta * beta_dev0[i]);
		z0[2] = a1[i];
		z0[1] = 1 - a1[i];
		for(j in 1:nsess){
			po[1,(j + 1),1] = 0;
			po[2,(j + 1),1]=(p ^ j) * (1 - p) ^ (nsess - j);
			}
			
		po[2,1,1] = (1 - p) ^ (nsess);
		for(k in 1:nsites){
			tr[2,2] = inv_logit(beta + sigma_u * alpha_dev0[i] + alpha1_mean + sigma_alpha1 * alpha_dev1[i] + sigma_alpha2 * alpha_dev2[hab[k],i]);
			tr[1,2] = inv_logit(beta + sigma_u * alpha_dev0[i] + sigma_alpha2 * alpha_dev2[hab[k],i]);
			tr[1,1] = 1 - tr[1,2];
			tr[2,1] = 1 - tr[2,2];
			pz[i,k,1,1] = (z0[1] * tr[1,1] + z0[2] * tr[2,1]) * po[1,sY[i,k,1],visit[k]];
			pz[i,k,1,2] = (z0[1] * tr[1,2] + z0[2] * tr[2,2]) * po[2,sY[i,k,1],visit[k]];
			Z[i,k,1] = pz[i,k,1,2] / (pz[i,k,1,1] + pz[i,k,1,2]);
			for(t in 1:(nyears - 1)){
				pz[i,k,(t + 1),1] = (pz[i,k,t,1] * tr[1,1] + pz[i,k,t,2] * tr[2,1]) * po[1,sY[i,k,(t + 1)],1];
				pz[i,k,(t + 1),2] = (pz[i,k,t,1] * tr[1,2] + pz[i,k,t,2] * tr[2,2]) * po[2,sY[i,k,(t + 1)],1];
				Z[i,k,(t + 1)] = pz[i,k,(t + 1),2] / (pz[i,k,(t + 1),1] + pz[i,k,(t + 1),2]);
				}
				}
				}
				}

model{
	alpha_dev1 ~ normal(0,1);
	alpha_dev0 ~ normal(0,1);
	for(h in 1:7){
	  alpha_dev2[h] ~ normal(0,1);
	}
	beta_dev0 ~ normal(0,1);
	for(i in 1:(nspp)){
		for(k in 1:nsites){
			target += log(sum(pz[i,k,nyears,]));
			}
			}
			}

generated quantities{
	int SR [nsites,nyears];
	int temp [nspp];  

	for(k in 1:nsites){  
		for(t in 1:nyears){
			for(i in 1:nspp){
			  temp[i] = bernoulli_rng(Z[i,k,t]);
			} 
				SR[k,t] = sum(temp);
    }
    }
    }


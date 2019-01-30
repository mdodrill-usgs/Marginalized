// Marginalized Integrated Population Model - Brown Trout Data  
    
data{
 int NAZsamps;                            // Number of samples AGFD
 int ts[NAZsamps];
 vector [NAZsamps] AZeff;
 int bAZ[NAZsamps, 3];
 int seasNO[23];
 int bNOc[23,3];
 int NOpasses[23];
 int NCH;                       // Number of capture histories
 int FR[NCH];
 int last[NCH];
 int bCH[NCH, 23];
 int sumf[NCH];
 int spawn[23];               // Indicator for spawning month
}

parameters{
  // matrix<lower = -2, upper = 2> [3, 4]lphi;
matrix<lower = 0, upper = 1>[3, 4]bphi;  // make this a parm -- between 0,1
  // real<lower = 0.01, upper = 4> sd_lphi;                  
  real<lower = 0, upper = 1> bpsi2;
  vector<lower = 0, upper = 1>[3] bpsi1;
  real<lower = 0.1, upper = 2> sd_blp;
  vector[4] mu_blp;
  matrix[23, 3]blp_pass;
  vector[3] AZadj;
  vector<lower = 0, upper = 1000>[2] IN;
  real<lower = 0.1, upper = 4> sd_beta;                  
  vector[17] beta_eps;
  real<lower = -6, upper = 0> lbeta_0;                                
  real<lower = 0, upper = 6> mu_I;                                    
  real<lower = 0.01, upper = 3> sd_I;
  vector[68] I;
  matrix[NAZsamps, 3]blpAZ;
}

transformed parameters{
  
  real btrans[4, 4, 4];                                                     
  matrix[23,3]bp_pass;
  real bp[23,4,4];                                                          
  vector[3] mu_AZ;
  matrix[69,3] bN;                                                          
  vector[17] wA;
  vector[17] Beta;
  
  // for(j in 1:3){  // remove for no prior version
  //   for(k in 1:4){
  //     bphi[j,k] = inv_logit(lphi[j,k]);
  //   }
  // }
  
  // define transition matrix that combines survival and growth parameters
  for(i in 1:4){
    btrans[1,3,i] = 0;
    btrans[1,4,i] = 1 - bphi[1,i];
    btrans[2,1,i] = 0;
    btrans[2,2,i] = bphi[2,i] * (1 - bpsi2);
    btrans[2,3,i] = bphi[2,i] * bpsi2;
    btrans[2,4,i] = 1 - bphi[2,i];
    btrans[3,1,i] = 0;
    btrans[3,2,i] = 0;
    btrans[3,3,i] = bphi[3,i];
    btrans[3,4,i] = 1 - bphi[3,i];
    btrans[4,1,i] = 0;
    btrans[4,2,i] = 0;
    btrans[4,3,i] = 0;
    btrans[4,4,i] = 1;
  }
  
  // size class one transitions are done separately because growth is allowed to vary between seasons
  for(i in 1:3){
    btrans[1,1,i] = bphi[1,i] * (1 - bpsi1[i]);
    btrans[1,2,i] = bphi[1,i] * bpsi1[i];
  }
  
  btrans[1,1,4] = 0;
  btrans[1,2,4] = bphi[1,4];
  
  // this loop calculates actual per pass pcaps for each trip and modifies based on # of passes
  for(j in 1:23){
    for(k in 1:3){
      bp_pass[j,k] = inv_logit(blp_pass[j,k]);
      bp[j,k,k] = 1 - pow((1 - bp_pass[j,k]), NOpasses[j]);                 
      bp[j,k,4] = 1 - bp[j,k,k];
    }
    
    bp[j,1,2] = 0;
    bp[j,1,3] = 0;
    bp[j,2,1] = 0;
    bp[j,2,3] = 0;
    bp[j,3,1] = 0;
    bp[j,3,2] = 0;
    bp[j,4,1] = 0;
    bp[j,4,2] = 0;
    bp[j,4,3] = 0;
    bp[j,4,4] = 1;
  }
  
  // calculate offset for each size classes of AZGF effort and calculate expected pcap
  mu_AZ[1] = mu_blp[1] + AZadj[1];
  mu_AZ[2] = mu_blp[2] + AZadj[2];
  mu_AZ[3] = mu_blp[3] + AZadj[3];
  
  bN[1,1] = 0;                                     
  bN[1,2] = IN[1];
  bN[1,3] = IN[2];
  
  //calculate latent abundance of brown trout from fall 2000 to end of 2017
  for(j in 1:17){
    for(k in 1:3){
      bN[((j-1) * 4 + k + 1),1] = btrans[1,1,k] * bN[((j-1) * 4 + k),1];
      bN[((j-1) * 4 + k + 1),2] = btrans[1,2,k] * bN[((j-1) * 4 + k),1] + btrans[2,2,k] * bN[((j-1) * 4 + k),2];
      bN[((j-1) * 4 + k + 1),3] = btrans[2,3,k] * bN[((j-1) * 4 + k),2] + btrans[3,3,k] * bN[((j-1) * 4 + k),3] + exp(I[((j-1) * 4 + k)]);
    }
    
    // BNT recruits produced in fall as a function weighted sum of adults (wA) and reprodutive rate (Beta) in winter
    wA[j] = (bN[((j - 1) * 4 + 2),2] + 4 * bN[((j - 1) * 4 + 2),3]);
    Beta[j] = exp(lbeta_0 + beta_eps[j]);
    
    // between summer and fall all bnt graduate to sz 2 & recruits show up
    bN[(1 + j * 4),1] = wA[j] * Beta[j];
    bN[(1 + j * 4),2] = btrans[1,2,4] * bN[(j * 4),1] + btrans[2,2,4] * bN[(j * 4),2];
    bN[(1 + j * 4),3] = btrans[2,3,4] * bN[(j * 4),2] + btrans[3,3,4] * bN[(j * 4),3] + exp(I[(j * 4)]);
  }
}
  
model{
  real pz[NCH, 23, 4];                                    
  matrix[NAZsamps, 3]bpAZ;
  matrix[NAZsamps, 3]blamAZ;
  matrix[23, 3]blamNO;
  vector[4] temp;
  
  ///////////////////////////
  // lphi is the logit of survival and is given a prior based on the Lorenzen
  // function and the average mass of fish in each size class during each season.
  // variation from the priod mode is determined by an estimated variance
  // parameter (sd_lphi)
  
  // see note above on removing this
  
  // lphi[1,1] ~ normal(1.08, sd_lphi);
  // lphi[1,2] ~ normal(1.14, sd_lphi);
  // lphi[1,3] ~ normal(1.26, sd_lphi);
  // lphi[1,4] ~ normal(1.38, sd_lphi);
  // lphi[2,1] ~ normal(1.96, sd_lphi);
  // lphi[2,2] ~ normal(1.96, sd_lphi);
  // lphi[2,3] ~ normal(2.02, sd_lphi);
  // lphi[2,4] ~ normal(2.02, sd_lphi);
  // lphi[3,1] ~ normal(2.29, sd_lphi);
  // lphi[3,2] ~ normal(2.29, sd_lphi);
  // lphi[3,3] ~ normal(2.29, sd_lphi);
  // lphi[3,4] ~ normal(2.29, sd_lphi);
  // for(i in 1:3){
  //   for(j in 1:4){
  //     bphi[i,j] ~ uniform(0.01, .99);
  //   }
  // }
  
  ///////////////////////////
  
  // mean pcaps per pass on a logit scale for three size classes, plus largest size class during spawning season
  for(i in 1:4){
    mu_blp[i] ~ normal(-3, 2);                                     // take this out of loop? 
  }
  
  // (done above in transformed) this loop calculates actual per pass pcaps for each trip and modifies based on # of passes
   for(j in 1:23){
    blp_pass[j,1] ~ normal(mu_blp[1], sd_blp);
    blp_pass[j,2] ~ normal(mu_blp[2], sd_blp);
    blp_pass[j,3] ~ normal(mu_blp[spawn[j]], sd_blp);
   }
   
  //  for(k in 1:NCH){
  //   pz[k,sumf[k],1] = (1 == bCH[k,sumf[k]]);
  //   pz[k,sumf[k],2] = (2 == bCH[k,sumf[k]]);
  //   pz[k,sumf[k],3] = (3 == bCH[k,sumf[k]]);
  //   pz[k,sumf[k],4] = 0;
  //   
  //   for(j in sumf[k]:(last[k] - 1)){
  //     for(i in 1:4){
  //       // pz[k,(j + 1),i] = inprod(pz[k,j,], btrans[,i,seasNO[(j + 1)]]) * bp[j,i,bCH[k,(j + 1)]]
  //       // pz[k,(j + 1),i] = dot_product(pz[k,j,], btrans[,i,seasNO[(j + 1)]]) * bp[j,i,bCH[k,(j + 1)]];  //write like HBC model, check this
  //     }
  //   }
  //   
  //   target += FR[k] * log(sum(pz[k,last[k],]));                            
  //   
  // }
  
  for(k in 1:NCH){
		pz[k,sumf[k],1] = (1 == bCH[k,sumf[k]]);
		pz[k,sumf[k],2] = (2 == bCH[k,sumf[k]]);
		pz[k,sumf[k],3] = (3 == bCH[k,sumf[k]]);
		pz[k,sumf[k],4] = 0;
		for(t in sumf[k]:(last[k] - 1)){
			for(i in 1:4){
				for(j in 1:4){
					temp[j] = pz[k,t,j] * btrans[j,i,seasNO[(t + 1)]] * bp[t,i,bCH[k,(t + 1)]];
					}
				pz[k,(t + 1),i] = sum(temp);
				}
			}
		target += FR[k] * log(sum(pz[k,last[k],]));                             
		}
  
  //////////////////////////
  // (done above in transformed) calculate offset for each size classes of AZGF effort and calculate expected pcap
  AZadj[1] ~ normal(0,1);                                                  // vectorize
  AZadj[2] ~ normal(0,1);
  AZadj[3] ~ normal(0,1);
  
   
  // calculate actual immigration in each interval on log scale
  for(j in 1:68){                                                  // vectorize here?
    I[j] ~ normal(mu_I, sd_I);
  }
  
  for(j in 1:17){                                                  // vectorize here?
    beta_eps[j] ~ normal(0, sd_beta);
  }
  
  //////////////////////////
  // 2000 - 2017 AZGF data
  for(j in 1:3){
    for(k in 2:3){
      blpAZ[j,k] ~ normal(mu_AZ[k], sd_blp);
      bpAZ[j,k] = inv_logit(blpAZ[j,k]);
      blamAZ[j,k] = bpAZ[j,k] * bN[ts[j],k] * AZeff[j] / 35;
      bAZ[j,k] ~ poisson(blamAZ[j,k]);
    }
  }
 
  
  for(j in 4:NAZsamps){
    for(k in 1:3){
      blpAZ[j,k] ~ normal(mu_AZ[k], sd_blp);
      bpAZ[j,k] = inv_logit(blpAZ[j,k]);
      blamAZ[j,k] = bpAZ[j,k] * bN[ts[j],k] * AZeff[j] / 35;
      bAZ[j,k] ~ poisson(blamAZ[j,k]);
    }
  }
  

  // 2012 - 2017 NO: starts in april 2012
  for(j in 1:23){
    for(k in 1:3){
      blamNO[j,k] = bp[j,k,k] * 0.08 * bN[(j + 46),k];
      bNOc[j,k] ~ poisson(blamNO[j,k]);
    }
  }
}

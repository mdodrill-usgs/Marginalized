data{
 int NAZsamps;                            // Number of samples AGFD
 int ts[NAZsamps];
 vector [NAZsamps] AZeff;
 int bAZ [NAZsamps, 3];
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
  matrix<lower = 0, upper = 1>[3, 4]bphi;  // make this a parm -- between 0,1
  real<lower = 0, upper = 1> bpsi2;
  vector<lower = 0, upper = 1>[3] bpsi1;
  vector<lower = 0, upper = 1>[4] mu_blp;
  vector<lower = 0, upper = 1>[3] mu_AZ;
  vector<lower = 0, upper = 1000>[2] IN;                              
  vector<lower = -4, upper = 8>[68] I;
  vector<lower = 0, upper = 4>[17] Beta;
}

transformed parameters{
  real btrans[4, 4, 4];                                                     
  matrix[23,3]bp_pass;
  real bp[23,4,4];                                                          
  matrix[69,3] bN;                                                          
  vector[17] wA;


for(j in 1:23){
    bp_pass[j,1]=mu_blp[1];
    bp_pass[j,2]=mu_blp[2];
    bp_pass[j,3]=mu_blp[spawn[j]];
   }
     
  
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
    
    // between summer and fall all bnt graduate to sz 2 & recruits show up
    bN[(1 + j * 4),1] = wA[j] * Beta[j];
    bN[(1 + j * 4),2] = btrans[1,2,4] * bN[(j * 4),1] + btrans[2,2,4] * bN[(j * 4),2];
    bN[(1 + j * 4),3] = btrans[2,3,4] * bN[(j * 4),2] + btrans[3,3,4] * bN[(j * 4),3] + exp(I[(j * 4)]);
  }
}
  
model{
  real pz[NCH, 23, 4];                                    
  matrix[NAZsamps, 3]blamAZ;
  matrix[23, 3]blamNO;
  vector[4] temp;
  
  
  // (done above in transformed) this loop calculates actual per pass pcaps for each trip and modifies based on # of passes

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
    
  // calculate actual immigration in each interval on log scale
 
  //////////////////////////
  // 2000 - 2017 AZGF data
  for(j in 1:3){
    for(k in 2:3){
      blamAZ[j,k] = mu_AZ[k] * bN[ts[j],k] * AZeff[j] / 35;
      bAZ[j,k] ~ poisson(blamAZ[j,k]);
    }
  }
 
  
  for(j in 4:NAZsamps){
    for(k in 1:3){
      blamAZ[j,k] = mu_AZ[k] * bN[ts[j],k] * AZeff[j] / 35;
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

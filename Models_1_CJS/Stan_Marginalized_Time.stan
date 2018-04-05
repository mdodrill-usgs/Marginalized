// CJS Model with time varying survival (s) and capture probability (p)
    
data{
  int<lower = 1> NsumCH;
  int<lower = 1> n_occasions;
  int<lower = 1, upper = 2> sumCH[NsumCH, n_occasions];
  int<lower = 1> sumf[NsumCH];
  int<lower = 1> sumFR[NsumCH];
}

parameters{
  real<lower = 0, upper = 1> s[n_occasions - 1];   // 3 month survivals 
  real<lower = 0, upper = 1> p[n_occasions - 1];   // capture probability
}

transformed parameters{
  simplex[2] tr[2,n_occasions - 1];
  simplex[2] pmat[2,n_occasions - 1];
  
  for(k in 1:n_occasions - 1){
    tr[1,k,1] = s[k];
    tr[1,k,2] = (1 - s[k]);
    tr[2,k,1] = 0;
    tr[2,k,2] = 1;
    
    pmat[1,k,1] = p[k];
    pmat[1,k,2] = (1 - p[k]);
    pmat[2,k,1] = 0;
    pmat[2,k,2] = 1;
  }
}
  
model{
    vector[2] pz[n_occasions];
    
    p ~ uniform(0,1);
    
    for(i in 1:NsumCH){  
      pz[sumf[i],1] = 1;
      pz[sumf[i],2] = 0;
      
      for(k in (sumf[i] + 1):n_occasions){ 
        pz[k,1] = pz[k-1,1] * tr[1,k-1,1] * pmat[1,k-1,sumCH[i,(k)]];
        pz[k,2] = (pz[k-1, 1] * tr[1,k-1,2] + pz[k-1, 2]) * pmat[2,k-1,sumCH[i,(k)]];
      }  
      
      target += sumFR[i] * log(sum(pz[n_occasions])); 
    }
  }
    


data {
int<lower=1> NsumCH;
int<lower=1,upper=5> sumCH[NsumCH, 27];
int<lower=1,upper=26> sumf[NsumCH];
int<lower=1,upper=3> season[26];
int<lower=1> sumFR[NsumCH];
int<lower=1> LCRs [22];
int<lower=1> LCRns [4];
int<lower=1> CRs [24];
int<lower=1> CRns [2];
}

parameters {
real<lower=0,upper=1> s_i[4]; // 3 month survivals for small and large chub in LCR and CR
real<lower=0,upper=1> tau; // proportion of CR adults residing in observable location
vector<lower=0,upper=1> [2] g [3]; // seasonal growth for two rivers
vector<lower=0,upper=1> [4] m [3]; // seasonal movement rates for 2 sizes and locations
vector<lower=0,upper=1> [2] p_lcr [22]; // recapture probability in LCR for two size classes
vector<lower=0,upper=1> [2] p_cr [24]; // recapture probability in CR for two size classes
}

transformed parameters {
vector<lower=0,upper=1> [4] s[3];
simplex[7] tr[7,3];
simplex[5] p[7,26];

for (j in 1:4){
s[1,j]=s_i[j]; // 3 month survival for spring to summer interval
s[2,j]=s_i[j]; // 3 month survival for summer to fall interval
s[3,j]=s_i[j]*s_i[j]; // 6 month survival for fall to spring interval
} 
for (i in 1:3){
tr[1,i,1]=s[i,1]*(1-g[i,1])*(1-m[i,1]);
tr[1,i,2]=s[i,1]*g[i,1]*(1-m[i,2]);
tr[1,i,3]=s[i,1]*(1-g[i,1])*m[i,1]*tau;
tr[1,i,4]=s[i,1]*g[i,1]*m[i,2]*tau;
tr[1,i,5]=s[i,1]*(1-g[i,1])*m[i,1]*(1-tau);
tr[1,i,6]=s[i,1]*g[i,1]*m[i,2]*(1-tau);
tr[1,i,7]=1-s[i,1];
tr[2,i,1]=0;
tr[2,i,2]=s[i,2]*(1-m[i,2]);
tr[2,i,3]=0;
tr[2,i,4]=s[i,2]*m[i,2]*tau;
tr[2,i,5]=0;
tr[2,i,6]=s[i,2]*m[i,2]*(1-tau);
tr[2,i,7]=1-s[i,2];
tr[3,i,1]=s[i,3]*(1-g[i,2])*m[i,3];
tr[3,i,2]=s[i,3]*g[i,2]*m[i,4];
tr[3,i,3]=s[i,3]*(1-g[i,2])*(1-m[i,3]);
tr[3,i,4]=s[i,3]*g[i,2]*(1-m[i,4]);
tr[3,i,5]=0;
tr[3,i,6]=0;
tr[3,i,7]=1-s[i,3];
tr[4,i,1]=0;
tr[4,i,2]=s[i,4]*m[i,4];
tr[4,i,3]=0;
tr[4,i,4]=s[i,4]*(1-m[i,4]);
tr[4,i,5]=0;
tr[4,i,6]=0;
tr[4,i,7]=1-s[i,4];
tr[5,i,1]=tr[3,i,1];
tr[5,i,2]=tr[3,i,2];
tr[5,i,3]=0;
tr[5,i,4]=0;
tr[5,i,5]=tr[3,i,3];
tr[5,i,6]=tr[3,i,4];
tr[5,i,7]=tr[3,i,7];
tr[6,i,1]=0;
tr[6,i,2]=tr[4,i,2];
tr[6,i,3]=0;
tr[6,i,4]=0;
tr[6,i,5]=0;
tr[6,i,6]=tr[4,i,4];
tr[6,i,7]=tr[4,i,7];
for (j in 1:6){
tr[7,i,j]=0;
}
tr[7,i,7]=1;
}

for (i in 1:22){
p[1,LCRs[i],1]=p_lcr[i,1];
p[2,LCRs[i],2]=p_lcr[i,2];
}
// No sampling in LCR during 1st, 4th, 7th, 10th recap periods (summers of 2009 - 2012)
for (i in 1:4){
p[1,LCRns[i],1]=0;
p[2,LCRns[i],2]=0;
}
for (i in 1:24){
p[3,CRs[i],3]=p_cr[i,1];
p[4,CRs[i],4]=p_cr[i,2];
}
// No sampling in CR during 3rd and 6th recap periods (spring of 2010 & 2011)
for (i in 1:2){
p[3,CRns[i],3]=0;
p[4,CRns[i],4]=0;
}
for (i in 1:26){
p[1,i,2]=0;
p[1,i,3]=0;
p[1,i,4]=0;
p[1,i,5]=1-p[1,i,1];
p[2,i,1]=0;
p[2,i,3]=0;
p[2,i,4]=0;
p[2,i,5]=1-p[2,i,2];
p[3,i,1]=0;
p[3,i,2]=0;
p[3,i,4]=0;
p[3,i,5]=1-p[3,i,3];
p[4,i,1]=0;
p[4,i,2]=0;
p[4,i,3]=0;
p[4,i,5]=1-p[4,i,4];
p[5,i,1]=0;
p[5,i,2]=0;
p[5,i,3]=0;
p[5,i,4]=0;
p[5,i,5]=1;
p[6,i,1]=0;
p[6,i,2]=0;
p[6,i,3]=0;
p[6,i,4]=0;
p[6,i,5]=1;
p[7,i,1]=0;
p[7,i,2]=0;
p[7,i,3]=0;
p[7,i,4]=0;
p[7,i,5]=1;
}
}

model {
   real temp[7]; 
   vector[7] pz[27]; 

for (k in 1:NsumCH) {  
       for (j in 1:7)
pz[sumf[k], j] = (j == sumCH[k, sumf[k]]);     
       for (t in (sumf[k] + 1):27) { 
for (i in 1:7) { 
for (j in 1:7) 
temp[j] = pz[t - 1, j] * tr[j, season[(t - 1)], i] * p[i, t - 1, sumCH[k, t]]; 
pz[t, i] = sum(temp); 
} 
} 
    target +=sumFR[k]*log(sum(pz[27])); 
}
}


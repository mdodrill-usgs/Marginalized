###############################################################################
#                                                                      March 18
#        Fitting an Integrated Population Model to Brown Trout Data
#        Marginalized Stan version 
#
#  Notes:
#  * The model runs from the fall of 2000 to fall of 2017 on a seasonal basis
#  * We define three size states based on total length in mm 
#    - 0 - 200; 200 - 350; 350 +
#
#  To do: 
#  * 
#
###############################################################################
setwd('C:\\Users\\mdodrill\\Desktop\\Fish_Git\\marginalized_2\\Models_3_IPM')
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

#-----------------------------------------------------------------------------#
# read in data
NO_catch <- read.csv(paste0(getwd(), "/Data/", "NO_catch.csv"), header = TRUE)
AZGF_catch <- read.csv(paste0(getwd(), "/Data/", "AZGF_catch.csv"), header = TRUE)
MR_data <- read.csv(paste0(getwd(), "/Data/", "bCH.csv"), header = FALSE)

# extract/reformat data
bNOc <- as.matrix(NO_catch[,1:3])
NOpasses <- NO_catch$NOpasses
seasNO <- NO_catch$seasNO
spawn <- ifelse(seasNO == 1, 4, 3)
bAZ <- as.matrix(AZGF_catch[,1:3])
ts <- AZGF_catch$ts
AZeff <- AZGF_catch$AZeff
NAZsamps <- length(AZeff)

findlast <- function(x){ifelse(x[24] == 1, 23, max(which(x[1:23] != 4)))}
last <- apply(MR_data, 1, findlast)

bCH <- as.matrix(MR_data[,1:23])
NCH <- length(last)
FR <- rep(1, NCH)
findfirst <- function(x){which(x != 4)[1]}
sumf <- apply(bCH, 1, findfirst)

#-----------------------------------------------------------------------------#
# working with the Stan model in a seperate tab... see "Stan_M...stan"

# sink("Stan_M...stan")
# cat("
#     

# Put the model code here...

#     ",fill=TRUE)
# sink()    


#-----------------------------------------------------------------------------#

# sm.params <- c("s", "p")
# 
# sm.data <- list(NsumCH = NsumCH, n_occasions = n.occasions, sumCH = sumCH,
#                 sumf = sumf, sumFR = sumFR)

sm.data <- list(NAZsamps = NAZsamps, ts = ts, AZeff = AZeff, bAZ = bAZ,
                   seasNO = seasNO, bNOc = bNOc, NOpasses = NOpasses, ones = FR,
                   FR = FR, last = last, bCH = bCH, NCH = NCH, sumf = sumf,
                   spawn = spawn)

# BM_JM.par <- c('beta.I', 'bphi', 'bpsi1', 'bpsi2', 'mu.blp', 'sd.blp', "lbeta.0",
#                "mu.I", "I", "Beta", "IN", "AZadj", "sd.I", "sd.lphi", "sd.blp",
#                'sd.beta', "bN")

sm.params = c('blp_pass')

# MCMC settings
ni = 10
nt = 1
nb = 5
nc = 1


# Call Stan from R 
SM.c <- stan("Stan_Marginalized.stan",
             data = sm.data,
             pars = sm.params,
             control = list(max_treedepth = 14, adapt_delta = .925),
             chains = nc, iter = ni, thin = nt, seed = 1) 

SM.c

#-----------------------------------------------------------------------------#
MST_FPCA_CAR <- function(FPCAout, # Output from function MST_FPCA_univariate_decomposition, including the follows:
                           # dat.s1: data.frame in first dimension in long format with four labeled columns (described below)
                            # DATA.FRAME COLUMNS: 
                            # id: region IDs (vector of length ngrid*nregion) 
                            # gridPoints: follow-up time (vector of length ngrid*nregion) 
                            # Y1: hospitalization rate data (vector of length ngrid*nregion) 
                            # dt: follow-up time points in order (vector of length ngrid*nregion) 
                           # dat.s2: data.frame in second dimension in long format with four labeled columns (as dat.s1)
                           # sig1Est: estimated measurement error variance in first dimension (scalar)
                           # sig2Est: estimated measurement error variance in second dimension (scalar)
                           # mu1Est: estimated mean function in first-dimension(vector of length ngrid)
                           # mu2Est: estimated mean function in second-dimension(vector of length ngrid)
                           # sig1Est: estimated error variance in first-dimension (scalar)
                           # sig2Est: estimated error variance in second-dimension (scalar)
                           # efun11: estimated univariate eigenfunctions in first dimension (matrix of dimension ngrid*M) 
                           # efun22: estimated univaraite eigenfunctions in second dimension (matrix of dimension ngrid*M) 
                           # e.scores11: estimated univariate PC eigenscores in first dimension (matrix of dimension nregion*2*M)
                           # e.scores22: estimated univariate PC eigenscores in second dimension (matrix of dimension nregion*2*M)
                         MCARout, # Output from MST_FPCA_MCAR including the following:
                           # psi11Est: estimated multivariate eigenfunctions in first dimension (matrix of dimension ngrid*M)
                           # psi22Est: estimated multivariate eigenfunctions in second dimension (matrix of dimension ngrid*M)
                         Adj.Mat # Adjacency matrix from the map (0-1 matrix of dimension nregion*nregion)
){
  
#############################################################################
## Description: Function for the second MCMC (Estimation Algorithm steps 4 in Section 2.2) 
##              including estimation of multivariate PC scores.
## Definition:  nregion: number of regions, 
##              ngrid: number of follow-up time points,
##              M: number of eigen components in each dimension.
## Args:        see above
## Returns:     list()
##              xi.mat: all posterior samples of region-specific PC scores (matrix of dimension (nregion*M)*2500)
##              xiEst: estimates of region-specific PC scores (matrix of dimension nregion*(2*M))
##              alphaEst: estimates of spatial variance parameters (vector of length M)
##              sigma2Est: estimate of measurement error variances (scalar)
##              nuEst: estimate of spatial smoothing parameter (scalar)
#############################################################################
  
# Install missing packages
list.of.packages <- c("refund", "fda", "mgcv", "MASS", "caTools", "locpol", 
                        "KernSmooth", "fANCOVA", "mgcv", "mvtnorm", "spdep", 
                        "fields", "R2WinBUGS", "pbugs", "MCMCpack", "mcmcplots")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) 
  
# Load packages
library(refund)
library(fda)
library(mgcv)
library(MASS)
library(caTools)
library(locpol)
library(KernSmooth)
library(fANCOVA)
library(mvtnorm)
library(spdep)
library(fields)
library(Matrix)
library(nimble)
library(tidyverse)
library(fields)
library(R2WinBUGS)
library(pbugs)
library(MCMCpack)
library(mcmcplots)
  
## Define time points, number of regions, time points and eigen components in each dimension
nregion <- data.T$nregion
ngrid <- data.T$ngrid
M <- data.T$M
gridPoints <- data.T$gridPoints

# Input the generated two dimensional data
# 1st dimension
Y1 <- Y1
# 2nd dimension
Y2 <- Y2

###########################################################
## Using MCMC (2nd MCMC) to estimate multivariate 
#  PC scores by WinBUGs
#########################################################

# WinBUGS model for the 2nd MCMC (CAR part)

CARmodel = function(){
  
  # Likelihood
  for (i in 1:nregion){
    for (j in 1:nmonth){
      Y[i, j] ~ dnorm(Mean1[i,j], invSigma2)
      Mean1[i,j] <- mu1[j] + xi1[i]*psi1[j] + xi2[i]*psi2[j] + xi3[i]*psi3[j] + 
                    xi4[i]*psi4[j] + xi5[i]*psi5[j] + xi6[i]*psi6[j]
      Y[i+nregion, j] ~ dnorm(Mean2[i,j], invSigma2)
      Mean2[i,j] <- mu2[j] + xi1[i]*psi1.2[j] + xi2[i]*psi2.2[j] + xi3[i]*psi3.2[j] + 
                    xi4[i]*psi4.2[j] + xi5[i]*psi5.2[j] + xi6[i]*psi6.2[j]
    }
  }
  
  # Impose CAR prior on mutivariate PC scores xi's
  xi1[1:nregion] ~ car.proper(ceros[], Cw[], adj[], num[], Mv[], invAlpha1, nu)
  xi2[1:nregion] ~ car.proper(ceros[], Cw[], adj[], num[], Mv[], invAlpha2, nu)
  xi3[1:nregion] ~ car.proper(ceros[], Cw[], adj[], num[], Mv[], invAlpha3, nu)
  xi4[1:nregion] ~ car.proper(ceros[], Cw[], adj[], num[], Mv[], invAlpha4, nu)
  xi5[1:nregion] ~ car.proper(ceros[], Cw[], adj[], num[], Mv[], invAlpha5, nu)
  xi6[1:nregion] ~ car.proper(ceros[], Cw[], adj[], num[], Mv[], invAlpha6, nu)
  
  # Impose Uniform prior on the spatial smoothing parameter
  nu ~ dunif(nu.inf, nu.sup)
  
  # Define constants
  for(i in 1:nregion){ceros[i] <- 0}
  nu.inf <- min.bound(Cw[], adj[], num[], Mv[])
  nu.sup <- max.bound(Cw[], adj[], num[], Mv[])   
  
  #other super priors
  invAlpha1 ~ dgamma(1, 1) 
  alpha1 <- 1/invAlpha1
  
  invAlpha2 ~ dgamma(1, 1) 
  alpha2 <- 1/invAlpha2
  
  invAlpha3 ~ dgamma(1, 1) 
  alpha3 <- 1/invAlpha3
  
  invAlpha4 ~ dgamma(1, 1) 
  alpha4 <- 1/invAlpha4
  
  invAlpha5 ~ dgamma(1, 1) 
  alpha5 <- 1/invAlpha5
  
  invAlpha6 ~ dgamma(1, 1) 
  alpha6 <- 1/invAlpha6
 
  invSigma2 ~ dgamma(2, sigmaEst)
  sigest <- 1/invSigma2
  
} 


# Collect data needed for WinBUGS
# Number of total follow-up months
nmonth <- ngrid

# Adjacency matrix related elements needed for the built-in `proper.car` function in WinBUGS
W <- Adj.Mat
adj.weights.num <- as.carAdjacency(W)
CM <- as.carCM(adj.weights.num$adj, adj.weights.num$weights, adj.weights.num$num)
Cw <- CM$C
Mv <- CM$M

# Construct D matrix from number of neighbors for each region
D <- rowSums(W) 
D <- diag(D)

# Constructed outcome matrix from elements in two dimensions
Y <- rbind(Y1, Y2)

# Estimated mean function
mu1 <- FPCAout$mu1Est
mu2 <- FPCAout$mu2Est

# Estimated overall measurement error
sigmaEst <- mean(c(FPCAout$sig1Est, FPCAout$sig2Est))

# Estimated multivariate eigenfunctions
# 1st dimension
psi1 <- MCARout$psi11Est[,1]
psi2 <- MCARout$psi11Est[,2]
psi3 <- MCARout$psi11Est[,3]
psi4 <- MCARout$psi11Est[,4]
psi5 <- MCARout$psi11Est[,5]
psi6 <- MCARout$psi11Est[,6]
# 2nd dimension
psi1.2 <- MCARout$psi22Est[,1]
psi2.2 <- MCARout$psi22Est[,2]
psi3.2 <- MCARout$psi22Est[,3]
psi4.2 <- MCARout$psi22Est[,4]
psi5.2 <- MCARout$psi22Est[,5]
psi6.2 <- MCARout$psi22Est[,6]


# WinBUGS data list
bug_data <- list(nmonth = nmonth, nregion = nregion, Y = Y, 
                 adj = adj.weights.num$adj, num = adj.weights.num$num, 
                 Cw = Cw, Mv = Mv, mu1 = mu1, mu2 = mu2, sigmaEst = sigmaEst, 
                 psi1 = psi1, psi2=psi2, psi3=psi3, psi4=psi4, psi5=psi5, 
                 psi6=psi6, 
                 psi1.2 = psi1.2, psi2.2=psi2.2, psi3.2=psi3.2, psi4.2=psi4.2, 
                 psi5.2=psi5.2, psi6.2=psi6.2)

# Initital values
bug_init <- function(){
  list(invAlpha1 = 1, invAlpha2 = 1, invAlpha3 = 1, invAlpha4 = 1, 
       invAlpha5 = 1, invAlpha6 = 1,
       invSigma2 = rnorm(1, 1/sigmaEst, .0001), nu = rnorm(1, .9, .0001), 
       xi1 = rnorm(nregion, 0, .0001), xi2 = rnorm(nregion, 0, .0001), 
       xi3 = rnorm(nregion, 0, .0001), xi4 = rnorm(nregion, 0, .0001), 
       xi5 = rnorm(nregion, 0, .0001), xi6 = rnorm(nregion, 0, .0001))
}

# Paramaters to save
bug_param <- c("alpha1","alpha2", "alpha3", "alpha4", "alpha5", "alpha6",
               "sigest", "nu", "xi1", "xi2", "xi3","xi4", "xi5", "xi6")

# Model fit paralelly using WinBUGS
out = pbugs(data = bug_data, inits = bug_init, par = bug_param, 
            model = CARmodel, n.iter = 3500, n.burnin = 1000, DIC = F, 
            bugs.seed = 1, debug = FALSE, n.thin = 1, 
            bugs.directory="E:/QQ/Winbug/winbugs14_full_patched/WinBUGS14/")

# Collect result from all MCMC runs regarding multivariate PC scores (prepared for inference)
xi.mat <- rbind(t(out$sims.list$xi1),t(out$sims.list$xi2), t(out$sims.list$xi3))

# Obtain MCMC sample means
outsum <- out$summary[,c('mean')]

# Estimate region-specific multivariate PC scores
xi1Est <- rep(0, nregion)
xi2Est <- rep(0, nregion)
xi3Est <- rep(0, nregion)
xi4Est <- rep(0, nregion)
xi5Est <- rep(0, nregion)
xi6Est <- rep(0, nregion)
for (i in 1:nregion){
  xi1Est[i] <- outsum[paste('xi1[',i,']', sep = '')]
  xi2Est[i] <- outsum[paste('xi2[',i,']', sep = '')]
  xi3Est[i] <- outsum[paste('xi3[',i,']', sep = '')]
  xi4Est[i] <- outsum[paste('xi4[',i,']', sep = '')]
  xi5Est[i] <- outsum[paste('xi5[',i,']', sep = '')]
  xi6Est[i] <- outsum[paste('xi6[',i,']', sep = '')]
}
xiEst <- cbind(xi1Est, xi2Est, xi3Est, xi4Est, xi5Est, xi6Est)

# Estimate spatial variance papameter \alpha_\ell's
alphaEst <- rep(0, (2*M))
for (i in 1:(2*M)){
  alphaEst[i] <- outsum[paste('alpha',i,sep = '')]
}

# Estimate measurement error variance 
sigma2Est <- outsum["sigest"]

# Estimat spatial smoothing parameter \nu
nuEst <- outsum["nu"]

## Construct output
out <- list(xi.mat = xi.mat, xiEst = xiEst, alphaEst = alphaEst, 
            sigma2Est = sigma2Est, nuEst = nuEst)
return(out)
}

CARout <- MST_FPCA_CAR(FPCAout, 
                       MCARout, 
                       Adj.Mat)


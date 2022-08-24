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
                           # efun11: estimated univariate eigenfunctions in first dimension (matrix of dimension ngrid*M1) 
                           # efun22: estimated univaraite eigenfunctions in second dimension (matrix of dimension ngrid*M2) 
                           # e.scores11: estimated univariate PC eigenscores in first dimension (matrix of dimension nregion*M1)
                           # e.scores22: estimated univariate PC eigenscores in second dimension (matrix of dimension nregion*M2)
                           # M1: number of eigencomponents chosen by FVE in first diemension (scalar)
                           # M2: number of eigencomponents chosen by FVE in second diemension (scalar)
                         MCARout, # Output from MST_FPCA_MCAR including the following:
                           # psi11Est: estimated multivariate eigenfunctions in first dimension (matrix of dimension ngrid*(M1+M2))
                           # psi22Est: estimated multivariate eigenfunctions in second dimension (matrix of dimension ngrid*(M1+M2))
                         L, # An integer between 1 and (M1 + M2) defining the number of eigencomponents included in the final CAR fitting.
                            # The estimated spatial variance parameters alphaEst are used as a guidance in choosing L in applications.
                            # In our implementation, we choose L = M1 + M2.
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
##              xi.mat: all posterior samples of region-specific PC scores (matrix slices of dimension (7500*(M1+M2)*ngregion))
##              xiEst: estimates of region-specific PC scores (matrix of dimension nregion*(M1+M2))
##              alphaEst: estimates of spatial variance parameters (vector of length (M1+M2))
##              sigma2Est: estimate of measurement error variances (scalar)
##              nuEst: estimate of spatial smoothing parameter (scalar)
##              L: final number of eigencomponents (interger between 1 and M1+M2)
#############################################################################
  
# Install missing packages
list.of.packages <- c("refund", "fda", "mgcv", "MASS", "caTools", "locpol", 
                        "KernSmooth", "fANCOVA", "mgcv", "mvtnorm", "spdep", 
                        "fields", "R2WinBUGS", "pbugs")
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
  
## Define time points, number of regions, time points and eigen components in each dimension
nregion <- data.T$nregion
ngrid <- data.T$ngrid
gridPoints <- data.T$gridPoints

# Input the generated two dimensional data
# 1st dimension
Y1 <- Y1
# 2nd dimension
Y2 <- Y2

# Eigencomponents chosen in each dimension
M1 <- FPCAout$M1
M2 <- FPCAout$M2

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
      Mean1[i,j] <- mu1[j] + inprod2(xi[1:noutcome,i], psi1.mat[j,1:noutcome])
      Y[i+nregion, j] ~ dnorm(Mean2[i,j], invSigma2)
      Mean2[i,j] <- mu2[j] + inprod2(xi[1:noutcome,i], psi2.mat[j,1:noutcome])
    }
  }
  
  # Impose CAR prior on mutivariate PC scores xi's
  for (k in 1:noutcome){
    xi[k, 1:nregion] ~ car.proper(ceros[], Cw[], adj[], num[], Mv[], invAlpha[k], nu)
    invAlpha[k] ~ dgamma(1, 1) 
    alpha[k] <- 1/invAlpha[k]
  }
  
  # Impose Uniform prior on the spatial smoothing parameter
  nu ~ dunif(nu.inf, nu.sup)
  
  # Define constants
  for(i in 1:nregion){ceros[i] <- 0}
  nu.inf <- min.bound(Cw[], adj[], num[], Mv[])
  nu.sup <- max.bound(Cw[], adj[], num[], Mv[])   
  
  # Other super priors
  invSigma2 ~ dgamma(2, sigmaEst)
  sigest <- 1/invSigma2
  
} 


# Collect data needed for WinBUGS
# Number of total follow-up months
nmonth <- ngrid
# Total number of eigencomponents from two dimensions
noutcome <- L

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
psi1.mat <- MCARout$psi11Est
# 2nd dimension
psi2.mat <- MCARout$psi22Est


# WinBUGS data list
bug_data <- list(nmonth = nmonth, nregion = nregion, noutcome = noutcome,
                 Y = Y, adj = adj.weights.num$adj, num = adj.weights.num$num, 
                 Cw = Cw, Mv = Mv, mu1 = mu1, mu2 = mu2, sigmaEst = sigmaEst, 
                 psi1.mat = psi1.mat, psi2.mat = psi2.mat)

# Initital values
bug_init <- function(){
  list(invAlpha = rep(1, noutcome),
       invSigma2 = rnorm(1, 1/sigmaEst, .0001), nu = rnorm(1, .9, .0001), 
       xi = matrix(rnorm(nregion * noutcome, 1, 0.0001), nrow = noutcome))
}

# Paramaters to save
bug_param <- c("alpha", "sigest", "nu", "xi")

# Model fit paralelly using WinBUGS
out = pbugs(data = bug_data, inits = bug_init, par = bug_param, 
            model = CARmodel, n.iter = 3500, n.burnin = 1000, DIC = F, 
            bugs.seed = 1, debug = FALSE, n.thin = 1, 
            bugs.directory="E:/QQ/Winbug/winbugs14_full_patched/WinBUGS14/")

# Obtain MCMC sample means
outsum <- out$summary[,c('mean')]

# Estimate region-specific multivariate PC scores
xiEst <- matrix(rep(NA, nregion*noutcome), nrow = noutcome)
for (i in 1:noutcome){
  for (j in 1:nregion){
    xiEst[i,j] <- outsum[paste('xi[',i,',',j,']', sep = '')]
  }
}
xiEst <- t(xiEst)

# Collect estimated multivariate PC scores from all MCMC samples (used for inference)
xi.mat <- out$sims.list$xi

# Estimate spatial variance papameter \alpha_\ell's
alphaEst <- rep(0, L)
for (i in 1:L){
  alphaEst[i] <- outsum[paste('alpha[',i,']', sep = '')]
}

# Estimate measurement error variance 
sigma2Est <- outsum["sigest"]

# Estimat spatial smoothing parameter \nu
nuEst <- outsum["nu"]


## Construct output
out <- list(xi.mat = xi.mat, xiEst = xiEst, alphaEst = alphaEst, 
            sigma2Est = sigma2Est, nuEst = nuEst, L = L)
return(out)
}




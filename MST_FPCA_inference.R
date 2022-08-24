MST_FPCA_inference <- function(FPCAout, # Output from function MST_FPCA_univariate_decomposition, including the follows:
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
                               CARout, # Output from MST_FPCA_CAR including the following: 
                                 # xi.mat: all posterior samples of region-specific PC scores (matrix slices of dimension (7500*(M1+M2)*ngregion))
                                 # xiEst: estimates of region-specific PC scores (matrix of dimension nregion*(M1+M2))
                                 # alphaEst: estimates of spatial variance parameters (vector of length (M1+M2))
                                 # sigma2Est: estimate of measurement error variances (scalar)
                                 # nuEst: estimate of spatial smoothing parameter (scalar)
                                 # L: final number of eigencomponents (interger between 1 and M1+M2)
                               Adj.Mat # Adjacency matrix from the map (0-1 matrix of dimension nregion*nregion)
){
  
#############################################################################
## Description: Function for obtaining prediction and inference for multivariate trajectories (Estimation Algorithm steps 5 in Section 2.2) 
## Definition:  nregion: number of regions, 
##              ngrid: number of follow-up time points,
##              M: number of eigen components in each dimension.
## Args: see above
## Returns:     list()
##              Y1.hat: Region-specific trajectory predictions in first dimension (matrix of dimension nregion*ngrid)
##              Y2.hat: Region-specific trajectory predictions in second dimension (matrix of dimension nregion*ngrid)
##              R.traj1.SD:  Pointwise standard errors of the region-specific trajectory predictions for first dimension (matrix of dimension nregion*ngrid)
##              R.traj2.SD:  Pointwise standard errors of the region-specific trajectory predictions for second dimension (matrix of dimension nregion*ngrid)
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

# Final number of eigencomponents
L <- CARout$L

## Reconstruct multivariate region-specific trajectories
Y1.hat <- Y2.hat <- Y1

for(i in 1:nregion){
  Y1.hat[i,] <- FPCAout$mu1Est + CARout$xiEst[i,(1:L)] %*% t(MCARout$psi11Est[,(1:L)])
  Y2.hat[i,] <- FPCAout$mu2Est + CARout$xiEst[i,(1:L)] %*% t(MCARout$psi22Est[,(1:L)])
}

## Obtain the standard deviation matrix for constructing 95% confidence intervals
R.traj1.SD <- R.traj2.SD <- matrix(0,nrow = nregion, ncol = ngrid)
for(i in 1:nregion){
  # Pointwise 95% confidence intervals for region-specific trajectory 
  Cov.Est <- cov(CARout$xi.mat[,1:L,i])
  R.traj1.SD[i,] <- sqrt(diag(MCARout$psi11Est[,1:L] %*% Cov.Est %*% t(MCARout$psi11Est[,1:L])))
  R.traj2.SD[i,] <- sqrt(diag(MCARout$psi22Est[,1:L] %*% Cov.Est %*% t(MCARout$psi22Est[,1:L])))
}

## Construct output
out <- list(Y1.hat = Y1.hat, Y2.hat = Y2.hat, 
            R.traj1.SD = R.traj1.SD, R.traj2.SD = R.traj2.SD)
return(out)
}
  

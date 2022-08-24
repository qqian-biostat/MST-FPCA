MST_FPCA_simulation <- function(numRegion=2, # Number of regions (scalar, input 1 if you want 367 regions, input 2 if you want 49 regions)
                                sigma=0.02, # Measurement error variance (scalar)
){
#############################################################################
## Description: Function for simulating one data set under the simulation design described
##              in Section 4.
## Args: see above
## Definition: nregion: number of regions; ngrid: number of time points; M: number of eigen components in each dimension.
## Returns: list()
#           Y1: first dimension outcome (matrix of dimension nregion*ngrid)
#           Y2: second dimension outcome (matrix dimension of nregion*ngrid)
#           Adj.Mat: adjacent matrix (matrix of dimension nregion*nregion)
#           dat.True: list of all true values, which contains:
#             nregion: number of regions (scalar)
#             ngrid: number of follow-up time points (scalar)
#             gridPoints: follow up time points (vector of length ngrid)
#             M: number of PC components in each dimension (scalar)
#             psi1.True: multivariate eigenfunction in the first dimension (matrix of dimension ngrid*M)
#             psi2.True: multivariate eigenfunction in the second dimension (matrix of dimension ngrid*M)
#             nu.True: follow-up time points (vector of length ngrid)
#             W: spatial correlation parameter (scalar)
#             alpha1.True, alpha2.True, alpha3.True: spatial variance components (scalar)
#             sigma.True: measurement error variance (scalar)
#             region.PC: region-specific PC scores (matrix of dimension nregion*M)
#             Y1.True:  outcome in the first dimension without measurement errors (matrix of dimension nregion*ngrid)
#             Y2.True: outcome in the second dimension without measurement errors (matrix of dimension nregion*ngrid)  
 
################################################################################## 

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
  
# Define the grid points used for eigenfunctions
ngrid <- 24 # 2 year follow up
M <- 3 # number of eigen components in each dimension

# Generate multivariate eigenfunctions in each dimension
fbs <- create.fourier.basis(rangeval = c(0, 2), nbasis = 10)
gridPoints <- seq(0, 1, length.out = ngrid)

# Multivariate eigenfunction in the 1st dimension 
psi1.True <- eval.basis(gridPoints, fbs)[,c(4, 5, 8)]

# Multivariate eigenfunction in the 2nd dimension 
gridPoints2 <- seq(1, 2, length.out = ngrid)
psi2.True <- psi22.True <- eval.basis(gridPoints2, fbs)[,c(4, 9, 8)]
for(e in 1:M){ 
  set.seed(123+e)
  psi2.True[,e] <- (-1)^sample(1:2, 1, replace = T) * psi22.True[,e]
}

# Load adjacency matrix from maps
if(numRegion==1){
  # Merged HSA map
  load("E:/QQ/two outcome simulation/fpca_mst/mst-fpca/simulation_standardization/final/AdjMat_367.RData")
  # Adjacency matrix from merged HSA map
  W <- NewMergedData$W
} else {
  # US states map 
  load("E:/QQ/two outcome simulation/fpca_mst/mst-fpca/simulation_standardization/final/49StateAdjMat.RData")
}

# Number of total regions
nregion <- nrow(W)

# Number of neighbors for each region
D <- rowSums(W) 


# Spatial correlaion parameter
nu.True <- 0.9

# Covariance matrix of region-specific multivariate PC scores
covmat <- solve(diag(D) - nu.True * W)

# Spatial variance component \alpha_l
alpha1.True <- 2.0  
alpha2.True <- 1.2
alpha3.True <- 0.4

# measurement error
sigma.True <- sigma

#######################################################
# Construct Data set
#######################################################
# Region-specific deviation 
# Generate region-specific PC scores from multivariate normal distribution
set.seed(12)
region.PC1 <- mvrnorm(1, rep(0, nregion), covmat * alpha1.True)  
region.PC2 <- mvrnorm(1, rep(0, nregion), covmat * alpha2.True)  
region.PC3 <- mvrnorm(1, rep(0, nregion), covmat * alpha3.True)
# Store the generated PC scores into `region.PC` matrix
region.PC <- array(rep(0, nregion * M), c(nregion, M))
region.PC[,1] <- region.PC1
region.PC[,2] <- region.PC2
region.PC[,3] <- region.PC3


# Generate outcome
# 1st dimension
Y1 <- Y1.2 <- matrix(0.0, nrow = nregion, ncol = ngrid)
# set.seed(243135+1)
eps1e <- matrix(rnorm(nregion * ngrid, 0, sqrt(sigma.True)), nrow = nregion)
for(i in 1:nregion){
  #i=1
  Y1[i,] <- region.PC[i,1] * psi1.True[,1] + region.PC[i,2] * psi1.True[,2] + 
    region.PC[i,3] * psi1.True[,3] + eps1e[i,]
  Y1.2[i,] <- region.PC[i,1] * psi1.True[,1] + region.PC[i,2] * psi1.True[,2] + 
    region.PC[i,3] * psi1.True[,3] # true outcome values
}


# 2nd dimension
Y2 <- Y2.2 <- matrix(0.0, nrow = nregion, ncol = ngrid)
# set.seed(281+1)
eps2e <- matrix(rnorm(nregion * ngrid, 0, sqrt(sigma.True)), nrow = nregion)
for(i in 1:nregion){
  #i=1
  Y2[i,] <- region.PC[i,1] * psi2.True[,1] + region.PC[i,2] * psi2.True[,2] + 
    region.PC[i,3] * psi2.True[,3] + eps2e[i,]
  Y2.2[i,] <- region.PC[i,1] * psi2.True[,1] + region.PC[i,2] * psi2.True[,2] + 
    region.PC[i,3] * psi2.True[,3] # true outcome values
}

# Store all the related true values into `data.True` list
data.True <- list(nregion = nregion, ngrid = ngrid, gridPoints = gridPoints, 
                  M = M, psi1.True = psi1.True, psi2.True = psi2.True, 
                  nu.True = nu.True, alpha1.True = alpha1.True, 
                  alpha2.True = alpha2.True, alpha3.True = alpha3.True, 
                  sigma.True = sigma.True, region.PC = region.PC, 
                  Y1.True = Y1.2, Y2.True = Y2.2)

# Generate output
out <- list(Y1 = Y1,Y2 = Y2, Adj.Mat = W, data.True = data.True)
return(out)

}













MST_FPCA_univariate_decomposition <- function(Y1 = Y1,          # Outcome in the 1st dimension, matrix of dimension nregion * ngrid
                                              Y2 = Y2,          # Outcome in the 2nd dimension, matrix of dimension nregion * ngrid
                                              FVEvalue = 0.99   # A scalar between 0 to 1 to choose number of eigencomponents in each dimension. 
                                                                # to retain enough information at the initial step of the algorithm, it is better 
                                                                # to use a large FVEvalue (close to 1). In our implementation, we use FVEvalue = 0.99. 
){
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

#############################################################################
## Description: Function for FPCA decomposition (Estimation Algorithm steps 1 in Section 2.2) of MST-FPCA, 
##              including estimation of univariate eigenfunctions and eigenvalues in each dimension. 
## Definition:  nregion: number of regions, 
##              ngrid: number of follow-up time points,
##              M: number of eigen components in each dimension,
##              gridPoints: follow-up time points.
## Args:        see above
## Returns:     list()
##              dat.s1: data.frame in the first dimension with columns c("id", "gridPoints", "Y1", "dt"), for further analysis using MST-FPCA
#                 DATA.FRAME COLUMNS (data is stored in long format): 
#                  id: region IDs (vector of length ngrid*nregion) 
#                  gridPoints: follow-up time (vector of length ngrid*nregion) 
#                  Y1: hospitalization rate data (vector of length ngrid*nregion) 
#                  dt: follow-up time points in order (vector of length ngrid*nregion) 
##              dat.s2: data.frame in the second dimension with columns c("id", "gridPoints", "Y2", "dt"), for further analysis using MST-FPCA
#                 DATA.FRAME COLUMNS (data is stored in long format) is the same as dat.s1. 
##              mu1Est: estimated mean function in first-dimension (vector of length ngrid)
##              mu2Est: estimated mean function in second-dimension (vector of length ngrid)
##              sig1Est: estimated error variance in first-dimension (scalar)
##              sig2Est: estimated error variance in second-dimension (scalar)
##              efun11: estimated univariate eigenfunctions in first dimension (matrix of dimension ngrid*M) 
##              efun22: estimated univaraite eigenfunctions in second dimension (matrix of dimension ngrid*M) 
##              e.scores11: estimated univariate PC eigenscores in first dimension (matrix of dimension nregion*M)
##              e.scores22: estimated univariate PC eigenscores in second dimension (matrix of dimension nregion*M)
##              M1: number of eigencomponents chosen by FVE in first diemension (scalar)
##              M2: number of eigencomponents chosen by FVE in second diemension (scalar)
#############################################################################

## Define number of regions, time points and eigen components in each dimension
nregion <- data.T$nregion
ngrid <- data.T$ngrid
gridPoints <- data.T$gridPoints

# Input the generated two dimensional data
# 1st dimension
Y1 <- Y1
# 2nd dimension
Y2 <- Y2

###########################################
## Univariate FPCA in each dimension
###########################################

## FPCA in the 1st dimension 
# Construct long-format dataframe 
dat.s1 <- data.frame(na.omit(cbind(rep(1:nregion, each = ngrid), 
                                   rep(gridPoints, nregion), 
                                   as.vector(t(Y1[,]))))) 
colnames(dat.s1) <- c("id", "gridPoints", "Y1")

# Calculate the mean
rr1 <- smooth.spline(dat.s1$gridPoints, dat.s1$Y1, cv = F)
mu1Est <- rr1$y

# Calculate the raw covariances
# Center the data
x0 <- rep(gridPoints, each = ngrid)
x1 <- rep(gridPoints, ngrid)
dat.s1$dt <- round(dat.s1$gridPoints * (ngrid - 1)) + 1 
X1.c <- dat.s1$Y1 - mu1Est[dat.s1$dt]
Xmat1 <- matrix(X1.c, nrow = ngrid) 

# Calculate the raw covariance
Xcov1 <- tcrossprod(Xmat1) / nregion

# 2-d smoothing of the raw covariance
diag(Xcov1) <- NA
cov.mat1 <- gam(as.vector(Xcov1) ~ te(x0, x1, k = 10, bs = "ps"), method = "REML")
grids2d <- data.frame(x0 = x0, x1 = x1)
covmat.c1 <- matrix(predict(cov.mat1, newdata = grids2d), nrow = ngrid) 
cov1 <- (covmat.c1 + t(covmat.c1)) / 2 #  Symmetrize covariance function

# Estimate measurement errors
Xcov1 <- tcrossprod(Xmat1) / nregion 
diag <- diag(Xcov1)
loess_diag <- suppressWarnings(loess.as(gridPoints, diag, degree = 1,
                                        criterion = "gcv", user.span = NULL, plot = F))  
sig1Est <- mean(loess_diag$fitted - diag(cov1))

# Eigen decomposition
eigen_temp11 <- eigen(cov1, symmetric = TRUE)
eigen_temp11$values <- eigen_temp11$values[which(eigen_temp11$values > 0)]  # Obtain positive eigenvalues
eigen_temp11$vectors <- eigen_temp11$vectors[, 1:length(eigen_temp11$values)] 

# Normalizing
for(e in 1:length(eigen_temp11$values)){ 
  #e=1
  normal.factor <- trapz(gridPoints, eigen_temp11$vectors[, e]^2)
  eigen_temp11$vectors[, e] <- eigen_temp11$vectors[, e] / sqrt(2*normal.factor)
  eigen_temp11$values[e] <- eigen_temp11$values[e] * 2*normal.factor
}

# Number of eigencomponents in the 1st dimension
M1 <- length(which(cumsum(eigen_temp11$values) / sum(eigen_temp11$values) < FVEvalue)) + 1 

# Store the estimated univariate eigenfunctions in the first dimension
efun11 <- eigen_temp11$vectors[,1:M1]

# Estimate univariate PC scores in the first dimension
Yc11  <- t(Xmat1) # Centered data
e.scores11 <- matrix(0, nrow = nregion, ncol = M1)
for(i in 1:nregion){
  e.scores11[i,] <- apply(matrix(rep(Yc11[i,], M1), nrow = ngrid) * efun11[,1:M1],2, mean) 
}

## FPCA in the 2nd dimension 
# Construct long-format dataframe 
dat.s2 <- data.frame(na.omit(cbind(rep(1:nregion, each = ngrid), 
                                   rep(gridPoints, nregion), 
                                   as.vector(t(Y2[,]))))) 
colnames(dat.s2) <- c("id", "gridPoints", "Y2")

# Calculate the mean
rr2 <- smooth.spline(dat.s2$gridPoints, dat.s2$Y2, cv = F)
mu2Est <- rr2$y

# Calculate the raw covariances
# Center the data
x0 <- rep(gridPoints, each = ngrid)
x1 <- rep(gridPoints, ngrid)
dat.s2$dt <- round(dat.s2$gridPoints * (ngrid - 1)) + 1 
X2.c <- dat.s2$Y2 - mu2Est[dat.s2$dt]
Xmat2 <- matrix(X2.c, nrow = ngrid)

# Calculate the raw covariance
Xcov2 <- tcrossprod(Xmat2) / nregion

# 2-d smoothing of the raw covariance
diag(Xcov2) <- NA
cov.mat2 <- gam(as.vector(Xcov2) ~ te(x0, x1, k = 10, bs = "ps"), method = "REML")
grids2d <- data.frame(x0 = x0, x1 = x1)
covmat.c2 <- matrix(predict(cov.mat2, newdata = grids2d), nrow = ngrid) 
cov2 <- (covmat.c2 + t(covmat.c2)) / 2 #  Symmetrize covariance function

# Estimate measurement errors
Xcov2 <- tcrossprod(Xmat2) / nregion  
diag2 <- diag(Xcov2)
loess_diag <- suppressWarnings(loess.as(gridPoints, diag2, degree = 1,
                                        criterion = "gcv", user.span = NULL, plot = F))  
sig2Est <- mean(loess_diag$fitted - diag(cov2))

# Eigen decomposition
eigen_temp22 <- eigen(cov2, symmetric = TRUE)
eigen_temp22$values <- eigen_temp22$values[which(eigen_temp22$values > 0)]  # Obtain positive eigenvalues
eigen_temp22$vectors <- eigen_temp22$vectors[, 1:length(eigen_temp22$values)] 

# Normalizing
for(e in 1:length(eigen_temp22$values)){ 
  normal.factor <- trapz(gridPoints, eigen_temp22$vectors[, e]^2)
  eigen_temp22$vectors[, e] <- eigen_temp22$vectors[, e] / sqrt(2*normal.factor)
  eigen_temp22$values[e] <- eigen_temp22$values[e] * 2*normal.factor
}

# Number of eigencomponents in the 1st dimension
M2 <- length(which(cumsum(eigen_temp22$values) / sum(eigen_temp22$values) < FVEvalue)) + 1 

# Store the estimated univariate eigen functions in the 2nd dimension
efun22 <- eigen_temp22$vectors[,1:M2]

# Estimate univariate PC scores in the 2nd dimension
Yc22       <- t(Xmat2) # Centered data
e.scores22 <- matrix(0, nrow = nregion, ncol = M2)
for(i in 1:nregion){
  e.scores22[i,] <- apply(matrix(rep(Yc22[i,], M2), nrow = ngrid) * efun22[,1:M2],2, mean) 
}


## Construct output
out <- list(dat.s1 = dat.s1, dat.s2 = dat.s2, 
            mu1Est = mu1Est, mu2Est = mu2Est,
            sig1Est = sig1Est, sig2Est = sig2Est,
            efun11 = efun11, efun22 = efun22, 
            e.scores11 = e.scores11, e.scores22 = e.scores22,
            M1 = M1, M2 = M2)
return(out)
}





  



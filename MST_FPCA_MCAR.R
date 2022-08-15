MST_FPCA_MCAR <- function(FPCAout, # Output from function MST_FPCA_univariate_decomposition, including the follows:
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
                           # e.scores11: estimated univariate PC eigenscores in first dimension (matrix of dimension nregion*M)
                           # e.scores22: estimated univariate PC eigenscores in second dimension (matrix of dimension nregion*M)
                          Y1,  # Generated outcome in first dimension (matrix of dimension nregion * ngrid)
                          Y2,  # Generated outcome in second dimension (matrix of dimension nregion * ngrid)
                          Adj.Mat # Adjacency matrix from the map (0-1 matrix of dimension nregion*nregion)
){
#############################################################################
## Description: Function for the first MCMC (Estimation Algorithm steps 2-3 in Section 2.2) 
##              including estimation the between eigencomponent matrix, as well as mutivariate eigenfucntions
## Definition:  nregion: number of regions, 
##              ngrid: number of follow-up time points,
##              M: number of eigen components in each dimension.
## Args:        see above
## Returns:     list()
##              psi11Est: estimated multivariate eigenfunctions in first dimension (matrix of dimension ngrid*2*M)
##              psi22Est: estimated multivariate eigenfunctions in second dimension (matrix of dimension ngrid*2*M)
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
  
  
## Define time points, number of regions, time points, and eigen components in each dimension
nregion <- data.T$nregion
ngrid <- data.T$ngrid
M <- data.T$M
gridPoints <- data.T$gridPoints
  
# Input the generated two dimensional data
# 1st dimension
Y1 <- Y1
# 2nd dimension
Y2 <- Y2

## WinBUGS model for the 1st MCMC (MCAR part)
MCARmodel = function(){
  # Likelihood
  for (i in 1:nregion){
    for (j in 1:noutcome){
      mat.Xi[i,j] ~ dnorm(Xi[i,j], invTao2)
    }
  }
  
  # Definition of the matrix Xi
  for (i in 1:nregion){
    for (j in 1:noutcome){
      Xi[i,j] <- inprod(tPsi[1:j,i], roDis[j, 1:j])
    }
  }
  
  # Psi matrix(spatially dependent random effects)
  # Impose each row of tPsi matrix a CAR prior
  for (j in 1:noutcome){
    tPsi[j,1:nregion] ~ car.proper(ceros[], Cw[], adj[], num[], Mv[], 1, nu)
  }
  
  # Impose Uniform prior for nu
  for(i in 1:nregion){ceros[i] <- 0}
  nu ~ dunif(nu.inf, nu.sup)
  nu.inf <- min.bound(Cw[], adj[], num[], Mv[])
  nu.sup <- max.bound(Cw[], adj[], num[], Mv[])
  
  
  # Prior distributions of correlation between diseases
  #(cell not used the Cholesky decompostion are fixed to 0)
  for (i in 1:(noutcome-1)){
    for (j in (i+1):noutcome){
      roDis[i,j] <- 0
    }
  }
  for (i in 1:noutcome){
      roDis[i,i] ~ dlnorm(-3, 1)
  }
  for (i in 2:noutcome){
    for (j in 1:(i-1)){
      roDis[i,j] ~ dnorm(0, 5)
    }
  }
  
  # Other priors
  invTao2 ~ dgamma(2, 0.02)
  tao2 <- 1/invTao2
} 


# Collect all the data needed for WinBUGS
# Define the score matrix
mat.Xi <- cbind(FPCAout$e.scores11, FPCAout$e.scores22)

# Number of total eigen components from two dimensions
noutcome <- M * 2

# Adjacency matrix related elements needed for the built-in `proper.car` function in WinBUGS
W <- Adj.Mat
adj.weights.num <- as.carAdjacency(W)
CM <- as.carCM(adj.weights.num$adj, adj.weights.num$weights, adj.weights.num$num)
Cw <- CM$C
Mv <- CM$M


# Construct D matrix from number of neighbors for each region
D <- rowSums(W) 
D <- diag(D)


# Data list for MCMC
MCARdata <- list(mat.Xi = mat.Xi, adj = adj.weights.num$adj, 
                 num = adj.weights.num$num, Cw = Cw, Mv = Mv, 
                 nregion = nregion, noutcome = noutcome)

# Define initial values 
MCARinits = function(){
  list(tPsi = matrix(rnorm(nregion * noutcome), nrow = noutcome), 
       nu = runif(1, 0, 1), 
       roDis = matrix(c(runif(1, 0, 2), rep(NA, 5), rnorm(1), 
                        runif(1, 0, 2), rep(NA, 4), rnorm(2), 
                        runif(1, 0, 2), rep(NA, 3), rnorm(3), 
                        runif(1, 0, 2), rep(NA, 2), rnorm(4), 
                        runif(1, 0, 2), rep(NA, 1), rnorm(5), 
                        runif(1, 0, 2)), ncol = noutcome, byrow = T), 
       invTao2 = rnorm(1,50, 0.0001))
}

# Paramaters to save
MCARparam = c("tao2","nu", "roDis")

# Model fit paralelly using WinBUGS
MCARmcmc = pbugs(data = MCARdata, inits = MCARinits, par = MCARparam, model = MCARmodel, 
              n.iter = 20000, n.burnin = 5000, 
              DIC = F, bugs.seed = 1, debug = FALSE, n.thin = 2, 
              bugs.directory="E:/QQ/Winbug/winbugs14_full_patched/WinBUGS14/")

# Collect results from MCMC
sum <- MCARmcmc$summary[,c('mean')]
roDis <- matrix(rep(NA, noutcome*noutcome), ncol = noutcome)
# roDis
for (i in 1:noutcome){
  for (j in 1:noutcome){
    if (i ==j){
      roDis[i,i] <- sum[paste('roDis[',i,',',i,']',sep = '')]
    } else{
      roDis[i,j] <- sum[paste('roDis[',i,',',j,']', sep = '')]
      roDis[j,i] <- 0
    }
  }
}

# Estimate between eigencomponents matrix Sigma_b
sigb <- roDis %*% t(roDis)

# Estimate spatial smoothing parameter nu
MCAR.nuEst <- sum['nu']

# Estimate the error variance
tao2Est <- sum['tao2']

## Using MFPCA algorithm to estimate multivariate eigenfunctions
# Eigen analysis of Sigma_b matrix
z.eg = svd(sigb)

# Estimate multivariate eigenfunctions from their univariate counterparts
psi11 <- matrix(0, nrow = ngrid, ncol = (2*M))
psi22 <- matrix(0, nrow = ngrid, ncol = (2*M))
for(j in 1:(2*M)){
  # j = 1
  psi11[,j] = FPCAout$efun11 %*% z.eg$u[1:M,j] 
  psi22[,j] = FPCAout$efun22 %*% z.eg$u[(M+1):(2*M),j]
}

# Normalizing
for(e in 1:(2*M)){ 
  normal.factor <- trapz(gridPoints, psi11[, e]^2)
  psi11[, e] <- psi11[, e] / sqrt(2*normal.factor)
}

for(e in 1:(2*M)){ 
  normal.factor <- trapz(gridPoints, psi22[, e]^2)
  psi22[, e] <- psi22[, e] / sqrt(2*normal.factor)
}

psi11Est <- psi11
psi22Est <- psi22


# Construct output
out <- list(psi11Est = psi11Est, psi22Est = psi22Est)

return(out)
}


MCARout <- MST_FPCA_MCAR(FPCAout,
                        Y1, 
                        Y2,
                        Adj.Mat)

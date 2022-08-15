## MST_FPCA_tutorial.R
#############################################################################
## Description: A step-by-step implementation of MST-FPCA and the associated  
## procedures described in "Multivariate spatiotemporal functional principal 
## component analysis for modeling hospitalization and mortality rates in the
## U.S. Dialysis Population".  
#############################################################################
## Functions implemented: 
## MST_FPCA_simulation.R, 
## MST_FPCA_univariate_decomposition.R, 
## MST_FPCA_MCAR.R, 
## MST_FPCA_CAR.R, 
## MST_FPCAR_inference.R.
#############################################################################
## Tutorial Outline:
## 1. Simulate two-dimensional outcome data (MST_FPCA_simulation.R)
## 2. Perform MST-FPCA estimation (MST_FPCA_univariate_decomposition.R, MST_FPCA_MCAR.R, MST_FPCA_CAR.R)
## 3. Prediction and inference on multivariate trajectories (MST_FPCA_inference.R)
## 4. Visualization of MST-FPCA results
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

#############################################################################
# 1. Simulate two-dimensional data
#############################################################################

# Simulate one dataset from the simulation design described 
data.G <- MST_FPCA_simulation(numRegion = 2, sigma = .02)  # MST_FPCA_simulation.R

# Data frame used for estimation and inference
Y1 <- data.G[[1]]
Y2 <- data.G[[2]]

# Adjacency matrix from the map
Adj.Mat <- data.G[[3]]

# Save the underlying true values when generating the data
data.T <- data.G[[4]]

#############################################################################
# 2. Perform MST-FPCA estimation
#############################################################################

# NOTE: Performing MST-FPCA estimation steps 1 (univariate FPCA decomposition) with 49 regions will take less than one minute.

FPCAout <- MST_FPCA_univariate_decomposition(Y1 = Y1, Y2 = Y2)  # MST_FPCA_univariate_decompostion.R

# NOTE: Performing MST-FPCA estimation steps 2-3 (1st MCMC and estimate multivariate eigenfunctions) with 49 regions will take less than one minute.

MCARout <- MST_FPCA_MCAR(FPCAout, Y1, Y2, Adj.Mat) # MST_FPCA_MCAR.R

# NOTE: Performing MST-FPCA estimation steps 4 (2nd MCMC to estimate multivariate PC scores) with 49 regions will take less than one minute.

CARout <- MST_FPCA_CAR(FPCAout, MCARout, Adj.Mat) # MST_FPCA_CAR.R

#############################################################################
# 3. Prediction and inference on two outcome trajectories
#############################################################################

PREDout <- MST_FPCA_inference(FPCAout, MCARout, CARout, Adj.Mat)

#############################################################################
# 4. Visualization of MST-FPCA results
#############################################################################  
# Define the grid points used for the multivariate eigenfunctions
ngrid <- data.T$ngrid # 2 year follow up
gridPoints <- data.T$gridPoints

# Define number of eigen components in each dimension
M <- data.T$M

# True multivariate eigenfunctions
psi1.True <- data.T$psi1.True # 1st dimension
psi2.True <- data.T$psi2.True # 2nd dimension

## Align eigenfunctions
# Note: the estimated eigenfunctions may flip signs. Here we match the sign of the
# estimated eigenfunctions to the true eigenfunctions for plotting  
MSDE.eifun <- function(x,y){
  err1 <- trapz(gridPoints, (x-y)^2)
  err2 <- trapz(gridPoints, (x+y)^2)
  return(c(min(err1, err2),which.min(c(err1, err2))))
}

# Estimated multivariate eigenfunctions
psi1Est <- MCARout$psi11Est[,1:M]
psi2Est <- MCARout$psi22Est[,1:M]

# Mean squred deviation errors
# 1st dimension
MSDE.psi11 <- MSDE.eifun(psi1Est[,1], psi1.True[,1])
MSDE.psi12 <- MSDE.eifun(psi1Est[,2], psi1.True[,2]) 
MSDE.psi13 <- MSDE.eifun(psi1Est[,3], psi1.True[,3])
# 2nd dimension
MSDE.psi21 <- MSDE.eifun(psi2Est[,1], psi2.True[,1])
MSDE.psi22 <- MSDE.eifun(psi2Est[,2], psi2.True[,2]) 
MSDE.psi23 <- MSDE.eifun(psi2Est[,3], psi2.True[,3])

# Align 
# 1st dimension
fun1.sign.1 <- MSDE.psi11[2] * (-2) + 3
psi1Est[,1] <- psi1Est[,1] * fun1.sign.1
fun1.sign.2 <- MSDE.psi12[2] * (-2) + 3
psi1Est[,2] <- psi1Est[,2] * fun1.sign.2
fun1.sign.3 <- MSDE.psi13[2] * (-2) + 3
psi1Est[,3] <- psi1Est[,3] * fun1.sign.3

# 2nd dimension
fun2.sign.1 <- MSDE.psi21[2] * (-2) + 3
psi2Est[,1] <- psi2Est[,1] * fun2.sign.1
fun2.sign.2 <- MSDE.psi22[2] * (-2) + 3
psi2Est[,2] <- psi2Est[,2] * fun2.sign.2
fun2.sign.3 <- MSDE.psi23[2] * (-2) + 3
psi2Est[,3] <- psi2Est[,3] * fun2.sign.3

# Plot estimates of multivariate eigenfunctions
# 1st dimension eigenfunctions
par(mfrow = c(2,2))
matplot(gridPoints, psi1.True, "l", ylim = c(-2,2), xaxs = "i", main = "(a) True: 1st dimension",
        cex.main=2,xlab = "", ylab = "", lwd = 2, col = c("black", "red", "blue"))
legend("bottomleft", legend = c("1", "2", "3"), col = c("black", "red", "blue"), 
       lty = 1:3, cex = 1)
title(xlab = "time",ylab=expression(paste(psi^(1),(t))), line=2, cex.lab=1.6)

matplot(gridPoints, psi1Est, "l", ylim = c(-2,2), 
        xaxs = "i", main = "(b) Estimated: 1st dimension",cex.main=2,xlab = "", ylab = "", lwd = 2,
        col = c("black", "red", "blue"))
title(xlab = "time",ylab=expression(paste(widehat(psi)^(1),(t))), line=2, cex.lab=1.6)

# second dimension eigenfunctions
matplot(gridPoints, psi2.True, "l", ylim = c(-2,2), xaxs = "i", main = "(c) True: 2nd dimension",
        cex.main=2,xlab = "", ylab = "", lwd = 2, col = c("black", "red", "blue"))
legend("bottomleft", legend = c("1", "2", "3"), col = c("black", "red", "blue"), 
       lty = 1:3, cex = 1)
title(xlab = "time",ylab=expression(paste(psi^(2),(t))), line=2, cex.lab=1.6)

matplot(gridPoints, psi2Est, "l", ylim = c(-2,2), 
        xaxs = "i", main = "(d) Estimated: 2nd dimension",cex.main=2,xlab = "", ylab = "", lwd = 2,
        col = c("black", "red", "blue"))
title(xlab = "time",ylab=expression(paste(widehat(psi)^(2),(t))), line=2, cex.lab=1.6)



# Region-specific predicted trajectory for the first region
# Region id for plots
R.ID <- 1
# True trajectory for the 1st region
R.traj1.T <- data.T$Y1.True[R.ID,] #1st dimension
R.traj2.T <- data.T$Y2.True[R.ID,] #2nd dimension
# Predicted trajectory for the 1st region
R.traj1.Est <- PREDout$Y1.hat[R.ID,] #1st dimension
R.traj2.Est <- PREDout$Y2.hat[R.ID,] #2nd dimension
# Pointwise 95% confidence intervals 
R.traj1.SD <- PREDout$R.traj1.SD[R.ID,] #1st dimension
R.traj1.upper <- R.traj1.Est + 1.96 * R.traj1.SD
R.traj1.lower <- R.traj1.Est - 1.96 * R.traj1.SD
R.traj2.SD <- PREDout$R.traj2.SD[R.ID,] #2nd dimension
R.traj2.upper <- R.traj2.Est + 1.96 * R.traj2.SD
R.traj2.lower <- R.traj2.Est - 1.96 * R.traj2.SD
# Visualization of trajectory prediction for the first dimension
par(mfrow = c(1,2))
ylim1 <- min(R.traj1.lower, R.traj1.T, R.traj1.Est)
ylim2 <- max(R.traj1.upper, R.traj1.T, R.traj2.Est)
plot(gridPoints, R.traj1.T, "l", lwd = 2, ylim = c(ylim1,ylim2), xaxs = "i", 
     main = "1st dimension",cex.main = 2, xlab = "", ylab = "", col = "black")
title(xlab = "Time",ylab = "Trajectory", line=2, cex.lab=1.6)
lines(gridPoints, R.traj1.Est, col = "blue", lwd = 2, lty = 1)
lines(gridPoints, R.traj1.upper, col = "blue", lwd = 1, lty = 2)
lines(gridPoints, R.traj1.lower, col = "blue", lwd = 1, lty = 2)

ylim1 <- min(R.traj2.lower, R.traj2.T, R.traj2.Est)
ylim2 <- max(R.traj2.upper, R.traj2.T, R.traj2.Est)
plot(gridPoints, R.traj2.T, "l", lwd = 2, ylim = c(ylim1,ylim2), xaxs = "i",
     main = "2nd dimension",cex.main = 2, xlab = "", ylab = "", col = "black")
title(xlab = "Time",ylab = "Trajectory", line=2, cex.lab=1.6)
lines(gridPoints, R.traj2.Est, col = "blue", lwd = 2, lty = 1)
lines(gridPoints, R.traj2.upper, col = "blue", lwd = 1, lty = 2)
lines(gridPoints, R.traj2.lower, col = "blue", lwd = 1, lty = 2)


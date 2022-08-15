CONTENTS OF THIS FOLDER ——————————————

MST_FPCA_tutorial.R : A step-by-step implementation of MST-FPCA and the associated procedures described in "Multivariate Spatiotemporal Functional Principal Component Analysis for Modeling Hospitalization and Mortality Rates in the Dialysis Population".

MST_FPPCA_simulation.R : Function for simulating one data set under the simulation design described in Section 4 of "Multivariate Spatiotemporal Functional Principal Component Analysis for Modeling Hospitalization and Mortality Rates in the Dialysis Population".

MST_FPCA_univariate_decomposition.R : Function for univariate FPCA decomposition (estimation steps 1 in the estimation algorithm) described in "Multivariate Spatiotemporal Functional Principal Component Analysis for Modeling Hospitalization and Mortality Rates in the Dialysis Population", including estimation of mean function, univariate eigenfunctions and eigenscores.

MST_FPCA_MCAR.R : Function for the 1st MCMC estimation (estimation step 2-3 in the estimation algorithm) of MST-FPCA model described in "Multivariate Spatiotemporal Functional Principal Component Analysis for Modeling Hospitalization and Mortality Rates in the Dialysis Population", including estimation of within eigencomponent variation, multivariate eigenfunctions.

MST_FPCA_CAR.R : Function for the 2nd MCMC estimation (estimation step 4 in the estimation algorithm) of MST-FPCA model described in "Multivariate Spatiotemporal Functional Principal Component Analysis for Modeling Hospitalization and Mortality Rates in the Dialysis Population", including estimation of spatial variance parameters, measurement error variance and multivariate PC scores.

MST_FPCA_inference.R : Function for obtaining prediction and inference for multivariate trajectories described in Section 2.2 of "Multivariate Spatiotemporal Functional Principal Component Analysis for Modeling Hospitalization and Mortality Rates in the Dialysis Population".


INTRODUCTION ——————————————

The contents of this folder allow for implementation of the MST-FPCA estimation and inference described in "Multivariate Spatiotemporal Functional Principal Component Analysis for Modeling Hospitalization and Mortality Rates in the Dialysis Population". Users can simulate a sample data frame (MST_FPCA_simulation.R) and apply the proposed estimation algorithm (MST_FPCA_univariate_decomposition.R, MST_FPCA_MCAR.R, MST_FPCA_CAR.R). Also, we include tools to perform prediction and inference on multvariate (hospitalization and mortality rate) trajectories (MST_FPCA_inference.R), allowing users to obtain region-specific predicted hospitalization and mortality rate trajectories as well as their pointwise confidence intervals. Detailed instructions on how to perform the aforementioned procedures, make predictions of region-specific hospitalization and mortality rate trajectories and visualize results are included in MST_FPCA_tutorial.R.

REQUIREMENTS ——————————————

The included R programs require R 3.5.3 (R Core Team, 2019) and the packages listed in MST_FPCA_tutorial.R, as well as preinstalled WinBUGS software.

INSTALLATION ——————————————

Load the R program files into the global environment and install required packages using commands in MST_FPCA_tutorial.R

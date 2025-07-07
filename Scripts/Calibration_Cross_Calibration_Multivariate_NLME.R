# -------------------------------------------------------------------------------------------------------------
# Title: Calibration and Cross-Model Calibration of Nonlinear Multivariate Mixed-Effects Models - R Script
# Author: Albert Ciceu
# Affiliation: Austrian Research Center For Forests (BFW)
# Corresponding Author: Albert Ciceu, Email: albert.ciceu@bfw.gv.at
#
# Description: This script accompanies the article "Nonlinear multilevel seemingly unrelated height-diameter 
# and crown length mixed-effects models for the southern Transylvanian forests, Romania" 
# https://doi.org/10.1016/j.fecs.2025.100322
# It implements calibration and cross-model calibration of multivariate mixed-effects models.
# For a detailed numeric example, refer to Appendix A of the article, and for further information, 
# see the main text.
#
# Date: 17/09/2024
# -------------------------------------------------------------------------------------------------------------

#----------------------------------------------- First Scenario ------------------------------------------------

# We aim to calibrate the BE model using a CU that includes two SP. 
# We will calibrate both models using a sample of four BE trees, with two trees selected from each SP.

# First, we load the required packages, two predefined functions and the model:
# - fun_HT_CL.R: the multivariate function based on Equation 17 from the main text.
# - function_Z.R: builds the Z matrix.
# - SUR_BEECH2.rds: is the model fit for BE species

# The second section defines a function to predict random effects 
# when observations for all model components are available. The function takes two arguments:
# a fitted model and a dataset (plot), and returns the predicted random effects.

getwd()
# Load required packages

require(Matrix)      # Matrix operations and sparse matrices
require(magic)       # Matrix manipulation functions
require(nlme)        # Nonlinear mixed-effects models
require(tidyverse)   # Data manipulation and visualization

# Function used to fit the model (check the Functions folder)
source('Functions/fun_HT_CL.R')
#Import function that builts the Z matrix (check the Functions folder)
source('Functions/function_Z.R')

# Load model
model <- readRDS('Output/SUR_BEECH2.rds')


#data used for calibration
data_plot<- data.frame(
  unique_ID = c("KII13", "KII13", "KII13", "KII13"),
  SSP = c(1, 1, 2, 2),
  tree_ID = c(6, 16, 7, 4),
  dbh = c(8.467, 19.735, 24.987, 39.502),
  HT_m = c(8.700, 9.100, 24.500, 23.900),
  CL = c(3.600, 6.400, 8.600, 9.300),
  QMD = c(22.591, 22.591, 30.524, 30.524),
  ddom = c(40.215, 40.215, 54.125, 54.125),
  G_species = c(21.307, 21.307, 40.980, 40.980),
  G = c(24.051, 24.051, 40.980, 40.980),
  N = c(600, 600, 560, 560),
  BALMOD_c = c(5.457, 4.400, 6.193, 4.623),
  gini = c(0.308, 0.308, 0.356, 0.356)
)

View(data_plot)

# Function to predict random effects for a multivariate model using observations from all components.
prediction_random_effects_sur <- function(model, data_plot) {
  
  # Extract fixed effects from the model
  beta <- fixed.effects(model) 
  
  # Prepare data for the specific plot, removing rows with NA in H or CL
  dfplot<-data_plot[!is.na(data_plot$HT_m)& !is.na(data_plot$CL),]
  nobs <- nrow(dfplot) # Number of observations (trees) in the plot
  
  # Extract correlation matrix of the model
  corepsilon <- corMatrix(model$modelStruct[2]$corStruct)[[1]] # Correlation matrix
  B.power <- coef(model$modelStruct$varStruct$B) # parameters of the power type function
  sdepsilon <- model$sigma * c(1, exp(coef(model$modelStruct$varStruct$A))) # Standard deviations of the two models
  
  # Compute variances
  v1 <- sdepsilon[1] ^ 2 * dfplot$dbh ^ (2 * B.power[1])
  v2 <- sdepsilon[2] ^ 2 * dfplot$dbh ^ (2 * B.power[2])
  
  V1 <- if (nobs == 1) {
    corepsilon[1, 2] * sqrt(v1 * v2)
  } else {
    corepsilon[1, 2] * diag(sqrt(v1 * v2))
  }
  
  # Construct the R matrix
  R <- rbind(cbind(matrix(0, nobs, nobs), V1),
             cbind(V1, matrix(0, nobs, nobs)))
  diag(R) <- c(v1, v2)
  
  # Determine the number of subplots
  nSSP <- length(unique(dfplot$SSP))
  
  # Extract and construct D matrix
  Dplot <- sigma(model)^2 * pdMatrix(model$modelStruct$reStruct)$plot
  Dsubplot <- sigma(model)^2 * pdMatrix(model$modelStruct$reStruct)$subplot
  
  # Create a list of D matrices for each subplot
  Dlist <- replicate(nSSP + 1, Dsubplot, simplify = FALSE)
  Dlist[[1]] <- Dplot
  D <- bdiag(Dlist) # Combine D matrices into a block-diagonal matrix
  
  # Determine parameter lengths
  length.param <- length(fixed.effects(model))
  length.param.rf <- ncol(ranef(model, level = 1))
  
  # Identify random effect parameters
  rf_position <- as.integer(names(fixed.effects(model)) %in% names(ranef(model, level = 1)))
  
  # Create A matrix for fixed effects
  A <- diag(rep(1, length.param))
  A <- do.call(rbind, replicate(nobs, A, simplify = FALSE))
  
  # Create B matrix for random effects
  num_ones <- sum(rf_position)
  num_elements <- length(rf_position)
  Bijk <- matrix(0, nrow = num_elements, ncol = num_ones)
  
  # Populate Bijk matrix
  col_index <- 1
  for (i in 1:num_elements) {
    if (rf_position[i] == 1) {
      Bijk[i, col_index] <- 1
      col_index <- col_index + 1
    }
  }
  
  # Create zero block matrix
  zeroblock <- matrix(0, ncol = length.param.rf, nrow = length.param)
  Bplot <- do.call(rbind, replicate(nobs, Bijk, simplify = FALSE))
  
  # Construct B matrix for plots
  B1 <- rbind(
    do.call(rbind, replicate(length(dfplot$tree_ID[dfplot$SSP == 1]), Bijk, simplify = FALSE)),
    do.call(rbind, replicate(length(dfplot$tree_ID[dfplot$SSP == 2]), zeroblock, simplify = FALSE))
  )
  
  # Handle cases for multiple subplots
  if (nSSP == 1) {
    B <- cbind(Bplot, B1)
  } else {
    B2 <- rbind(
      do.call(rbind, replicate(length(dfplot$tree_ID[dfplot$SSP == 1]), zeroblock, simplify = FALSE)),
      do.call(rbind, replicate(length(dfplot$tree_ID[dfplot$SSP == 2]), Bijk, simplify = FALSE))
    )
    B <- cbind(Bplot, B1, B2)
  }
  
  # Initialize random effects
  b0 <- rep(0, length.param.rf * (1 + nSSP))
  bOld <- b0 - 1
  history <- matrix(ncol = length(b0), nrow = 0)
  
  # Iteratively update random effects
  while (sum((b0 - bOld) ^ 2) > 1e-20) {
    y <- c(dfplot$HT_m, dfplot$CL) # Response vector
    cf <- A %*% beta + B %*% b0 # Compute coefficients
    
    for (i in 1:length.param) {
      assign(paste0('cf', i), cf[seq(i, nrow(cf), by = length.param)])
    }
    
    # Model formulation (Check the main manuscript)
    a1 <- cf1 + cf2 * log(dfplot$QMD) + cf3 * log(dfplot$G_species) + cf4 * dfplot$ddom
    a2 <- cf5 + cf6 * log(dfplot$QMD) + cf7 * sqrt(dfplot$N)
    b1 <- cf8 + cf9 * log(dfplot$QMD) + cf10 * dfplot$gini
    b2 <- cf11 + cf12 * log(dfplot$QMD) + cf13 * dfplot$G_species + cf14 * dfplot$BALMOD_c
    
    # Predict H and CL
    predHT <- fun_HT_CL(dbh = dfplot$dbh, which = 1, a1, a2, b1, b2)
    predCL <- fun_HT_CL(dbh = dfplot$dbh, which = 0, a1, a2, b1, b2)
    pred <- c(predHT, predCL)
    
    # Compute gradients for H and CL predictions
    X_HT <- attributes(predHT)$gradient
    X_CL <- attributes(predCL)$gradient
    
    # Build Z matrix
    Z <- function_Z(X_HT, X_CL, dfplot, nSSP)
    
    # Estimate  random effects
    bOld <- b0
    b0 <- D %*% t(Z) %*% solve(Z %*% D %*% t(Z) + R) %*% (y - pred + Z %*% b0)
    
    # Record history of random effects
    history <- rbind(history, t(b0))
    cat(as.numeric(b0), "\n")
  }
  
  return(b0) # Return final estimates of random effects
}

# Perform b prediction and display results for data_plot. 
(b<-prediction_random_effects_sur(model = model, data_plot = data_plot))
# The first two random effects are plot level random efects for H and the next two for CL
# The 5th and 6th  random effects are for SP1 for H and the next two for CL.
# The 9th and 10th  random effects are for SP2 for H and the next two for CL
round(b,3)



#--------------------------------------------Second scenario------------------------------------------------------
# With the same dataset we cross-calibrate the model using only the H observations
# The following function takes two arguments:
# a fitted model and a dataset (plot), and returns the predicted random effects.


cross_prediction_random_effects_sur <- function(model, data_plot) { 
  
  # Extract the H model fixed-effects parameters
  beta <- fixed.effects(model)[str_detect(names(fixed.effects(model)), pattern = 'a1.|a2')]
  
  # Number of observations and unique subplots
  nobs <- nrow(data_plot)
  nSSP <- length(unique(data_plot$SSP))
  
  # Extract correlation and variance matrices from the model
  corepsilon <- corMatrix(model$modelStruct[2]$corStruct)[[1]]
  B_power <- coef(model$modelStruct$varStruct$B)
  sdepsilon <- model$sigma * c(1, exp(coef(model$modelStruct$varStruct$A)))
  
  # Compute variance components for the H model
  v1 <- sdepsilon[1]^2 * data_plot$dbh^(2 * B_power[1])
  R0 <- if (nobs == 1) v1 else diag(v1)
  
  # Construct D and C matrices for random effects
  Dplot <- sigma(model)^2 * pdMatrix(model$modelStruct$reStruct)$plot
  Dsubplot <- sigma(model)^2 * pdMatrix(model$modelStruct$reStruct)$subplot
  
  # Extract and combine diagonal blocks for D and C matrices
  D0plot <- Dplot[1:2, 1:2]
  D0subplot <- Dsubplot[1:2, 1:2]
  
  Dlist <- c(list(D0plot), rep(list(D0subplot), nSSP))
  D0 <- as.matrix(bdiag(Dlist))
  
  Clist <- c(list(Dplot[, 1:2]), rep(list(Dsubplot[, 1:2]), nSSP))
  C <- as.matrix(bdiag(Clist))
  
  # Determine the number of fixed and random parameters
  model_fe <- names(beta)
  model_re <- names(ranef(model, level = 1))[str_detect(names(ranef(model, level = 1)), pattern = 'a1.|a2')]
  
  length_param <- length(model_fe)
  length_param_rf <- length(model_re)
  rf_position <- ifelse(model_fe %in% model_re, 1, 0)
  
  # Construct A matrix
  A <- do.call(rbind, replicate(nobs, diag(rep(1, length_param)), simplify = FALSE))
  
  # Construct B matrix based on random effect positions
  num_ones <- sum(rf_position)
  Bijk <- matrix(0, nrow = length_param, ncol = num_ones)
  col_index <- 1
  for (i in 1:length_param) {
    if (rf_position[i] == 1) {
      Bijk[i, col_index] <- 1
      col_index <- col_index + 1
    }
  }
  
  zeroblock <- matrix(0, ncol = length_param_rf, nrow = length_param)
  Bplot <- do.call(rbind, replicate(nobs, Bijk, simplify = FALSE))
  
  # Combine blocks for different subplots (SP)
  B1 <- rbind(
    do.call(rbind, replicate(length(data_plot$tree_ID[data_plot$SSP == 1]), Bijk, simplify = FALSE)),
    do.call(rbind, replicate(length(data_plot$tree_ID[data_plot$SSP == 2]), zeroblock, simplify = FALSE))
  )
  
  if (nSSP == 1) {
    B <- cbind(Bplot, B1)
  } else {
    B2 <- rbind(
      do.call(rbind, replicate(length(data_plot$tree_ID[data_plot$SSP == 1]), zeroblock, simplify = FALSE)),
      do.call(rbind, replicate(length(data_plot$tree_ID[data_plot$SSP == 2]), Bijk, simplify = FALSE))
    )
    B <- cbind(Bplot, B1, B2)
  }
  
  # Initialize random effects and set up for iterative estimation
  b0 <- rep(0, length_param_rf * (1 + nSSP))
  bOld <- b0 - 1
  history <- matrix(ncol = length(b0), nrow = 0)
  
  # Iteratively estimate random effects
  while (sum((b0 - bOld)^2) > 1e-20) {
    y <- c(data_plot$HT_m) # Response vector
    cf <- A %*% beta + B %*% b0 # Compute coefficients
    
    for (i in 1:length_param) {
      assign(paste0('cf', i), cf[seq(i, nrow(cf), by = length_param)])
    }
    
    # Model formulation for H predictions (Check the main manuscript)
    a1 <- cf1 + cf2 * log(data_plot$QMD) + cf3 * log(data_plot$G_species) + cf4 * data_plot$ddom
    a2 <- cf5 + cf6 * log(data_plot$QMD) + cf7 * sqrt(data_plot$N)
    b1 <- 0
    b2 <- 0
    
    # Predict H using the defined function
    predHT <- fun_HT_CL(dbh = data_plot$dbh, which = 1, a1, a2, b1, b2)
    pred <- c(predHT) 
    
    # Construct Z0 matrix for random effects
    X_HT <- matrix(attributes(predHT)$gradient[, 1:2], ncol = 2)
    
    if (nSSP == 1) {
      Z0 <- cbind(matrix(X_HT, ncol = 2), matrix(X_HT, ncol = 2))
    } else {
      group1_rows <- nrow(data_plot[data_plot$SSP == 1,])
      mat_group1 <- matrix(X_HT[1:group1_rows,], ncol = 2)
      mat_group2 <- matrix(X_HT[(group1_rows + 1):nrow(data_plot),], ncol = 2)
      Z0 <- cbind(X_HT, bdiag(matrix(mat_group1, ncol = 2), matrix(mat_group2, ncol=2)))
    }
    
    # Pool the random effects associated with the H model
    if(nSSP==1){index<-c(1,2,5,6)} else{index<-c(1,2,5,6,9,10)}
    
    # Estimate  random effects
    bOld <- b0
    b <- C %*% t(Z0) %*% solve(Z0 %*% D0 %*% t(Z0) + R0) %*% (y - pred + Z0 %*% b0)
    b0 <- b[index]
    
    # Record history of random effects
    history <- rbind(history, t(b0))
    cat(as.numeric(b), "\n")
  }
  
  return(b) # Return final estimates of random effects
}

# Perform b prediction and display results for data_plot. 
(b0<-cross_prediction_random_effects_sur(model = model, data_plot = data_plot))
# The first two random effects are plot level random efects for H and the next two for CL
# The 5th and 6th  random effects are for SP1 for H and the next two for CL.
# The 9th and 10th  random effects are for SP2 for H and the next two for CL
round(b0,3)
#------------------------------------------------------------------------------------------------------------
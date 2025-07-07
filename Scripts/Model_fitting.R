# -------------------------------------------------------------------------------------------------------------
# Title: Data Formatting, Model Formulation, and Fitting Procedures for Nonlinear Multivariate Mixed-Effects Models
# Author: Albert Ciceu,  Lauri Meht?talo
# Affiliation: Austrian Research Center for Forests (BFW)
# Corresponding Author: Albert Ciceu, Email: albert.ciceu@bfw.gv.at
#
# Description: This script accompanies the article "Nonlinear multilevel seemingly unrelated height-diameter 
# and crown length mixed-effects models for the southern Transylvanian forests, Romania" 
# https://doi.org/10.1016/j.fecs.2025.100322
# It covers the procedures for data formatting, model formulation, and fitting of nonlinear multivariate 
# mixed-effects models. For more details, refer to the main text of the article.
#
# Date: 17/09/2024
# -------------------------------------------------------------------------------------------------------------

# Load necessary libraries
require(tidyverse) # For data manipulation and visualization
require(nlme)      # For fitting nonlinear mixed-effects models

#------------------------------------------DATA FORMATTING-------------------------------------------------

# Load and filter data from CSV file
inv_reg <- read.csv('Output/Regression_data.csv')  # Read data from CSV file
inv_reg <- inv_reg %>% filter(Species == 'BE')     # Filter data for species 'BE'

# Combine data into a single dataframe for model fitting
data_nlme <- as.data.frame(rbind(
  cbind(
    plot = inv_reg$unique_ID,         # Unique ID for production unit and cluster unit
    subplot = inv_reg$SSP,            # Subplot identifier
    tree = inv_reg$tree_ID,           # Tree identifier
    index = 1,                        # Model indicator for Height (H) measurements
    y = inv_reg$HT_m,                 # Height measurements
    conHT = 1,                        # Dummy variable indicating H model
    dbh = inv_reg$dbh,                # Diameter at breast height
    qmd = inv_reg$QMD,                # Mean quadratic diameter
    ddom = inv_reg$ddom,              # Dominant diameter
    g_species = inv_reg$G_species,    # Basal area of the species per hectare
    N = inv_reg$N,                    # Number of trees per hectare
    balmod = inv_reg$BALMOD_c,        # BALMOD variable (see article)
    gini = inv_reg$gini               # Gini coefficient
  ),
  cbind(
    plot = inv_reg$unique_ID,         # Unique ID for production unit and cluster unit
    subplot = inv_reg$SSP,            # Subplot identifier
    tree = inv_reg$tree_ID,           # Tree identifier
    index = 2,                        # Model indicator for Crown Length (CL) measurements
    y = inv_reg$CL,                   # Crown Length measurements
    conHT = 0,                        # Dummy variable indicating CL model
    dbh = inv_reg$dbh,                # Diameter at breast height
    qmd = inv_reg$QMD,                # Mean quadratic diameter
    ddom = inv_reg$ddom,              # Dominant diameter
    g_species = inv_reg$G_species,    # Basal area of the species per hectare
    N = inv_reg$N,                    # Number of trees per hectare
    balmod = inv_reg$BALMOD_c,        # BALMOD variable (see article)
    gini = inv_reg$gini               # Gini coefficient
  )
))

# Convert all character columns (except 'plot') to numeric
data_nlme <- data_nlme %>% mutate(across(c(where(is.character), -c(plot)), as.numeric))

# Function to compute H or CL and its gradient
fun_HT_CL <- function(dbh, which, a1, a2, b1, b2) {
  # Calculate the H or CL based on the model
  value <- which * (1.3 + a1 * (dbh / (1 + dbh)) ^ a2) +
    (1 - which) * ((dbh / (b1 + exp(b2) * dbh)) ^ 2)
  
  # Compute the gradient of the function
  grad <- cbind(
    a1 = which * (dbh / (dbh + 1)) ^ a2,
    a2 = a1 * which * ((dbh / (dbh + 1)) ^ a2) * log(dbh / (dbh + 1)),
    b1 = ((2 * dbh^2) * (which - 1)) / (exp(b2) * dbh + b1)^3,
    b2 = ((2 * dbh^3) * (which - 1) * exp(b2)) / (b1 + dbh * exp(b2))^3
  )
  
  # Attach gradient information to the value
  attr(value, "gradient") <- grad
  return(value)
}

# Fit the nonlinear multivariate multilevel mixed-effects model
SUR_BEECH2 <- nlme(
  y ~ fun_HT_CL(dbh, conHT, a1, a2, b1, b2),  # Model formula specifying the function to use
  fixed = list(
    a1 ~ log(qmd) + log(g_species) + ddom,     # Fixed effects for parameter a1
    a2 ~ log(qmd) + sqrt(N),                    # Fixed effects for parameter a2
    b1 ~ log(qmd) + gini,                       # Fixed effects for parameter b1
    b2 ~ log(qmd) + g_species + balmod          # Fixed effects for parameter b2
  ),
  random = a1 + a2 + b1 + b2 ~ 1 | plot / subplot,  # Random effects structure at multiple levels
  corr = corSymm(form = ~ index | plot / subplot / tree),  # Correlation structure within the data
  weights = varComb(varIdent(form = ~ 1 | index),  # Combine variance functions; Estimate variance for each component
                    varPower(-0.5, form = ~ dbh | index)), # Variance function to handle heteroscedasticity
  data = data_nlme,  # Data used for fitting the model
  start = c(-41, 23, 1.6, -0.19, -3.8, 6, -0.15, -0.01, 0.6, -0.8, -0.5, -0.3, 0.001, 0.01),  # Initial values for parameters
  method = "REML",  # Restricted Maximum Likelihood method for model fitting
  control = lmeControl(
    returnObject = TRUE,  # Return the fitted object for further inspection
    msVerbose = TRUE,     # Print detailed information during the fitting process
    maxIter = 1000,       # Maximum number of iterations for the optimization algorithm
    msMaxIter = 10000,    # Maximum number of iterations for the marginal likelihood optimization
    niterEM = 1000,       # Number of EM algorithm iterations for refining variance-covariance estimates
    opt = "nlminb"        # Optimization algorithm used for parameter estimation
  )
)

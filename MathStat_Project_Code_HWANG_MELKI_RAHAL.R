# -- MATHEMATICAL STATISTICS PROJECT -- #
# -- PROJECT 1 - CCA, STOCHASTIC, EM -- #
# ------------------------------------- #
# HWANG Sungyun, MELKI Paul, RAHAL Mira #
# ------------------------------------- #
# ------- R version 3.6.1 ------------- #
# ------------------------------------- #


# PRELIMINARY LIBRARY IMPORTS ------------------------
# ----------------------------------------------------
library(mice)
library(ggplot2)
library(gridExtra)
library(TestDataImputation)
library(Amelia) # ! Please make sure to use Amelia version 1.7.3 or above !
# ----------------------------------------------------


# PRELIMINARY MANIPULATIONS
# BEFORE ANSWERING QUESTIONS
# ----------------------------------------------------
# STEP 0.1. ------------------------------------------
# Subset the dataset, keeping only the needed columns:
# salary (waget)
# education (educ) measured by the number of school years
# indicator of experience (exper)
# a variable (year88) indicating (label 1) whether the indvidual has been surveyed in 1988
# create a new variable (lwaget) the log(waget)
# ----------------------------------------------------
dataset <- test[test$year88 == 1, c("waget", "educ", "exper", "year88")]
dataset$lwaget <- log(dataset$waget)

# STEP 0.2. -------------------------------------------
# Run the linear regression model as specified, on the
# original dataset and keep the results.
# The regression to run is log(waget) on 'educ' and 'exper'
# -----------------------------------------------------
original_regression <- lm(lwaget ~ educ + exper, data = dataset)

# STEP 0.3. -------------------------------------------
# Create M amputed datasets using the 'mice' package 
# with missing values on the dependent variable "lwaget"
# with a proportion of missing values of 40%
# with MCAR missing mechanism
# -----------------------------------------------------
M <- 30    # number of amputed datasets to create
# Iterate M times, each time creating a new amputed dataset with the same specifications.
# Store all the information about the amputations in the same array called "amputations".
amputations <- list()
for (i in 1:M) {
  amputations[[i]] <-
    as.data.frame(ampute(
      data = dataset,
      prop = 0.4,
      mech = "MCAR",
      patterns = c(1,1,1,1,0)
    )$amp)
  amputations[[i]]$waget <- exp(amputations[[i]]$lwaget)
}
# We now have all the imputed datasets in one list

# -----------------------------------------------------
# -----------------------------------------------------


# QUESTION 1 --------- COMPLETE CASE ANALYSIS ---------
# -----------------------------------------------------
# STEP 1.1. -------------------------------------------
# Perform M simple single "imputations" (not really imputations)  
# on the M datasets with missing data using Complete Case Analysis
# this method is easy to apply as it works on removing each row that contains
# at least one missing (NA) value
# ----------------------------------------------------

# STEP 1.2. ------------------------------------------
# Run M linear regressions of log(waget) on 'educ' and 'exper'
# using the datasets obtained after application of CCA
# We save all the results of the all the regressions in the below
# list called 'q1_regression_results'.
# Furthermore, we create an array called 'q1_estimation_biases_avg' of length 3
# which will contain, respectively, the estimation bias for beta0 (intercept),
# beta1 (coefficient of 'educ'), beta2 (coefficient of 'exper').
# In addition, we create an array called 'q1_estimated_variances_avg' of length 3
# which will contain, respectively, the estimation bias for beta0 (intercept),
# beta1 (coefficient of 'educ'), beta2 (coefficient of 'exper').
# ----------------------------------------------------
q1_regression_results <- list()
q1_estimated_biases <- matrix(0, nrow = M, ncol = 3)
q1_estimated_variances <- matrix(0, nrow = M, ncol = 3)
q1_imputed_datasets <- list()

# BEGIN LOOP
for (i in 1:M) {
  
  # STEP 1.1
  q1_imputed_datasets[[i]] <- as.data.frame(na.omit(amputations[[i]]))
  # STEP 1.2
  q1_regression_results[[i]] <- lm(lwaget ~ educ + exper, data = q1_imputed_datasets[[i]])
  for (j in 1:3) {
    # Fill in the biases and variances
    q1_estimated_biases[j] <- q1_estimated_biases[j] + q1_regression_results[[i]]$coefficients[j] - original_regression$coefficients[j]
    q1_estimated_biases[i, j] <- q1_regression_results[[i]]$coefficients[j] - original_regression$coefficients[j]
    q1_estimated_variances[i, j] <- vcov(q1_regression_results[[i]])[j,j]
  }
}

q1_estimated_biases_avg <- colSums(q1_estimated_biases) / M
q1_estimated_variances_avg <- colSums(q1_estimated_variances) / M

# STEP 1.3. ------------------------------------------
# Plot the difference "Beta_educ_m" - "Beta_educ_original"
# in order to show how close the coefficients of 
# education estimated from the datasets on which CCA has
# been applied, to the coefficient of education estimated 
# using the complete dataset.
# Before plotting, we gather all education coefficients
# from the M datasets in one array.'
# ----------------------------------------------------
q1_educ_coefficients = c(1:M)
for (i in 1:M) {
  q1_educ_coefficients[i] <- q1_regression_results[[i]]$coefficients[2]
}

q1_differences_forPlotting <- data.frame(
  "dataset" = c(1:M),
  "educ_bias" = q1_estimated_biases[,2]
)

q1_educ_biases_plot <- ggplot(data = q1_differences_forPlotting, aes(x, y)) +
  geom_point(data = q1_differences_forPlotting, aes(x = dataset, y = educ_bias), shape = 23, fill = "darkblue") +
  geom_hline(yintercept = 0, color = "red") + 
  ylim(-0.05, 0.05) + 
  xlab(expression('m (Dataset number)')) + 
  ylab(expression(hat(beta)['educ']^'m'~'-'~hat(beta)['educ']^'C')) + 
  ggtitle(expression("Complete Case Analysis:"~hat(beta)['educ']^'m'~'-'~hat(beta)['educ']^'C'~"for the M (=30) datasets."))
q1_educ_biases_plot

# -----------------------------------------------------
# -----------------------------------------------------


# QUESTION 2 -------- STOCHASTIC REGRESSION -----------
# -----------------------------------------------------
# STEP 2.1. -------------------------------------------
# Perform B multiple imputations on each of the M 
# amputed datasets using the Stochastic Regression 
# imputation method. 
# The parameter 'm' of the 'mice' function allows for 
# the specification of the number of multiple imputations
# to perform. We set this parameter equal to B, in this
# project it is set to 5. However, we wrote the code 
# in a flexible way which allows us to perform any positive
# number of multiple imputations by giving a new value to B.
# -----------------------------------------------------

# STEP 2.2. ------------------------------------------
# Run M linear regressions of log(waget) on 'educ' and 'exper'
# using the datasets imputed by Stochastic Regression Imputation
# We save all the results of the all the regressions in the below
# list called 'q2_regression_results_stochastic'.
# Furthermore, we create an array called 'q2_estimated_biases_stochastic_avg' of length 3
# which will contain, respectively, the estimation bias for beta0 (intercept),
# beta1 (coefficient of 'educ'), beta2 (coefficient of 'exper').
# In addition, we create an array called 'q2_estimated_variances_stochastic_avg' of length 3
# which will contain, respectively, the estimation bias for beta0 (intercept),
# beta1 (coefficient of 'educ'), beta2 (coefficient of 'exper').
# To finish, we plot a graph showing the values imputed by Deterministic 
# Regression and Stochastic Regression, in order to show the importance of 
# Stochastic Regression Imputation in preserving the correlation between Y and X.
# ----------------------------------------------------
B <- 5
q2_imputed_datasets_stochastic <- list()
q2_regression_results_stochastic <- list()
q2_coefficients_stochastic <- matrix(0, M, 3)
q2_estimated_biases_stochastic <- matrix(0, nrow = M, ncol = 3)
q2_estimated_variances_stochastic <- matrix(0, nrow = M, ncol = 3)

# BEGIN LOOP
for (i in 1:M) {
  # STEP 2.1
  q2_imputed_datasets_stochastic[[i]] <-
    as.data.frame(complete(mice(
      data = amputations[[i]],
      method = "norm.nob",
      m = B,
      maxit = 1,
      seed = 1
    )))
  q2_imputed_datasets_stochastic[[i]]$waget <- exp(q2_imputed_datasets_stochastic[[i]]$waget)
  # STEP 2.2
  q2_regression_results_stochastic[[i]] <- lm(lwaget ~ educ + exper, data = q2_imputed_datasets_stochastic[[i]])
  for (j in 1:3) {
    # Fill in the biases and variances
    q2_coefficients_stochastic[i, j] <- q2_regression_results_stochastic[[i]]$coefficients[j]
    q2_estimated_biases_stochastic[i, j] <- q2_regression_results_stochastic[[i]]$coefficients[j] - original_regression$coefficients[j]
    q2_estimated_variances_stochastic[i, j] <- vcov(q2_regression_results_stochastic[[i]])[j,j]
  }
}

# Calculate the averages of the estimated biases and the estimated variances
q2_estimated_biases_stochastic_avg <- colSums(q2_estimated_biases_stochastic) / M
q2_estimated_variances_stochastic_avg <- colSums(q2_estimated_variances_stochastic) / M

# STEP 2.3. ------------------------------------------
# Plot the difference Beta_educ_m - Beta_educ_original 
# in order to show how close the coefficients of 
# education estimated from the datasets imputed by Stochastic Regression, 
# to the coefficient of education estimated using the complete dataset.
# Before plotting, we gather all education coefficients
# from the M datasets in one array.'
# ----------------------------------------------------
q2_educ_coefficients = c(1:M)
for (i in 1:M) {
  q2_educ_coefficients[i] <- q2_regression_results_stochastic[[i]]$coefficients[2]
}

q2_differences_forPlotting <- data.frame(
  "dataset" = c(1:M),
  "educ_bias" = q2_estimated_biases_stochastic[,2]
)

q2_educ_biases_plot <- ggplot(data = q2_differences_forPlotting, aes(x, y)) +
  geom_point(data = q2_differences_forPlotting, aes(x = dataset, y = educ_bias), shape = 21, fill = "green") +
  geom_hline(yintercept = 0, color = "red") + 
  ylim(-0.05, 0.05) + 
  xlab(expression('m (Dataset number)')) + 
  ylab(expression(hat(beta)['educ']^'m'~'-'~hat(beta)['educ']^'C')) + 
  ggtitle(expression("Stochastic Regression:"~hat(beta)['educ']^'m'~'-'~hat(beta)['educ']^'C'~"for the M (=30) datasets."))
q2_educ_biases_plot

# STEP 2.4. -------------------------------------------
# In order to show the importance of Stochastic Regression Imputation,
# and why it was created and is used, we impute one of the amputed datasets
# using both Deterministic Regression Imputation and Stochastic Regression
# Imputation. We then show a scatter of 'educ' data against 'lwaget' data. 
# -----------------------------------------------------
stochastic_imputed <- q2_imputed_datasets_stochastic[[1]][is.na(amputations[[1]]$lwaget), ]

deterministic_imputed <-
  as.data.frame(complete(mice(
    data = amputations[[1]],
    method = "norm.predict",
    m = B,
    maxit = 1,
    seed = 1
  )))[is.na(amputations[[1]]$lwaget), ]

ggplot(data = dataset, aes(x, y)) +
  geom_point(data = dataset, aes(x = educ, y = lwaget, color = "Original data")) + 
  geom_point(data = stochastic_imputed, aes(x = educ, y = lwaget, color = "Data Imputed by Stochastic Regression")) + 
  geom_point(data = deterministic_imputed, aes(x = educ, y = lwaget, color = "Data Imputed by Deterministic Regression")) +
  ggtitle("Comparison between Stochastic and Deterministic Regression") + 
  xlab("Education (Years)") + 
  ylab(expression(log(waget)))

# -----------------------------------------------------
# -----------------------------------------------------


# QUESTION 3 ------------- EM ALGORITHM ---------------
# -----------------------------------------------------
# In this question we focus on the EM Algorithm 
# (Expectation-Maximisation) which imputes the data numerically.
# We implement two versions of this 
# algorithm using the Amelia package developed at 
# Harvard University: simple imputation using EM Algorithm
# and bootstrap imputation using EM algorithm.
# -----------------------------------------------------

# In the following loop running on all amputed 
# datasets with missing data, we do the following tasks: 
# STEP 3.1 --------------------------------------------
# Impute the data using the EM single imputation
# STEP 3.2 --------------------------------------------
# Impute the data using the Bootstrap EM imputation
# STEP 3.3 --------------------------------------------
# Run the linear regression model on each newly imputed dataset
# STEP 3.4 --------------------------------------------
# Calculate and save the biases and variances
# of estimated coefficients obtained in each regression.
# ----------------------------------------------------- 
q3_imputed_datasets_single <- list()
q3_imputed_datasets_bootstrap <- list()
q3_regression_results_single <- list()
q3_regression_results_bootstrap <- list()
q3_coefficients_single <- matrix(0, nrow = M, ncol = 3)
q3_coefficients_bootstrap <- matrix(0, nrow = M, ncol = 3)
q3_estimated_biases_single <- matrix(0, nrow = M, ncol = 3)
q3_estimated_biases_single <- matrix(0, nrow = M, ncol = 3)
q3_estimated_biases_bootstrap <- matrix(0, nrow = M, ncol = 3)
q3_estimated_variances_single <- matrix(0, nrow = M, ncol = 3)
q3_estimated_variances_bootstrap <- matrix(0, nrow = M, ncol = 3)

# BEGIN LOOP 
for (i in 1:M) {
  # STEP 3.1
  q3_imputed_datasets_single[[i]] <- amelia(amputations[[i]][, c("waget", "exper", "educ", "lwaget")], boot.type = "none")$imputations$imp5
  # STEP 3.2
  q3_imputed_datasets_bootstrap[[i]] <- amelia(amputations[[i]][, c("waget", "exper", "educ", "lwaget")])$imputations$imp5
  # STEP 3.3
  q3_regression_results_single[[i]] <- lm(
    lwaget ~ educ + exper, data = q3_imputed_datasets_single[[i]]
  )
  q3_regression_results_bootstrap[[i]] <- lm(
    lwaget ~ educ + exper, data = q3_imputed_datasets_bootstrap[[i]]
  )
  # STEP 3.4
  for (j in 1:3) {
    q3_coefficients_single[i, j] <- q3_regression_results_single[[i]]$coefficients[j]
    q3_estimated_biases_single[i, j] <- q3_regression_results_single[[i]]$coefficients[j] - original_regression$coefficients[j]
    q3_estimated_variances_single[i, j] <- vcov(q3_regression_results_single[[i]])[j,j]
    q3_coefficients_bootstrap[i, j] <- q3_regression_results_bootstrap[[i]]$coefficients[j]
    q3_estimated_biases_bootstrap[i, j] <- q3_regression_results_bootstrap[[i]]$coefficients[j] - original_regression$coefficients[j]
    q3_estimated_variances_bootstrap[i, j] <- vcov(q3_regression_results_bootstrap[[i]])[j,j]
  }
}

q3_estimated_bias_single_avg <- colSums(q3_estimated_biases_single) / M
q3_estimated_variance_single_avg <- colSums(q3_estimated_variances_single) / M
q3_estimated_bias_bootstrap_avg <- colSums(q3_estimated_biases_bootstrap) / M
q3_estimated_variance_bootstrap_avg <- colSums(q3_estimated_variances_bootstrap) / M

q3_estimated_bias_single_avg
q3_estimated_bias_bootstrap_avg
q3_estimated_variance_single_avg
q3_estimated_variance_bootstrap_avg

# STEP 3.5. ------------------------------------------
# Plot the difference Beta_educ_m - Beta_educ_original 
# in order to show how close the coefficients of 
# education estimated from the datasets imputed by EM Algorithm,
# simple and bootstrap versions, to the coefficient of education 
# estimated using the complete dataset.
# Before plotting, we gather all education coefficients
# from the M datasets in one array, for each type of EM.
# ----------------------------------------------------
q3_educ_coefficients_bootstrap <- c(1:M)
q3_educ_coefficients_single <- c(1:M)
for (i in 1:M) {
  q3_educ_coefficients_single[i] <- q3_regression_results_single[[i]]$coefficients[2]
  q3_educ_coefficients_bootstrap[i] <- q3_regression_results_bootstrap[[i]]$coefficients[2]
}

q3_differences_forPlotting <- data.frame(
  "dataset" = c(1:M),
  "educ_bias_single" = abs(q3_estimated_biases_single[,2]),
  "educ_bias_bootstrap" = abs(q3_estimated_biases_bootstrap[,2])
)

q3_educ_biases_plot <- ggplot(data = q3_differences_forPlotting, aes(x, y)) +
  geom_point(data = q3_differences_forPlotting, aes(x = dataset, y = educ_bias_single, fill="Single EM"), shape = 21) +
  geom_point(data = q3_differences_forPlotting, aes(x = dataset, y = educ_bias_bootstrap, fill="Bootstrap EM"), shape = 21) +
  geom_hline(yintercept = 0, color = "red") + 
  ylim(-0.05, 0.05) + 
  xlab(expression('m (Dataset number)')) + 
  ylab(expression(hat(beta)['educ']^'m'~'-'~hat(beta)['educ']^'C')) + 
  ggtitle(expression("EM Algorithm Imputation:"~hat(beta)['educ']^'m'~'-'~hat(beta)['educ']^'C'~"for the M (=30) datasets."))
q3_educ_biases_plot

q3_educ_biases_plot_zoomIN <- ggplot(data = q3_differences_forPlotting, aes(x, y)) +
  geom_point(data = q3_differences_forPlotting, aes(x = dataset, y = educ_bias_single, fill="Single EM"), shape = 21, color="gray") +
  geom_point(data = q3_differences_forPlotting, aes(x = dataset, y = educ_bias_bootstrap, fill="Bootstrap EM"), shape = 21) +
  geom_hline(yintercept = 0, color = "red") + 
  ylim(0, 0.002) + 
  xlab(expression('m (Dataset number)')) + 
  ylab(expression(hat(beta)['educ']^'m'~'-'~hat(beta)['educ']^'C')) + 
  ggtitle(expression("Stochastic Regression:"~hat(beta)['educ']^'m'~'-'~hat(beta)['educ']^'C'~"for the M (=30) datasets."))
q3_educ_biases_plot_zoomIN

# -----------------------------------------------------
# -----------------------------------------------------

# END.
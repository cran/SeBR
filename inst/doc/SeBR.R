## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6, fig.height = 4
)

## ----install------------------------------------------------------------------
# install.packages("devtools")
# devtools::install_github("drkowal/SeBR")
library(SeBR) 

## ----sims---------------------------------------------------------------------
set.seed(123) # for reproducibility

# Simulate data from a transformed linear model:
dat = simulate_tlm(n = 200,  # number of observations
                   p = 10,   # number of covariates 
                   g_type = 'step' # type of transformation (here, positive data)
                   )
# Training data:
y = dat$y; X = dat$X 

# Testing data:
y_test = dat$y_test; X_test = dat$X_test 

## ----sblm---------------------------------------------------------------------
# Fit the semiparametric Bayesian linear model:
fit = sblm(y = y, 
           X = X, 
           X_test = X_test)

names(fit) # what is returned

## ----ppd, echo = FALSE--------------------------------------------------------
# Posterior predictive checks on testing data: empirical CDF

# Construct a plot w/ the unique testing values:
y0 = sort(unique(y_test))
plot(y0, y0, type='n', ylim = c(0,1),
     xlab='y', ylab='F_y', main = 'Posterior predictive ECDF: testing data')

# Add the ECDF of each simulated predictive dataset
temp = sapply(1:nrow(fit$post_ypred), function(s)
  lines(y0, ecdf(fit$post_ypred[s,])(y0), 
        col='gray', type ='s'))

# Add the ECDF for the observed testing data
lines(y0, ecdf(y_test)(y0),  # ECDF of testing data
     col='black', type = 's', lwd = 3)

## ----pi-----------------------------------------------------------------------
# Evaluate posterior predictive means and intervals on the testing data:
plot_pptest(fit$post_ypred, 
            y_test, 
            alpha_level = 0.10) # coverage should be >= 90% 

## ----comp, echo = FALSE-------------------------------------------------------
# Summarize the transformation:

# post_g contains draws of g evaluated at the sorted and unique y-values:
y0 = sort(unique(y)) 

# Construct a plot:
plot(y0, fit$post_g[1,], type='n', ylim = range(fit$post_g),
     xlab = 'y', ylab = 'g(y)', main = "Posterior draws of the transformation")

# Add the posterior draws of g:
temp = sapply(1:nrow(fit$post_g), function(s)
  lines(y0, fit$post_g[s,], col='gray')) # posterior draws

# Add the posterior mean of g:
lines(y0, colMeans(fit$post_g), lwd = 3) # posterior mean

# Add the true transformation, rescaled for easier comparisons:
lines(y, 
      scale(dat$g_true)*sd(colMeans(fit$post_g)) + mean(colMeans(fit$post_g)), 
      type='p', pch=2)

legend('bottomright', c('Truth'), pch = 2) # annotate the true transformation

## ----comp-theta---------------------------------------------------------------

# Summarize the parameters (regression coefficients):

# Posterior means:
coef(fit)

# Check: correlation with true coefficients
cor(dat$beta_true[-1],
    coef(fit)[-1]) # excluding the intercept

# 95% credible intervals:
theta_ci = t(apply(fit$post_theta, 2, quantile, c(.025, 0.975)))

# Check: agreement on nonzero coefficients?
which(theta_ci[,1] >= 0 | theta_ci[,2] <=0) # 95% CI excludes zero
which(dat$beta_true != 0) # truly nonzero

## ----sims-qr------------------------------------------------------------------
# Simulate data from a heteroskedastic linear model (no transformation):
dat = simulate_tlm(n = 200,  # number of observations
                   p = 10,   # number of covariates 
                   g_type = 'box-cox', lambda = 1, # no transformation
                   heterosked = TRUE # heteroskedastic errors
                   )
# Training data:
y = dat$y; X = dat$X 

# Testing data:
y_test = dat$y_test; X_test = dat$X_test 

## ----fit-qr-------------------------------------------------------------------
# Quantile to target:
tau = 0.05

# (Traditional) Bayesian quantile regression:
fit_bqr = bqr(y = y, 
           X = X, 
           tau = tau, 
           X_test = X_test,
           verbose = FALSE  # omit printout
)

# Semiparametric Bayesian quantile regression:
fit = sbqr(y = y, 
           X = X, 
           tau = tau, 
           X_test = X_test,
           verbose = FALSE # omit printout
)
      
names(fit) # what is returned

## ----ppd-bqr, echo = FALSE----------------------------------------------------
# Posterior predictive checks on testing data: empirical CDF

# Construct a plot w/ the unique testing values:
y0 = sort(unique(y_test))
plot(y0, y0, type='n', ylim = c(0,1),
     xlab='y', ylab='F_y', 
     main = 'Posterior predictive ECDF: testing data (black) \n bqr (red) vs. sbqr (gray)')

# Add the ECDF of each simulated predictive (bqr) dataset
temp = sapply(1:nrow(fit_bqr$post_ypred), function(s)
  lines(y0, ecdf(fit_bqr$post_ypred[s,])(y0), 
        col='red', type ='s'))

# Same, for sbqr:
temp = sapply(1:nrow(fit$post_ypred), function(s)
  lines(y0, ecdf(fit$post_ypred[s,])(y0), 
        col='gray', type ='s'))

# Add the ECDF for the observed testing data
lines(y0, ecdf(y_test)(y0),  # ECDF of testing data
     col='black', type = 's', lwd = 3)

## ----bqr-test-----------------------------------------------------------------
# Quantile point estimates:
q_hat_bqr = fitted(fit_bqr) 

# Empirical quantiles on testing data:
(emp_quant_bqr = mean(q_hat_bqr >= y_test))

# Evaluate posterior predictive means and intervals on the testing data:
(emp_cov_bqr = plot_pptest(fit_bqr$post_ypred, 
                           y_test, 
                           alpha_level = 0.10))

## ----sbqr-test----------------------------------------------------------------
# Quantile point estimates:
q_hat = fitted(fit) 

# Empirical quantiles on testing data:
(emp_quant_sbqr = mean(q_hat >= y_test))

# Evaluate posterior predictive means and intervals on the testing data:
(emp_cov_sbqr = plot_pptest(fit$post_ypred, 
                            y_test, 
                            alpha_level = 0.10))

## ----sim-gp-------------------------------------------------------------------
# Training data:
n = 200 # sample size
x = seq(0, 1, length = n) # observation points

# Testing data:
n_test = 1000 
x_test = seq(0, 1, length = n_test) 

# True inverse transformation:
g_inv_true = function(z) 
  qbeta(pnorm(z), 
        shape1 = 0.5, 
        shape2 = 0.1) # approx Beta(0.5, 0.1) marginals

# Training observations:
y = g_inv_true(
  sin(2*pi*x) + sin(4*pi*x) + .25*rnorm(n)
             ) 

# Testing observations:
y_test = g_inv_true(
  sin(2*pi*x_test) + sin(4*pi*x_test) + .25*rnorm(n)
             ) 

plot(x_test, y_test, 
     xlab = 'x', ylab = 'y',
     main = "Training (gray) and testing (black) data")
lines(x, y, type='p', col='gray', pch = 2)

## ----fit-bc-------------------------------------------------------------------
# Fit the Box-Cox Gaussian process model:
fit_bc = bgp_bc(y = y, 
           locs = x,
           locs_test = x_test)

# Fitted values on the testing data:
y_hat_bc = fitted(fit_bc)

# 90% prediction intervals on the testing data:
pi_y_bc = t(apply(fit_bc$post_ypred, 2, quantile, c(0.05, .95))) 

# Average PI width:
(width_bc = mean(pi_y_bc[,2] - pi_y_bc[,1]))

# Empirical PI coverage:
(emp_cov_bc = mean((pi_y_bc[,1] <= y_test)*(pi_y_bc[,2] >= y_test)))

# Plot these together with the actual testing points:
plot(x_test, y_test, type='n', 
     ylim = range(pi_y_bc, y_test), xlab = 'x', ylab = 'y', 
     main = paste('Fitted values and prediction intervals: \n Box-Cox Gaussian process'))

# Add the intervals:
polygon(c(x_test, rev(x_test)),
        c(pi_y_bc[,2], rev(pi_y_bc[,1])),
        col='gray', border=NA)
lines(x_test, y_test, type='p') # actual values
lines(x_test, y_hat_bc, lwd = 3) # fitted values

## ----fit-gp-------------------------------------------------------------------
# Fit the semiparametric Gaussian process model:
fit = sbgp(y = y, 
           locs = x,
           locs_test = x_test)

names(fit) # what is returned

coef(fit) # estimated regression coefficients (here, just an intercept)

## ----oos-gp-------------------------------------------------------------------
# Fitted values on the testing data:
y_hat = fitted(fit)

# 90% prediction intervals on the testing data:
pi_y = t(apply(fit$post_ypred, 2, quantile, c(0.05, .95))) 

# Average PI width:
(width = mean(pi_y[,2] - pi_y[,1]))

# Empirical PI coverage:
(emp_cov = mean((pi_y[,1] <= y_test)*(pi_y[,2] >= y_test)))

# Plot these together with the actual testing points:
plot(x_test, y_test, type='n', 
     ylim = range(pi_y, y_test), xlab = 'x', ylab = 'y', 
     main = paste('Fitted values and prediction intervals: \n semiparametric Gaussian process'))

# Add the intervals:
polygon(c(x_test, rev(x_test)),
        c(pi_y[,2], rev(pi_y[,1])),
        col='gray', border=NA)
lines(x_test, y_test, type='p') # actual values
lines(x_test, y_hat, lwd = 3) # fitted values


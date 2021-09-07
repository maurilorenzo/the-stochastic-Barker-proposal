# stan model Logistic Regression
source("logistic_regression_setup.R")

library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = 1)

betaPrior = function(d) {
  return(diag(rep(1, d)))
}

# simulated data
d = 100
N = 1000  
set.seed(123)
dataLR_stan = genDataLR(d, N + 250, algorithm = 'STAN')
names(dataLR_stan)
dataLR_stan$theta[1:20]
dataset_stan = dataLR_stan$data

# Get prior variance for beta parameters
betaVar = betaPrior(d)

dataset = list("N" = N, "X" = dataset_stan$X[1:N,], "y" = dataset_stan$y[1:N], "d" = d, "Sigma0" = betaVar)
stan_output_new = stan("C:\\Users\\Lorenzo\\Desktop\\DSBA\\tesi\\logistic_regression_model.stan", data = dataset, iter = 100000, chains = 1)
list_of_draws_new <- extract(stan_output_new)
stan_means_new = colMeans(list_of_draws_new$beta)
stan_vars_new = colSds(list_of_draws_new$beta)^2



# arrythmia - keeping 100 covariates
y_stan = train_set$y
x_stan = train_set[,-1]
d = dim(arrhythmia)[2] - 1
N = dim(train_set)[1]
betaVar = betaPrior(d)
arrhythmia_stan_new = list("N" = N, "X" = x_stan, "y" = y_stan, "d" = d, "Sigma0" = betaVar)

stan_arrhythmia_output_new = stan("C:\\Users\\Lorenzo\\Desktop\\DSBA\\tesi\\logistic_regression_model.stan", data = arrhythmia_stan_new, iter = 100000, chains = 1)

list_of_draws_arrythmia_new <- extract(stan_arrhythmia_output_new)
print(names(list_of_draws_arrythmia_new))
stan_means_arrhythmia_new = colMeans(list_of_draws_arrythmia_new$beta)
stan_vars_arrhythmia_new = colSds(list_of_draws_arrythmia_new$beta)^2


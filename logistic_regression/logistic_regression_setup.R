##### 6.1 Logistic Regression Model
source("../utils.R")
source("../ksd.R")

library(coda)
library(mvtnorm)
library(sBIC)
library(ISLR)
library(plot3D)




# grad log - likelihood
gradLogL = function(data, theta){
  # it returns the gradient of the log likelihood of the logistic regression model 
  # args:
  # theta (vector d * 1) = parameter vector of the model
  # data = dataset: data$y (vector n * 1) = response variable; data$X (matrix n * d) = covariates
  # return: 
  # res (vector d * 1) = gradient of the log likelihood
  res = t(as.matrix(data[-1])) %*% ((as.matrix(data[1])- exp(as.matrix(data[-1]) %*% theta) 
                                    / (1 + exp(as.matrix(data[-1]) %*% theta))))
  return(res)
}




# grad log - prior
gradLogP = function(theta) {
  # it return the gradient log prior (mutlivariate normal with diagonal covariance)
  # PARAMS:
  # theta (vector d * 1) = parameter vector
  # OUTPUT:
  # gradient of the log prior (vector d * 1)
  return(- theta)
}

# grad log - posterior
gradLogPost = function(data, theta) {
  # it returns the gradient of the log posterior of the logistic regression model with normal prior
  # PARAMS:
  # theta (vector d * 1) = parameter vector of the model
  # data = dataset: data$y (vector n * 1) = response variable; data$X (matrix n * d) = covariates
  # OUTPUT: 
  # res (vector d * 1) = gradient of the log posterior
  res = gradLogP(theta) + gradLogL(data, theta)
  return(res)
}


genDataLR = function (d, N, algorithm = 'SGMCMC', reg_type = 'logistic'){
  # it generates the data set used to fit the logistic regression
  # PARAMS:
  # d (scalar) = dimensionality of the parameter
  # N (scalar) = size of the dataset created
  # OUTPUT:
  # list: data (list): data$X (matrix N * d), data$y (vector N * 1)
  #       theta (vector d * 1) = vector of parameter used to generate the data
  #beta_sds = c(0.1, 1, 10)
  var = 100/d
  theta = rnorm(d, 0, sqrt(var))
  Sigma0 = genCov(d, 0.4)
  X = rmvnorm(N, rep(0, d), Sigma0)
  if (reg_type == 'logistic'){
    proba = apply(X, 1, function(x_i) 1 / (1+ exp(-sum(x_i * theta))))
  }
  else{
    proba = apply(X, 1, function(x_i) (pnorm(sum(x_i * theta))))
  }
  y = rbinom(N, 1, proba)
  if (algorithm == 'SGMCMC'){data = data.frame("y" = y, "X" = X)}
  else{data = list("X" = X, "y" = y)}
  res = list("data" = data, "theta" = theta, 'Sigma' = Sigma0)
  return(res)
}



genCov = function(d, rho){
  # it generates the covariance matrix
  # Of the form diag(Sigma0) = 1 and Sigma0[i,j] ~ Unif(-rho, rho)^(i-j)
  # The aim is to emulate standardised covariates with a small amount of correlation
  # PARAMS:
  # d (scalar) = size of the matrix
  # rho (scalar) = parameter controlling the correlation among the covariates
  # OUTPUT:
  # res (matrix d * d) = covariance matrix
  res = diag(d)
  for (i in 1:d-1){
    for (j in (i+1):d){
      res[i,j] = runif(1, -rho, rho) ^ (j - i)
      res[j,i] = res[i,j]
    }
  }
  return(res)
}


logLossLR = function(data, theta){
  # this function compute the log loss of the Logistic Regression
  # PARAMS:
  # data (data.frame N by (d+1)) = held out dataset
  # theta (d-dimesional vector) = parameters
  # OUTPUT:
  # logLoss (scalar)
  N = dim(data)[1]
  x_theta = rowSums(t(theta * t(data[names(data)!='y'])))
  p = (1 + exp(- x_theta)) ^ (-1)
  ll = prod(p ^ (data$y) * (1 - p) ^ (1 - data$y))
  logLoss = - 1 / N * log(ll)
  return(logLoss)
}

x_theta = function(df, x){
  return(rowSums(t(x * t(df[names(df)!='y']))))
}

logLossLR_alt = function(data, samples){
  N = dim(data)[1]
  x_t = apply(samples, 1, function(x) x_theta(data, x))
  ps = (1 + exp(- x_t)) ^ (-1)
  print(dim(ps))
  p = rowMeans(ps)
  ll = prod((p ^ (data$y)) * ((1 - p) ^ (1 - data$y)))
  logLoss = - 1 / N * log(ll, exp(1))
  return(logLoss)
}



avgLogLossLR = function(samples, data){
  # this function compute the average log loss using the last 1000 MCMC samples
  # PARAMS:
  # data (data.frame N by (d+1)) = held out dataset
  # samples (matrix N by d) = MCMC samples of the parameters
  # OUTPUT:
  # res (scalar)
  N = dim(samples)[1]
  res = apply(samples, 1, function(x) logLossLR(data, x))
  res = sum(res) / N
  return(res)
}


sim_logistic_regression = function(h, sigma_2s, sigma_2s_alternative, S, N, n, d, ksd = TRUE, ksd_thinned = TRUE, moments = FALSE, stan_means=0, stan_vars=0){
  # helper function used to simulate the sgld and stochastic barker with 
  # different value of the step size and compute the ksd, minimum ess and 
  # the log loss on a held out set
  set.seed(123)
  
  dataLR = genDataLR(d, N + 250)
  grad_log_l = gradLogL
 
  print(stan_means[1])
  dataset = dataLR$data[1:N,]
  held_out_test = dataLR$data[c((N+1):(N+250)),]
  theta = dataLR$theta
  print(theta[1:10])
  Sigma = dataLR$Sigma
  
  # starting point of the algorithm
  theta0 = rnorm(d, 0, 1)
  res_sgld = data.frame(t(apply(matrix(h), 1, function(x) (mcmc_logistic_regression(x, S, n, theta0, dataset, held_out_test, grad_log_l, gradLogP, SGLD_naive, ksd, ksd_thinned, moments, type=type, stan_means = stan_means, stan_vars = stan_vars)))))
  res_barker = data.frame(t(apply(matrix(sigma_2s), 1, function(x) (mcmc_logistic_regression(x, S, n, theta0, dataset, held_out_test, grad_log_l, gradLogP, barker, ksd, ksd_thinned, moments, type=type, stan_means = stan_means, stan_vars = stan_vars)))))
  res_barker_alternative = data.frame(t(apply(matrix(sigma_2s_alternative), 1, function(x) (mcmc_logistic_regression(x, S, n, theta0, dataset, held_out_test, grad_log_l, gradLogP, barker, ksd, ksd_thinned, moments, version='alternative', type=type, stan_means = stan_means, stan_vars = stan_vars)))))
  
  if (moments){
      col_names = c('step_size', 'ksd', 'log_loss_test', 'log_loss_train', 'min_ESS', 'median_ESS', 'bias_mean', 'bias_mean_max','bias_mean_relative' ,'bias_var', 'bais_var_max')
    }
  else {
      col_names = c('step_size', 'ksd', 'log_loss', 'log_loss_train','min_ESS', 'median_ESS')
  }
  
  colnames(res_sgld) = col_names
  colnames(res_barker) = col_names
  colnames(res_barker_alternative) = col_names
  
  plot_logistic_regression(res_barker, res_barker_alternative, res_sgld, n, ksd, moments)
  return(list("res_sgld" = res_sgld, "res_barker" = res_barker, "res_barker_alternative" = res_barker_alternative))
}


mcmc_logistic_regression = function(step_size, S, n, theta0, dataset, held_out_test, gradLogL, gradLogP, algorithm, ksd = TRUE, ksd_thinned = TRUE, moments = TRUE, version = NULL, stan_means = 0, stan_vars = 0, type='logistic'){
  options(digits=8)
  print(step_size)
  set.seed(123)
  if (is.null(version)){
    samples = algorithm(dataset, S, n, theta0, step_size, gradLogL, gradLogP)
  }
  else{
    samples = algorithm(dataset, S, n, theta0, step_size, gradLogL, gradLogP, version)
  }
  if (any(is.na(samples))){
    return(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
  }
  samples = samples[(S/2+1):S,]
  options(digits=15)
  if (type=='logistic'){
    LL_test = logLossLR_alt(held_out_test, samples)
    LL_train = logLossLR_alt(dataset, samples)
  }
  else{
    LL_test = rmse(held_out_test, samples)
    LL_train = rmse(dataset, samples)
  }
  ESS_samples = ESS(samples)
  min_ESS = min(ESS_samples)
  median_ESS = median(ESS_samples)
  
  if (moments){
    bias_mean_samples= mean(abs(colMeans(samples) - stan_means))
    bias_mean_max_samples = max(abs(colMeans(samples) - stan_means))
    bias_mean_relative_samples= mean(abs(colMeans(samples) - stan_means)/sqrt(stan_vars))
    bias_var_samples = mean(abs(colSds(samples)^2 - stan_vars))
    bias_var_max_samples = max(abs(colSds(samples)^2 - stan_vars))
  }
  options(digits = 15)
  if (ksd){
    if (ksd_thinned){
      samples = samples[seq(1, S/2+1, length.out = 10001),]
    }
    grad = t(apply(samples, 1, function(x) gradLogPost(dataset, x)))
    KSD_samples = imqKSD(samples, grad, c = 1, beta = 0.5)
  }
  else {KSD_samples = 0}
  options(digits=8)
  if(moments){
    print(c(step_size, KSD_samples, LL_test, LL_train, min_ESS, median_ESS, bias_mean_samples, bias_mean_max_samples, bias_mean_relative_samples, bias_var_samples, bias_var_max_samples))
    return(c(step_size, KSD_samples, LL_test, LL_train, min_ESS, median_ESS, bias_mean_samples, bias_mean_max_samples, bias_mean_relative_samples, bias_var_samples, bias_var_max_samples))
  }
  return(c(step_size, KSD_samples, LL_test, LL_train, min_ESS, median_ESS))
}


plot_logistic_regression = function(res_barker, res_barker_alternative, res_sgld, n, ksd, moments, log_scale=FALSE, max_=FALSE){
  par(mfrow = c(1,1))
  if (log_scale){
    res_barker[ c('bias_mean', 'bias_mean_relative','bias_var')] = log(1 + res_barker[, c('bias_mean', 'bias_mean_relative','bias_var')])
    res_barker_alternative[ c('bias_mean', 'bias_mean_relative','bias_var')] = log(1 + res_barker_alternative[, c('bias_mean', 'bias_mean_relative','bias_var')])
    res_sgld[ c('bias_mean', 'bias_mean_relative','bias_var')] = log(1 + res_sgld[, c('bias_mean', 'bias_mean_relative','bias_var')])
  }
  if (ksd){
    x_max_1 = max(c(max(res_barker['min_ESS']), max(res_barker_alternative['min_ESS']), max(res_sgld['min_ESS'])))
    x_max_2 = max(c(max(res_barker['median_ESS']), max(res_barker_alternative['median_ESS']), max(res_sgld['median_ESS'])))
    y_min = min(c(min(res_barker['ksd']), min(res_barker_alternative['ksd']), min(res_sgld['ksd'])))
    y_max = max(c(max(res_barker['ksd']), max(res_barker_alternative['ksd']), max(res_sgld['ksd'])))
    plot(res_barker[,'min_ESS'], res_barker[,'ksd'], main = paste('Speed of convergence vs Accuracy, n:', n),
         xlab = 'min ESS', ylab = 'ksd', xlim = c(0, x_max_1), ylim = c(y_min, y_max), type = 'o')
    points(res_barker_alternative[,'min_ESS'], res_barker_alternative[,'ksd'], col = 'blue', type = 'o')
    points(res_sgld[,'min_ESS'], res_sgld[,'ksd'], col = 'red', type = 'o')
    plot(res_barker[,'median_ESS'], res_barker[,'ksd'], main = paste('Speed of convergence vs Accuracy, n:', n),
         xlab = 'median ESS', ylab = 'ksd', xlim = c(0, x_max_1), ylim = c(y_min, y_max), type = 'o')
    points(res_barker_alternative[,'median_ESS'], res_barker_alternative[,'ksd'], col = 'blue', type = 'o')
    points(res_sgld[,'median_ESS'], res_sgld[,'ksd'], col = 'red', type = 'o')
  }
  if (moments){
    x_max_1 = max(c(max(res_barker['min_ESS']), max(res_barker_alternative['min_ESS']), max(res_sgld['min_ESS'])))
    x_max_2 = max(c(max(res_barker['median_ESS']), max(res_barker_alternative['median_ESS']), max(res_sgld['median_ESS'])))
    y_min_1 = min(c(min(res_barker['bias_mean']), min(res_barker_alternative['bias_mean']), min(res_sgld['bias_mean'])))
    y_max_1 = max(c(max(res_barker['bias_mean']), max(res_barker_alternative['bias_mean']), max(res_sgld['bias_mean'])))
    y_min_2 = min(c(min(res_barker['bias_var']), min(res_barker_alternative['bias_var']), min(res_sgld['bias_var'])))
    y_max_2 = max(c(max(res_barker['bias_var']), max(res_barker_alternative['bias_var']), max(res_sgld['bias_var'])))

    plot(res_barker[,'min_ESS'], res_barker[,'bias_mean'],
         xlab = 'Min ESS', ylab = 'Bias mean', xlim = c(0, x_max_1), ylim = c(y_min_1, y_max_1), type = 'o')
    points(res_barker_alternative[,'min_ESS'], res_barker_alternative[,'bias_mean'], col = 'blue', type = 'o')
    points(res_sgld[,'min_ESS'], res_sgld[,'bias_mean'], col = 'red', type = 'o')
    plot(res_barker[,'median_ESS'], res_barker[,'bias_mean'], 
         xlab = 'Median ESS', ylab = 'Bias mean', xlim = c(0, x_max_2), ylim = c(y_min_1, y_max_1), type = 'o')
    points(res_barker_alternative[,'median_ESS'], res_barker_alternative[,'bias_mean'], col = 'blue', type = 'o')
    points(res_sgld[,'median_ESS'], res_sgld[,'bias_mean'], col = 'red', type = 'o')
    
    plot(res_barker[,'min_ESS'], res_barker[,'bias_var'], 
         xlab = 'Min ESS', ylab = 'Bias variance', xlim = c(0, x_max_1), ylim = c(y_min_2, y_max_2), type = 'o')
    points(res_barker_alternative[,'min_ESS'], res_barker_alternative[,'bias_var'], col = 'blue', type = 'o')
    points(res_sgld[,'min_ESS'], res_sgld[,'bias_var'], col = 'red', type = 'o')
    plot(res_barker[,'median_ESS'], res_barker[,'bias_var'], 
         xlab = 'Median ESS', ylab = 'Bias var', xlim = c(0, x_max_2), ylim = c(y_min_2, y_max_2), type = 'o')
    points(res_barker_alternative[,'median_ESS'], res_barker_alternative[,'bias_var'], col = 'blue', type = 'o')
    points(res_sgld[,'median_ESS'], res_sgld[,'bias_var'], col = 'red', type = 'o')
    if(max_ == TRUE){
      y_max_3 = max(c(max(res_barker['bais_var_max']), max(res_barker_alternative['bais_var_max']), max(res_sgld['bais_var_max'])))
      plot(res_barker[,'min_ESS'], res_barker[,'bias_mean_max'], 
           xlab = 'Min ESS', ylab = 'Bias mean (max)', xlim = c(0, x_max_1), ylim = c(y_min_2, y_max_2), type = 'o')
      points(res_barker_alternative[,'min_ESS'], res_barker_alternative[,'bias_mean_max'], col = 'blue', type = 'o')
      points(res_sgld[,'min_ESS'], res_sgld[,'bias_mean_max'], col = 'red', type = 'o')
      plot(res_barker[,'median_ESS'], res_barker[,'bias_mean_max'], 
           xlab = 'Median ESS', ylab = 'Bias mean (max)', xlim = c(0, x_max_2), ylim = c(y_min_2, y_max_2), type = 'o')
      points(res_barker_alternative[,'median_ESS'], res_barker_alternative[,'bias_mean_max'], col = 'blue', type = 'o')
      points(res_sgld[,'median_ESS'], res_sgld[,'bias_mean_max'], col = 'red', type = 'o')
      plot(res_barker[,'min_ESS'], res_barker[,'bias_var_max'], 
          xlab = 'Min ESS', ylab = 'Bias var (max)', xlim = c(0, x_max_1), ylim = c(y_min_2, y_max_3), type = 'o')
      points(res_barker_alternative[,'min_ESS'], res_barker_alternative[,'bias_var_max'], col = 'blue', type = 'o')
      points(res_sgld[,'min_ESS'], res_sgld[,'bias_var_max'], col = 'red', type = 'o')
      plot(res_barker[,'median_ESS'], res_barker[,'bias_var_max'], 
           xlab = 'Median ESS', ylab = 'Bias var (max)', xlim = c(0, x_max_2), ylim = c(y_min_2, y_max_3), type = 'o')
      points(res_barker_alternative[,'median_ESS'], res_barker_alternative[,'bias_var_max'], col = 'blue', type = 'o')
      points(res_sgld[,'median_ESS'], res_sgld[,'bias_var_max'], col = 'red', type = 'o')
    }
  }
  
}

plot_trace_LR = function(idx, sgldLR, barkerLR, barkerLR_alternative, list_of_draws, stan_means, ylim = c(-10, 10), legend=TRUE){
  par(mfrow = c(1,1))
  plot(barkerLR[,idx], xlab = 'Iteration', ylab = paste(bquote(theta), idx), type = 'l', ylim = ylim)
  lines(barkerLR_alternative[,idx], col='blue')
  lines(sgldLR[,idx], col='red')
  lines(seq(50001, 100000), list_of_draws$beta[,idx], col=rgb(0, 0.39, 0, alpha = 0.5))
  abline(h=stan_means[idx], col = "dark green", lwd = 5)
  if (legend){
    legend('bottomright', legend=c("barker", "barker.2","sgld", "stan"),
         col=c("black", "blue","red", "dark green"), lty=19:19:19:19, cex=0.8)
  }
}

plot_3d = function(samples, idx_1=1, idx_2=2){
  x_c <- cut(samples[,idx_1], 50)
  y_c <- cut(samples[,idx_2], 50)
  z <- table(x_c, y_c)
  image2D(z=z, border="black")
}

histogram_combined_LR = function(idx, stan, sgld, barker, barker_alt, xlim, ylim=c(0,1), S=100000){
  dev.new()
  par(mfrow = c(1,3))
  hist(sgld[(S/2+1):S,idx], breaks=40, probability=TRUE, main = 'SGLD', xlab = 'theta', xlim = xlim, ylim=ylim)
  lines(density(stan$beta[,idx]), col='green')
  
  hist(barker[(S/2+1):S,idx], breaks=40, probability=TRUE, main = 'Barker', xlab = 'theta',  xlim = xlim, ylim=ylim)
  lines(density(stan$beta[,idx]), col='green')
  
  hist(barker_alt[(S/2+1):S,idx], breaks=40, probability=TRUE, main = 'Barker 2', xlab = 'theta',  xlim = xlim, ylim=ylim)
  lines(density(stan$beta[,idx]), col='green')
  par(mfrow = c(1,1))
}

sim_logistic_regression_arrhytmya = function(h, sigma_2s, sigma_2s_alternative, S, n, train_set, test_set, stan_means=0, stan_vars=0){
  # helper function used to simulate the sgld and stochastic barker with 
  # different value of the step size and compute the ksd, minimum ess and 
  # the log loss on a held out set
  set.seed(123)

  print(stan_means[1])
  # starting point of the algorithm
  d = ncol(train_set) - 1
  theta0 = rnorm(d, 0, 1)
  res_sgld = data.frame(t(apply(matrix(h), 1, function(x) (mcmc_logistic_regression(x, S, n, theta0, train_set, test_set, gradLogL, gradLogP, SGLD_naive, FALSE, FALSE, TRUE, type='logistic', stan_means = stan_means, stan_vars = stan_vars)))))
  res_barker = data.frame(t(apply(matrix(sigma_2s), 1, function(x) (mcmc_logistic_regression(x, S, n, theta0, train_set, test_set, gradLogL, gradLogP, barker, FALSE, FALSE, TRUE, type='logistic', stan_means = stan_means, stan_vars = stan_vars)))))
  res_barker_alternative = data.frame(t(apply(matrix(sigma_2s_alternative), 1, function(x) (mcmc_logistic_regression(x, S, n, theta0, train_set, test_set, gradLogL, gradLogP, barker, FALSE, FALSE, TRUE, version='alternative', type='logistic', stan_means = stan_means, stan_vars = stan_vars)))))
  
  col_names = c('step_size', 'ksd', 'log_loss_test', 'log_loss_train', 'min_ESS', 'median_ESS', 'bias_mean', 'bias_mean_max','bias_mean_relative' ,'bias_var', 'bias_var_max')
  col_names = c('step_size', 'ksd', 'log_loss_test', 'log_loss_train', 'min_ESS', 'median_ESS', 'bias_mean', 'bias_mean_max', 'bias_mean_relative','bias_var', 'bais_var_max')
  
  print(res_sgld)
  print(res_barker)
  print(res_barker_alternative)
  
  colnames(res_sgld) = col_names
  colnames(res_barker) = col_names
  colnames(res_barker_alternative) = col_names
  
  print(res_sgld)
  print(res_barker)
  print(res_barker_alternative)
  #plot_logistic_regression(res_barker, res_barker_alternative, res_sgld, n, ksd, moments)
  return(list("res_sgld" = res_sgld, "res_barker" = res_barker, "res_barker_alternative" = res_barker_alternative))
}

histogram_LR = function(idx, stan, sgld, barker, barker_alt, xlim, S=100000){
  hist(stan$beta[,idx], breaks=40, probability=TRUE, xlim=xlim, main = ' ', xlab=paste('theta ', idx))
  lines(density(stan$beta[,idx]), col='green')
  lines(density(sgld[(S/2+1):S,idx]), col='red')
  lines(density(barker_alt[(S/2+1):S,idx]), col='blue')
  lines(density(barker[(S/2+1):S,idx]))
}


p = function(z, grad_log_pi){
  return(1/(1+exp(-z*grad_log_pi)))
}

sim_p_LR = function(i, theta,  data){
  grad = gradLogL(data, theta) + gradLogP(theta)
  grad = grad[i]
  print(grad)
  z = seq(-5, 5, length.out=10000)
  p_s = p(z, grad)
  return(p_s)
}


sim_p_hat_LR = function(S, mini_batch_size, j, theta, data){
  grads = rep(0, S)
  p_matrix = matrix(0, ncol=10000, nrow =S)
  z = seq(-5, 5, length.out=10000)
  N = dim(data)[1]
  for(i in seq(1, S)){
    mini_batch = sample(N, mini_batch_size)
    mini_batch = data[mini_batch,]
    grad = N/mini_batch_size * gradLogL(mini_batch, theta) + gradLogP(theta)
    grad = grad[j]
    grads[i] = grad
    p_s = p(z, grad)
    p_matrix[i, ] = p_s
  }
  p_hat = colMeans(p_matrix)
  p_q_1 = colQuantiles(p_matrix, probs=0.25)
  p_q_3 = colQuantiles(p_matrix, probs=0.75)
  print(mean(grads))
  return(list('p_hat'= p_hat, 'q_1' = p_q_1, 'q_3' = p_q_3))
}

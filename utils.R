library(coda)
library(mvtnorm)
library(pracma)
library(matrixStats)
library(LaplacesDemon)
library(kernlab)
library(KSD)
library(sn)

# one-step euler discretization
one_step_langevin = function(theta, h, grad_P) {
  # PARAMS:
  # theta = the current value of the parameter
  # h = stepsize
  # grad_P = gradient of the log(p(theta))
  # it return the next state of the parameter (one step discretization of langevin dynamics)
  d = length(theta)
  theta_next = theta + h * grad_P + sqrt(2 * h) * rnorm(d, 0, 1)
  return(theta_next)
}

# Stochastic ULA
SGLD_naive = function(data, S, n, theta_0, h, grad_L, grad_p) {
  # simulate SGLD (stochastic variation of ULA)
  # theta evolves approximately through the one step-discretization of the langevin dynamics
  # at each step the gradient is estimated using a minibath
  # instead of the burn in, the first half of iterations is used to find a local mode via SGD 
  # that is used a initializer of the SGLD
  
  # PARAMS:
  # S = number of steps
  # n = mini-bacth size
  # theta_0 = starting point
  # h = stepsize
  # eta = initial learning rate of SGD
  # grad_L = gradient of the log - likelihood
  # grad_p = gradient of the log - prior
  
  
  # OUTPUT
  # (theta_1, ..., theta_S) = the entire path 
  
  # size of the dataset
  data = data.frame(data)
  N = dim(data)[1]
  # dimensionality of the parameter
  d = length(theta_0)
  thetas = matrix(NA, nrow = S, ncol = d)
  # initialize starting point
  thetas[1,] = theta_0
  # theta_new = is updated at each iter and represent the most recent value of the parameter
  
  for (i in 2:S) { 
    # compute unbiased estimate of the gradient
    # sample without replacement n out N units
    mini_batch = sample(N, n)
    mini_batch = data[mini_batch,]
    # etimate of the gradient as sum of the gradient of potential of the prior and of the unbiased estimated of gradient of the potential of likelihood
    grad_U1 = grad_p(theta_0) + N/n * grad_L(mini_batch, theta_0)

    theta_new = one_step_langevin(theta_0, h, grad_U1)
    
    theta_0 = theta_new
    if (any(is.na(theta_new))){
      print("NA value")
      break
    }
    thetas[i,] = theta_0
  }
  return(thetas)
}



oneStepBarker = function(theta, sigma_2, grad_P, version = 'standard') {
  # this function computes the next state where the process follows the barker dynamics
  # PARAMS:
  # theta = the current value of the parameter
  # sigma_2 = variance of the normal distribution
  # grad_P =  gradient of the log(p(theta))
  # it return the next state of the parameter
  d = length(theta)
  sigma = sqrt(sigma_2)
  if (version == 'standard') {z = sigma * rnorm(d)}
  else  {z = rnorm(d, sigma, 0.1 * sigma)}
  b = 1 / (1 + exp(-z * grad_P))
  u = runif(d)
  theta_next = theta + z * (2 * as.integer(u < b) -1)
  return(theta_next)
}

# stochastic barker dinamics
barker = function(data, S, n, theta_0, sigma, grad_L, grad_p, version='standard'){
  
  data = data.frame(data)
  N = dim(data)[1]
  # dimensionality of the parameter
  d = length(theta_0)
  thetas = matrix(NA, nrow = S, ncol = d)
  # initialize starting point
  thetas[1,] = theta_0
  # theta_new = is updated at each iter and represent the most recent value of the parameter
  
  for (i in 2:S) { 
    mini_batch = sample(N, n)
    mini_batch = data[mini_batch,]
    grad_U = grad_p(theta_0) + N/n * grad_L(mini_batch, theta_0)
    theta_new = oneStepBarker(theta_0, sigma, grad_U, version)
    
    theta_0 = theta_new
    if (any(is.na(theta_new))){
      print("NA value")
      break
    }
    thetas[i,] = theta_0
  }
  return(thetas)
}



# 1d simulation helper functions
grad_p1 = function(theta) {
  # compute the gradient of the prior 
  return(- (theta - mu0)  /  tau0 **2)
}  


grad_U1 = function(data, theta) {
  # compute estimate of the gradient of the log posterior (normal model with conjugate prior)
  return(sum((data - theta) / (sigma ** 2)))
}  


run_sim_1d = function(h, sigmas, sigmas_alternative, S, N, n, theta_0, tau_0, sigma, grad_L, grad_P, ksd = FALSE) {
  
  set.seed(123)
  generated_data = gen_data_1d(N, 0, tau_0, sigma)
  sigma_n = generated_data$var
  mu_n = generated_data$mean
  data = generated_data$data
  
  res_sgld = data.frame(t(apply(matrix(h), 1, function(x) (sim_normal_model(x, S, n, theta0, data, grad_L, grad_P, mu_n, sigma_n, SGLD_naive, ksd=ksd)))))
  res_barker = data.frame(t(apply(matrix(sigmas), 1, function(x) (sim_normal_model(x, S, n, theta0, data, grad_L, grad_P, mu_n, sigma_n, barker, ksd=ksd)))))
  res_barker_alternative = data.frame(t(apply(matrix(sigmas_alternative), 1, function(x) (sim_normal_model(x, S, n, theta0, data, grad_L, grad_P, mu_n, sigma_n, barker, version='alternative', ksd=ksd)))))
                                                                              
  col_names = c("step_size", "bias_mean", "bias_mean_relative", "bias_var","ESS", "ksd")
  colnames(res_sgld) = col_names
  colnames(res_barker) = col_names
  colnames(res_barker_alternative) = col_names
  
  print(res_sgld)
  print(res_barker)
  print(res_barker_alternative)
  plot_simulation_2(res_barker, res_barker_alternative, res_sgld, n, type = 'normal model', legend=FALSE)
  
  return(list('res_sgld' = res_sgld, 'res_barker' = res_barker, 'res_barker_alternative' = res_barker_alternative))
}

sim_normal_model = function(step_size, S, n, theta_0, data, grad_L, grad_P, mu_n, sigma_n, algorithm, version=NULL, ksd=FALSE){
  print(step_size)
  if (is.null(version)){
    sim = algorithm(data, S, n, theta_0, step_size, grad_L, grad_P)
  }
  else{
    sim = algorithm(data, S, n, theta_0, step_size, grad_L, grad_P, version)
  }
  
  
  if (ksd == FALSE){
    ksd_sim = 0
  }
  else{
    grad = data.frame(gradPost(sim[seq(floor(S/2), S, length.out = 10001)]))      
    ksd_sim = imqKSD(data.frame(sim[seq(floor(S/2), S, length.out = 10001)]), grad, 1, 0.5)
  }
  
  sim = sim[(floor(S/2)+1):S]
  ESS_sim= ESS(sim)
  ESS_sim_2 = effectiveSize(sim)
  
  return(c(step_size, abs(mean(sim)-mu_n), abs(mean(sim)-mu_n)/sqrt(sigma_n), abs(std(sim)**2-sigma_n), ESS_sim, ESS_sim_2, ksd_sim))
}



gen_data_1d = function(N, mu.0, tau.0, sigma) {
  # generate data from normal model with conjugate prior on the mean parameter (variances are known)
  # PARAMS:
  # mu.0 = prior mean
  # tau.0 = prior sd
  # sigma = model sd
  # OUTPUT:
  # returns a list with posterior params, the value of theta and data generated
  theta_0 = rnorm(1, mu.0, tau.0) 
  data = rnorm(N, theta_0, sigma) 
  x_bar = mean(data)
  sigma_n = 1/(1/(tau.0)**2 + N/(sigma**2)) # posterior variance
  mu_n = x_bar * (tau.0)**2 /((sigma**2) /N + tau.0**2) + mu.0 * sigma**2 /(N*((sigma**2)/N + tau.0**2)) # posterior mean
  return(list('mean' = mu_n, 'var' = sigma_n, "theta" = theta_0,"data" = data))
}



artificial_noise = function(algorithm, gradLogP, h, sigma, S, theta_0, version = NaN){
  d = length(theta_0)
  samples = matrix(NA, ncol=d, nrow=S)
  samples[1,] = theta_0
  for (i in seq(2, S)){
    grad = gradLogP(theta_0) + sigma * rnorm(d)
    if (is.na(version)==TRUE){ 
      theta_0 = algorithm(theta_0, h, grad)
      }
    else {theta_0 = algorithm(theta_0, h, grad, version = version)}
    samples[i,] = theta_0
  }
  return(samples)
}

sim_artificial_noise = function(sigma_2s, sigma_2s_alternative, hs, grad_log_p, sigma, S, theta_0){
  
  res_sgld = data.frame(t(apply(matrix(h), 1, function(x) (mcmc_artificial_noise(SGLD_naive, grad_log_p, x, sigma, S, theta_0)))))
  res_barker = data.frame(t(apply(matrix(sigma_2s), 1, function(x) (mcmc_artificial_noise(barker, grad_log_p, x, sigma, S, theta_0)))))
  res_barker_alternative = data.frame(t(apply(matrix(sigma_2s_alternative), 1, function(x) (mcmc_artificial_noise(barker, grad_log_p, x, sigma, S, theta_0, version='alternative')))))
  
  col_names = c("step_size", "mean", "var", "ESS")
  colnames(res_sgld) = col_names
  colnames(res_barker) = col_names
  colnames(res_barker_alternative) = col_names
  plot_simulation_1(res_barker, res_barker_alternative, res_sgld, sigma)
  
  return(list('res_sgld' = res_sgld, 'res_barker' = res_barker, 'res_barker_alternative' = res_barker_alternative))
}


mcmc_artificial_noise = function(algorithm, grad_log_p, h, sigma, S, theta_0, version=NULL, type='Normal', alpha=1){
  if (is.null(version)){
    chain = artificial_noise(algorithm, grad_log_p, h, sigma, S, theta_0)
  }
  else{
    chain = artificial_noise(algorithm, grad_log_p, h, sigma, S, theta_0, version='alternative')
  }
  chain = chain[(floor(S/2)+1):S,]
  
  mean_samples = mean(chain)
  if (algorithm == SGLD_naive){
    std_samples = sqrt((sigma^2 * h + 2)/(2 - h))
  }
  std_samples = std(chain)
  ess_samples = ESS(chain)
  if (type == 'Normal'){
    return(c(h, mean_samples, std_samples**2, ess_samples))
  }
  else{
    mean_skew = 1 * alpha / (sqrt(1 + alpha**2)) * sqrt(2 / pi)
    var_skew = 1 - 2*alpha**2/((1+alpha**2)*pi)
    bias_mean = abs(mean_samples - mean_skew)
    bias_var = abs(std_samples**2 - var_skew)
  }
  return(c(h, bias_mean, bias_var, ess_samples))
}

sim_artificial_noise_2 = function(sigma_2s, sigma_2s_alternative, hs, alpha, sigma, S, theta_0){
  
  res_sgld = data.frame(t(apply(matrix(h), 1, function(x) (mcmc_artificial_noise(SGLD_naive, grad_skew_normal, x, sigma, S, theta_0, type='skew', alpha=alpha)))))
  res_barker = data.frame(t(apply(matrix(sigma_2s), 1, function(x) (mcmc_artificial_noise(barker, grad_skew_normal, x, sigma, S, theta_0, type='skew', alpha=alpha)))))
  res_barker_alternative = data.frame(t(apply(matrix(sigma_2s_alternative), 1, function(x) (mcmc_artificial_noise(barker, grad_skew_normal, x, sigma, S, theta_0, version='alternative', type='skew', alpha=alpha)))))
  
  col_names = c("step_size", "bias_mean", "bias_var", "ESS")
  colnames(res_sgld) = col_names
  colnames(res_barker) = col_names
  colnames(res_barker_alternative) = col_names
  plot_simulation_2(res_barker, res_barker_alternative, res_sgld, sigma, legend=FALSE)
  
  return(list('res_sgld' = res_sgld, 'res_barker' = res_barker, 'res_barker_alternative' = res_barker_alternative))
}

grad_normal_1 = function(theta){
  return(-theta)
}

grad_skew_normal = function(theta, alpha){
  return(-theta + alpha * dnorm(alpha * theta) / pnorm(alpha * theta))
}



plot_simulation_1 = function(df.barker, df.barker.alt, df.sgld, sigma, xlim=15000, ylim=1){
  par(mfrow = c(1,1))
  print(sigma)
  x_min = min(min(df.barker[,'ESS']), min(df.sgld[,'ESS']))
  x_max = max(max(df.barker[,'ESS']), max(df.sgld[,'ESS']))
  y_min = min(min(df.barker[,'var']), min(df.sgld[,'theoretical_var'])) -1
  y_max = max(max(df.barker[,'var']), max(df.sgld[,'theoretical_var']))
  plot(df.barker[,'ESS'], df.barker[,'var'] - 1, #main = bquote(sigma == .(sigma)) ,
       xlab = 'ESS', ylab = 'Bias var', pch = 19, xlim = c(0, xlim), ylim = c(0, ylim))
  points(df.barker.alt[,'ESS'], df.barker.alt[,'var'] - 1,  pch = 19, col = 'blue')
  points(df.sgld[,'ESS'], df.sgld[,'theoretical_var'] - 1, pch = 19, col = 'red')
  
}

plot_simulation_2 = function(df_barker, df_barker_alternative, df_sgld, sigma, type = 'skew', legend = TRUE, xlim=15000, ylim1=0.42, ylim2=0.42){
  par(mfrow = c(1,1))
  x_min = min(min(df_barker[,'ESS']), min(df_sgld[,'ESS']))
  x_max = max(max(df_barker[,'ESS']), max(df_sgld[,'ESS']))
  y_min_1 = min(min(df_barker[,'bias_mean']), min(df_sgld[,'bias_mean']))
  y_max_1 = max(max(df_barker[,'bias_mean']), max(df_sgld[,'bias_mean']))
  y_min_2 = min(min(df_barker[,'bias_var']), min(df_sgld[,'bias_var']))
  y_max_2 = max(max(df_barker[,'bias_var']), max(df_sgld[,'bias_var']))
  plot(df_barker[,'ESS'], df_barker[,'bias_mean'], xlab = 'ESS', ylab = 'Bias mean', pch = 19, xlim = c(0, xlim), ylim = c(0, ylim1))
  points(df_barker_alternative[,'ESS'], df_barker_alternative[,'bias_mean'], pch = 19, col = 'blue')
  points(df_sgld[,'ESS'], df_sgld[,'bias_mean'], pch = 19, col = 'red')
  if (legend) {
  legend('topleft', legend=c("barker", "barker.2", "sgld"),
         col=c("black", "blue","red"), lty=19:19:19, cex=0.8)
  }
  plot(df_barker[,'ESS'], df_barker[,'bias_var'], xlab = 'ESS', ylab = 'Bias variance', pch = 19, xlim = c(0, xlim), ylim = c(0, ylim2))
  points(df_barker_alternative[,'ESS'], df_barker_alternative[,'bias_var'] , pch = 19, col = 'blue')
  points(df_sgld[,'ESS'], df_sgld[,'bias_var'] , pch = 19, col = 'red')
  if (legend){
  legend('topleft', legend=c("barker", "barker.2","sgld"),
         col=c("black", "blue","red"), lty=19:19:19, cex=0.8)
  }
}


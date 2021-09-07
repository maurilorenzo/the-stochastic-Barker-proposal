####### 6.3 Bayesian probabilistic matrix factorization
# this script contains helper function for the 6.3 section of the review paper
source("../utils.R")

library(mvtnorm)
library(zoo)

# init parameters 
initParams_bpfm = function(N, M, d, mu0, a0, b0){
  #sigma_U = rgamma(d, a0, b0)
  #sigma_V = rgamma(d, a0, b0)
  sigma_U = rnorm(d, -.1, 1)
  sigma_V = rnorm(d, -.1, .1)
  mu_U = rnorm(d, mu0, 1)
  mu_V = rnorm(d, mu0,1)
  U = t(rmvnorm(N, mu_U, 1 * diag(d)))
  V = t(rmvnorm(M, mu_V, 1 * diag(d)))
  return(list('sigma_U' = sigma_U, 'sigma_V' = sigma_V, 'mu_U' = mu_U, 'mu_V' = mu_V, 'U' = U, 'V' = V))
}
  


plot_gradient=function(samples, window=50, S=5000, div=100){
  plot(seq(window, S)/div, log(rollmean(samples$grads_U, k=window)), ylim=c(0, 15), ylab='log gradient', xlab='Epoch', type='l')
  points(seq(window,S)/div, log(rollmean(samples$grads_V, k=window)), type='l')
  points(seq(window,S)/div, log(rollmean(samples$grads_mu_U, k=window)), col='blue', type='l')
  points(seq(window,S)/div, log(rollmean(samples$grads_mu_V, k=window)), type='l', col='blue')
  points(seq(window,S)/div, log(rollmean(samples$grads_lambda_U, k=window)), col='red', type='l')
  points(seq(window,S)/div, log(rollmean(samples$grads_lambda_V, k=window)), type='l', col='red')
  plot(seq((S-99), S)/div, log(samples$grads_U[(S-99): S]), ylim=c(0, 15), ylab='log gradient', xlab='Epoch')
  points(seq((S-99),S)/div, log(samples$grads_V[(S-99): S]), pch=2)
  points(seq((S-99),S)/div, log(samples$grads_mu_U[(S-99): S]), col='blue')
  points(seq((S-99),S)/div, log(samples$grads_mu_V[(S-99): S]), pch=2, col='blue')
  points(seq((S-99),S)/div, log(samples$grads_lambda_U[(S-99): S]), col='red')
  points(seq((S-99),S)/div, log(samples$grads_lambda_V[(S-99): S]), pch=2, col='red')
}



plot_rmse=function(mf_sgld, mf_barker, mf_barker_alternative, div, train=FALSE){
  plot(seq((S/2+1),S)/div, mf_sgld$rMSEs_avg, col='red', xlab = 'Epoch', ylab = 'RMSE', ylim = c(0.8, 1.6), type='l')
  lines(seq((S/2+1),S)/div, mf_barker_alternative$rMSEs_avg, col='blue')
  lines(seq((S/2+1),S)/div, mf_barker$rMSEs)
  plot(seq(1,S)/div, mf_sgld$rMSEs, col='red', xlab = 'Epoch', ylab = 'RMSE', ylim = c(0.8, 1.6))
  lines(seq(1,S)/div, mf_barker_alternative$rMSEs, col='blue')
  lines(seq(1,S)/div, mf_barker$rMSEs)
}


# define helper function used to compute the gradient
f1 = function(x, U, V) {return(return(x['R'] - sum(U[,x['uid']] * V[, x['iid']])))}


f2_c = function(x, df,  U, V, mu_U, sigma_U, alpha, N_data, N_batch){
  index = df['iid'][df['uid'] == x]
  if (length(index) == 1){
    return(N_data / N_batch * alpha * t(df['res'][df['uid'] == x] * t(V[,index])))
  }
  return(rowSums(N_data / N_batch * alpha * t(df['res'][df['uid'] == x] * t(V[,index]))))
}

f2_b = function(grad_matrix, U, mu_U, sigma_U){
  return(grad_matrix - exp(sigma_U) * (U - mu_U))
}


f3_c = function(x, df, U, V, mu_V, sigma_V, alpha, N_data, N_batch){
  index = df['uid'][df['iid'] == x]
  n_j = length(index)
  if (length(index) == 1){
    return(N_data / N_batch * alpha * t(df['res'][df['iid'] == x] *  t(U[,index])))
  }
  return(rowSums(N_data / N_batch * alpha * t(df['res'][df['iid'] == x] *  t(U[,index]))))
}

f3_b = function(grad_matrix, V, mu_V, sigma_V) {
  return(grad_matrix - exp(sigma_V) * (V - mu_V))
}


grad_lambda_U_f = function(sigma_U, Us, mu_U, N,alpha0, beta){
  rs1.u = rowSums((Us - mu_U)**2)
  return(((N + 1)/2 + alpha0)  - ((mu_U**2 + rs1.u) / 2 + beta) * exp(sigma_U))
}

grad_lambda_V_f = function(sigma_V, Vs, mu_V, M, alpha0, beta){
  rs1.v = rowSums((Vs - mu_V)**2)
  return(((M + 1)/2 + alpha0)  - ((mu_V**2 + rs1.v ) / 2 + beta) * exp(sigma_V))
}

grad_mu_U_f = function(sigma_U, mu_U, Us, N){
  rs2.u = rowSums(Us)
  return(exp(sigma_U) * (rs2.u - (N + 1) * mu_U))
}

grad_mu_V_f = function(sigma_V, mu_V, Vs, M){
  rs2.v = rowSums(Vs)
  return(exp(sigma_V) * (rs2.v - (M + 1) * mu_V))
}

p = function(z, grad_log_pi){
  return(1/(1+exp(-z*grad_log_pi)))
}


BPMF = function(h, df, test, S, d, miniBatchSize = 800, algo = 'sgld',
                version = 'standard', burn_in = 0, grad = TRUE, 
                train_accuracy = FALSE, return_samples=FALSE){
  
  h_1 = h[1]
  h_2 = h[2]
  h_3 = h[3]
  
  
  alpha = 3
  alpha0 = 1
  beta = 5
  
  N = length(unique(df$uid))
  M = length(unique(df$iid))
  N_u = max(df$uid)
  N_i = max(df$iid)
  
  d = 20
  N_data = dim(df)[1]
  
  print(d)
  set.seed(123)
  params = initParams_bpfm(N_u, N_i, d, 0, 1, 5)
  U = params$U
  V = params$V
  mu_U = params$mu_U
  mu_V = params$mu_V
  sigma_U = params$sigma_U
  sigma_V = params$sigma_V
  
  rMSEs_avg = rep(NA, S - burn_in)
  rMSEs = rep(NA, S)
  rMSEs_avg_train = rep(NA, S - burn_in)
  rMSEs_train = rep(NA, S)
  
  
  if (grad) {
    grads_U = rep(NA, S)
    grads_V = rep(NA, S)
    grads_mu_U = rep(NA, S)
    grads_mu_V = rep(NA, S)
    grads_lambda_U = rep(NA, S)
    grads_lambda_V = rep(NA, S)
  }
  else{
    grads_U = 0
    grads_V = 0
    grads_mu_U = 0
    grads_mu_V = 0
    grads_lambda_U = 0
    grads_lambda_V = 0
  }
  
  L = dim(df)[1]
  # miniBatchSize = floor(0.1 * L)
  # standardize the rating
  meanR = mean(df$R)
  df$R = df$R - meanR
  
  
  for (i in seq(1, S)){
    minibatch = sample(L, miniBatchSize)
    minibatch = df[minibatch,]
    uis = sort(unique(minibatch$uid))
    vis = sort(unique(minibatch$iid))
    n = length(uis)
    m = length(vis)
    Us = U[, uis]
    Vs = V[, vis]
    uis.df = data.frame(uis)
    vis.df = data.frame(vis)
    # compute R - UTV
    minibatch['res'] = apply(minibatch, 1, function(x) (f1(x, U, V)))
    
    grad_U = matrix(0, ncol = N_u, nrow = d)
    grad_V = matrix(0, ncol = N_i, nrow = d)
    
    grad_sigma_U = grad_lambda_U_f(sigma_U, U, mu_U, N, alpha0, beta)
    grad_sigma_V = grad_lambda_V_f(sigma_V, V, mu_V, M, alpha0, beta)
    
    grad_mu_U = grad_mu_U_f(sigma_U, mu_U, U, N)
    grad_mu_V =  grad_mu_V_f(sigma_V, mu_V, V, M)
    

    grad_Us = apply(uis.df, 1, function(x1) (f2_c(x1, minibatch, U, V, mu_U, sigma_U, alpha, N_data, miniBatchSize)))
    grad_Vs = apply(vis.df, 1, function(x1) (f3_c(x1, minibatch, U, V, mu_V, sigma_V, alpha, N_data, miniBatchSize)))
    
    grad_U[, uis] = grad_Us
    grad_U = f2_b(grad_U, U, mu_U, sigma_U)

    grad_V[, vis] = grad_Vs
    grad_V = f3_b(grad_V, V, mu_V, sigma_V)
    
    
    if (grad){
      grads_lambda_U[i] = max(abs(grad_sigma_U))
      grads_lambda_V[i] = max(abs(grad_sigma_V))
      grads_mu_U[i] = max(abs(grad_mu_U))
      grads_mu_V[i] = max(abs(grad_mu_V))
      grads_U[i] = max(abs(grad_Us))
      grads_V[i] = max(abs(grad_Vs))
    }
    
    if (algo == 'sgld') {
      
      # one step langevin dynamics
      sigma_U = sigma_U + h_3/2 * grad_sigma_U +  sqrt(h_3) * rnorm(d)
      sigma_V = sigma_V + h_3/2 * grad_sigma_V +  sqrt(h_3) * rnorm(d)
      mu_U = mu_U + h_2/2 * grad_mu_U +  sqrt(h_2) * rnorm(d)
      mu_V = mu_V + h_2/2 * grad_mu_V +  sqrt(h_2) * rnorm(d)
      U = U + sqrt(h_1) * rnorm(N_u * d) + h_1/2 * grad_U
      #U[, uis] = U[,uis] + h_1/2 * grad_Us
      V = V + sqrt(h_1) * rnorm(N_i * d) + h_1/2 * grad_V
      #V[, vis] = V[,vis] + h_1/2 * grad_Vs
    }
    
    else {
      # one step barker
      sigma_U = oneStepBarker(sigma_U, h_3, grad_sigma_U, version)
      sigma_V = oneStepBarker(sigma_V, h_3, grad_sigma_V, version)
      mu_U = oneStepBarker(mu_U, h_2, grad_mu_U, version)
      mu_V = oneStepBarker(mu_V, h_2, grad_mu_V, version)
      U = matrix(oneStepBarker(as.vector(U), h_1, as.vector(grad_U), version), ncol = N_u)
      V = matrix(oneStepBarker(as.vector(V), h_1, as.vector(grad_V), version), ncol = N_i)
    }
    
    # predict new data and compute MSE
    current_pred = apply(test, 1, function(x) f4(x, U, V, meanR))
    rMSE = sqrt(mean((test$R - current_pred)**2))
    rMSEs[i] = rMSE
    if(train_accuracy){
      current_pred_train = apply(df, 1, function(x) f4(x, U, V, meanR))
      rMSE_train = sqrt(mean((df$R - current_pred_train)**2))
      rMSEs_train[i] = rMSE_train
    }
    if (i%%100 == 0) {
      cat("iter: ", i, "rMSE: ",rMSE)
      print("")
    }
    if (i > burn_in){
      if (i == burn_in +1) {
        pred = current_pred
        
      }
      else {
        pred = (current_pred + (i - burn_in - 1)*old_pred) / (i - burn_in)
      }
      rMSE = sqrt(mean((test$R - pred)**2))
      #print(rMSE)
      if (i%%100 == 0) {
        cat("iter: ", i, "rMSE avg: ",rMSE)
        print("")
      }
      old_pred = pred
      #print(current_pred)
      rMSEs_avg[i - burn_in] = rMSE
      
      
      if (is.na(rMSE)) {
        print("NaN values")
        
        return(rMSEs)
      }
    }
  }
  if(return_samples){
    samples = list('U'=U, 'V'=V, 'mu_U'=mu_U, 'mu_V'=mu_V, 'lambda_U'=sigma_U, 'lambda_V'=sigma_V)
  }
  else{
    samples = 0
  }
  
  return(list('rMSEs' = rMSEs, 'rMSEs_avg' = rMSEs_avg, 'rMSEs_train' = rMSEs_train, 'grads_U' = grads_U, 'grads_V' = grads_V, 'grads_mu_U' = grads_mu_U, 'grads_mu_V' = grads_mu_V, 'grads_lambda_U' = grads_lambda_U, 'grads_lambda_V' = grads_lambda_V, 'samples'=samples))
}

sim_p = function(u_i, v_j, df, U, V, mu_U, sigma_U, alpha){
  meanR = mean(df$R)
  df$R = df$R - meanR
  # compute residuals  
  df['res'] = apply(df, 1, function(x) (f1(x, U, V)))
  grad_U_i = f2_c(u_i, df, U, V, mu_U, sigma_U, alpha, 80000, 80000)
  grad_prior = f2_b(0, U, mu_U, sigma_U)
  grad_U_i = grad_U_i + grad_prior[,u_i]
  print(length(grad_U_i))
  print(grad_U_i)
  grad_U_ij = grad_U_i[v_j]
  print(grad_U_ij)
  z = seq(-5, 5, length.out=10000)
  p_s = p(z, grad_U_ij)
  return(p_s)
}


sim_p_hat = function(S, mini_batch_size, u_i, v_j, df, U, V, mu_U, sigma_U, alpha){
  L = dim(df)[1]
  meanR = mean(df$R)
  df$R = df$R - meanR
  p_matrix = matrix(0, ncol=10000, nrow=S)
  grad = rep(0, S)
  z = seq(-5, 5, length.out=10000)
  
  for(i in seq(1, S)){
    mini_batch = sample(L, mini_batch_size)
    mini_batch = df[mini_batch,]
    u_is = sort(unique(mini_batch$uid))
    n = length(uis)
    m = length(vis)
    #Us = U[, uis]
    #Vs = V[, vis]
    uis_df = data.frame(u_is)
    # compute R - UTV
    mini_batch['res'] = apply(mini_batch, 1, function(x) (f1(x, U, V)))
    if (u_i %in% u_is){
      grad_U_i = f2_c(u_i, mini_batch, U, V, mu_U, sigma_U, alpha, 80000, mini_batch_size)
    }
    else{
      grad_U_i = rep(0, d)
    }
    
    grad_prior = f2_b(0, U, mu_U, sigma_U)
    grad_U_i = grad_U_i + grad_prior[,u_i]
    grad_U_ij = grad_U_i[v_j]
    grad[i] = grad_U_ij
    p_s = p(z, grad_U_ij)
    p_matrix[i, ] = p_s
  }
  p_hat = colMeans(p_matrix)
  sd = colSds(p_matrix)
  quantiles = colQuantiles(p_matrix)
  print(mean(grad))
  return(list('p_hat'=p_hat, 'min'=quantiles[,1], 'q_1'=quantiles[,2], 'median'=quantiles[,3],'q_3'=quantiles[,4], 'max'=quantiles[,5], 'sd'=sd, 'p_matrix'=p_matrix))
}


plot_BPMF = function(barker, barker_alt, sgld, S, div, train=FALSE){
  plot(seq(1, S)/div, barker$rMSEs, ylab='RMSE', xlab='Epoch', type='l', ylim = c(0.8, 2.2))
  points(seq(1,  S)/div, barker_alt$rMSEs, pch=0, col='blue', type='l')
  points(seq(1, S)/div, sgld$rMSEs, pch=0, col='red', type='l')
  plot(seq((S/2+1), S)/div, barker$rMSEs_avg, ylab='RMSE', xlab='Epoch', type='l', ylim = c(0.8, 2.2))
  points(seq((S/2+1), S)/div, barker_alt$rMSEs_avg, type='l', col='blue')
  points(seq((S/2+1), S)/div, sgld$rMSEs_avg, type='l', col='red')
  if(train){
    plot(seq(1,S)/div, sgld$rMSEs_train, col='red', xlab = 'Epoch', ylab = 'RMSE', ylim = c(0.8, 1.6))
    lines(seq(1,S)/div, barker_alt$rMSEs_train, col='blue')
    lines(seq(1,S)/div, barker$rMSEs_train)
  }
}


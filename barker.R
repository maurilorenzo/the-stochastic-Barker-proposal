# barker proposal
library(sn)

x = seq(-5, 5, by=0.001)

sigma = 1

y = dnorm(x, 0, sigma)
#y_1 = dsn(x, alpha=1)
#y_2 = dsn(x, alpha=3)
#y_3 = dsn(x, alpha=5)
#y_4 = dsn(x, alpha=10)



z_base = dnorm(x, sigma, 0.1*sigma)
z_base_1 = dnorm(x, -sigma, 0.1*sigma)
z = 0.5 * z_base + 0.5 * z_base_1

p_1 = 1/(1+exp(-1*x))
p_2 = 1/(1+exp(-2*x))
p_3 = 1/(1+exp(-5*x))
p_4 = 1/(1+exp(-10*x))

y_1 = 2*y*p_1
y_2 = 2*y*p_2
y_3 = 2*y*p_3
y_4 = 2*y*p_4
z_1 = (1- rev(p_1))*z_base_1 + p_1*z_base
z_2 = (1- rev(p_2))*z_base_1 + p_2*z_base
z_3 = (1- rev(p_3))*z_base_1 + p_3*z_base
z_4 = (1- rev(p_4))*z_base_1 + p_4*z_base

z_1 = 2*z*p_1
z_2 = 2*z*p_2
z_3 = 2*z*p_3



plot(x, y, type='l', ylim=c(0, 1), xlab=expression(theta), ylab='pdf of proposal')
lines(x, y_1, col='red', lwd=0.5)
lines(x, y_2, col='red', lwd=1.5)
lines(x, y_3, col='red', lwd=2)
#lines(x, y_4, col='red', lwd=2.5)
lines(x, y)
plot(x, z, type='l', ylim=c(0, 4.3), xlim=c(-2,2), xlab=expression(theta), ylab='pdf of proposal')
lines(x, z_1, col='red', lwd=0.5)
lines(x, z_2, col='red', lwd=1.5)
lines(x, z_3, col='red', lwd=2)
lines(x, z_4, col='red', lwd=2.3)
lines(x, z)



dev.new()
par(mfrow=c(1,2))
plot(x, dnorm(x, 0, 0.1))
lines(x, dnorm(x, 0, 0.2))
lines(x, dnorm(x, 0, 0.3))
lines(x, dnorm(x, 0, 0.5))
lines(x, dnorm(x, 0, 1))

sigma = 0.3
x = seq(-1, 5, by=0.0001)
plot(x, dnorm(x, sigma, 0.1*sigma), type='l', xlim=c(0, 2), ylim=c(0, 14), lwd=0.3)
sigma = 0.4
lines(x, dnorm(x, sigma, 0.1*sigma), lwd=0.5)
sigma = 0.5
lines(x, dnorm(x, sigma, 0.1*sigma), lwd=1)
sigma = 0.75
lines(x, dnorm(x, sigma, 0.1*sigma), lwd=1.5)
sigma = 1
lines(x, dnorm(x, sigma, 0.1*sigma), lwd=2)


# simulation p_hat vs p
p = function(z, grad_log_pi){
  return(1/(1+exp(-z*grad_log_pi)))
}

p_hat = function(z, grad_log_pi, eta){
  return(1/(1-exp(-z*(grad_log_pi + eta))))
}


estimate_p_hat= function(z, grad_log_pi, eta_sd){
  etas = rnorm(1000, 0, eta_sd)
  p_hat_s = p_hat(z, grad_log_pi, etas)
  e_p_hat = mean(p_hat_s)
  sd_p_hat = std(p_hat_s)
  return(c(e_p_hat, sd_p_hat))
}

sim_p_hat = function(z_s_1, z_s_2, grad_log_pi, eta_sd){
  p_s = p(z_s_1, grad_log_p)
  e_p_hat_s = rep(0, length(z_s_2))
  sd_p_hat_s = rep(0, length(z_s_2))
  for(i in seq(1, length(z_s_2))){
    z = z_s_2[i]
    p_hat = estimate_p_hat(z, grad_log_pi, eta_sd)
    e_p_hat_s[i] = p_hat[1]
    sd_p_hat_s[i] = p_hat[2]
  }
  return(list('p'=p_s, 'p_hat'=e_p_hat_s, 'sd_p'=sd_p_hat_s))
}

sim_p_hat_sigmas =  function(z_s_1, z_s_2, grad_log_pi, eta_sd_s){
  
  p_s = matrix(0, ncol=length(z_s_1), nrow=length(eta_sd_s))
  p_hat_s = matrix(0, ncol=length(z_s_2), nrow=length(eta_sd_s))
  sd_p_hat_s = matrix(0, ncol=length(z_s_2), nrow=length(eta_sd_s))
  
  for(i in seq(1, legnth(eta_sd_s))){
    sigma = sta_sd_s[i]
    sim = sim_p_hat(z_s_1, z_s_2, grad_log_pi, sigma)
    p_s[i, ] = sim$p
    p_hat_s[i, ] = sim$p_hat
    sd_p_hat_s[i, ] = sim$sd_p_hat
  }
}


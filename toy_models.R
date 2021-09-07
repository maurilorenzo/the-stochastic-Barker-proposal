source("utils.R")


# set seed
set.seed(123)

# 1st sim (ch. 5.1.1 in thesis) 
# N(0, 1) with artificial noise in the gradient estimation (N(0, sigma^2))
# for different configuration of the hyperparameters, analyze the behaviour of mixing
# and accuracy (ESS vs Bias of Variance) 
S = 200000
theta_0 = 0

sigma = 0.02
hs = seq(0.001, 0.26, length.out=30)
sigma_2s = seq(0.003, 0.8, length.out=30)
sigma_2s_alternative = seq(0.003, 0.62, length.out=30)
sim_1_1 =  sim_artificial_noise(sigma_2s, sigma_2s_alternative,  hs, grad_normal_1, sigma, S, theta_0)


sigma = 0.1
hs = seq(0.001, 0.26, length.out=30)
sigma_2s = seq(0.003, 0.8, length.out=30)
sigma_2s_alternative = seq(0.003, 0.7, length.out=30)
sim_2_1 =  sim_artificial_noise(sigma_2s, sigma_2s_alternative, hs, grad_normal_1, sigma, S, theta_0)


sigma = 1
hs = seq(0.001, 0.26, length.out=30)
sigma_2s = seq(0.003, 0.95, length.out=30)
sigma_2s_alternative = seq(0.003, 0.72, length.out=30)
sim_3_1 =  sim_artificial_noise(sigma_2s,sigma_2s_alternative, hs, grad_normal_1, sigma, S, theta_0)


sigma = 5
hs = seq(0.001, 0.28, length.out=30)
sigma_2s = exp(seq(log(0.003, exp(1)), log(3.2, exp(1)), length.out = 40))
sigma_2s_alternative = exp(seq(log(0.003, exp(1)), log(3, exp(1)), length.out = 40))
sim_4_1 =  sim_artificial_noise(sigma_2s, sigma_2s_alternative, hs, grad_normal_1, sigma, S, theta_0)


sigma = 10
hs = seq(0.001, 0.27, length.out=30)
sigma_2s = exp(seq(log(0.003, exp(1)), log(11, exp(1)), length.out = 40))
sigma_2s_alternative = exp(seq(log(0.003, exp(1)), log(10, exp(1)), length.out = 40))
sim_5_1 =  sim_artificial_noise(sigma_2s, sigma_2s_alternative, hs, grad_normal_1, sigma, S, theta_0)






#######################################################################################
# 2nd sim (ch. 5.1.2 in the thesis)
# skew normal f(x) = 2*phi(x)*Phi(alpha*x)
# grad log(f(x)) = grad_log(phi(x)) + alpha * 1 / (Phi(alpha*x)) * phi(alpha * x)

# 2.a
# alpha = 2 (intermediate skewness)
alpha = 2
mean_skew = 1 * alpha / (sqrt(1 + alpha**2)) * sqrt(2 / pi)
var_skew = 1 - 2*alpha**2/((1+alpha**2)*pi)
mean_skew
var_skew

S = 200000
theta_0 = 0

hs = seq(0.005, 0.13, length.out=40)
sigma_2s = seq(0.011, 0.37, length.out = 40)
sigma_2s_alternative = seq(0.01, 0.34, length.out = 40)
sigma = 0.02
sim_skew_1_2 =  sim_artificial_noise_2(sigma_2s, sigma_2s_alternative, hs, alpha, sigma, S, theta_0)


hs = seq(0.002, 0.13, length.out=50)
sigma_2s = seq(0.018, 0.38, length.out = 50)
sigma_2s_alternative = seq(0.016, 0.32, length.out = 50)
sigma = 0.1
sim_skew_2_2 =  sim_artificial_noise_2(sigma_2s, sigma_2s_alternative, hs, alpha, sigma, S, theta_0)


hs = seq(0.0015, 0.13, length.out=50)
sigma_2s = seq(0.018, 0.39, length.out = 50)
sigma_2s_alternative = seq(0.018, 0.34, length.out = 50)
sigma = 1
sim_skew_3_2 =  sim_artificial_noise_2(sigma_2s,sigma_2s,  hs, alpha, sigma, S, theta_0)



sigma = 5
hs = seq(0.005, 0.145, length.out=50)
sigma_2s = exp(seq(log(0.018, exp(1)), log(1.4, exp(1)), length.out = 50))
sigma_2s_alternative = exp(seq(log(0.008, exp(1)), log(1.2, exp(1)), length.out = 50))
sim_skew_5_2 = sim_artificial_noise_2(sigma_2s, sigma_2s_alternative,  hs, alpha, sigma, S, theta_0)


sigma = 10
hs = seq(0.005, 0.15, length.out=50)
sigma_2s = exp(seq(log(0.005, exp(1)), log(6, exp(1)), length.out = 50))
sigma_2s_alternative = exp(seq(log(0.005, exp(1)), log(4.5, exp(1)), length.out = 50))
sim_skew_6_2 =  sim_artificial_noise_2(sigma_2s, sigma_2s_alternative, hs, alpha, sigma, S, theta_0)




# 2.b
# alpha = 10 (extreme skewness)

alpha = 10
mean_skew = 1 * alpha / (sqrt(1 + alpha**2)) * sqrt(2 / pi)
mean_skew
var_skew = 1 - 2*alpha**2/((1+alpha**2)*pi)
var_skew
sqrt(var_skew)

S = 200000
theta_0 = 0

hs = seq(0.001, 0.13, length.out=50)
sigma_2s = seq(0.011, 0.30, length.out = 50)
sigma_2s_alternative = seq(0.01, 0.27, length.out = 50)
sigma = 0.02
sim_1_10 =  sim_artificial_noise_2(sigma_2s, sigma_2s_alternative,hs, alpha, sigma, S, theta_0)


sigma_2s = seq(0.018, 0.30, length.out = 50)
sigma_2s_alternative = seq(0.016, 0.27, length.out = 50)
sigma = 0.1
sim_2_10 =  sim_artificial_noise_2(sigma_2s, sigma_2s_alternative, hs, alpha, sigma, S, theta_0)


sigma = 1
hs = seq(0.001, 0.13, length.out=50)
sigma_2s = seq(0.018, 0.35, length.out = 50)
sigma_2s_alternative = seq(0.018, 0.3, length.out = 50)
sim_3_10 =  sim_artificial_noise_2(sigma_2s, sigma_2s_alternative, hs, alpha, sigma, S, theta_0)


sigma = 5
hs = seq(0.001, 0.1, length.out=50)
sigma_2s = exp(seq(log(0.018, exp(1)), log(0.9, exp(1)), length.out = 50))
sigma_2s_alternative = exp(seq(log(0.008, exp(1)), log(1.2, exp(1)), length.out = 50))
sim_4_10 = sim_artificial_noise_2(sigma_2s, sigma_2s_alternative, hs, alpha, sigma, S, theta_0)


sigma = 10
hs = seq(0.001, 0.09, length.out=50)
sigma_2s = exp(seq(log(0.005, exp(1)), log(0.9, exp(1)), length.out = 40))
sigma_2s_alternative = exp(seq(log(0.005, exp(1)), log(4.5, exp(1)), length.out = 60))
sim_5_10 =  sim_artificial_noise_2(sigma_2s, sigma_2s_alternative, hs, alpha, sigma, S, theta_0)


#######################################################################################
# 3rd sim (ch. 5.2 in the thesis)
# univariate gaussian model

# generate data
# prior params
options(digits = 8)
mu_0 = 0
tau_0 = 8
# data sd
sigma = 8
N = 1000
set.seed(123)
data1d = gen_data_1d(N, mu0, tau0, sigma)
data1 = data1d$data
mu_n = data1d$mean
sigma_n = data1d$var
print(mu_n)
print(sigma_n)
print(sqrt(sigma_n))

## 500 k
S = 500000
# size of minibatch
n = 10
# define vector of step sizes
h = seq(0.0001, 0.004, length.out = 30)
sigmas = exp(seq(log(0.0001, exp(1)), log(0.045, exp(1)), length.out = 30)) #check
sigmas_alternative = exp(seq(log(0.0001, exp(1)), log(0.0275, exp(1)), length.out = 30))

set.seed(123)
sim_normal_10 = run_sim_1d_2(h, sigmas, sigmas_alternative, S, N, n, theta_0, tau_0, sigma, grad_U1, grad_p1, ksd = FALSE)
plot_simulation_2(sim_normal_10$res_barker, sim_normal_10$res_barker_alternative, sim_normal_10$res_sgld, n, type = 'normal model', legend=FALSE, xlim=5000, ylim2=0.4)


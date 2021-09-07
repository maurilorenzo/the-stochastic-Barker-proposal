# Logistic regression simulations (section 6.1)
source("logistic_regression_setup.R")


# 4th sim (ch. 5.3.1 in the thesis)
d = 100# dimensionality )
N = 1000 # Number of data points

options(digits=8)

set.seed(123)
dataLR = genDataLR(d, N+250)
dataset = dataLR$data[1:N,]
heldOutData = dataLR$data[c((N+1):(N+250)),]

theta = dataLR$theta
options(digits=3)
theta[1:20]
options(digits=8)
Sigma = dataLR$Sigma

# LR (glm)
LR.model = glm(y ~ . -1, data = dataset, family = binomial)
summary(LR.model)
options(digits=22)
logLossLR(heldOutData, LR.model$coefficients)
options(digits=8)
logLossLR(heldOutData, theta)


d = 100 # dimensionality 
N = 1000 # Number of data points
n=10 # minibathch size
S = 200000 # number of iter
h = c(0.0001, 0.0005, 0.001, 0.005, 0.0075)
sigma_2s = c(0.001, 0.0025, 0.05, 0.1, 0.2, 0.5, 0.75)
sigma_2s_alt= c(0.001, 0.0025,0.05, 0.1, 0.2, 1)
sim_logistic_regression_2d = sim_logistic_regression(h, sigma_2s, sigma_2s_alt, 
                                                     S, N, n, d, ksd = FALSE, 
                                                     ksd_thinned = FALSE, 
                                                     moments = TRUE, 
                                                     stan_means = stan_means_new,
                                                     stan_vars = stan_vars_new)

plot_logistic_regression(sim_logistic_regression_2d$res_barker,
                         sim_logistic_regression_2d$res_barker_alternative,
                         sim_logistic_regression_2d$res_sgld, n, FALSE, TRUE, max_=TRUE)





S =  10^5 # Number of SGMCMC iterations
theta0 = rnorm(d, 0, 1) # initialization


## median ESS = 500

# SGLD
h = 0.000865 
set.seed(123)
sgldLR_10 = SGLD_naive(dataset, S, n, theta0, h, gradLogL, gradLogP)

# Brker with normal noise
h = 0.045  
set.seed(123)
barkerLR_10 = barker(dataset, S, n, theta0, h, gradLogL, gradLogP)

# barker with bimodal noise
h = 0.0155  
set.seed(123)
barkerLR_alternative_10 = barker(dataset, S, n, theta0, h, gradLogL, gradLogP, version = 'alternative')



## median ESS = 1000

# SGLD
h = 0.003 # 690-997
set.seed(123)
sgldLR_10_2 = SGLD_naive(dataset, S, n, theta0, h, gradLogL, gradLogP)

# Barker with Normal noise
h = 0.9   # 702-1003
set.seed(123)
barkerLR_10_2 = barker(dataset, S, n, theta0, h, gradLogL, gradLogP)

# Barker with bimodal noise
h = 0.13  # 805-989
set.seed(123)
barkerLR_alternative_10_2 = barker(dataset, S, n, theta0, h, gradLogL, gradLogP, version = 'alternative')




# plots
idx_1 = 72
idx_2 = 98
stan_means_new[idx_1]
stan_vars_new[idx_1]
ylim = c(-10, 5)
plot_trace_LR(idx_1, sgldLR_10, barkerLR_10, barkerLR_alternative_10, list_of_draws_new, stan_means = stan_means_new, ylim = ylim, legend = FALSE)
xlim = ylim
histogram_combined_LR(idx_1, list_of_draws_new, sgldLR_10, barkerLR_10, barkerLR_alternative_10, xlim)
dev.off()
histogram_LR(idx_1, list_of_draws_new, sgldLR_10, barkerLR_10, barkerLR_alternative_10, xlim)

stan_means_new[idx_2]
stan_vars_new[idx_2]
ylim= c(-5, 10)
plot_trace_LR(idx_2, sgldLR_10, barkerLR_10, barkerLR_alternative_10, list_of_draws_new, stan_means = stan_means_new, ylim = ylim, legend = FALSE)
xlim = ylim
histogram_combined_LR(idx_2, list_of_draws_new, sgldLR_10, barkerLR_10, barkerLR_alternative_10, xlim)
dev.off()
histogram_LR(idx_2, list_of_draws_new, sgldLR_10, barkerLR_10, barkerLR_alternative_10, xlim)


stan_vars_new[idx_2]
ylim = c(-20, 5)
plot_trace_LR(idx_1, sgldLR_10_2, barkerLR_10_2, barkerLR_alternative_10_2, list_of_draws_new, stan_means = stan_means_new, ylim = ylim, legend = FALSE)
xlim = ylim
histogram_combined_LR(idx_1, list_of_draws_new, sgldLR_10_2, barkerLR_10_2, barkerLR_alternative_10_2, xlim)
histogram_LR(idx_1, list_of_draws_new, sgldLR_10_2, barkerLR_10_2, barkerLR_alternative_10_2, xlim)

stan_means_new[idx_2]
stan_vars_new[idx_2]
ylim = c(-5, 20)
plot_trace_LR(idx_2, sgldLR_10_2, barkerLR_10_2, barkerLR_alternative_10_2, list_of_draws_new, stan_means = stan_means_new, ylim = ylim, legend = FALSE)
xlim = ylim
histogram_combined_LR(idx_2, list_of_draws_new, sgldLR_10_2, barkerLR_10_2, barkerLR_alternative_10_2, xlim)
histogram_LR(idx_2, list_of_draws_new, sgldLR_10_2, barkerLR_10_2, barkerLR_alternative_10_2, xlim)





idx_3 = 2
idx_4 = 3

par(mfrow=c(2,2))
plot_3d(list_of_draws_new$beta, idx_3, idx_4)
plot_3d(sgldLR_10, idx_3, idx_4)
plot_3d(barkerLR_10, idx_3, idx_4)
plot_3d(barkerLR_alternative_10, idx_3, idx_4)

par(mfrow=c(2,2))
plot_3d(list_of_draws_new$beta, idx_3, idx_4)
plot_3d(sgldLR_10_2, idx_3, idx_4)
plot_3d(barkerLR_10_2, idx_3, idx_4)
plot_3d(barkerLR_alternative_10_2, idx_3, idx_4)



#######################################################################################
# 5th suìim (ch. 5.3.2 in the thesis)
######## Arrhythmia dataset 100 covariates
arrhythmia = read.csv("data/arrhythmia.data", header=FALSE, na.strings="?")
arrhythmia = arrhythmia[arrhythmia['V280'] < 16, ]
arrhythmia$y = 1
arrhythmia$y[arrhythmia['V280']==1] = 0
arrhythmia = arrhythmia[names(arrhythmia) != 'V280']
sum(arrhythmia$y==1)
sum(arrhythmia$y==0)
colSums(is.na(arrhythmia))
# eliminate columns with only zeros
arrhythmia = arrhythmia[, ! names(arrhythmia) %in% c('V11', 'V12', 'V13', 'V14', 'V15')]
arrhythmia = arrhythmia[, ! names(arrhythmia) %in% c('V20', 'V68', 'V70','V84', 'V205', 'V165', 'V157', 'V158', 'V146', 'V152', 'V140', 'V142', 'V144', 'V132', 'V133', 'V275', 'V265')]
arrhythmia = arrhythmia[,c(ncol(arrhythmia), 1:100)]
dim(arrhythmia)
head(arrhythmia)

# standard scaling
arrhythmia[,-1] <- scale(arrhythmia[,-1])

# train-test split
set.seed(123)
train = sample(1:nrow(arrhythmia), nrow(arrhythmia)*0.8) 
train = as.numeric(rownames(train_set))

train_set = arrhythmia[train, ]
test_set = arrhythmia[-train, ]


# simulation with SGLD and stochastic barker
S = 200000
n = 10
h = c(0.001, 0.005, 0.01, 0.015, 0.019)
sigma_2s = c(0.001, 0.01, 0.1, 1)
sigma_2s_alternative = c(0.001, 0.01, 0.1, 0.5, 1)

sim_arrhytmya_1 = sim_logistic_regression_arrhytmya(h, sigma_2s, sigma_2s_alternative, S, n, train_set, test_set, stan_means=stan_means_arrhythmia_new, stan_vars=stan_vars_arrhythmia_new)
plot_logistic_regression(sim_arrhytmya_1$res_barker, 
                         sim_arrhytmya_1$res_barker_alternative,
                         sim_arrhytmya_1$res_sgld, 1, TRUE, TRUE)


set.seed(123)
S = 100000
theta0 = rnorm(d)

# median ESS = 100

# SGLD
h = 0.0103 
set.seed(123)
sgldLR_arrythmia_1 = SGLD_naive(train_set, S, n, theta0, h, gradLogL, gradLogP)

# Barker gaussian noise
h = 0.165 
set.seed(123)
barker_arrythmia_1 = barker(train_set, S, n, theta0, h, gradLogL, gradLogP)

# Barker bimodal noise
h = 0.09 
set.seed(123)
barker_alternative_arrythmia_1 = barker(train_set, S, n, theta0, h, gradLogL, gradLogP, version='alternative')




bias_mean_sgld = abs(colMeans(sgldLR_arrythmia_1[(S/2+1):S,])-stan_means_arrhythmia_new)
idx_1 = which(bias_mean_sgld==max(bias_mean_sgld)) # 16
idx_1
bias_mean_barker_alt = abs(colMeans(barker_alternative_arrythmia_1[(S/2+1):S,])-stan_means_arrhythmia_new)
idx_2 = which(bias_mean_barker_alt==max(bias_mean_barker_alt)) #81 - 4.56
idx_2
set.seed(2)
idx_3 = sample(100, 1)
idx_4 = sample(100, 1)

ylim = c(-10, 20)
plot_trace_LR(idx_1, sgldLR_arrythmia_1, barker_arrythmia_1, barker_alternative_arrythmia_1, list_of_draws_arrythmia_new, stan_means_arrhythmia_new, ylim=ylim, legend=FALSE)
xlim = c(-5, 10)
histogram_LR(idx_1, list_of_draws_arrythmia_new, sgldLR_arrythmia_1, barker_arrythmia_1, barker_alternative_arrythmia_1, xlim)
histogram_combined_LR(idx_1, list_of_draws_arrythmia_new, sgldLR_arrythmia_1, barker_arrythmia_1, barker_alternative_arrythmia_1, xlim)

plot_trace_LR(idx_2, sgldLR_arrythmia_1, barker_arrythmia_1, barker_alternative_arrythmia_1, list_of_draws_arrythmia_new, stan_means_arrhythmia_new, ylim=ylim, legend=FALSE)
xlim = c(-5, 10)
histogram_LR(idx_2, list_of_draws_arrythmia_new, sgldLR_arrythmia_1, barker_arrythmia_1, barker_alternative_arrythmia_1, xlim)
histogram_combined_LR(idx_2, list_of_draws_arrythmia_new, sgldLR_arrythmia_1, barker_arrythmia_1, barker_alternative_arrythmia_1, xlim)

ylim = c(-15, 15)
plot_trace_LR(idx_3, sgldLR_arrythmia_1, barker_arrythmia_1, barker_alternative_arrythmia_1, list_of_draws_arrythmia_new, stan_means_arrhythmia_new, ylim=ylim, legend=FALSE)
xlim = c(-5, 5)
histogram_LR(idx_3, list_of_draws_arrythmia_new, sgldLR_arrythmia_1, barker_arrythmia_1, barker_alternative_arrythmia_1, xlim)
histogram_combined_LR(idx_3, list_of_draws_arrythmia_new, sgldLR_arrythmia_1, barker_arrythmia_1, barker_alternative_arrythmia_1, xlim)

ylim = c(-10, 20)
plot_trace_LR(idx_4, sgldLR_arrythmia_1, barker_arrythmia_1, barker_alternative_arrythmia_1, list_of_draws_arrythmia_new, stan_means_arrhythmia_new, ylim=ylim, legend=FALSE)
xlim = c(-5, 10)
histogram_LR(idx_4, list_of_draws_arrythmia_new, sgldLR_arrythmia_1, barker_arrythmia_1, barker_alternative_arrythmia_1, xlim)
histogram_combined_LR(idx_4, list_of_draws_arrythmia_new, sgldLR_arrythmia_1, barker_arrythmia_1, barker_alternative_arrythmia_1, xlim)



par(mfrow=c(2,2))
plot_3d(list_of_draws_arrythmia_new$beta, idx_1=idx_3, idx_2=idx_4)
plot_3d(sgldLR_arrythmia_1, idx_1=idx_3, idx_2=idx_4)
plot_3d(barker_arrythmia_1, idx_1=idx_3, idx_2=idx_4)
plot_3d(barker_alternative_arrythmia_1, idx_1=idx_3, idx_2=idx_4)


par(mfrow=c(2,2))
plot_3d(list_of_draws_arrythmia_new$beta, idx_1=idx_1, idx_2=idx_2)
plot_3d(sgldLR_arrythmia_1, idx_1=idx_1, idx_2=idx_2)
plot_3d(barker_arrythmia_1, idx_1=idx_1, idx_2=idx_2)
plot_3d(barker_alternative_arrythmia_1, idx_1=idx_5, idx_2=idx_2)




# simulation p_hat
theta = barker_alternative_arrythmia_1[100000,]
set.seed(123)
j = sample(100, 1)
j
ps = sim_p_LR(j, theta,  train_set)
p_hat = sim_p_hat_LR(10000, 10, j, theta, train_set)
z = seq(-5, 5, length.out=10000)
plot(z, ps, type='l', ylab='p')
lines(z, p_hat$p_hat, col='red')


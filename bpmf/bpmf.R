source('bpmf_setup.R')


#import data
df1 = read.delim("ml-100k/u1.base", header=FALSE, col.names = c("uid", "iid", "R", "info"))[c("uid", "iid", "R")]
test1 = read.delim("ml-100k/u1.test", header=FALSE, col.names = c("uid", "iid", "R", "info"))[c("uid", "iid", "R")]


# 6th sim (ch. 5.4 in the thesis)


########### 80
options(digits=8)
mini_batch = 80
S = 50000 # 50 epochs

# SGLD
h = c(0.00001,  0.00001, 0.00001)
mf_sgld_80 = BPMF(h, df1, test1, S, d, miniBatchSize = mini_batch, algo='sgld', burn_in=S/2, grad=FALSE, train_accuracy=FALSE)

# Barker with normal noise
h = c(0.0001, 0.0001, 0.0001)
mf_barker_80 = BPMF(h, df1, test1, S, d, miniBatchSize = mini_batch, algo = 'barker', burn_in = S/2, grad=FALSE, train_accuracy = FALSE)

# Barker with bimodal noise
mf_barker_alternative_80= BPMF(h, df1, test1, S, d, miniBatchSize = mini_batch, algo = 'barker', version='alternative', burn_in = S/2, train_accuracy = FALSE)


# plot rmse
plot_BPMF(mf_barker_80, mf_barker_alternative_80, mf_sgld_80, S=50000, div=1000)




########## 800
mini_batch = 800
S = 5000 # 50 ephocs

# SGLD
h = c(0.00005, 0.00005, 0.00005)
mf_sgld_800 = BPMF(h, df1, test1, S, d, miniBatchSize = mini_batch, algo = 'sgld', burn_in = S/2, grad = TRUE, train_accuracy = TRUE)

# Barker with normal noise
h = c(0.0005, 0.0005, 0.0005)
mf_barker_800 = BPMF(h, df1, test1, S, d, miniBatchSize = mini_batch, algo = 'barker', burn_in = S/2, grad=FALSE, train_accuracy = TRUE)

# Barker with bimodal noise
mf_barker_alternative_800 = BPMF(h, df1, test1, S, d, miniBatchSize = mini_batch, algo = 'barker', version='alternative', burn_in = S/2, return_samples=TRUE, train_accuracy = TRUE)


# plot rmse
plot_BPMF(mf_barker_800, mf_barker_alternative_800, mf_sgld_800, S=5000, div=100)

# plot SGLD grad
plot_gradient(mf_sgld_800)


######### 8000
mini_batch = 8000
S = 500 # 50 ephocs
h = c(0.0001, 0.0001, 0.0001)
mf_sgld_8000 = BPMF(h, df1, test1, S, d, miniBatchSize = mini_batch, algo = 'sgld', burn_in = S/2, grad = FALSE)


h = c(0.005, 0.005, 0.005)
mf_barker_8000 = BPMF(h, df1, test1, S, d, miniBatchSize = mini_batch, algo = 'barker', burn_in = S/2, grad=FALSE)

mf_barker_alternative_8000 = BPMF(h, df1, test1, S, d, miniBatchSize = mini_batch, algo = 'barker', version='alternative', burn_in = S/2, return_samples = TRUE)

# plot rmse
plot_BPMF(mf_barker_8000, mf_barker_alternative_8000, mf_sgld_8000, S=500, div=10)


# type='else'
mini_batch = 8000
S = 500 # 50 ephocs

# SGLD
h = c(0.0001, 0.0001, 0.0001)
mf_sgld_8000 = BPMF(h, df1, test1, S, d, miniBatchSize = mini_batch, algo = 'sgld', burn_in = S/2, grad = FALSE)

# Barker with normal noise
h = c(0.005, 0.005, 0.005)
mf_barker_8000 = BPMF(h, df1, test1, S, d, miniBatchSize = mini_batch, algo = 'barker', burn_in = S/2, grad=FALSE)

# Barker with bimodal noise
mf_barker_alternative_8000 = BPMF(h, df1, test1, S, d, miniBatchSize = mini_batch, algo = 'barker', version='alternative', burn_in = S/2, return_samples = TRUE)


# plot rmse
plot_BPMF(mf_barker_8000_e, mf_barker_alternative_8000_e, mf_sgld_8000_e, S=500, div=10)


# varying the step-size across the params
mini_batch = 800
S = 5000 # 50 ephocs

# SGLD
h = c(0.000075, 0.0005, 0.001)
mf_sgld_800_v_2 = BPMF(h, df1, test1, S, d, miniBatchSize = mini_batch, algo = 'sgld', burn_in = S/2, grad = FALSE, type='else')

# Barker with gaussian noise
h = c(0.0005, 0.005, 0.005)
mf_barker_800_v = BPMF(h, df1, test1, S, d, miniBatchSize = mini_batch, algo = 'barker', burn_in = S/2, grad=FALSE)

# Barker with isotropic noise
mf_barker_alternative_800_v = BPMF(h, df1, test1, S, d, miniBatchSize = mini_batch, algo = 'barker', version='alternative', burn_in = S/2)


plot_BPMF(mf_barker_800_v, mf_barker_alternative_800_v, mf_sgld_800_v_2, 
          S=5000, div=100, train=TRUE)


# p_hat simulation
samples = mf_barker_alternative_800_e$samples 
set.seed(123)
u_i = sample(dim(samples$U)[2], 1) #415
set.seed(123)
v_j = sample(dim(samples$U)[1], 1) # 15
u_i
v_j

p_s = sim_p(u_i, v_j, df1, samples$U, samples$V, samples$mu_U, samples$lambda_U, alpha=3)
set.seed(123)
sim_p_hat_1 = sim_p_hat(10000, 800, u_i, v_j, df1, samples$U, samples$V, samples$mu_U, samples$lambda_U, alpha=3)

x = seq(-5, 5, length.out=10000)

plot(x, p_s, type='l', xlab='z', ylab='p')
points(x, sim_p_hat_1$p_hat, col='red', type='l')

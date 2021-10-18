# The Stochastic Barker Proposal


## Introduction
Supporting R code of my master thesis.  

The *Stochastic Barker Proposal* (SBP) is the stochastic gradient version of the *Barker scheme*[[1]](#1), a recently introduced gradient based MCMC. The Barker scheme has been shown to outperform standard gradient based MCMC, as MALA and HMC, in terms of robustness to hyperparameter tuning.  
The peculiarity of Barker lies in how the gradient is used. In particular, Barker uses the jump kernel of a Markov jump process that approximately preserves the target distribution as proposal. The proposal is a skew-symmetric distribution where the degree of skewnwess depends on the magnitude of the gradient.  
The Stochastic Barker Proposal is obtained substituting the gradient of the potential function with an unbiased and computationally cheaper estimate based on a mini-batch of the original dataset.  

This represents the first study of the subject. The SBP is compared to the stochastic gradient Langevin dyanmics (SGLD) algorithm[[2]](#2) through simulations involving both simulated and real data.  
Firstly, I considered two toy models with respectively a standard Gaussian and a skew Normal distribution; in these examples, no data is involved but at each iteration isotropic Gaussian noise is added to the true gradient of the target to study how different noise variances affect the relation between mixing and accuracy.  
Secondly, I studied a Bayesian Normal model with simulated data. Thirdly, a Bayesian logistic regression for binary data classification applied to simulated data and the Arrhythmia dataset. Finally, a Bayesian probabilistic factorization model for predicting movie ratings based on the MovieLens data set. I ran all simulation keeping the hyperparameters fixed as this is what is usually done in practice.   

The SBP shows greater **robustness to hyperparameter tuning** when the posterior is irregular and a promising **predictive accuracy** on unseen data.


## References
<a id="1">[1]</a> 
Samuel Livingstone, Giacomo Zanella (2020),
The Barker proposal: combining robustness and efficiency in gradient-based MCMC,
arXiv:1908.11812 [stat.CO]

<a id="2">[2]</a> 
Max Welling, Yee Whye Teh (2011),
Bayesian learning via stochastic gradient langevin dynamics,
ICML'11

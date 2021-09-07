# the-stochastic-Barker-proposal
Supporting R code of my master thesis.  

The *stochastic Barker proposal* (SBP) is the stochastic gradient version of the *Barker scheme*, a recently introduced gradient based MCMC. The Barker scheme has been shown to outperform standard gradient based MCMC, as MALA and HMC, in terms of robustness to hyperparameter tuning.  
This represent the first study of the subject. SBP is compared to the stochastic gradient Langevin dyanmics (SGLD) algorithm through simulations involving both simulated and real data.  
Firstly, I considered two toy models with respectively a standard Gaussian and a skew Normal distribution; in these examples, no data were involved but at each iteration isotropic Gaussian noise was added to the true gradient of the target to study how different noise variances affect the relation between mixing and accuracy. 
Secondly, I studied a Bayesian Normal model with simulated data. Thirdly, a Bayesian logistic regression for binary data classification applied to simulated data and the Arrhythmia dataset. Finally, a Bayesian probabilistic factorization model for predicting movie ratings based on the MovieLens data set. I ran all simulation keeping the hyperparameters fixed as this is what is usually done in practice.

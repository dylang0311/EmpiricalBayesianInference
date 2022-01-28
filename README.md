# EmpiricalBayesianInference

This program seeks to use an empirical Bayesian approach to reconstruct signals from multiple measurement vectors (MMVs) of Fourier data assuming joint sparsity in the MMVs. The original work by [Zhang et al.](https://arxiv.org/abs/2103.15618) performed this process for real signals using variance-based joint sparsity (VBJS) to identify the support of the sparse domain to form the prior and the Metropolis Hastings MCMC algorithm to sample from the posterior of the signal. The signals considered here were either 1D and sparse in the spatial or edge domain or were 2D and sparse in the spatial domain.

This work seeks to expand on this approach by considering both real and complex signals, 2D signals that are sparse in the edge domain, forming the posterior in the spatial and frequency domains, and sampling using Hamiltonian MC. The software currently does not efficiently support the sampling of signals with length much larger than 500.

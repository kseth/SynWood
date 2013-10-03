SynWood
=======
In addition to the following, theses functionality are progressively transitioned to the synlik package https://bitbucket.org/cbarbu/synlik

LDsynLikSpat.R contains the likelihood functions and simulation with a model using the kernel

noKernelMCMC.R contains the likelihood functions and simulation with a model that doesn't use the kernel

functions_migration.R contains all the interfacing between R and C for the gillespie and probability matrix generation

functions_migration.c contains all the c code for running multiple gillespies, for computing the statistics, and for computing the probability matrix

functions_sampling.R contains the sampling code written in R

param.r contains the parameters being used for the gillespie simulations -> this file is an important file (it needs to be written in some better manner later)


To Do
------

Must convert MCMC code to a function

Must unify LDsynLikSpat.R and noKernelMCMC.R

Must implement likelihood functions for priors

Must bring in real data to test

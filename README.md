SynWood
=======

LDsynLikSpat.R contains the likelihood functions and MCMC simulation as written with using the OmniSample function in functions_sampling.R

functions_migration.R contains all the interfacing between R and C for the gillespie and probability matrix generation

functions_migration.c contains all the c code for running multiple gillespies, for computing the statistics, and for computing the probability matrix

functions_sampling.R contains the MCMC sampling code written in R 

param.r contains the parameters being used for the gillespie simulations -> this file is an important file (it needs to be written in some better manner later)

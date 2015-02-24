#!/bin/bash - 
#===============================================================================
#
#          FILE: SynWood_compile.sh
# 
#         USAGE: ./SynWood_compile.sh 
# 
#   DESCRIPTION: Allow to compile in one line everything needed for SynWood
#               Nota: with the MAKEFLAGS allows to analyse core-dumps:
#                   ulimit -c 2000
#                   R -d gdb
# 		    core-file core 
#               or directly navigage a segfault:
#                   ulimit -c 2000
#                   R -d gdb
#                   R 
#                   source("daisyChainMCMC.R")
#              if segfault will stop saying a lot and
#                   p myvariable 
#              will print the state of my variable!
#
#       OPTIONS: ---
#  REQUIREMENTS: gsl (on ubuntu libgsl0-dev)
#          BUGS: ---
#         NOTES: ---
#        AUTHORS: KS, CB
#  ORGANIZATION: 
#       CREATED: 05/09/2013 03:31:25 PM EDT
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error
rm *.o
R CMD SHLIB kmeans.c
MAKEFLAGS="CFLAGS=-g" R CMD SHLIB -lgsl -lgslcblas -lm functions_migration.c samlmu.f 
# R CMD SHLIB -lgsl -lgslcblas -lm functions_migration.c samlmu.f 

#!/bin/bash - 
#===============================================================================
#
#          FILE: SynWood_compile.sh
# 
#         USAGE: ./SynWood_compile.sh 
# 
#   DESCRIPTION: Allow to compile in one line everything needed for SynWood
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
R CMD SHLIB kmeans.c
R CMD SHLIB -lgsl -lgslcblas -lm functions_migration.c 

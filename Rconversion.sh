#!/bin/bash
# Script to convert genomicSimulation (C version) files
# to the library files used by genomicSimulation (R version)
#    -> uncomments //RPACKINSERT lines
#    -> removes C random generators
#    -> replaces C standard libraries with R library equivalents

awk 'NR<181 || NR>248' ./sim-operations.c | sed -e 's+//RPACKINSERT ++g' -e 's/fprintf(stderr,\(.*\)); exit([0-9]);/error(\1);/g' -e 's/fprintf(stderr,\(.*\));/warning(\1);/g' -e 's/randpoi/Rf_rpois/g' -e 's/randb()/(unif_rand() > 0.5)/g' -e 's/randlim\(.*\);/round(unif_rand() * \1);/g' -e 's/printf("/Rprintf("/g' -e 's+srand+//srand+g' -e 's+((double)rand() / (double)RAND_MAX)+unif_rand()+g' > ../sim-operations-for-R.c

awk 'NR<268 || NR>274' ./sim-operations.h | sed -e 's/<stdlib.h>/<R.h>/' -e 's/<stdio.h>/<Rinternals.h>/' -e 's/<math.h>/<Rmath.h>/' -e 's+<time.h>+<R_ext/Utils.h>+' -e 's+#define PI+//#define PI+' > ../sim-operations-for-R.h
#!/bin/bash
# Script to convert genomicSimulation (C version) files
# to the library files used by genomicSimulation (R version)
#    -> uncomments //RPACKINSERT lines
#    -> removes C random generators
#    -> replaces C standard libraries with R library equivalents

echo ">>> Deleting lines from sim-operations.c:"
awk 'NR>=179 && NR<=213' ./sim-operations.c

awk 'NR<179 || NR>213' ./sim-operations.c | sed \
	-e 's+//RPACKINSERT ++g' \
	-e 's/fprintf(stderr,\(.*\)); exit([0-9]);/error(\1);/g' \
	-e 's/fprintf(stderr,/warning(/g' \
	-e 's/rnd_pcg_seed( &d->rng, RNGseed )//g' \
	-e 's/RND_U32 RNGseed//g' \
	-e 's/gsc_shuffle_up_to([^,]*,/shuffle_up_to(/g' \
	-e 's/gsc_randpoi(.*,/Rf_rpois(/g' \
	-e 's/rnd_pcg_range([^,]*,0,1)/(unif_rand() > 0.5)/g' \
	-e 's/rnd_pcg_range([^,]*,0,\(.*\));/round(unif_rand() * (\1));/g' \
	-e 's/rnd_pcg_nextf([^,]*);/unif_rand();/g' \
	-e 's/printf("/Rprintf("/g' -e 's+srand+//srand+g' \
	-e 's+((double)rand() / (double)RAND_MAX)+unif_rand()+g' \
	-e 's/size_t/unsigned int/g' \
	-e 's/SIZE_MAX/UINT_MAX/g' \
	> ../sim-operations-for-R.c
	
echo ">>> Deleting lines from sim-operations.h:"
awk 'NR>=1097 && NR<=1118' ./sim-operations.h

awk 'NR<1097 || NR>1118' ./sim-operations.h | sed \
	-e 's+#include "lib/rnd.h"++g' \
	-e 's/RND_U32 RNGseed//g' \
	-e 's+rnd_pcg_t rng;+//CRANDOMGENERATOR+g' \
	-e 's/shuffle_up_to([^,]*,/shuffle_up_to(/g' \
	-e 's/<stdlib.h>/<R.h>/' \
	-e 's/<stdio.h>/<Rinternals.h>/' \
	-e 's/<math.h>/<Rmath.h>/' \
	-e 's+<time.h>+<R_ext/Utils.h>+' \
	-e 's+#define PI+//#define PI+' \
	-e 's/size_t/unsigned int/g' \
	> ../sim-operations-for-R.h
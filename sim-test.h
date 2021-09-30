#ifndef SIM_OPERATIONS_H
#include "sim-operations.h"
#endif

#include <assert.h>

#define TOL 0.00001

const char HELPER_GENOTYPES[] = "name\tG01\tG02\tG03\tG04\tG05\tG06\n"
	"m1\tTT\tTT\tTT\tTA\tTT\tAT\n"
	"m3\tTT\tTT\tTA\tTA\tTT\tTT\n"
	"m2\tAA\tAA\tAA\tAA\tTT\tAA";
const char HELPER_MAP[] = "marker chr pos\n" "m3 3 15\n"
	"m2 1 8.3\n" "m1 1 5.2";
const char HELPER_EFF[] = "m1 A -0.8\n" "m2 A -0.1\n"
	"m3 A 0.1\n" "m1 T 0.9\n"
	"m3 T -0.1\n" "m2 T -0.5";
const char HELPER_PLAN[] = "G01\tG02\tG03\tG05\n"
	"G01\tG03\tG05\tG06\n" "G05\tG06\tG01\tG04";

float calculate_heterozygosity(SimData* d, int group_number);

int test_loaders(SimData* d);
int test_effect_calculators(SimData *d, int g0);
int test_optimal_calculators(SimData *d);
int test_crossing(SimData *d, int g0);
int test_crossing_unidirectional(SimData *d, int g0);
int test_crossing_from_file(SimData *d, char* fname);
int test_crossing_selfing(SimData *d, int g1);
int test_deletors(SimData *d, int g0);
int test_block_generator(SimData *d);
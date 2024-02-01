#include "sim-operations.h"

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
const char HELPER_EFF2[] = "m1 A 1\n";
const char HELPER_PLAN[] = "G01\tG02\tG03\tG05\n"
	"G01\tG03\tG05\tG06\n" "G05\tG06\tG01\tG04";

float calculate_heterozygosity(SimData* d, GroupNum group_number);
int compareFiles(char* f1, char* f2); // returns 0 if matching

GroupNum test_loaders(SimData* d);

GroupNum test_grouping(SimData *d, GroupNum g0);
GroupNum test_specific_splits(SimData *d, GroupNum g0);
GroupNum test_random_splits(SimData *d, GroupNum g0);
GroupNum test_labels(SimData *d, GroupNum g0);

int test_effect_calculators(SimData *d, GroupNum g0);
int test_optimal_calculators(SimData *d, GroupNum g0);

int test_data_access(SimData* d, GroupNum gp);
int test_iterators(SimData* d, GroupNum gp);
int test_getters(SimData* d, GroupNum gp);

int test_crossing(SimData *d, GroupNum g0);
GroupNum test_crossing_unidirectional(SimData *d, GroupNum g0);
GroupNum test_crossing_from_file(SimData *d, char* fname);
GroupNum test_crossing_selfing(SimData *d, GroupNum g1);
int test_crossing_randomly(SimData *d, GroupNum g1);

int test_deletors(SimData *d, GroupNum g0);

int test_block_generator(SimData *d);

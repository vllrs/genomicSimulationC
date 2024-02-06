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

const char TEST1_TRUTH_save_group_alleles[] = "5	m1	m2	m3\n"
        "G01	TT	AA	TT\n"
        "G02	TT	AA	TT\n"
        "F17	TT	AT	TT\n"
        "F18	TT	AA	AT\n"
        "s12	TT	AT	TT\n"
        "s16	AT	AA	TT\n"
        "12	TT	AA	TT\n"
        "13	TT	AA	TT\n";
const char TEST1_TRUTH_save_allele_matrix[] = "	m1	m2	m3\n"
        "G01	TT	AA	TT\n"
        "G02	TT	AA	TT\n"
        "G03	TT	AA	TA\n"
        "G04	TA	AA	TA\n"
        "G05	TT	TT	TT\n"
        "G06	AT	AA	TT\n"
        "F17	TT	AT	TT\n"
        "F18	TT	AA	AT\n"
        "F19	TT	AA	TA\n"
        "F110	TT	TA	TT\n"
        "F111	TA	AA	TT\n"
        "s12	TT	AT	TT\n"
        "s13	TT	AA	TT\n"
        "s14	TT	AA	TT\n"
        "s15	TT	TA	TT\n"
        "s16	AT	AA	TT\n"
        "12	TT	AA	TT\n"
        "13	TT	AA	TT\n"
        "14	TT	AA	TT\n"
        "15	TT	TT	TT\n"
        "16	TT	AA	TT\n";
const char TEST1_TRUTH_save_transposed_allele_matrix[] =
        "	G01	G02	G03	G04	G05	G06	F17	F18	F19	F110	F111	s12	s13	s14	s15	s16	12	13	14	15	16\n"
        "m1	TT	TT	TT	TA	TT	AT	TT	TT	TT	TT	TA	TT	TT	TT	TT	AT	TT	TT	TT	TT	TT\n"
        "m2	AA	AA	AA	AA	TT	AA	AT	AA	AA	TA	AA	AT	AA	AA	TA	AA	AA	AA	AA	TT	AA\n"
        "m3	TT	TT	TA	TA	TT	TT	TT	AT	TA	TT	TT	TT	TT	TT	TT	TT	TT	TT	TT	TT	TT\n";
const char TEST1_TRUTH_save_transposed_group_alleles[] =
        "5	G01	G02	F17	F18	s12	s16	12	13\n"
        "m1	TT	TT	TT	TT	TT	AT	TT	TT\n"
        "m2	AA	AA	AT	AA	AT	AA	AA	AA\n"
        "m3	TT	TT	TT	AT	TT	TT	TT	TT\n";
const char TEST1_TRUTH_save_count_matrix[] =
        "	G01	G02	G03	G04	G05	G06	F17	F18	F19	F110	F111	s12	s13	s14	s15	s16	12	13	14	15	16\n"
        "m1	0	0	0	1	0	1	0	0	0	0	1	0	0	0	0	1	0	0	0	0	0\n"
        "m2	2	2	2	2	0	2	1	2	2	1	2	1	2	2	1	2	2	2	2	0	2\n"
        "m3	0	0	1	1	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0\n";
const char TEST1_TRUTH_save_count_matrix_of_group[] =
        "5	G01	G02	F17	F18	s12	s16	12	13\n"
        "m1	2	2	2	2	2	1	2	2\n"
        "m2	0	0	1	0	1	0	0	0\n"
        "m3	2	2	2	1	2	2	2	2\n";
const char TEST1_TRUTH_save_bvs[] = "1	G01	1.400000\n"
        "2	G02	1.400000\n"
        "3	G03	1.600000\n"
        "4	G04	-0.100000\n"
        "5	G05	0.600000\n"
        "6	G06	-0.300000\n"
        "7	F17	1.000000\n"
        "8	F18	1.600000\n"
        "9	F19	1.600000\n"
        "10	F110	1.000000\n"
        "11	F111	-0.300000\n"
        "0	s12	1.000000\n"
        "0	s13	1.400000\n"
        "0	s14	1.400000\n"
        "0	s15	1.000000\n"
        "0	s16	-0.300000\n"
        "12		1.400000\n"
        "13		1.400000\n"
        "14		1.400000\n"
        "15		0.600000\n"
        "16		1.400000\n";
const char TEST1_TRUTH_save_group_bvs[] = "1	G01	1.400000\n"
        "2	G02	1.400000\n"
        "7	F17	1.000000\n"
        "8	F18	1.600000\n"
        "0	s12	1.000000\n"
        "0	s16	-0.300000\n"
        "12		1.400000\n"
        "13		1.400000\n";
const char TEST1_TRUTH_save_marker_blocks[] = "Chrom	Pos	Name	Class	Markers\n"
        "0	0	b0	b	m1;m2;\n"
        "0	0	b0	b	m3;\n";
const char TEST1_TRUTH_save_local_bvs[] = "G01_1 0.800000 -0.100000\n"
        "G01_2 0.800000 -0.100000\n"
        "G02_1 0.800000 -0.100000\n"
        "G02_2 0.800000 -0.100000\n"
        "G03_1 0.800000 -0.100000\n"
        "G03_2 0.800000 0.100000\n"
        "G04_1 0.800000 -0.100000\n"
        "G04_2 -0.900000 0.100000\n"
        "G05_1 0.400000 -0.100000\n"
        "G05_2 0.400000 -0.100000\n"
        "G06_1 -0.900000 -0.100000\n"
        "G06_2 0.800000 -0.100000\n"
        "F17_1 0.800000 -0.100000\n"
        "F17_2 0.400000 -0.100000\n"
        "F18_1 0.800000 0.100000\n"
        "F18_2 0.800000 -0.100000\n"
        "F19_1 0.800000 -0.100000\n"
        "F19_2 0.800000 0.100000\n"
        "F110_1 0.400000 -0.100000\n"
        "F110_2 0.800000 -0.100000\n"
        "F111_1 0.800000 -0.100000\n"
        "F111_2 -0.900000 -0.100000\n"
        "s12_1 0.800000 -0.100000\n"
        "s12_2 0.400000 -0.100000\n"
        "s13_1 0.800000 -0.100000\n"
        "s13_2 0.800000 -0.100000\n"
        "s14_1 0.800000 -0.100000\n"
        "s14_2 0.800000 -0.100000\n"
        "s15_1 0.400000 -0.100000\n"
        "s15_2 0.800000 -0.100000\n"
        "s16_1 -0.900000 -0.100000\n"
        "s16_2 0.800000 -0.100000\n"
        "12_1 0.800000 -0.100000\n"
        "12_2 0.800000 -0.100000\n"
        "13_1 0.800000 -0.100000\n"
        "13_2 0.800000 -0.100000\n"
        "14_1 0.800000 -0.100000\n"
        "14_2 0.800000 -0.100000\n"
        "15_1 0.400000 -0.100000\n"
        "15_2 0.400000 -0.100000\n"
        "16_1 0.800000 -0.100000\n"
        "16_2 0.800000 -0.100000\n";
const char TEST1_TRUTH_save_group_local_bvs[] = "G01_1 0.800000 -0.100000\n"
        "G01_2 0.800000 -0.100000\n"
        "G02_1 0.800000 -0.100000\n"
        "G02_2 0.800000 -0.100000\n"
        "F17_1 0.800000 -0.100000\n"
        "F17_2 0.400000 -0.100000\n"
        "F18_1 0.800000 0.100000\n"
        "F18_2 0.800000 -0.100000\n"
        "s12_1 0.800000 -0.100000\n"
        "s12_2 0.400000 -0.100000\n"
        "s16_1 -0.900000 -0.100000\n"
        "s16_2 0.800000 -0.100000\n"
        "12_1 0.800000 -0.100000\n"
        "12_2 0.800000 -0.100000\n"
        "13_1 0.800000 -0.100000\n"
        "13_2 0.800000 -0.100000\n";
const char TEST1_TRUTH_save_full_pedigrees[] = "1	G01\n"
        "2	G02\n"
        "3	G03\n"
        "4	G04\n"
        "5	G05\n"
        "6	G06\n"
        "7	F17=(G03,G05)\n"
        "8	F18=(G04,G02)\n"
        "9	F19=(G02,G04)\n"
        "10	F110=(G05,G03)\n"
        "11	F111=(G01,G06)\n"
        "0	s12=(F17=(G03,G05))\n"
        "0	s13=(F18=(G04,G02))\n"
        "0	s14=(F19=(G02,G04))\n"
        "0	s15=(F110=(G05,G03))\n"
        "0	s16=(F111=(G01,G06))\n"
        "12	\n"
        "13	\n"
        "14	\n"
        "15	\n"
        "16	\n"
        "17	F217=(F111=(G01,G06),F110=(G05,G03))\n";
const char TEST1_TRUTH_save_group_full_pedigrees[] = "1	G01\n"
        "2	G02\n"
        "7	F17=(G03,G05)\n"
        "8	F18=(G04,G02)\n"
        "0	s12=(F17=(G03,G05))\n"
        "0	s16=(F111=(G01,G06))\n"
        "12	\n"
        "13	\n";
const char TEST1_TRUTH_save_one_step_pedigrees[] = "G01		\n"
        "G02		\n"
        "G03		\n"
        "G04		\n"
        "G05		\n"
        "G06		\n"
        "F17	G03	G05\n"
        "F18	G04	G02\n"
        "F19	G02	G04\n"
        "F110	G05	G03\n"
        "F111	G01	G06\n"
        "s12	F17	F17\n"
        "s13	F18	F18\n"
        "s14	F19	F19\n"
        "s15	F110	F110\n"
        "s16	F111	F111\n"
        "12		\n"
        "13		\n"
        "14		\n"
        "15		\n"
        "16		\n"
        "F217	F111	F110\n";
const char TEST1_TRUTH_save_group_one_step_pedigrees[] = "G01		\n"
        "G02		\n"
        "F17	G03	G05\n"
        "F18	G04	G02\n"
        "s12	F17	F17\n"
        "s16	F111	F111\n"
        "12		\n"
        "13		\n";
const char TEST1_TRUTH_sayg_genotype_header[] = "	m1	m2	m3";
const char TEST1_TRUTH_sayg_genotype_bodyrow[] = "0	TT	AA	TT";
const char TEST1_TRUTH_sayg_bv_bodyrow[] = "0		1.400000";
const char TEST1_TRUTH_sayg_pedigree_bodyrow[] = "0	=(16)";

float calculate_heterozygosity(SimData* d, GroupNum group_number);
int compareFiles(char* f1, char* f2); // returns 0 if matching
int compareFileToString(char* filename, const char* target);
int compareRepeatingFileToTable(char* filename, unsigned int expectedNRows, const char* header, const char* body);

GroupNum just_load(SimData* d);
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

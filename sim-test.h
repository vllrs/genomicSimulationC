#include "sim-operations.h"

//#define GSC_DEPRECATED_VERBOSE
//#include "names_stopgap.h"

#include <assert.h>

#define TOL 0.00001

const char HELPER_GENOTYPES[] = "name\tG01\tG02\tG03\tG04\tG05\tG06\n"
	"m1\tTT\tTT\tTT\tTA\tTT\tAT\n"
    "am3\tTT\tTT\tTA\tTA\tTT\tTT\n"
	"m2\tAA\tAA\tAA\tAA\tTT\tAA";
const char HELPER_MAP[] = "marker chr pos\n" "am3 3 15\n"
	"m2 1 8.3\n" "m1 1 5.2";
const char HELPER_EFF[] = "m1 A -0.8\n" "m2 A -0.1\n"
    "am3 A 0.1\n" "m1 T 0.9\n"
    "am3 T -0.1\n" "m2 T -0.5";
const char HELPER_EFF2[] = "m1 A 1\n";
const char HELPER_PLAN[] = "G01\tG02\tG03\tG05\n"
	"G01\tG03\tG05\tG06\n" "G05\tG06\tG01\tG04";

const char TEST1_TRUTH_save_group_alleles[] = "5	m1	m2	am3\n"
        "G01	TT	AA	TT\n"
        "G02	TT	AA	TT\n"
        "F17	TT	AT	AT\n"
        "F18	AT	AA	TT\n"
        "s12	TT	AT	TT\n"
        "s16	TT	AA	TT\n"
        "12	TT	AA	TT\n"
        "13	TT	AA	TT\n";
const char TEST1_TRUTH_save_allele_matrix[] = "	m1	m2	am3\n"
        "G01	TT	AA	TT\n"
        "G02	TT	AA	TT\n"
        "G03	TT	AA	TA\n"
        "G04	TA	AA	TA\n"
        "G05	TT	TT	TT\n"
        "G06	AT	AA	TT\n"
        "F17	TT	AT	AT\n"
        "F18	AT	AA	TT\n"
        "F19	TT	AA	TT\n"
        "F110	TT	TA	TT\n"
        "F111	TT	AA	TT\n"
        "s12	TT	AT	TT\n"
        "s13	TT	AA	TT\n"
        "s14	TT	AA	TT\n"
        "s15	TT	TA	TT\n"
        "s16	TT	AA	TT\n"
        "12	TT	AA	TT\n"
        "13	TT	AA	TT\n"
        "14	TT	AA	TT\n"
        "15	TT	TT	TT\n"
        "16	TT	AA	TT\n";
const char TEST1_TRUTH_save_transposed_allele_matrix[] =
        "	G01	G02	G03	G04	G05	G06	F17	F18	F19	F110	F111	s12	s13	s14	s15	s16	12	13	14	15	16\n"
        "m1	TT	TT	TT	TA	TT	AT	TT	AT	TT	TT	TT	TT	TT	TT	TT	TT	TT	TT	TT	TT	TT\n"
		"m2	AA	AA	AA	AA	TT	AA	AT	AA	AA	TA	AA	AT	AA	AA	TA	AA	AA	AA	AA	TT	AA\n"
		"am3	TT	TT	TA	TA	TT	TT	AT	TT	TT	TT	TT	TT	TT	TT	TT	TT	TT	TT	TT	TT	TT\n";
const char TEST1_TRUTH_save_transposed_group_alleles[] =
        "5	G01	G02	F17	F18	s12	s16	12	13\n"
		"m1	TT	TT	TT	AT	TT	TT	TT	TT\n"
		"m2	AA	AA	AT	AA	AT	AA	AA	AA\n"
		"am3	TT	TT	AT	TT	TT	TT	TT	TT\n";
const char TEST1_TRUTH_save_count_matrix[] =
        "	G01	G02	G03	G04	G05	G06	F17	F18	F19	F110	F111	s12	s13	s14	s15	s16	12	13	14	15	16\n"
        "m1	0	0	0	1	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0\n"
        "m2	2	2	2	2	0	2	1	2	2	1	2	1	2	2	1	2	2	2	2	0	2\n"
        "am3	0	0	1	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0\n";
const char TEST1_TRUTH_save_count_matrix_of_group[] =
        "5	G01	G02	F17	F18	s12	s16	12	13\n"
        "m1	2	2	2	1	2	2	2	2\n"
        "m2	0	0	1	0	1	0	0	0\n"
        "am3	2	2	1	2	2	2	2	2\n";
const char TEST1_TRUTH_save_bvs[] = "1	G01	1.400000\n"
        "2	G02	1.400000\n"
        "3	G03	1.600000\n"
        "4	G04	-0.100000\n"
        "5	G05	0.600000\n"
        "6	G06	-0.300000\n"
        "7	F17	1.200000\n"
        "8	F18	-0.300000\n"
        "9	F19	1.400000\n"
        "10	F110	1.000000\n"
        "11	F111	1.400000\n"
        "0	s12	1.000000\n"
        "0	s13	1.400000\n"
        "0	s14	1.400000\n"
        "0	s15	1.000000\n"
        "0	s16	1.400000\n"
        "12		1.400000\n"
        "13		1.400000\n"
        "14		1.400000\n"
        "15		0.600000\n"
        "16		1.400000\n";
const char TEST1_TRUTH_save_group_bvs[] = "1	G01	1.400000\n"
        "2	G02	1.400000\n"
        "7	F17	1.200000\n"
        "8	F18	-0.300000\n"
        "0	s12	1.000000\n"
        "0	s16	1.400000\n"
        "12		1.400000\n"
        "13		1.400000\n";
/*const char TEST1_TRUTH_save_marker_blocks[] = "Chrom	Pos	Name	Class	Markers\n"
        "0	0	b0	b	m1;m2;\n"
        "0	0	b0	b	am3;\n";*/
const char TEST1_TRUTH_save_marker_blocks_blocksonly[] = "m1;m2;\nam3;\n";
const char TEST1_TRUTH_save_marker_blocks_chrinfo[] = "Chrom	Len	Markers\n"
        "0	3.100000	m1;m2;\n"
        "1	0.000000	am3;\n";
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
        "F17_1 0.800000 0.100000\n"
        "F17_2 0.400000 -0.100000\n"
        "F18_1 -0.900000 -0.100000\n"
        "F18_2 0.800000 -0.100000\n"
        "F19_1 0.800000 -0.100000\n"
        "F19_2 0.800000 -0.100000\n"
        "F110_1 0.400000 -0.100000\n"
        "F110_2 0.800000 -0.100000\n"
        "F111_1 0.800000 -0.100000\n"
        "F111_2 0.800000 -0.100000\n"
        "s12_1 0.800000 -0.100000\n"
        "s12_2 0.400000 -0.100000\n"
        "s13_1 0.800000 -0.100000\n"
        "s13_2 0.800000 -0.100000\n"
        "s14_1 0.800000 -0.100000\n"
        "s14_2 0.800000 -0.100000\n"
        "s15_1 0.400000 -0.100000\n"
        "s15_2 0.800000 -0.100000\n"
        "s16_1 0.800000 -0.100000\n"
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
        "F17_1 0.800000 0.100000\n"
        "F17_2 0.400000 -0.100000\n"
        "F18_1 -0.900000 -0.100000\n"
        "F18_2 0.800000 -0.100000\n"
        "s12_1 0.800000 -0.100000\n"
        "s12_2 0.400000 -0.100000\n"
        "s16_1 0.800000 -0.100000\n"
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
        "17	F217=(F111=(G01,G06),F19=(G02,G04))\n";
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
        "F217	F111	F19\n";
const char TEST1_TRUTH_save_group_one_step_pedigrees[] = "G01		\n"
        "G02		\n"
        "F17	G03	G05\n"
        "F18	G04	G02\n"
        "s12	F17	F17\n"
        "s16	F111	F111\n"
        "12		\n"
        "13		\n";
const char TEST1_TRUTH_sayg_genotype_header[] = "	m1	m2	am3";
const char TEST1_TRUTH_sayg_genotype_bodyrow[] = "0	TT	AA	TT";
const char TEST1_TRUTH_sayg_bv_bodyrow[] = "0		1.400000";
const char TEST1_TRUTH_sayg_pedigree_bodyrow[] = "0	=(16)";


const char* const TEST1_MAPLOADERS[] = {
	// Regular most-expected file
	"marker chr pos\n"
	"first 3 3.2\n"
	"second 3 43.2\n"
	"3rd 1A 1e2\n"
	"fourth 1A 108.3\n"
	"5 1A 110\n"
	,
	// Mixed spacing 
	"marker\tchr pos\n"
	"first\t\t3\t3.2\n"
	"second 3\t43.2\n"
	"3rd 1A 108.3\n"
	"fourth 1A 110\n"
	"5 1A 1e2\n"
	,
	// No final line 
	"marker chr pos\n"
	"first 3 3.2\n"
	"second 3 43.2\n"
	"3rd 1A 108.3\n"
	"fourth 1A 110\n"
	"5 1A 1e2"
	,
	// Rearranged column order 
	"chr marker pos\n"
	"3 first 3.2\n"
	"3 second 43.2\n"
    "1A 3rd 108.3\n"
    "1A fourth 110\n"
    "1A 5 1e2\n"
	,
	// No header line 
	"first 3 3.2\n"
	"second 3 43.2\n"
	"3rd 1A 108.3\n"
	"fourth 1A 110\n"
    "5 1A 1e2"
	,
	// Mix order
	 "first 3 3.2\n"
	"5 1A 1e2\n"
	"3rd 1A 108.3\n"
	"second 3 43.2\n"
	"fourth 1A 110\n"
	,
	// Leaves out a marker, includes a fake marker that should be ignored, and gives them different chrs and positions
	"pos marker chr\n"
	"43.2 first 1B\n"
	"3.2 second 1B\n"
	"1.083e2 3rd 3\n"
	"100 fourth 3\n"
    "110 fake5 3\n"
	,
	// Too many columns on one row
	"first 3 3.2\n"
	"5 1A 1e2\n"
	"3rd 1A 108.3 //commint_oops\n"
	"second 3 43.2\n"
	"fourth 1A 110\n"
	,
	// Single line with header 
	"chr pos marker\n"
    "1A 3700 second\n"
	,
	// Single line without header
	"second 1A 3700"
	,
	// Single line and no valid markers
	"pos marker chr\n"
    "110 fake5 3\n"
	,
	// Too many columns on first row
	"pos marker chr pos\n"
    "110 fake5 3 109.9\n"
	,
	// Too few columns on first row
	"chr pos\n"
	"3N 110"
};

const char* const TEST1_GENOMATRIX_LOADERS[] = {
	// markers cols, tabs, with corner, with endline. "B"
	"corner\tm1\tm2\n"
	"g1\tAA\tTA\n"
	"g2\tAT\tAA\n"
	"oog3\tAA\tTT\n"
	,
	// markers rows, tabs, with corner, with endline. "C"
	"corner\tg1\tg2\toog3\n"
	"m1\tAA\tAT\tAA\n"
	"m2\tTA\tAA\tTT\n"
	,
	// markers rows, spaces, with corner, with endline. "D"
	"coner g1 g2 oog3\n"
	"m1 AA AT AA\n"
	"m2 TA AA TT\n"
	,
	// markers rows, Spaces, with corner, with \r\n endlines. "E"
	"coner g1 g2 oog3\r\n"
	"m1 AA AT AA\r\n"
    "m2 TA AA TT\r\n"
	,
	// markers rows, Mixed line endings and spacings, with corner. "F" 
	"coner\tg1\tg2\toog3\r\n"
	"m1,AA AT AA\n"
	"m2 TA AA,  TT\r\n"
	,
	// Markers rows, tabs, with corner, no endline. "G"
	"corner\tg1\tg2\toog3\n"
	"m1\tAA\tAT\tAA\n"
	"m2\tTA\tAA\tTT"
	,
	// Markers rows, spaces, with corner, no endline. "H"
	"1 g1 g2 oog3\n"
	"m1 AA AT AA\r\n"
	"m2 TA AA TT"
	,
	// Markers rows, spaces, no corner, no endline. "I"
	" g1 g2 oog3\n"
	"m1 AA AT AA\n"
	"m2 TA AA TT"
	,
	// Markers rows, tabs, no corner, no endline. "J" 
	"\tg1\tg2\toog3\n"
	"m1\tAA\tAT\tAA\n"
	"m2\tTA\tAA\tTT"
	,
	// Markers rows, tabs, no corner or corner space, no endline "K"
	"g1\tg2\toog3\n"
	"m1\tAA\tAT\tAA\n"
    "m2\tTA\tAA\tTT"
	,
	// Markers rows, tabs, no corner, with endline. "L"
	"\tg1\tg2\toog3\n"
	"m1\tAA\tAT\tAA\n"
	"m2\tTA\tAA\tTT\n"
	,
	// Markers columns, tabs, corner, no end line. "M"
	"corner\tm1\tm2\n"
	"g1\tAA\tTA\n"
	"g2\tAT\tAA\n"
	"oog3\tAA\tTT"
	,
	// Markers columns, tabs, no corner or corner space, no endline. "N"
	"m1\tm2\n"
	"g1\tAA\tTA\n"
	"g2\tAT\tAA\n"
	"oog3\tAA\tTT"
	,
	// Markers columns, tabs, no corner, with endline. "O"
	"\tm1\tm2\n"
	"g1\tAA\tTA\n"
	"g2\tAT\tAA\n"
	"oog3\tAA\tTT\n"
	,
	// Markers columns, commas, no corner or corner space, with endline. "P"
	"m1,m2\n"
	"g1,AA,TA\n"
	"g2,AT,AA\n"
	"oog3,AA,TT\n"
	,
	// Markers columns, commas, no corner, with endline, ignore empty extra cell in row. "Q"
	",m1,m2\n"
	"g1,AA,TA\n"
	"g2,AT,AA,\n"
	"oog3,AA,TT\n"
	,
	// Markers columns, commas, no corner, no endline, ignore extra cell in row. "R"
	" ,m1,m2\n"
	"g1,AA,TA\n"
	"g2,AT,AA,BB\n"
	"oog3,AA, TT"
	,
	// Markers columns, commas, no corner, with endline, slashpairs. "S"
	",m1,m2\n"
	"g1,A/A,T/A\n"
	"g2,A/T,A/A,\n"
	"oog3,A/A,T/T\n"
	,
	// Markers rows, commas, no corner, with endline, slashpairs. "T"
	"\tg1\tg2\toog3\n"
	"m1\tA/A\tA/T\tA/A\n"
	"m2\tT/A\tA/A\tT/T\n"
	,
	// Markers rows, commas, no corner, with endline, slashpairs, rearranged markers. "U"
	"\tg1\tg2\toog3\n"
	"m2\tT/A\tA/A\tT/T\n"
	"m1\tA/A\tA/T\tA/A\n"
	,
	// Markers columns, commas, no corner, with endline, slashpairs, rearranged markers "V"
	" ,m2,m1\n"
	"g1,TA,AA\n"
	"g2,AA,AT,\n"
	"oog3,TT,AA\n"
	,
	// Markers columns, commas, no corner, rearranged markers, fake marker that should be ignored. "W"
	" ,m2,m1,omg\n"
	"g1,TA,AA,BB\n"
	"g2,AA,AT,BA\n"
	"oog3,TT,AA,CA\n"
	,
	// Markers columns, commas, no corner, rearranged markers, central fake marker that should be ignored. "X"
	" ,m2,omg,m1\n"
	"g1,TA,BB,AA\n"
	"g2,AA,BA,AT\n"
	"oog3,TT,CA,AA\n"
	,
	// Markers rows, commas, no corner or corner space, fake marker that should be ignored. "Y"
	"g1\tg2\toog3\n"
	"m1\tAA\tAT\tAA\n"
	"notthisone\t:)\tAA\tTT\n"
	"m2\tTA\tAA\tTT\n"
	,
	// Markers columns, blank genotype names, commas, rearranged markers. "Z"
	" m2,m1\n"
	" TA,AA\n"
	" AA,AT\n"
	" TT,AA\n"
	,
    // Markers as rows, no genotype names "ZA"
    "\n"
    "m2\tT/A\tA/A\tT/T\n"
    "m1\tA/A\tA/T\tA/A\n"
    ,
    // Markers as rows, no genotype names v2 "ZB"
    "m2\tT/A\tA/A\tT/T\n"
    "m1\tA/A\tA/T\tA/A\n"
    ,
    // Markers as columns, no genotype names "ZC"
    ",m1,m2\n"
    ",A/A,T/A\n"
    ",A/T,A/A,\n"
    ",A/A,T/T\n"
    ,
    // Markers as columns, no genotype names v2 "ZD"
    "corner,m1,m2\n"
    ",A/A,T/A\n"
    ",A/T,A/A,\n"
    ",A/A,T/T\n"
    ,
	// Markers rows, single marker with header, with corner.
	"1 g1 g2 oog3\n"
    "m1 AA AT AA\n"
	,
	// Markers rows, single row no header
	"m2 AA TA AA\n"
	,
	// Markers rows, single body column with header, no corner or corner space
	"g1\n"
	"m2\tT/A\n"
	"m1\tA/A"
	,
	// Markers rows, single body column no header
	"m2\tT/A\n"
    "m1\tA/A"
	,
	// Markers columns, single row with header, with corner
	"20 ,m1,m2\n"
	"g1,AA,TA\n"
	//,
	// markers columns single marker noheader -> fails, can't guarantee marker order.
	// "g1,AA,TA\n"
	,
	// Markers columns, single column and single row with header 
	"m2\n"
	"g1,AA\n"
	,
    // Markers rows, single row no header
    "m1 0 - 2"
    ,
    // Markers as rows, single row, no body
    " g1 g2 oog3"
    ,
    // Markers as columns, single row, no body
    " m1 m2"
    ,
    // Markers as rows, single column, no body
    "\nm1\nm2"
    ,
    // Markers as rows, no marker name matches
    "\tg1\tg2\toog3\n"
    "m24\tT/A\tA/A\tT/T\n"
    "m14\tA/A\tA/T\tA/A\n"
    ,
    // Markers as columns, no marker name matches
    ",m14,m24\n"
    "g1,A/A,T/A\n"
    "g2,A/T,A/A,\n"
    "oog3,A/A,T/T\n"
    ,
	// Normal: counts as body
	" ,m2,m1\n"
	"g1,1,1\n"
	"g2,0,1\n"
	"oog3,1,2\n"
	,
	// Normal: IUPAC encodings as body
	" ,m1,m2\n"
	"g1,G,A\n"
	"g2,T,C\n"
	"oog3,R,R\n"
	"oog4,R,R\n"
	"oog5,Y,Y\n"
	"oog6,Y,Y\n"
	"oog7,M,M\n"
	"oog8,M,M\n"
	"oog9,K,K\n"
	"oog10,K,K\n"
	"oogA,S,S\n"
	"oogB,S,S\n"
	"oogC,W,W\n"
	"oogD,W,W\n"
	"oogBb,S,S\n"
};
const char TEST1_2MARKER_MAP[] = "m1 1 10\n" "m2 1 20";

const char TEST1_EFFLOADERS_MAP[] = "marker chr pos\n" "first 3 3.2\n" "second 3 43.2\n";
const char* const TEST1_EFFLOADERS[] = {
	// Regular most-expected file 
	"marker allele eff\n"
	"first A 0.5\n"
	"first T -0.5\n"
	"second A 1e2\n"
	"second T 1e-2\n"
	,
	// Tabs, no final newline
	"marker\tallele\teff\n"
	"first\tA\t0.5\n"
	"first\tT\t-0.5\n"
	"second\tA\t1e2\n"
	"second\tT\t1e-2"
	,
	// Mixed spacing 
	"marker, allele, eff\n"
	"first \tA 0.5\n"
	"first \tT -0.5\n"
	"second \tA 1e2\n"
	"second T,,,1e-2\n"
	,
	// Rearranged rows 
	"marker allele eff\n"
	"second A 1e2\n"
	"first A 0.5\n"
	"first T -0.5\n"
	"second T 1e-2"
	,
	// No header
	"second A 1e2\n"
	"first A 0.5\n"
	"first T -0.5\n"
	"second T 1e-2\n"
	,
	// Rearranged columns
	"marker eff allele \n"
	"second 1e2 A\n"
	"first 0.5 A\n"
	"second 1e-2 T\n"
    "first -0.5 T\n"
	,
	// assorted alleles
	"marker allele eff\n"
	"second 8 1e-2\n"
	"first A 0.5\n"
	"first T -0.5\n"
	"second T 1e2\n"
	,
	// only one row
	"marker allele eff\n"
    "first A 0.5\n"
	,
	// only one row no header
	"first a 0.5\n"
	,
	// Discard markers not tracked by the simulation
	"marker eff allele \n"
	"second 1e2 A\n"
	"first 0.5 A\n"
	"3rd 1e-2 B\n"
	"first -0.5 T\n"
	,
	// Too many columns on one row 
	"marker eff allele\n"
	"second 1e2 A ?\n"
	"first 0.5 A //comment forgotten here\n"
	"second 1e-2 T\n"
	"first -0.5 T\n"
	,
	// Duplicated marker/allele pair
	"marker eff allele\n"
	"second 1e2 A\n"
	"first 0.5 T\n"
	"second 1e-2 T\n"
	"first -0.5 T\n"
	,
	// File with no valid lines
	"marker eff allele\n"
	"first 0.5 A //comment forgotten here\n"
	"fifth 1e-2 T\n"
};

float calculate_heterozygosity(SimData* d, GroupNum group_number);
int compareFiles(char* f1, char* f2); // returns 0 if matching
int compareFileToString(char* filename, const char* target);
int compareRepeatingFileToTable(char* filename, unsigned int expectedNRows, const char* header, const char* body);

struct MultiIDSet just_load(SimData* d);
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
int test_targeted_crossing(SimData* d, GroupNum g1);

int test_deletors(SimData *d, GroupNum g0);

int test_block_generator(SimData *d);

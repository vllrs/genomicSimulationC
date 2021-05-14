#include "sim-test.h"
//#include <string.h>

float calculate_heterozygosity(SimData* d, int group_number) {
	int hetcount = 0;
	int gn = get_group_size(d, group_number);
	char** galleles = get_group_genes(d, group_number, gn);
	
	// uses subjects as the first index
	for (int i = 0; i < gn; i++) {
		for (int j = 0; j < d->n_markers; j += 2) {
			if (galleles[i][j] != galleles[i][j + 1]) {
				hetcount += 1;
			}
		}
	}
	
	free(galleles);
	return (float) hetcount / (gn * d->n_markers);
}

int test_loaders(SimData* d) {
	FILE* fp; 
	if ((fp = fopen("a-test.txt", "w")) == NULL) {
		fprintf(stderr, "Failed to create file.\n");
		exit(1);
	}
	fwrite(HELPER_GENOTYPES, sizeof(char), strlen(HELPER_GENOTYPES), fp);
	fclose(fp);
	if ((fp = fopen("a-test-map.txt", "w")) == NULL) {
		fprintf(stderr, "Failed to create file.\n");
		exit(1);
	}
	fwrite(HELPER_MAP, sizeof(char), strlen(HELPER_MAP), fp);
	fclose(fp);
	if ((fp = fopen("a-test-eff.txt", "w")) == NULL) {
		fprintf(stderr, "Failed to create file.\n");
		exit(1);
	}
	fwrite(HELPER_EFF, sizeof(char), strlen(HELPER_EFF), fp);
	fclose(fp);
	
	int g0 = load_all_simdata(d, "a-test.txt", "a-test-map.txt", "a-test-eff.txt");
	
	remove("a-test.txt");
	remove("a-test-map.txt");
	remove("a-test-eff.txt");
	
	assert(d->n_markers == 3); // all markers loaded
	assert(strcmp(d->markers[0], "m1") == 0); // all markers ordered right
	assert(strcmp(d->markers[1], "m2") == 0);
	assert(strcmp(d->markers[2], "m3") == 0);
	printf("...SNP marker names loaded correctly\n");

	assert(d->map.n_chr == 3); // correct number of chromosomes
	assert(d->map.chr_ends[0] == 0);
	assert(d->map.chr_ends[1] == 2);
	assert(d->map.chr_ends[2] == 2);
	assert(abs(d->map.chr_lengths[0] - 3.1) < TOL);
	// @other lengths don't matter?
	assert(abs(d->map.positions[0].position - 5.2) < TOL);
	assert(d->map.positions[0].chromosome == 1);
	assert(abs(d->map.positions[1].position - 8.3) < TOL);
	assert(d->map.positions[1].chromosome == 1);
	assert(abs(d->map.positions[2].position - 15) < TOL);
	assert(d->map.positions[2].chromosome == 3);
	printf("...genome map loaded correctly\n");
	
	assert(d->e.effects.rows == 2);
	assert(d->e.effects.cols == 3);
	// don't mind which order the effects are in so just figure it out, then check based on that.
	int apos = 0;
	if (d->e.effect_names[0] == 'A') {
		assert(d->e.effect_names[0] == 'A');
		assert(d->e.effect_names[1] == 'T');
	} else {
		apos = 1;
		assert(d->e.effect_names[0] == 'T');
		assert(d->e.effect_names[1] == 'A');
	}
	assert(abs(d->e.effects.matrix[apos][0] - (-0.8)) < TOL);
	assert(abs(d->e.effects.matrix[apos][1] - (-0.1)) < TOL);
	assert(abs(d->e.effects.matrix[apos][2] - (0.1)) < TOL);
	int tpos = 1 - apos;
	assert(abs(d->e.effects.matrix[tpos][0] - (0.9)) < TOL);
	assert(abs(d->e.effects.matrix[tpos][1] - (-0.5)) < TOL);
	assert(abs(d->e.effects.matrix[tpos][2] - (-0.1)) < TOL);
	printf("...effect values loaded correctly\n");
	
	assert(d->m != NULL);
	assert(d->current_id == 6);
	assert(d->m->n_markers == 3);
	assert(d->m->n_subjects == 6);
	assert(d->m->ids[0] == 1);
	assert(d->m->ids[1] == 2);
	assert(d->m->ids[2] == 3);
	assert(d->m->ids[3] == 4);
	assert(d->m->ids[4] == 5);
	assert(d->m->ids[5] == 6);
	assert(d->m->pedigrees[0][0] == 0 && d->m->pedigrees[1][0] == 0);
	assert(d->m->pedigrees[0][1] == 0 && d->m->pedigrees[1][1] == 0);
	assert(d->m->pedigrees[0][2] == 0 && d->m->pedigrees[1][2] == 0);
	assert(d->m->pedigrees[0][3] == 0 && d->m->pedigrees[1][3] == 0);
	assert(d->m->pedigrees[0][4] == 0 && d->m->pedigrees[1][4] == 0);
	assert(d->m->pedigrees[0][5] == 0 && d->m->pedigrees[1][5] == 0);
	assert(d->m->groups[0] == g0);
	assert(d->m->groups[1] == g0);
	assert(d->m->groups[2] == g0);
	assert(d->m->groups[3] == g0);
	assert(d->m->groups[4] == g0);
	assert(d->m->groups[5] == g0);
	assert(g0 > 0);

	// might not be important that thele load in this order but we'll ask for it anyway. 
	assert(strcmp(d->m->subject_names[0], "G01") == 0);
	assert(strcmp(d->m->subject_names[1], "G02") == 0);
	assert(strcmp(d->m->subject_names[2], "G03") == 0);
	assert(strcmp(d->m->subject_names[3], "G04") == 0);
	assert(strcmp(d->m->subject_names[4], "G05") == 0);
	assert(strcmp(d->m->subject_names[5], "G06") == 0);
	assert(strncmp(d->m->alleles[0],"TTAATT", 6) == 0); // G01
	assert(strncmp(d->m->alleles[1],"TTAATT", 6) == 0); // G02
	assert(strncmp(d->m->alleles[2],"TTAATA", 6) == 0); // G03
	assert(strncmp(d->m->alleles[3],"TAAATA", 6) == 0); // G04
	assert(strncmp(d->m->alleles[4],"TTTTTT", 6) == 0); // G05
	assert(strncmp(d->m->alleles[5],"ATAATT", 6) == 0); // G06
	printf("...genotypes loaded correctly\n");
	
	return g0;
}

int test_effect_calculators(SimData *d, int g0) {
	DecimalMatrix dec = calculate_fitness_metric_of_group(d, g0);
	
	assert(dec.rows == 1);
	assert(dec.cols == 6);
	assert(abs(dec.matrix[0][0] - 1.4) < TOL);
	assert(abs(dec.matrix[0][1] - 1.4) < TOL);
	assert(abs(dec.matrix[0][2] - 1.6) < TOL);
	assert(abs(dec.matrix[0][3] - (-0.1)) < TOL);
	assert(abs(dec.matrix[0][4] - 0.6) < TOL);
	assert(abs(dec.matrix[0][5] - (-0.3)) < TOL);
	printf("...GEBVs calculated correctly\n");
	
	delete_dmatrix(&dec);
	return 0;
}

int test_optimal_calculators(SimData *d) {
	char* ig = calculate_ideal_genotype(d);
	assert(ig[0] == 'T');
	assert(ig[1] == 'A');
	assert(ig[2] == 'A');
	free(ig);
	
	double optimal = calculate_optimal_gebv(d);
	assert(abs(optimal - 1.8) < TOL);
	printf("...Optimal genotype and GEBV calculated correctly\n");
	
	return 0;
}

int test_crossing(SimData *d, int g0) {
	int gall = test_crossing_unidirectional(d, g0);
	
	FILE* fp;
	if ((fp = fopen("a-test-plan.txt", "w")) == NULL) {
		fprintf(stderr, "Failed to create file.\n");
		exit(1);
	}
	fwrite(HELPER_PLAN, sizeof(char), strlen(HELPER_PLAN), fp);
	fclose(fp);
	int gfile = test_crossing_from_file(d, "a-test-plan.txt");
	remove("a-test-plan.txt");
	
	int gselfed = test_crossing_selfing(d, gall);
	
	assert(gselfed != g0 && gfile != gall && gfile != g0 && gfile != gselfed);
	printf("...group number allocations are correct\n");
	
	return 0;
}

int test_crossing_unidirectional(SimData *d, int g0) {
	GenOptions g = {.will_name_subjects=TRUE, .subject_prefix="F1", .family_size=1,
		.will_track_pedigree=TRUE, .will_allocate_ids=TRUE,
		.filename_prefix="atestF1", .will_save_pedigree_to_file=TRUE,
		.will_save_effects_to_file=TRUE, .will_save_genes_to_file=TRUE,
		.will_save_to_simdata=TRUE};
	//AlleleMatrix* a = make_all_unidirectional_crosses(&sd, 0, g);
	//sd.m->next_gen = a;
	int g1 = make_all_unidirectional_crosses(d, g0, g);
	
	assert(g1 != g0);
	assert(d->m->n_subjects == 21);
	assert(d->m->n_markers == 3);
	//assert(strcmp(sd.m->name, "F1g") == 0);
	assert(strcmp(d->m->subject_names[6], "F107") == 0);
	assert(strcmp(d->m->subject_names[7], "F108") == 0);
	assert(strcmp(d->m->subject_names[8], "F109") == 0);
	assert(strcmp(d->m->subject_names[9], "F110") == 0);
	assert(strcmp(d->m->subject_names[10], "F111") == 0);
	assert(strcmp(d->m->subject_names[11], "F112") == 0);
	assert(strcmp(d->m->subject_names[12], "F113") == 0);
	assert(strcmp(d->m->subject_names[13], "F114") == 0);
	assert(strcmp(d->m->subject_names[14], "F115") == 0);
	assert(strcmp(d->m->subject_names[15], "F116") == 0);
	assert(strcmp(d->m->subject_names[16], "F117") == 0);
	assert(strcmp(d->m->subject_names[17], "F118") == 0);
	assert(strcmp(d->m->subject_names[18], "F119") == 0);
	assert(strcmp(d->m->subject_names[19], "F120") == 0);
	assert(strcmp(d->m->subject_names[20], "F121") == 0);
	assert(d->m->pedigrees[0][6] == 1 && d->m->pedigrees[1][6] == 2); //assumes you're doing the top triangle of crosses
	assert(d->m->pedigrees[0][7] == 1 && d->m->pedigrees[1][7] == 3);
	assert(d->m->pedigrees[0][8] == 1 && d->m->pedigrees[1][8] == 4);
	assert(d->m->pedigrees[0][9] == 1 && d->m->pedigrees[1][9] == 5);
	assert(d->m->pedigrees[0][10] == 1 && d->m->pedigrees[1][10] == 6);
	assert(d->m->pedigrees[0][11] == 2 && d->m->pedigrees[1][11] == 3);
	assert(d->m->pedigrees[0][12] == 2 && d->m->pedigrees[1][12] == 4);
	assert(d->m->pedigrees[0][13] == 2 && d->m->pedigrees[1][13] == 5);
	assert(d->m->pedigrees[0][14] == 2 && d->m->pedigrees[1][14] == 6);
	assert(d->m->pedigrees[0][15] == 3 && d->m->pedigrees[1][15] == 4);
	assert(d->m->pedigrees[0][16] == 3 && d->m->pedigrees[1][16] == 5);
	assert(d->m->pedigrees[0][17] == 3 && d->m->pedigrees[1][17] == 6);
	assert(d->m->pedigrees[0][18] == 4 && d->m->pedigrees[1][18] == 5);
	assert(d->m->pedigrees[0][19] == 4 && d->m->pedigrees[1][19] == 6);
	assert(d->m->pedigrees[0][20] == 5 && d->m->pedigrees[1][20] == 6);
	assert(strncmp(d->m->alleles[6], "TTAATT", 6) == 0); // G01 x G02
	assert(strncmp(d->m->alleles[9], "TTATTT", 6) == 0); // G01 x G05
	assert(strncmp(d->m->alleles[11], "TTAATT", 6) == 0 || strncmp(d->m->alleles[11], "TTAATA", 6) == 0); 
	printf("...crossed all pairs correctly\n");
	// does not check that crossing over is happening.
	
	return g1;
}

int test_crossing_from_file(SimData *d, char* fname) {
	// check we can load a plan from a file.
	GenOptions g2 = BASIC_OPT;
	g2.will_track_pedigree = TRUE;
	g2.family_size = 2;
	//g2.will_save_pedigree_to_file = TRUE;
	//2.filename_prefix = "atest-dc";
	int bp = make_double_crosses_from_file(d, fname, g2);
	assert(d->m->pedigrees[0][21] == d->m->ids[6] && d->m->pedigrees[1][21] == 17);
	assert(d->m->pedigrees[0][23] == d->m->ids[7] && d->m->pedigrees[1][23] == 21);
	assert(d->m->pedigrees[0][25] == d->m->ids[20] && d->m->pedigrees[1][25] == 9);
	assert(d->m->pedigrees[0][22] == 7 && d->m->pedigrees[1][22] == 17);
	assert(d->m->pedigrees[0][24] == 8 && d->m->pedigrees[1][24] == 21);
	assert(d->m->pedigrees[0][26] == 21 && d->m->pedigrees[1][26] == 9);
	assert(d->m->n_subjects == 27);
	assert(d->m->n_markers == 3);
	printf("...crossed combinations from file correctly\n");
	
	return bp;
}

int test_crossing_selfing(SimData *d, int g1) {
	float h1 = calculate_heterozygosity(d, g1);
	int g1selfed = self_n_times(d, 5, g1, BASIC_OPT);
	float h2 = calculate_heterozygosity(d,  g1selfed);
	
	assert(g1selfed != g1);
	//printf("Heterozygousity reduction from selfing: %f %f\n", h2, h1);
	assert(h1 - h2 > 0);
	assert(d->m->n_subjects == 42);
	assert(d->m->n_markers == 3);
	printf("...selfing function correctly reduced heterozygosity by %f%%\n", (h2-h1)*100);
	
	return g1selfed;
}

int test_deletors(SimData *d, int g0) {
	int ngroups1;
	int* groups1 = get_existing_groups(d, &ngroups1);
	
	delete_group(d, g0);
	int ngroups2;
	int* groups2 = get_existing_groups(d, &ngroups2);
	
	assert(ngroups1 - ngroups2 == 1);
	for (int i = 0; i < ngroups1; ++i) {
		if (i == ngroups1 - 1 || groups1[i] != groups2[i]) {
			assert(groups1[i] == g0);
			break;
		}
	}
	printf("...group of genotypes cleared correctly\n");
	
	delete_simdata(d);
	printf("...SimData cleared correctly\n");
	
	return 0;
}

/* main, for testing. Only uses a small dataset. */
int main(int argc, char* argv[]) {
	printf("Testing functionality ...");

	// test random number generators

	// test matrix operations

	// test SimData loaders
	printf("\nNow testing loader functions:\n");
	SimData sd = EMPTY_SIMDATA;
	int g0 = test_loaders(&sd);
	printf("\t\t-> Loader functions all clear\n");
	
	// test effect calculators
	printf("\nNow testing GEBV calculator:\n");
	test_effect_calculators(&sd, g0);
	test_optimal_calculators(&sd);
	printf("\t\t-> GEBV calculators all clear\n");
	
	// test crossers
	printf("\nNow testing crossing functions:\n");
	test_crossing(&sd, g0);
	printf("\t\t-> Crossing functions all clear\n");
	
	// test SimData deletors.
	printf("\nNow testing deletor functions:\n");
	test_deletors(&sd, g0);
	printf("\t\t-> Deletor functions all clear\n");
	
	printf("\n------- All tests passed. -------\n");
	return 0;
}
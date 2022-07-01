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
    assert(fabs(d->map.chr_lengths[0] - 3.1) < TOL);
	// @other lengths don't matter?
    assert(fabs(d->map.positions[0].position - 5.2) < TOL);
	assert(d->map.positions[0].chromosome == 1);
    assert(fabs(d->map.positions[1].position - 8.3) < TOL);
	assert(d->map.positions[1].chromosome == 1);
    assert(fabs(d->map.positions[2].position - 15) < TOL);
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
    assert(fabs(d->e.effects.matrix[apos][0] - (-0.8)) < TOL);
    assert(fabs(d->e.effects.matrix[apos][1] - (-0.1)) < TOL);
    assert(fabs(d->e.effects.matrix[apos][2] - (0.1)) < TOL);
	int tpos = 1 - apos;
    assert(fabs(d->e.effects.matrix[tpos][0] - (0.9)) < TOL);
    assert(fabs(d->e.effects.matrix[tpos][1] - (-0.5)) < TOL);
    assert(fabs(d->e.effects.matrix[tpos][2] - (-0.1)) < TOL);
	printf("...effect values loaded correctly\n");

	assert(d->current_id == 6);
    assert(d->m); // != NULL
	assert(d->m->n_markers == 3);
	assert(d->m->n_genotypes == 6);
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
	assert(strcmp(d->m->names[0], "G01") == 0);
	assert(strcmp(d->m->names[1], "G02") == 0);
	assert(strcmp(d->m->names[2], "G03") == 0);
	assert(strcmp(d->m->names[3], "G04") == 0);
	assert(strcmp(d->m->names[4], "G05") == 0);
	assert(strcmp(d->m->names[5], "G06") == 0);
	assert(strncmp(d->m->alleles[0],"TTAATT", 6) == 0); // G01
	assert(strncmp(d->m->alleles[1],"TTAATT", 6) == 0); // G02
	assert(strncmp(d->m->alleles[2],"TTAATA", 6) == 0); // G03
	assert(strncmp(d->m->alleles[3],"TAAATA", 6) == 0); // G04
	assert(strncmp(d->m->alleles[4],"TTTTTT", 6) == 0); // G05
	assert(strncmp(d->m->alleles[5],"ATAATT", 6) == 0); // G06
	printf("...genotypes loaded correctly\n");

	/*unsigned int g1 = load_more_transposed_genes_to_simdata(d, "a-test.txt");

	assert(d->m != NULL);
	assert(d->current_id == 12);
	assert(d->m->n_markers == 3);
	assert(d->m->n_genotypes == 12);
	assert(d->m->groups[0] == g0);
	assert(d->m->groups[1] == g0);
	assert(d->m->groups[2] == g0);
	assert(d->m->groups[3] == g0);
	assert(d->m->groups[4] == g0);
	assert(d->m->groups[5] == g0);
	assert(g0 > 0);
	assert(d->m->groups[6] == g1);
	assert(d->m->groups[7] == g1);
	assert(d->m->groups[8] == g1);
	assert(d->m->groups[9] == g1);
	assert(d->m->groups[10] == g1);
	assert(d->m->groups[11] == g1);
	assert(g1 > 0);
	assert(g1 != g0);

	delete_group(d, g1);*/

	return g0;
}

int test_effect_calculators(SimData *d, int g0) {
	DecimalMatrix dec = calculate_group_bvs(d, g0);

	assert(dec.rows == 1);
	assert(dec.cols == 6);
    assert(fabs(dec.matrix[0][0] - 1.4) < TOL);
    assert(fabs(dec.matrix[0][1] - 1.4) < TOL);
    assert(fabs(dec.matrix[0][2] - 1.6) < TOL);
    assert(fabs(dec.matrix[0][3] - (-0.1)) < TOL);
    assert(fabs(dec.matrix[0][4] - 0.6) < TOL);
    assert(fabs(dec.matrix[0][5] - (-0.3)) < TOL);
	printf("...GEBVs calculated correctly\n");

	delete_dmatrix(&dec);
	return 0;
}

int test_optimal_calculators(SimData *d) {
	char* ig = calculate_optimal_alleles(d);
	assert(ig[0] == 'T');
	assert(ig[1] == 'A');
	assert(ig[2] == 'A');
	free(ig);

	double optimal = calculate_optimum_bv(d);
    assert(fabs(optimal - 1.8) < TOL);

	double unoptimal = calculate_minimum_bv(d);
    assert(fabs(unoptimal + 2.8) < TOL);

	printf("...Optimal genotype and GEBV calculated correctly\n");

	return 0;
}

int test_crossing(SimData *d, int g0) {
	int gall = test_crossing_unidirectional(d, g0);

    test_crossing_randomly(d, g0);

	FILE* fp;
	if ((fp = fopen("a-test-plan.txt", "w")) == NULL) {
		fprintf(stderr, "Failed to create file.\n");
		exit(1);
	}
	fwrite(HELPER_PLAN, sizeof(char), strlen(HELPER_PLAN), fp);
	fclose(fp);
	int gfile = test_crossing_from_file(d, "a-test-plan.txt");
	remove("a-test-plan.txt");

    int gselfed = test_crossing_selfing(d, g0);

	assert(gselfed != g0 && gfile != gall && gfile != g0 && gfile != gselfed);
	printf("...group number allocations are correct\n");

	return 0;
}

int test_crossing_unidirectional(SimData *d, int g0) {
	GenOptions g = {.will_name_offspring=TRUE, .offspring_name_prefix="F1", .family_size=1,
		.will_track_pedigree=TRUE, .will_allocate_ids=TRUE,
		.filename_prefix="atestF1", .will_save_pedigree_to_file=TRUE,
		.will_save_bvs_to_file=TRUE, .will_save_alleles_to_file=TRUE,
		.will_save_to_simdata=TRUE};
	//AlleleMatrix* a = make_all_unidirectional_crosses(&sd, 0, g);
	//sd.m->next_gen = a;
	int g1 = make_all_unidirectional_crosses(d, g0, g);

	assert(g1 != g0);
	assert(d->m->n_genotypes == 21);
	assert(d->m->n_markers == 3);
	//assert(strcmp(sd.m->name, "F1g") == 0);
	assert(strcmp(d->m->names[6], "F107") == 0);
	assert(strcmp(d->m->names[7], "F108") == 0);
	assert(strcmp(d->m->names[8], "F109") == 0);
	assert(strcmp(d->m->names[9], "F110") == 0);
	assert(strcmp(d->m->names[10], "F111") == 0);
	assert(strcmp(d->m->names[11], "F112") == 0);
	assert(strcmp(d->m->names[12], "F113") == 0);
	assert(strcmp(d->m->names[13], "F114") == 0);
	assert(strcmp(d->m->names[14], "F115") == 0);
	assert(strcmp(d->m->names[15], "F116") == 0);
	assert(strcmp(d->m->names[16], "F117") == 0);
	assert(strcmp(d->m->names[17], "F118") == 0);
	assert(strcmp(d->m->names[18], "F119") == 0);
	assert(strcmp(d->m->names[19], "F120") == 0);
	assert(strcmp(d->m->names[20], "F121") == 0);
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


	// @should check that the saved files are correct too before deleting them
    remove("atestF1-bv.txt");
	remove("atestF1-genotype.txt");
	remove("atestF1-pedigree.txt");

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
	assert(d->m->n_genotypes == 27);
	assert(d->m->n_markers == 3);
	printf("...crossed combinations from file correctly\n");

	return bp;
}

int test_crossing_selfing(SimData *d, int g1) {
    int oldsize = d->m->n_genotypes;
    GenOptions opt = BASIC_OPT;
    opt.will_track_pedigree = TRUE;
	float h1 = calculate_heterozygosity(d, g1);
    int g1selfed = self_n_times(d, 5, g1, opt);
	float h2 = calculate_heterozygosity(d,  g1selfed);

	assert(g1selfed != g1);
	//printf("Heterozygousity reduction from selfing: %f %f\n", h2, h1);
	assert(h1 - h2 > 0);
    assert(d->m->n_genotypes == oldsize + 6);
	assert(d->m->n_markers == 3);
    assert(d->m->groups[oldsize] == g1selfed && d->m->groups[oldsize + 5] == g1selfed && d->m->groups[oldsize + 6] != g1selfed);
    assert(d->m->pedigrees[0][oldsize + 0] == d->m->ids[0] && d->m->pedigrees[1][oldsize + 0] == d->m->ids[0]);
    assert(d->m->pedigrees[0][oldsize + 4] == d->m->ids[4] && d->m->pedigrees[1][oldsize + 4] == d->m->ids[4]);
	printf("...selfing function correctly reduced heterozygosity by %f%%\n", (h2-h1)*100);

    // test doubled haploids
    int g1dhap = make_doubled_haploids(d, g1, opt);
    assert(g1dhap != g1);
    assert(d->m->n_genotypes == oldsize + 2*6);
    assert(d->m->n_markers == 3);
    assert(calculate_heterozygosity(d,  g1dhap) == 0);
    assert(d->m->groups[oldsize + 6] == g1dhap && d->m->groups[oldsize + 11] == g1dhap && d->m->groups[oldsize + 12] != g1dhap);
    assert(d->m->pedigrees[0][oldsize + 6] == d->m->ids[0] && d->m->pedigrees[1][oldsize + 6] == d->m->ids[0]);
    assert(d->m->pedigrees[0][oldsize + 10] == d->m->ids[4] && d->m->pedigrees[1][oldsize + 10] == d->m->ids[4]);
    assert(strcmp(d->m->alleles[oldsize + 6],"TTAATT") == 0);
    assert(strcmp(d->m->alleles[oldsize + 11],"AAAATT") == 0 || strcmp(d->m->alleles[oldsize + 11],"TTAATT") == 0);

    delete_group(d, g1dhap);

    // test cloning
    int g1clones = make_clones(d, g1, TRUE, opt);
    assert(g1clones != g1);
    assert(d->m->n_genotypes == oldsize + 2*6);
    assert(d->m->n_markers == 3);
    assert(fabsf(calculate_heterozygosity(d,  g1clones) - h1) < TOL);
    assert(d->m->groups[oldsize + 6] == g1clones && d->m->groups[oldsize + 11] == g1clones && d->m->groups[oldsize + 12] != g1clones);
    assert(d->m->pedigrees[0][oldsize + 6] == d->m->ids[0] && d->m->pedigrees[1][oldsize + 6] == d->m->ids[0]);
    assert(d->m->pedigrees[0][oldsize + 10] == d->m->ids[4] && d->m->pedigrees[1][oldsize + 10] == d->m->ids[4]);
    assert(strcmp(d->m->alleles[oldsize + 6],"TTAATT") == 0);
    assert(strcmp(d->m->alleles[oldsize + 11],"ATAATT") == 0);
    assert(strcmp(d->m->names[1],d->m->names[oldsize + 7]) == 0);
    assert(strcmp(d->m->names[3],d->m->names[oldsize + 9]) == 0);

    delete_group(d, g1clones);

	return g1selfed;
}

void test_crossing_randomly(SimData *d, int g1) {
    // we test it correctly does its crossing randomly (requiring a bit of human input)
    // we do not test that the genes were correctly crossed
    // created 7 Apr 2022

    GenOptions gopt = BASIC_OPT;
    gopt.will_track_pedigree = TRUE;
    // Test random crossing seems about right
    int g2 = cross_random_individuals( d , g1, 4, 0, gopt);
    int* g2ixs = get_group_indexes(d, g2, -1);

    assert(get_group_size(d, g2) == 4);
    int g2minid = d->m->ids[g2ixs[0]];
    int g2maxid = d->m->ids[g2ixs[3]];
    fprintf(stdout, "Should be random parents: %d, %d, %d, %d\n",
            d->m->pedigrees[0][g2ixs[0]], d->m->pedigrees[0][g2ixs[1]],
            d->m->pedigrees[0][g2ixs[2]], d->m->pedigrees[0][g2ixs[3]]);
    fprintf(stdout, "Should be random parents: %d, %d, %d, %d\n\n",
            d->m->pedigrees[1][g2ixs[0]], d->m->pedigrees[1][g2ixs[1]],
            d->m->pedigrees[1][g2ixs[2]], d->m->pedigrees[1][g2ixs[3]]);

    printf("...crossed randomly within a group\n");

    // Test random crossing between two groups seems about right.
    gopt.family_size = 2;
    int g3 = cross_randomly_between( d, g1, g2, 3, 0, 0, gopt);

    assert(get_group_size(d, g3) == 6);
    assert(d->m->pedigrees[0][g2ixs[3] + 1] == d->m->pedigrees[0][g2ixs[3] + 2]); // family size works
    assert(d->m->pedigrees[1][g2ixs[3] + 1] == d->m->pedigrees[1][g2ixs[3] + 2]); // family size works
    assert(d->m->pedigrees[1][g2ixs[3] + 3] >= g2minid && d->m->pedigrees[1][g2ixs[3] + 3] <= g2maxid ); //right parent groupings
    assert(d->m->pedigrees[0][g2ixs[3] + 3] < g2minid); //right parent groupings
    fprintf(stdout, "Should be random parents: %d, %d, %d\n",
            d->m->pedigrees[0][g2ixs[3] + 1], d->m->pedigrees[0][g2ixs[3] + 3],
            d->m->pedigrees[0][g2ixs[3] + 5]);
    fprintf(stdout, "Should be random parents: %d, %d, %d\n",
            d->m->pedigrees[1][g2ixs[3] + 1], d->m->pedigrees[1][g2ixs[3] + 3],
            d->m->pedigrees[1][g2ixs[3] + 5]);

    delete_group(d, g3);

    gopt.family_size = 1;
    int tempg = split_from_group( d, 1, g2ixs + 1);
    assert(get_group_size(d, tempg) == 1);
    int g4sizetobe = 6;
    int g4 = cross_randomly_between( d, g1, tempg, g4sizetobe, 1, 0, gopt );
    assert(get_group_size(d, g4) == g4sizetobe);
    int combine[2] = {tempg, g2};
    g2 = combine_groups( d, 2, combine );

    assert(get_group_size(d, g4) == g4sizetobe);
    for (int i = 1; i < g4sizetobe; ++i) {
        for (int j = i + 1; j <= g4sizetobe; ++j) {
            assert(d->m->pedigrees[1][g2ixs[3] + i] == d->m->pedigrees[1][g2ixs[3] + j]);// right parent repetition
            assert(d->m->pedigrees[0][g2ixs[3] + i] != d->m->pedigrees[0][g2ixs[3] + j]); // should be different first parents
        }
    }
    delete_group( d, g4 );

    //int g5 = cross_randomly_between( d, g1, g2, 2000, 1, 1, gopt); // expect error
    //delete_group( d, g5 );
    printf("...crossed randomly between two groups\n");
    delete_group( d, g2 );
    free( g2ixs );
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

int test_block_generator(SimData *d) {
	MarkerBlocks b = create_n_blocks_by_chr(d, 2);

	assert(b.num_blocks == 4);
	assert(b.num_markers_in_block[0] == 1);
	assert(b.num_markers_in_block[1] == 1);
	assert(b.num_markers_in_block[2] == 1);
	assert(b.num_markers_in_block[3] == 0);
	assert(b.markers_in_block[0][0] == 0);
	assert(b.markers_in_block[1][0] == 1);
	assert(b.markers_in_block[2][0] == 2);
	assert(b.markers_in_block[3] == NULL);

	printf("...chr slicer correctly deals with multi-marker and single-marker chrs\n");

	/*printf("\nNum blocks: %d", b.num_blocks);
	for (int i = 0; i < b.num_blocks; ++i) {
		printf("\n%d: ", b.num_markers_in_block[i]);
		for (int j = 0; j < b.num_markers_in_block[i]; ++j) {
			printf("%d ", b.markers_in_block[i][j]);
		}
	}
	fflush(stdout);*/

	delete_markerblocks(&b);

	printf("...MarkerBlocks deletor works\n");

	b = create_n_blocks_by_chr(d, 4);

	assert(b.num_blocks == 8);
	assert(b.num_markers_in_block[0] == 1);
	assert(b.num_markers_in_block[1] == 0);
	assert(b.num_markers_in_block[2] == 0);
	assert(b.num_markers_in_block[3] == 1);
	assert(b.num_markers_in_block[4] == 1);
	assert(b.num_markers_in_block[5] == 0);
	assert(b.num_markers_in_block[6] == 0);
	assert(b.num_markers_in_block[7] == 0);
	assert(b.markers_in_block[0][0] == 0);
	assert(b.markers_in_block[1] == NULL);
	assert(b.markers_in_block[2] == NULL);
	assert(b.markers_in_block[3][0] == 1);
	assert(b.markers_in_block[4][0] == 2);
	assert(b.markers_in_block[5] == NULL);
	assert(b.markers_in_block[6] == NULL);
	assert(b.markers_in_block[7] == NULL);

	printf("...chr slicer correctly deals with empty blocks\n");

	delete_markerblocks(&b);

	return 0;
}

int test_data_access(SimData* d, int gp) {
    assert(get_group_size(d, gp) == 6);
    char** alleles = get_group_genes(d, gp, 6);
	assert(strncmp(alleles[0],"TTAATT", 6) == 0); // G01
	assert(strncmp(alleles[1],"TTAATT", 6) == 0); // G02
	assert(strncmp(alleles[2],"TTAATA", 6) == 0); // G03
	assert(strncmp(alleles[3],"TAAATA", 6) == 0); // G04
	assert(strncmp(alleles[4],"TTTTTT", 6) == 0); // G05
	assert(strncmp(alleles[5],"ATAATT", 6) == 0); // G06

    char** names = get_group_names(d, gp, -1);
    assert(strcmp(names[0], "G01") == 0);
	assert(strcmp(names[1], "G02") == 0);
	assert(strcmp(names[2], "G03") == 0);
	assert(strcmp(names[3], "G04") == 0);
	assert(strcmp(names[4], "G05") == 0);
	assert(strcmp(names[5], "G06") == 0);

    unsigned int* ids = get_group_ids(d, gp, -1);
    for (int i = 0; i < 6; ++i) {
        assert(ids[i] == i+1);
    }

    int* indexes = get_group_indexes(d, gp, 6);
    for (int i = 0; i < 6; ++i) {
        assert(indexes[i] == i);
    }

    double* bvs = get_group_bvs(d, gp, 6);
    assert(fabs(bvs[0] - 1.4) < TOL);
    assert(fabs(bvs[1] - 1.4) < TOL);
    assert(fabs(bvs[2] - 1.6) < TOL);
    assert(fabs(bvs[3] - (-0.1)) < TOL);
    assert(fabs(bvs[4] - 0.6) < TOL);
    assert(fabs(bvs[5] - (-0.3)) < TOL);

	// missing the parent and pedigree checks but this group doesn't have info for that anyway

    return 0;
}

// Returns 0 if matching.
int compareFiles(char* f1, char* f2) {
	FILE* fp1 = fopen(f1, "r");
	FILE* fp2 = fopen(f2, "r");

	char c1, c2;

	do {
		c1 = fgetc(fp1);
		c2 = fgetc(fp2);

		if (c1 != c2) return -1;

	} while (c1 != EOF && c2 != EOF);

	if (c1 == EOF && c2 == EOF) return 0;
	else return -1;
}


/* main, for testing. Only uses a small dataset. */
int main(int argc, char* argv[]) {
	printf("Testing functionality ...");

	// test random number generators

	// test matrix operations

	// test SimData loaders
	printf("\nNow testing loader functions:\n");
	SimData* d = create_empty_simdata();
	int g0 = test_loaders(d);
	printf("\t\t-> Loader functions all clear\n");

	// test effect calculators
	printf("\nNow testing GEBV calculator:\n");
	test_effect_calculators(d, g0);
	test_optimal_calculators(d);
	printf("\t\t-> GEBV calculators all clear\n");

	// test crossers
	printf("\nNow testing crossing functions:\n");
	test_crossing(d, g0);
	printf("\t\t-> Crossing functions all clear\n");

	//test blocking
	printf("\nNow testing blocking functions:\n");
	test_block_generator(d);
	printf("\t\t-> Blocking functions all clear\n");

	//test file savers
	printf("\nNow testing saver functions:\n");
	printf("TODO Saver tests not implemented yet\n");

	//test data access functions
	printf("\nNow testing data access functions:\n");
	test_data_access(d, g0);
	printf("TODO Data access tests not fully implemented yet\n");

	// test SimData deletors.
	printf("\nNow testing deletor functions:\n");
	test_deletors(d, g0);
	printf("\t\t-> Deletor functions all clear\n");

	printf("\n------- All tests passed. -------\n");


	//testing new grouping functions
	d = create_empty_simdata();
	g0 = load_all_simdata(d, "./gt_parents_mr2_50-trimto-5000.txt",
			 "./genetic-map_5112-trimto5000.txt",
			 "./qtl_mr2.eff-processed.txt");

	/*printf("\nNow testing split into individuals\n");
	assert(get_group_size(d, g0) == 50);
	int results[50];
	split_into_individuals(d, g0, &results);
	assert(results[0] == 2);
	assert(results[49] == 51);
	assert(get_group_size(d, results[23]) == 1);*/

	// split into families
	/*int nc = 23;
	int crosses[2][nc];
	for (int i = 0; i < nc; ++i) {
        crosses[0][i] = i % 50;
        crosses[1][i] = 3*(i+1) % 50;
	}
	GenOptions gens =  {.will_name_offspring=FALSE, .offspring_name_prefix=NULL, .family_size=243,
		.will_track_pedigree=TRUE, .will_allocate_ids=TRUE,
		.filename_prefix="testcross", .will_save_pedigree_to_file=FALSE,
		.will_save_bvs_to_file=FALSE, .will_save_alleles_to_file=FALSE,
		.will_save_to_simdata=TRUE};
	int g1 = cross_these_combinations(d,nc,crosses,gens);

    assert(get_group_size(d, g1) == nc*243);
	int results[nc];
	//split_into_families(d, g1, results);
	split_into_halfsib_families(d,g1,2, results);
	int r1;
	int** r2 = get_existing_group_counts(d, &r1);

    assert(results[0] == 3);
	assert(get_group_size(d, results[r1 - 2]) == 243);
	//assert(get_group_size(d, results[0]) == 2*243);*/

    //int gns[6];
    //get_n_new_group_nums(d,6,gns);

    /*int sequence[100];
    for (int i = 0; i < 100; ++i) {
        sequence[i] = i;
    }
    shuffle_up_to(sequence,100,20);*/

    //int g0b = split_randomly_into_two(d, g0);
    //int g0b = split_evenly_into_two(d, g0);
    //int g0bs[50];
    //split_randomly_into_n(d,g0,5,g0bs);
    //split_evenly_into_n(d,g0,3,g0bs);
    /*double probs[3];
    probs[0] = 0.1;
    probs[1] = 0.95;
    split_by_probabilities_into_n(d,g0,3,probs,g0bs);*/

    /*// there's too little for my buckets~
    int buckets[3];
    buckets[0] = 46;
    buckets[1] = 13;
    split_by_specific_counts_into_n(d,g0,3,buckets,g0bs);*/

    //int a = get_group_size(d, g0);
    //a = get_group_size(d,g0b);
    /*// all genotypes transferred check
    int a = 0;
    for (int i = 0; i < 5; ++i) {
        a += get_group_size(d, g0bs[i]);
    }*/

    // This should not segfault
    int m = split_by_bv(d, g0, 55, 1);
    assert(get_group_size(d, m) == 50);


	printf("\nAll done\n");
	delete_simdata(d);

	/*clock_t c;

	printf("\n--------Timing tests--------------\n");
	c = clock();
	SimData* sd = create_empty_simdata();
	int fg0 = load_all_simdata(sd, "./gt_parents_mr2_50-trimto-5000.txt",
			 "./genetic-map_5112-trimto5000.txt",
			 "./qtl_mr2.eff-processed.txt");
    c = clock() - c;
	printf("Loading took %f seconds to run\n", (double)c / CLOCKS_PER_SEC);

	c = clock();
	GenOptions g = {.will_name_offspring=FALSE, .offspring_name_prefix=NULL, .family_size=1,
		.will_track_pedigree=TRUE, .will_allocate_ids=TRUE,
		.filename_prefix="testcross", .will_save_pedigree_to_file=FALSE,
		.will_save_bvs_to_file=FALSE, .will_save_alleles_to_file=FALSE,
		.will_save_to_simdata=TRUE};
    int f = cross_random_individuals(sd, fg0, 100000, g);
	c = clock() - c;
	printf("Random crossing took %f seconds to run\n", (double)c / CLOCKS_PER_SEC);

	//int ngroups = 0;
	//int** egroups = get_existing_group_counts(sd, &ngroups);
	//for (int i = 0; i < ngroups; ++i) {
	//	printf("Group %d has %d members\n", egroups[0][i], egroups[1][i]);
	//}
	//free(egroups[0]);
	//free(egroups[1]);
	//free(egroups);

	c = clock();
	FILE* fp = fopen("testcross-alleles.txt","w");
	save_group_alleles(fp, sd, f);
	fclose(fp);
	c = clock() - c;
	printf("Saving genotypes to file took %f seconds to run\n", (double)c / CLOCKS_PER_SEC);

	c = clock();
	fp = fopen("testcross-bvs.txt","w");
	save_group_bvs(fp, sd, f);
	fclose(fp);
	c = clock() - c;
	printf("Saving GEBVs to file took %f seconds to run\n", (double)c / CLOCKS_PER_SEC);*/

	//delete_simdata(sd);


	return 0;
}

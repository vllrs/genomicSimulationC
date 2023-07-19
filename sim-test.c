#include "sim-test.h"

float calculate_heterozygosity(SimData* d, int group_number) {
	int hetcount = 0;
	int gn = get_group_size(d, group_number);
    char* galleles[gn];
    get_group_genes(d, group_number, gn, galleles);

	// uses subjects as the first index
	for (int i = 0; i < gn; i++) {
		for (int j = 0; j < d->n_markers; j += 2) {
			if (galleles[i][j] != galleles[i][j + 1]) {
				hetcount += 1;
			}
		}
	}

	return (float) hetcount / (gn * d->n_markers);
}

int test_loaders2(SimData* d) {
    // Test genotype matrix loading
    /*
name g1 g2 g3
m1 AA AA AA
m2 AT AA TT
    */
    FILE* fp;
    const char HEADER_TAB[] = "\tg1\tg2\tg3";
    const char HEADER_SP[] = " g1 g2 g3";
    const char BODY1_TAB[] = "m1\tAA\tAA\tAA";
    const char BODY1_SP[] = "m1 AA AA AA";
    const char BODY2_TAB[] = "m2\tAT\tAA\tTT";
    const char BODY2_SP[] = "m2 AT AA TT";
    if ((fp = fopen("1test.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    fwrite("name", sizeof(char), strlen("name"), fp);
    fwrite(HEADER_TAB, sizeof(char), strlen(HEADER_TAB), fp);
    fputc('\n', fp);
    fwrite(BODY1_TAB, sizeof(char), strlen(BODY1_TAB), fp);
    fputc('\n', fp);
    fwrite(BODY2_TAB, sizeof(char), strlen(BODY2_TAB), fp);
    fclose(fp);
    //fputc('\n', fp);

    // @@

    // Test replacing tabs with spaces
    if ((fp = fopen("1test.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    fwrite("name", sizeof(char), strlen("name"), fp);
    fwrite(HEADER_SP, sizeof(char), strlen(HEADER_SP), fp);
    fputc('\n', fp);
    fwrite(BODY1_SP, sizeof(char), strlen(BODY1_SP), fp);
    fputc('\n', fp);
    fwrite(BODY2_SP, sizeof(char), strlen(BODY2_SP), fp);
    fclose(fp);

    // @@

    // test extra free line at the end
    if ((fp = fopen("1test.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    fwrite("1", sizeof(char), strlen("name"), fp);
    fwrite(HEADER_SP, sizeof(char), strlen(HEADER_SP), fp);
    fputc('\n', fp);
    fwrite(BODY1_SP, sizeof(char), strlen(BODY1_SP), fp);
    fputc('\n', fp);
    fwrite(BODY2_SP, sizeof(char), strlen(BODY2_SP), fp);
    fclose(fp);
    fputc('\n', fp);

    // @@


    return 0;
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

int test_labels(SimData *d, int g0) {
    // g0: contents of "a-test.txt", contains 6 genotypes

    // Can create a first label.
    assert(d->n_labels == 0);
    const int label1value = 0;
    const int label1 = create_new_label(d, label1value);
    const int label1index = get_index_of_label(d, label1);
    assert(d->n_labels == 1);
    assert(label1 == 1); // expected
    assert(label1index == 0); // expected
    assert(d->label_ids != NULL);
    assert(d->label_ids[label1index] == label1);
    assert(d->label_defaults != NULL);
    assert(d->label_defaults[label1index] == label1value);
    assert(d->m->n_labels == 1);
    assert(d->m->labels != NULL);
    assert(d->m->labels[label1index] != NULL);
    assert(d->m->labels[label1index][0] == label1value);
    assert(d->m->labels[label1index][1] == label1value);
    assert(d->m->labels[label1index][2] == label1value);
    assert(d->m->labels[label1index][3] == label1value);
    assert(d->m->labels[label1index][4] == label1value);
    assert(d->m->labels[label1index][5] == label1value);

    // Can create a second label
    const int label2value = 10;
    const int label2 = create_new_label(d, label2value);
    const int label2index = get_index_of_label(d, label2);
    assert(d->n_labels == 2);
    assert(label2 == 2); // expected
    assert(label2index == 1); // expected
    assert(d->label_ids != NULL);
    assert(d->label_ids[label1index] == label1);
    assert(d->label_ids[label2index] == label2);
    assert(d->label_defaults != NULL);
    assert(d->label_defaults[label1index] == label1value);
    assert(d->label_defaults[label2index] == label2value);
    assert(d->m->n_labels == 2);
    assert(d->m->labels != NULL);
    assert(d->m->labels[label1index] != NULL);
    assert(d->m->labels[label2index] != NULL);
    assert(d->m->labels[label1index][0] == label1value);
    assert(d->m->labels[label1index][1] == label1value);
    assert(d->m->labels[label1index][2] == label1value);
    assert(d->m->labels[label1index][3] == label1value);
    assert(d->m->labels[label1index][4] == label1value);
    assert(d->m->labels[label1index][5] == label1value);

    assert(d->m->labels[label2index][0] == label2value);
    assert(d->m->labels[label2index][1] == label2value);
    assert(d->m->labels[label2index][2] == label2value);
    assert(d->m->labels[label2index][3] == label2value);
    assert(d->m->labels[label2index][4] == label2value);
    assert(d->m->labels[label2index][5] == label2value);

    // Can change a default
    const int newlabel2default = 8;
    set_label_default(d,label2,newlabel2default);
    assert(d->n_labels == 2);
    assert(d->label_ids != NULL);
    assert(d->label_ids[label1index] == label1);
    assert(d->label_ids[label2index] == label2);
    assert(d->label_defaults != NULL);
    assert(d->label_defaults[label1index] == label1value);
    assert(d->label_defaults[label2index] == newlabel2default);
    assert(d->m->n_labels == 2);
    assert(d->m->labels != NULL);
    assert(d->m->labels[label1index] != NULL);
    assert(d->m->labels[label2index] != NULL);
    assert(d->m->labels[label2index][0] == label2value);
    assert(d->m->labels[label2index][1] == label2value);
    assert(d->m->labels[label2index][2] == label2value);
    assert(d->m->labels[label2index][3] == label2value);
    assert(d->m->labels[label2index][4] == label2value);
    assert(d->m->labels[label2index][5] == label2value);


    // Can set values of a label
    const int newlabel1value = 2;
    set_labels_to_const(d,g0,label1,newlabel1value);
    assert(d->n_labels == 2);
    assert(d->m->n_labels == 2);
    assert(d->m->labels != NULL);
    assert(d->m->labels[label1index] != NULL);
    assert(d->m->labels[label2index] != NULL);
    assert(d->m->labels[label1index][0] == newlabel1value);
    assert(d->m->labels[label1index][1] == newlabel1value);
    assert(d->m->labels[label1index][2] == newlabel1value);
    assert(d->m->labels[label1index][3] == newlabel1value);
    assert(d->m->labels[label1index][4] == newlabel1value);
    assert(d->m->labels[label1index][5] == newlabel1value);

    // Can set values of a label (more labels than needed)
    const unsigned int newlabel2values[] = {11, 12, 13, 14, 15, 16, 17};
    set_labels_to_values(d, g0, 0, label2, 7, newlabel2values);
    assert(d->n_labels == 2);
    assert(d->m->n_labels == 2);
    assert(d->m->labels != NULL);
    assert(d->m->labels[label1index] != NULL);
    assert(d->m->labels[label2index] != NULL);
    assert(d->m->labels[label2index][0] == newlabel2values[0]);
    assert(d->m->labels[label2index][1] == newlabel2values[1]);
    assert(d->m->labels[label2index][2] == newlabel2values[2]);
    assert(d->m->labels[label2index][3] == newlabel2values[3]);
    assert(d->m->labels[label2index][4] == newlabel2values[4]);
    assert(d->m->labels[label2index][5] == newlabel2values[5]);

    // Can set value of a label (only some values set)
    const unsigned int newerlabel1values[] = {0,1,-1,100};
    set_labels_to_values(d, g0, 1, label1, 4, newerlabel1values);
    assert(d->n_labels == 2);
    assert(d->m->n_labels == 2);
    assert(d->m->labels != NULL);
    assert(d->m->labels[label1index] != NULL);
    assert(d->m->labels[label2index] != NULL);
    assert(d->m->labels[label1index][0] == newlabel1value);
    assert(d->m->labels[label1index][1] == newerlabel1values[0]);
    assert(d->m->labels[label1index][2] == newerlabel1values[1]);
    assert(d->m->labels[label1index][3] == newerlabel1values[2]);
    assert(d->m->labels[label1index][4] == newerlabel1values[3]);
    assert(d->m->labels[label1index][5] == newlabel1value);

    // Test name-setters here too, they use the same procedure as set_labels_to_values
    // (Tests when setting whole body not just group)
    const char* newNames[] = {"name0", "aname", "two", "THR33", "7", "name$$"};
    set_names_to_values(d,0,0,6,newNames);
    assert(strncmp(d->m->names[0],newNames[0],strlen(newNames[0])) == 0);
    assert(strncmp(d->m->names[1],newNames[1],strlen(newNames[1])) == 0);
    assert(strncmp(d->m->names[2],newNames[2],strlen(newNames[2])) == 0);
    assert(strncmp(d->m->names[3],newNames[3],strlen(newNames[3])) == 0);
    assert(strncmp(d->m->names[4],newNames[4],strlen(newNames[4])) == 0);
    assert(strncmp(d->m->names[5],newNames[5],strlen(newNames[5])) == 0);
    const char* betterNames[] = {"G01", "G02", "G03", "G04", "G05", "G06"};
    set_names_to_values(d,g0,0,6,betterNames);
    assert(strncmp(d->m->names[0],betterNames[0],strlen(betterNames[0])) == 0);
    assert(strncmp(d->m->names[1],betterNames[1],strlen(betterNames[1])) == 0);
    assert(strncmp(d->m->names[2],betterNames[2],strlen(betterNames[2])) == 0);
    assert(strncmp(d->m->names[3],betterNames[3],strlen(betterNames[3])) == 0);
    assert(strncmp(d->m->names[4],betterNames[4],strlen(betterNames[4])) == 0);
    assert(strncmp(d->m->names[5],betterNames[5],strlen(betterNames[5])) == 0);

    // Can increment a label
    const int increment = -4+1;
    increment_labels(d, g0, label1, 1);
    increment_labels(d, g0, label1, -4);
    assert(d->n_labels == 2);
    assert(d->m->n_labels == 2);
    assert(d->m->labels != NULL);
    assert(d->m->labels[label1index] != NULL);
    assert(d->m->labels[label2index] != NULL);
    assert(d->m->labels[label1index][0] == newlabel1value + increment);
    assert(d->m->labels[label1index][1] == newerlabel1values[0] + increment);
    assert(d->m->labels[label1index][2] == newerlabel1values[1] + increment);
    assert(d->m->labels[label1index][3] == newerlabel1values[2] + increment);
    assert(d->m->labels[label1index][4] == newerlabel1values[3] + increment);
    assert(d->m->labels[label1index][5] == newlabel1value + increment);

    // Need a second group just for testing purposes
    GenOptions g = {.family_size = 1,
                    .will_track_pedigree = FALSE,
                    .will_name_offspring = FALSE,
                    .offspring_name_prefix = NULL,
                    .will_allocate_ids = FALSE,
                    .filename_prefix = NULL,
                    .will_save_alleles_to_file = FALSE,
                    .will_save_bvs_to_file = FALSE,
                    .will_save_pedigree_to_file = FALSE,
                    .will_save_to_simdata = TRUE};
    int f1 = cross_random_individuals(d,g0,4,0,g);
    // Check that the labels are set to their defaults
    assert(d->m->labels[label1index][6+0] == label1value);
    assert(d->m->labels[label1index][6+1] == label1value);
    assert(d->m->labels[label1index][6+2] == label1value);
    assert(d->m->labels[label1index][6+3] == label1value);
    assert(d->m->labels[label2index][6+0] == newlabel2default);
    assert(d->m->labels[label2index][6+1] == newlabel2default);
    assert(d->m->labels[label2index][6+2] == newlabel2default);
    assert(d->m->labels[label2index][6+3] == newlabel2default);

    set_labels_to_const(d,f1,label1,newlabel1value + increment);
    set_labels_to_values(d,f1,0,label2,4,newerlabel1values);

    // test can split by label (across groups)
    int groupB = split_by_label_value(d, 0, label1, newlabel1value + increment);
    assert(g0 != groupB && f1 != groupB);
    assert(get_group_size(d, groupB) == 2+4);
    int Bindexes[2+4];
    get_group_indexes(d, groupB, 2+4, Bindexes);
    assert(Bindexes[0] == 0);
    assert(Bindexes[1] == 5);
    assert(Bindexes[2] == 6 && Bindexes[3] == 7 && Bindexes[4] == 8 && Bindexes[5] == 9);

    int outtakes[4] = {6,7,8,9};
    int f1outtakes = split_from_group(d,4,outtakes);
    int toCombine[2] = {f1, f1outtakes};
    f1 = combine_groups(d, 2, toCombine);
    assert(get_group_size(d, f1) == 4);
    toCombine[0] = g0;
    toCombine[1] = groupB;
    g0 = combine_groups(d, 2, toCombine);
    assert(get_group_size(d, g0) == 6);

    // test can split by label range (across groups)
    groupB = split_by_label_range(d, 0, label2, 1, 15);
    assert(g0 != groupB && f1 != groupB);
    assert(get_group_size(d, groupB) == 6);
    get_group_indexes(d, groupB, 6, Bindexes);
    assert(Bindexes[0] == 0);
    assert(Bindexes[1] == 1);
    assert(Bindexes[2] == 2);
    assert(Bindexes[3] == 3);
    assert(Bindexes[4] == 4);
    assert(Bindexes[5] == 7);

    outtakes[0] = 7;
    f1outtakes = split_from_group(d,1,outtakes);
    toCombine[0] = f1;
    toCombine[1] = f1outtakes;
    f1 = combine_groups(d, 2, toCombine);
    assert(get_group_size(d, f1) == 4);
    toCombine[0] = g0;
    toCombine[1] = groupB;
    g0 = combine_groups(d, 2, toCombine);
    assert(get_group_size(d, g0) == 6);

    // and can split from group (within group)
    groupB = split_by_label_value(d, g0, label1, newlabel1value + increment);
    assert(g0 != groupB && f1 != groupB);
    assert(get_group_size(d, groupB) == 2);
    get_group_indexes(d, groupB, 2, Bindexes);
    assert(Bindexes[0] == 0);
    assert(Bindexes[1] == 5);

    toCombine[0] = g0;
    toCombine[1] = groupB;
    g0 = combine_groups(d, 2, toCombine);
    assert(get_group_size(d, g0) == 6);
    assert(get_group_size(d, f1) == 4);


    // and can split from group by label range (within group)
    groupB = split_by_label_range(d, g0, label2, 1, 15);
    assert(g0 != groupB && f1 != groupB);
    assert(get_group_size(d, groupB) == 5);
    get_group_indexes(d, groupB, 5, Bindexes);
    assert(Bindexes[0] == 0);
    assert(Bindexes[1] == 1);
    assert(Bindexes[2] == 2);
    assert(Bindexes[3] == 3);
    assert(Bindexes[4] == 4);

    toCombine[0] = g0;
    toCombine[1] = groupB;
    g0 = combine_groups(d, 2, toCombine);
    assert(get_group_size(d, g0) == 6);
    assert(get_group_size(d, f1) == 4);
    delete_group(d,f1);

    // Check label deletion works
    delete_label(d,label1);
    assert(d->n_labels == 1);
    assert(d->m->n_labels == 1);
    assert(d->label_ids[0] == label2);
    assert(d->label_defaults[0] == newlabel2default);
    assert(d->m->labels[0][0] == 11);

    printf("...Label manipulation runs correctly\n");

    return g0;
}

int test_random_splits(SimData *d, int g0) {
    // g0: contents of "a-test.txt", contains 6 genotypes

    // Can split into 2; repeat 10x for confidence
    for (int i = 0; i < 10; ++i) {
        int grpB = split_evenly_into_two(d, g0);
        assert(get_group_size(d,g0) == 3);
        assert(get_group_size(d,grpB) == 3);

        int toMerge[2] = {g0, grpB};
        g0 = combine_groups(d,2,toMerge);
        assert(get_group_size(d,g0) == 6);
    }

    // Can split into 3; repeat 10x for confidence
    for (int i = 0; i < 10; ++i) {
        int grpB[3];
        split_evenly_into_n(d,g0,3,grpB);
        assert(get_group_size(d,grpB[0]) == 2);
        assert(get_group_size(d,grpB[1]) == 2);
        assert(get_group_size(d,grpB[2]) == 2);

        g0 = combine_groups(d,3,grpB);
        assert(get_group_size(d,g0) == 6);
    }

    // Can split into bins; repeat 10x for confidence
    for (int i = 0; i < 10; ++i) {
        int grpB[3];
        int grpBsizes[3] = {1,3,2};
        split_by_specific_counts_into_n(d,g0,3,grpBsizes,grpB);
        assert(get_group_size(d,grpB[0]) == 1);
        assert(get_group_size(d,grpB[1]) == 3);
        assert(get_group_size(d,grpB[2]) == 2);

        g0 = combine_groups(d,3,grpB);
        assert(get_group_size(d,g0) == 6);
    }

    // In future replace these with something that uses the random generator seed
    // Can split by coinflip
    int grpBsizes = 0;
    for (int i = 0; i < 10; ++i) {
        int grpB = split_randomly_into_two(d, g0);
        int grpBsize = get_group_size(d,grpB);
        assert(get_group_size(d,g0) + grpBsize == 6);
        grpBsizes += grpBsize;

        int toMerge[2] = {g0, grpB};
        g0 = combine_groups(d,2,toMerge);
        assert(get_group_size(d,g0) == 6);
    }
    assert(grpBsizes > 20 && grpBsizes < 40); // about 1.3% chance of failure by random chance

    // Can split by thirds
    grpBsizes = 0;
    int grpBsizes2 = 0;
    for (int i = 0; i < 10; ++i) {
        int grpB[3];
        split_randomly_into_n(d,g0,3,grpB);
        int grpBsize = get_group_size(d,grpB[0]);
        int grpBsize2 = get_group_size(d,grpB[1]);
        assert(get_group_size(d,grpB[2]) + grpBsize + grpBsize2 == 6);
        grpBsizes += grpBsize;
        grpBsizes2 += grpBsize2;

        g0 = combine_groups(d,3,grpB);
        assert(get_group_size(d,g0) == 6);
    }
    assert(grpBsizes > 9 && grpBsizes < 34); // about 1.3% chance of failure by random chance
    assert(grpBsizes2 > 9 && grpBsizes2 < 34);

    // Can split by any probability
    grpBsizes = 0;
    for (int i = 0; i < 10; ++i) {
        int grpB[2];
        double grpBprobs = 3.0/7.0;
        split_by_probabilities_into_n(d,g0,2,&grpBprobs,grpB);
        int grpBsize = get_group_size(d,grpB[0]);
        assert(get_group_size(d,grpB[1]) + grpBsize == 6);
        grpBsizes += grpBsize;

        g0 = combine_groups(d,2,grpB);
        assert(get_group_size(d,g0) == 6);
    }
    assert(grpBsizes > 16 && grpBsizes < 39); // about 0.75% chance of failure by random chance

    grpBsizes = 0;
    for (int i = 0; i < 10; ++i) {
        int grpB[2];
        double grpBprobs[2] = {1./12.,11./12.};
        split_by_probabilities_into_n(d,g0,2,grpBprobs,grpB);
        int grpBsize = get_group_size(d,grpB[0]);
        assert(get_group_size(d,grpB[1]) + grpBsize == 6);
        grpBsizes += grpBsize;

        g0 = combine_groups(d,2,grpB);
        assert(get_group_size(d,g0) == 6);
    }
    assert(grpBsizes < 11); // about 1% chance of failure by random chance

    printf("...Random group splitting runs correctly\n");

    return g0;
}

int test_specific_splits(SimData *d, int g0) {
    // g0: contents of "a-test.txt", contains 6 genotypes
    // and there are few enough genotypes they all fit into the first AlleleMatrix

    // Can split into individuals and recombine
    const int g0size = get_group_size(d, g0);
    assert(g0size == 6); // pre-test test validity requirement
    int g0indivs[8];
    for (int i = 0; i < 8; ++i) {
        g0indivs[i] = 0;
    }
    split_into_individuals(d, g0, g0indivs);
    assert(g0indivs[0] > 0);
    assert(g0indivs[1] > 0);
    assert(g0indivs[1] != g0indivs[0]);
    assert(g0indivs[2] > 0);
    assert(g0indivs[2] != g0indivs[0] && g0indivs[2] != g0indivs[1]);
    assert(g0indivs[3] > 0);
    assert(g0indivs[3] != g0indivs[0] && g0indivs[3] != g0indivs[1] && g0indivs[3] != g0indivs[2]);
    assert(g0indivs[4] > 0);
    assert(g0indivs[4] != g0indivs[0] && g0indivs[4] != g0indivs[1]
            && g0indivs[4] != g0indivs[2] && g0indivs[4] != g0indivs[3]);
    assert(g0indivs[5] > 0);
    assert(g0indivs[5] != g0indivs[0] && g0indivs[5] != g0indivs[1] && g0indivs[5] != g0indivs[2]
            && g0indivs[5] != g0indivs[3] && g0indivs[5] != g0indivs[4]);
    assert(g0indivs[6] == 0);
    assert(g0indivs[7] == 0);

    g0 = combine_groups(d,6,g0indivs);
    assert(get_group_size(d,g0) == 6);

    // Can split into families
    GenOptions g = {.family_size = 5,
                    .will_track_pedigree = TRUE,
                    .will_name_offspring = FALSE,
                    .offspring_name_prefix = NULL,
                    .will_allocate_ids = FALSE,
                    .filename_prefix = NULL,
                    .will_save_alleles_to_file = FALSE,
                    .will_save_bvs_to_file = FALSE,
                    .will_save_pedigree_to_file = FALSE,
                    .will_save_to_simdata = TRUE};
    int f1 = cross_random_individuals(d, g0, 3, 1, g);
    int f1size = 5*3;
    int families[f1size];
    for (int i = 0; i < f1size; ++i) {
        families[i] = 0;
    }
    split_into_families(d,f1,families);
    assert(families[0] > 0);
    assert(families[1] > 0);
    assert(families[2] > 0);
    for (int i = 3; i < f1size; ++i) {
        assert(families[i] == 0);
    }
    assert(families[0] != families[1] && families[1] != families[2] && families[0] != families[2]);
    assert(get_group_size(d,families[0]) == 5);
    assert(get_group_size(d,families[1]) == 5);
    assert(get_group_size(d,families[2]) == 5);
    // Check all groups contain families
    int f1family1[5];
    assert(get_group_indexes(d,families[0],5,f1family1) == 5);
    assert(d->m->pedigrees[0][f1family1[0]] == d->m->pedigrees[0][f1family1[1]] &&
            d->m->pedigrees[0][f1family1[0]] == d->m->pedigrees[0][f1family1[2]] &&
            d->m->pedigrees[0][f1family1[0]] == d->m->pedigrees[0][f1family1[3]] &&
            d->m->pedigrees[0][f1family1[0]] == d->m->pedigrees[0][f1family1[4]]);
    assert(d->m->pedigrees[1][f1family1[0]] == d->m->pedigrees[1][f1family1[1]] &&
            d->m->pedigrees[1][f1family1[0]] == d->m->pedigrees[1][f1family1[2]] &&
            d->m->pedigrees[1][f1family1[0]] == d->m->pedigrees[1][f1family1[3]] &&
            d->m->pedigrees[1][f1family1[0]] == d->m->pedigrees[1][f1family1[4]]);

    int f1family2[5];
    assert(get_group_indexes(d,families[1],5,f1family2) == 5);
    assert(d->m->pedigrees[0][f1family1[0]] == d->m->pedigrees[0][f1family1[1]] &&
            d->m->pedigrees[0][f1family1[0]] == d->m->pedigrees[0][f1family1[2]] &&
            d->m->pedigrees[0][f1family1[0]] == d->m->pedigrees[0][f1family1[3]] &&
            d->m->pedigrees[0][f1family1[0]] == d->m->pedigrees[0][f1family1[4]]);
    assert(d->m->pedigrees[1][f1family1[0]] == d->m->pedigrees[1][f1family1[1]] &&
            d->m->pedigrees[1][f1family1[0]] == d->m->pedigrees[1][f1family1[2]] &&
            d->m->pedigrees[1][f1family1[0]] == d->m->pedigrees[1][f1family1[3]] &&
            d->m->pedigrees[1][f1family1[0]] == d->m->pedigrees[1][f1family1[4]]);

    int f1family3[5];
    assert(get_group_indexes(d,families[2],5,f1family3) == 5);
    assert(d->m->pedigrees[0][f1family1[0]] == d->m->pedigrees[0][f1family1[1]] &&
            d->m->pedigrees[0][f1family1[0]] == d->m->pedigrees[0][f1family1[2]] &&
            d->m->pedigrees[0][f1family1[0]] == d->m->pedigrees[0][f1family1[3]] &&
            d->m->pedigrees[0][f1family1[0]] == d->m->pedigrees[0][f1family1[4]]);
    assert(d->m->pedigrees[1][f1family1[0]] == d->m->pedigrees[1][f1family1[1]] &&
            d->m->pedigrees[1][f1family1[0]] == d->m->pedigrees[1][f1family1[2]] &&
            d->m->pedigrees[1][f1family1[0]] == d->m->pedigrees[1][f1family1[3]] &&
            d->m->pedigrees[1][f1family1[0]] == d->m->pedigrees[1][f1family1[4]]);

    // Check separate groups don't represent the same family
    assert(d->m->pedigrees[0][f1family1[0]] != d->m->pedigrees[0][f1family2[0]] &&
            d->m->pedigrees[1][f1family1[0]] != d->m->pedigrees[1][f1family2[0]]);
    assert(d->m->pedigrees[0][f1family1[0]] != d->m->pedigrees[0][f1family3[0]] &&
            d->m->pedigrees[1][f1family1[0]] != d->m->pedigrees[1][f1family3[0]]);
    assert(d->m->pedigrees[0][f1family2[0]] != d->m->pedigrees[0][f1family3[0]] &&
            d->m->pedigrees[1][f1family2[0]] != d->m->pedigrees[1][f1family3[0]]);

    delete_group(d, families[0]);
    delete_group(d, families[1]);
    delete_group(d, families[2]);

    // Can split into halfsib families
    int combinations[2][4];
    combinations[0][0] = 0; combinations[1][0] = 1;
    combinations[0][1] = 0; combinations[1][1] = 2;
    combinations[0][2] = 3; combinations[1][2] = 1;
    combinations[0][3] = 4; combinations[1][3] = 5;
    int fhs = cross_these_combinations(d,4,combinations[0], combinations[1],g);

    int halfsibfamilies[5];
    for (int i = 0; i < f1size; ++i) {
        halfsibfamilies[i] = 0;
    }
    split_into_halfsib_families(d,fhs,1,halfsibfamilies);
    assert(halfsibfamilies[0] > 0);
    assert(halfsibfamilies[1] > 0);
    assert(halfsibfamilies[2] > 0);
    assert(halfsibfamilies[3] == 0);
    assert(halfsibfamilies[4] == 0);

    // Check all halfsib families share the same parent 1
    int fhsfamily1size = get_group_size(d,halfsibfamilies[0]); // assumes the largest hsfamily will be the first one seen
    assert(fhsfamily1size == 10);
    int fhsfamily1[10];
    assert(get_group_indexes(d,halfsibfamilies[0],10,fhsfamily1) == fhsfamily1size);
    assert(d->m->pedigrees[0][fhsfamily1[0]] == d->m->pedigrees[0][fhsfamily1[1]] &&
            d->m->pedigrees[0][fhsfamily1[0]] == d->m->pedigrees[0][fhsfamily1[2]] &&
            d->m->pedigrees[0][fhsfamily1[0]] == d->m->pedigrees[0][fhsfamily1[3]] &&
            d->m->pedigrees[0][fhsfamily1[0]] == d->m->pedigrees[0][fhsfamily1[4]] &&
            d->m->pedigrees[0][fhsfamily1[0]] == d->m->pedigrees[0][fhsfamily1[5]] &&
            d->m->pedigrees[0][fhsfamily1[0]] == d->m->pedigrees[0][fhsfamily1[6]] &&
            d->m->pedigrees[0][fhsfamily1[0]] == d->m->pedigrees[0][fhsfamily1[7]] &&
            d->m->pedigrees[0][fhsfamily1[0]] == d->m->pedigrees[0][fhsfamily1[8]] &&
            d->m->pedigrees[0][fhsfamily1[0]] == d->m->pedigrees[0][fhsfamily1[9]]);

    int fhsfamily2size = get_group_size(d,halfsibfamilies[1]);
    assert(fhsfamily2size == 5);
    int fhsfamily2[5];
    assert(get_group_indexes(d,halfsibfamilies[1],5,fhsfamily2) == fhsfamily2size);
    assert(d->m->pedigrees[0][fhsfamily2[0]] == d->m->pedigrees[0][fhsfamily2[1]] &&
            d->m->pedigrees[0][fhsfamily2[0]] == d->m->pedigrees[0][fhsfamily2[2]] &&
            d->m->pedigrees[0][fhsfamily2[0]] == d->m->pedigrees[0][fhsfamily2[3]] &&
            d->m->pedigrees[0][fhsfamily2[0]] == d->m->pedigrees[0][fhsfamily2[4]]);

    int fhsfamily3size = get_group_size(d,halfsibfamilies[2]);
    assert(fhsfamily3size == 5);
    int fhsfamily3[fhsfamily3size];
    assert(get_group_indexes(d,halfsibfamilies[2],5,fhsfamily3) == fhsfamily3size);
    assert(d->m->pedigrees[0][fhsfamily3[0]] == d->m->pedigrees[0][fhsfamily3[1]] &&
            d->m->pedigrees[0][fhsfamily3[0]] == d->m->pedigrees[0][fhsfamily3[2]] &&
            d->m->pedigrees[0][fhsfamily3[0]] == d->m->pedigrees[0][fhsfamily3[3]] &&
            d->m->pedigrees[0][fhsfamily3[0]] == d->m->pedigrees[0][fhsfamily3[4]]);

    // and that none of the across-groups have the same parent 1
    assert(d->m->pedigrees[0][fhsfamily1[0]] != d->m->pedigrees[0][fhsfamily2[0]]);
    assert(d->m->pedigrees[0][fhsfamily1[0]] != d->m->pedigrees[0][fhsfamily3[0]]);
    assert(d->m->pedigrees[0][fhsfamily2[0]] != d->m->pedigrees[0][fhsfamily3[0]]);

    // Then check halfsibs in the other direction
    fhs = combine_groups(d,3,halfsibfamilies);
    for (int i = 0; i < f1size; ++i) {
        halfsibfamilies[i] = 0;
    }
    split_into_halfsib_families(d,fhs,2,halfsibfamilies);
    assert(halfsibfamilies[0] > 0);
    assert(halfsibfamilies[1] > 0);
    assert(halfsibfamilies[2] > 0);
    assert(halfsibfamilies[3] == 0);
    assert(halfsibfamilies[4] == 0);

    // Check all halfsib families share the same parent 1
    fhsfamily1size = get_group_size(d,halfsibfamilies[0]); // assumes the largest hsfamily will be the first one seen
    assert(fhsfamily1size == 10);
    assert(get_group_indexes(d,halfsibfamilies[0],10,fhsfamily1) == fhsfamily1size);
    assert(d->m->pedigrees[1][fhsfamily1[0]] == d->m->pedigrees[1][fhsfamily1[1]] &&
            d->m->pedigrees[1][fhsfamily1[0]] == d->m->pedigrees[1][fhsfamily1[2]] &&
            d->m->pedigrees[1][fhsfamily1[0]] == d->m->pedigrees[1][fhsfamily1[3]] &&
            d->m->pedigrees[1][fhsfamily1[0]] == d->m->pedigrees[1][fhsfamily1[4]] &&
            d->m->pedigrees[1][fhsfamily1[0]] == d->m->pedigrees[1][fhsfamily1[5]] &&
            d->m->pedigrees[1][fhsfamily1[0]] == d->m->pedigrees[1][fhsfamily1[6]] &&
            d->m->pedigrees[1][fhsfamily1[0]] == d->m->pedigrees[1][fhsfamily1[7]] &&
            d->m->pedigrees[1][fhsfamily1[0]] == d->m->pedigrees[1][fhsfamily1[8]] &&
            d->m->pedigrees[1][fhsfamily1[0]] == d->m->pedigrees[1][fhsfamily1[9]]);

    fhsfamily2size = get_group_size(d,halfsibfamilies[1]);
    assert(fhsfamily2size == 5);
    assert(get_group_indexes(d,halfsibfamilies[1],5,fhsfamily2) == fhsfamily2size);
    assert(d->m->pedigrees[1][fhsfamily2[0]] == d->m->pedigrees[1][fhsfamily2[1]] &&
            d->m->pedigrees[1][fhsfamily2[0]] == d->m->pedigrees[1][fhsfamily2[2]] &&
            d->m->pedigrees[1][fhsfamily2[0]] == d->m->pedigrees[1][fhsfamily2[3]] &&
            d->m->pedigrees[1][fhsfamily2[0]] == d->m->pedigrees[1][fhsfamily2[4]]);

    fhsfamily3size = get_group_size(d,halfsibfamilies[2]);
    assert(fhsfamily3size == 5);
    assert(get_group_indexes(d,halfsibfamilies[2],5,fhsfamily3) == fhsfamily3size);
    assert(d->m->pedigrees[1][fhsfamily3[0]] == d->m->pedigrees[1][fhsfamily3[1]] &&
            d->m->pedigrees[1][fhsfamily3[0]] == d->m->pedigrees[1][fhsfamily3[2]] &&
            d->m->pedigrees[1][fhsfamily3[0]] == d->m->pedigrees[1][fhsfamily3[3]] &&
            d->m->pedigrees[1][fhsfamily3[0]] == d->m->pedigrees[1][fhsfamily3[4]]);

    // and that none of the across-groups have the same parent 1
    assert(d->m->pedigrees[1][fhsfamily1[0]] != d->m->pedigrees[1][fhsfamily2[0]]);
    assert(d->m->pedigrees[1][fhsfamily1[0]] != d->m->pedigrees[1][fhsfamily3[0]]);
    assert(d->m->pedigrees[1][fhsfamily2[0]] != d->m->pedigrees[1][fhsfamily3[0]]);

    fhs = combine_groups(d,3,halfsibfamilies);

    // Can create a group from indexes
    int splitters[4] = {0, 2, 5, 6};
    int g0b = split_from_group(d, 4, splitters);
    assert(get_group_size(d,g0) == 3 && get_group_size(d,g0b) == 4);
    int g0bi[4];
    assert(get_group_indexes(d,g0b,4,g0bi) == 4);
    assert(g0bi[0] == 0 && g0bi[1] == 2 && g0bi[2] == 5 && g0bi[3] == 6);

    int n6 = split_from_group(d, 1, splitters+3);
    int combiners[2] = {fhs, n6};
    fhs = combine_groups(d, 2, combiners);
    delete_group(d, fhs);

    int g0members[6] = {0,1,2,3,4,5};
    g0 = split_from_group(d,6,g0members);

    printf("...Specified group splitting runs correctly\n");

    return g0;
}


int test_grouping(SimData *d, int g0) {
    g0 = test_labels(d, g0);
    g0 = test_random_splits(d, g0);
    g0 = test_specific_splits(d, g0);

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

int test_optimal_calculators(SimData *d, int g0) {
    char* ig = calculate_optimal_alleles(d);
    assert(ig[0] == 'T');
    assert(ig[1] == 'A');
    assert(ig[2] == 'A');
    free(ig);

    double optimal = calculate_optimum_bv(d);
    assert(fabs(optimal - 1.8) < TOL);

    double unoptimal = calculate_minimum_bv(d);
    assert(fabs(unoptimal + 2.8) < TOL);

    char* founderhaplo = calculate_optimal_available_alleles(d, g0);
    assert(founderhaplo[0] == 'T');
    assert(founderhaplo[1] == 'A');
    assert(founderhaplo[2] == 'A');
    free(founderhaplo);

    double founderoptimal = calculate_optimal_available_bv(d, g0);
    assert(fabs(founderoptimal - 1.8) < TOL);

    int factorout[2] = {4,5};
    int g0partial = split_from_group(d,2, factorout);
    char* founderhaplo2 = calculate_optimal_available_alleles(d, g0partial);
    assert(founderhaplo2[0] == 'T');
    assert(founderhaplo2[1] == 'A');
    assert(founderhaplo2[2] == 'T');
    free(founderhaplo2);

    double founderoptimal2 = calculate_optimal_available_bv(d, g0partial);
    assert(fabs(founderoptimal2 - 1.4) < TOL);

    int recombine[2] = {g0,g0partial};
    int newg0 = combine_groups(d,2,recombine);
    assert(g0 == newg0); // for validity of further tests

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
    GenOptions g = {.will_name_offspring=TRUE,
                    .offspring_name_prefix="F1",
                    .family_size=1,.will_track_pedigree=TRUE,
                    .will_allocate_ids=TRUE,
                    .filename_prefix="atestF1",
                    .will_save_pedigree_to_file=TRUE,
                    .will_save_bvs_to_file=TRUE,
                    .will_save_alleles_to_file=TRUE,
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
    //printf("%s -> %s\n", d->m->alleles[0], d->m->alleles[oldsize + 6]);
    assert(strncmp(d->m->alleles[oldsize + 6],d->m->alleles[0],sizeof(char)*6) == 0);
    assert(strncmp(d->m->alleles[oldsize + 11],"AAAATT",sizeof(char)*6) == 0 ||
            strncmp(d->m->alleles[oldsize + 11],"TTAATT",sizeof(char)*6) == 0);
    printf("...doubled haploid function works correctly\n");

    delete_group(d, g1dhap);

    // test cloning
    opt.will_allocate_ids = FALSE;
    int g1clones = make_clones(d, g1, TRUE, opt);
    assert(g1clones != g1);
    assert(d->m->n_genotypes == oldsize + 2*6);
    assert(d->m->n_markers == 3);
    assert(fabsf(calculate_heterozygosity(d,  g1clones) - h1) < TOL);
    assert(d->m->groups[oldsize + 6] == g1clones && d->m->groups[oldsize + 11] == g1clones && d->m->groups[oldsize + 12] != g1clones);
    assert(d->m->pedigrees[0][oldsize + 6] == d->m->ids[0] && d->m->pedigrees[1][oldsize + 6] == d->m->ids[0]);
    assert(d->m->pedigrees[0][oldsize + 10] == d->m->ids[4] && d->m->pedigrees[1][oldsize + 10] == d->m->ids[4]);
    assert(strncmp(d->m->alleles[oldsize + 6],d->m->alleles[0],sizeof(char)*6) == 0);
    assert(strncmp(d->m->alleles[oldsize + 11],"ATAATT",sizeof(char)*6) == 0);
    assert(strcmp(d->m->names[1],d->m->names[oldsize + 7]) == 0);
    assert(strcmp(d->m->names[3],d->m->names[oldsize + 9]) == 0);
    printf("...cloning function works correctly\n");

    opt.will_allocate_ids = TRUE;
    int gap = cross_random_individuals(d , g1, 5, 0, opt);
    int newparents[2];
    newparents[0] = g1clones;
    newparents[1] = cross_random_individuals(d , g1, 5, 0, opt);
    int cloneparents2 = combine_groups(d, 2, newparents);
    int g2clones = make_clones(d, cloneparents2, FALSE, opt);

    delete_group(d, cloneparents2);
    delete_group(d, gap);
    delete_group(d, g2clones);

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
    int g2ixs[4];
    assert(get_group_indexes(d, g2, -1, g2ixs) == 4);

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
}

int test_deletors(SimData *d, int g0) {
    int groups1[50];
    int ngroups1 = get_existing_groups(d, 50, groups1);

    int groups1b[50];
    int groupcounts1b[50];
    int ngroups1b = get_existing_group_counts(d, 40, groups1b, groupcounts1b);
    assert(ngroups1b == ngroups1);
    for (int i = 0; i < ngroups1; ++i) {
        assert(groups1[i] == groups1b[i]);
    }

	delete_group(d, g0);
    int groups2[50];
    int ngroups2 = get_existing_groups(d, 50, groups2);

	assert(ngroups1 - ngroups2 == 1);
	for (int i = 0; i < ngroups1; ++i) {
		if (i == ngroups1 - 1 || groups1[i] != groups2[i]) {
			assert(groups1[i] == g0);
            assert(groups1b[i] == g0);
            assert(groupcounts1b[i] == 6);
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

int test_iterators(SimData* d, int gp) {
    GenOptions g = {
        .will_name_offspring = FALSE,
        .offspring_name_prefix = NULL,
        .family_size = 1,
        .will_track_pedigree = TRUE,
        .will_allocate_ids = FALSE,
        .filename_prefix = NULL,
        .will_save_pedigree_to_file = FALSE,
        .will_save_bvs_to_file = FALSE,
        .will_save_alleles_to_file = FALSE,
        .will_save_to_simdata = TRUE
    };
    int combos[2][3];
    combos[0][0] = 0; combos[1][0] = 0;
    combos[0][1] = 1; combos[1][1] = 2;
    combos[0][2] = 1; combos[1][2] = 5;
    int f1 = cross_these_combinations(d, 3, combos[0], combos[1], g);
    const int label1 = create_new_label(d, 0);
    increment_labels(d,f1,label1,1);

    // Bidirectional global iterator; starting forwards
    BidirectionalIterator it1 = create_bidirectional_iter(d, 0);
    int i = 0;
    for (; i < 6; ++i) {
        GenoLocation gl = next_forwards(&it1);
        assert(isValidLocation(gl));
        assert(get_id(gl) == i+1);
        assert(get_label_value(gl,get_index_of_label(d,label1)) == 0);
        assert(get_group(gl) == gp);
    }
    for (; i < 6+3; ++i) {
        GenoLocation gl = next_forwards(&it1);
        assert(isValidLocation(gl));
        assert(get_id(gl) == 0);
        assert(get_label_value(gl,get_index_of_label(d,label1)) == 1);
        assert(get_group(gl) == f1);
    }
    assert(!isValidLocation(next_forwards(&it1)));
    for (int j = 0; j < 3-1; ++j) {
        assert(get_id(next_backwards(&it1)) == 0);
    }
    assert(get_id(next_backwards(&it1)) == 6);

    delete_bidirectional_iter(&it1);

    // Bidirectional group iterator; starting backwards
    BidirectionalIterator it2 = create_bidirectional_iter(d, gp);
    assert(get_id(next_backwards(&it2)) == 6);
    assert(!isValidLocation(next_forwards(&it2)));
    assert(get_id(next_backwards(&it2)) == 5);
    assert(get_id(next_backwards(&it2)) == 4);
    assert(get_id(next_backwards(&it2)) == 3);
    assert(get_id(next_backwards(&it2)) == 2);
    assert(get_id(next_forwards(&it2)) == 3);
    assert(get_id(next_backwards(&it2)) == 2);
    assert(get_id(next_backwards(&it2)) == 1);
    assert(!isValidLocation(next_backwards(&it2)));

    // RandomAccess global iterator
    RandomAccessIterator it3 = create_randomaccess_iter(d, 0);
    assert(get_group(next_get_nth(&it3, 6)) == f1);
    assert(get_group(next_get_nth(&it3, 2)) == gp);
    assert(get_id(next_get_nth(&it3, 2)) == 3);
    assert(get_id(next_get_nth(&it3, 5)) == 6);
    assert(!isValidLocation(next_get_nth(&it3,-1)));
    assert(!isValidLocation(next_get_nth(&it3,10)));

    delete_randomaccess_iter(&it3);

    // RandomAccess group iterator
    char* newnames[4] = {"f1a", "f1b", "f1c", "f1d"};
    set_names_to_values(d,f1,0,4,newnames);
    RandomAccessIterator it4 = create_randomaccess_iter(d, f1);
    assert(get_group(next_get_nth(&it4, 0)) == f1);
    assert(strncmp(get_name(next_get_nth(&it4, 0)),"f1a",5)==0);
    assert(get_first_parent(next_get_nth(&it4, 0)) == 1);
    assert(get_first_parent(next_get_nth(&it4, 2)) == 2);
    assert(get_second_parent(next_get_nth(&it4, 1)) == 3);

    delete_randomaccess_iter(&it4);

    // Empty group iterator
    RandomAccessIterator it5 = create_randomaccess_iter(d,4);
    assert(!isValidLocation(next_get_nth(&it5,0)));

    delete_group(d, f1);
    delete_label(d, get_index_of_label(d, label1));

    printf("...iterators work correctly.\n");

    return 0;
}

int test_getters(SimData* d, int gp) {
    assert(get_group_size(d, gp) == 6);
    char* alleles[6];
    assert(get_group_genes(d, gp, 6, alleles) == 6);
    assert(strncmp(alleles[0],"TTAATT", 6) == 0); // G01
    assert(strncmp(alleles[1],"TTAATT", 6) == 0); // G02
    assert(strncmp(alleles[2],"TTAATA", 6) == 0); // G03
    assert(strncmp(alleles[3],"TAAATA", 6) == 0); // G04
    assert(strncmp(alleles[4],"TTTTTT", 6) == 0); // G05
    assert(strncmp(alleles[5],"ATAATT", 6) == 0); // G06

    char* names[6];
    assert(get_group_names(d, gp, -1, names) == 6);
    assert(strcmp(names[0], "G01") == 0);
    assert(strcmp(names[1], "G02") == 0);
    assert(strcmp(names[2], "G03") == 0);
    assert(strcmp(names[3], "G04") == 0);
    assert(strcmp(names[4], "G05") == 0);
    assert(strcmp(names[5], "G06") == 0);

    unsigned int ids[10];
    assert(get_group_ids(d, gp, -1, ids) == 6);
    for (int i = 0; i < 6; ++i) {
        assert(ids[i] == i+1);
    }

    int indexes[6];
    assert(get_group_indexes(d, gp, 6, indexes) == 6);
    for (int i = 0; i < 6; ++i) {
        assert(indexes[i] == i);
    }

    double bvs[6];
    assert(get_group_bvs(d, gp, 6, bvs) == 6);
    assert(fabs(bvs[0] - 1.4) < TOL);
    assert(fabs(bvs[1] - 1.4) < TOL);
    assert(fabs(bvs[2] - 1.6) < TOL);
    assert(fabs(bvs[3] - (-0.1)) < TOL);
    assert(fabs(bvs[4] - 0.6) < TOL);
    assert(fabs(bvs[5] - (-0.3)) < TOL);

    int combos[2][3];
    combos[0][0] = 0; combos[1][0] = 0;
    combos[0][1] = 1; combos[1][1] = 2;
    combos[0][2] = 1; combos[1][2] = 5;
    GenOptions g = {
        .will_name_offspring = FALSE,
        .offspring_name_prefix = NULL,
        .family_size = 1,
        .will_track_pedigree = TRUE,
        .will_allocate_ids = FALSE,
        .filename_prefix = NULL,
        .will_save_pedigree_to_file = FALSE,
        .will_save_bvs_to_file = FALSE,
        .will_save_alleles_to_file = FALSE,
        .will_save_to_simdata = TRUE
    };
    int f1 = cross_these_combinations(d,3,combos[0],combos[1],g);

    unsigned int p1s[3]; unsigned int p2s[6];
    assert(get_group_parent_ids(d, f1, UNINITIALISED, 1, p1s) == 3);
    assert(get_group_parent_ids(d, f1, 3, 2, p2s) == 3);
    assert(p1s[0] == 1);
    assert(p1s[1] == 2);
    assert(p1s[2] == 2);
    assert(p2s[0] == 1);
    assert(p2s[1] == 3);
    assert(p2s[2] == 6);

    char* p1ns[3]; char* p2ns[3];
    assert(get_group_parent_names(d, f1, 3, 1, p1ns) == 3);
    assert(get_group_parent_names(d, f1, UNINITIALISED, 2, p2ns) == 3);
    assert(strncmp(p1ns[0],"G01",sizeof(char)*5)==0);
    assert(strncmp(p1ns[1],"G02",sizeof(char)*5)==0);
    assert(strncmp(p1ns[2],"G02",sizeof(char)*5)==0);
    assert(strncmp(p2ns[0],"G01",sizeof(char)*5)==0);
    assert(strncmp(p2ns[1],"G03",sizeof(char)*5)==0);
    assert(strncmp(p2ns[2],"G06",sizeof(char)*5)==0);

    delete_group(d,f1);

    // Missing full pedigree getter check.

    printf("...legacy getters work correctly.\n");

    return 0;
}

int test_data_access(SimData* d, int gp) {
    test_iterators(d,gp);

    test_getters(d,gp);

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
    unsigned int randomSeed = time(NULL);
    if (argc > 1) {
        randomSeed = strtol(argv[1], 0, 0);
    }

	printf("Testing functionality ...");

	// test random number generators

	// test matrix operations

	// test SimData loaders
	printf("\nNow testing loader functions:\n");
    SimData* d = create_empty_simdata(randomSeed);
	int g0 = test_loaders(d);
	printf("\t\t-> Loader functions all clear\n");

    printf("\nNow testing group manipulation functions:\n");
    g0 = test_grouping(d, g0);
    printf("\t\t-> Group manipulation functions all clear\n");

    //test data access functions
    printf("\nNow testing data access functions:\n");
    test_data_access(d, g0);

	// test effect calculators
	printf("\nNow testing GEBV calculator:\n");
	test_effect_calculators(d, g0);
    test_optimal_calculators(d, g0);
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

	// test SimData deletors.
	printf("\nNow testing deletor functions:\n");
	test_deletors(d, g0);
	printf("\t\t-> Deletor functions all clear\n");

	printf("\n------- All tests passed. -------\n");


	//testing new grouping functions
    d = create_empty_simdata(randomSeed);
    g0 = load_all_simdata(d, "./gt_parents_mr2_50-trimto-5000.txt",
			 "./genetic-map_5112-trimto5000.txt",
             "./qtl_mr2.eff-processed.txt");

    /*g0 = load_all_simdata(d, "./gt_parents_mr2_3000x30000.txt",
                          "./genetic-map_huge30000.txt",
                          "./eff_huge30000.txt");*/

    /*g0 = load_all_simdata(d, "../../HealthAndHorns/AACo_HD_geno_101122/AACOg2.txt",
                              "../../HealthAndHorns/AACo_HD_geno_101122/AACo_HD_geno_101122.map",
                              "./eff_huge30000.txt");
    fprintf(stdout, "REALLY HERE\n");*/

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
    int g1 = cross_these_combinations(d,nc,crosses[0],crosses[1],gens);

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
    GenOptions gens =  {.will_name_offspring=TRUE, .offspring_name_prefix="cr", .family_size=243,
        .will_track_pedigree=TRUE, .will_allocate_ids=TRUE,
        .filename_prefix="testcross", .will_save_pedigree_to_file=FALSE,
        .will_save_bvs_to_file=FALSE, .will_save_alleles_to_file=FALSE,
        .will_save_to_simdata=TRUE};
    cross_random_individuals(d, g0, 10, 0, gens);
    cross_random_individuals(d, g0, 10, 0, gens);
    cross_random_individuals(d, g0, 10, 0, gens);
    cross_random_individuals(d, g0, 10, 0, gens);
    cross_random_individuals(d, g0, 10, 0, gens);
    int m = split_by_bv(d, g0, 55, 1);
    assert(get_group_size(d, m) == 50);

    delete_simdata(d);

    // Testing 5 allele effects bug? Replicability.
    d = create_empty_simdata(randomSeed);
    g0 = load_all_simdata(d, "./gt_parents_mr2_50-trimto-5000.txt",
             "./genetic-map_5112-trimto5000.txt",
             "./qtl_5test.txt");


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

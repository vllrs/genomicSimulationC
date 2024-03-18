#include "sim-test.h"

float calculate_heterozygosity(SimData* d, GroupNum group_number) {
	int hetcount = 0;
    int gn = get_group_size(d, group_number);
    assert(gn < 1000);
    char* galleles[1000];
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

int test_savers(unsigned int rseed) {
    SimData* d = create_empty_simdata(rseed);
    GroupNum g0 = just_load(d);
    assert(d->n_markers == 3);
    assert(d->n_eff_sets == 1);
    assert(d->m->n_genotypes == 6);

    // create some interesting groups
    GenOptions g = {.will_name_offspring=GSC_TRUE,
            .offspring_name_prefix="F1",
            .family_size=1,.will_track_pedigree=GSC_TRUE,
            .will_allocate_ids=GSC_TRUE,
            .filename_prefix="",
            .will_save_pedigree_to_file=GSC_FALSE,
            .will_save_bvs_to_file=0,
            .will_save_alleles_to_file=GSC_FALSE,
            .will_save_to_simdata=GSC_TRUE};
    GroupNum f1 = make_random_crosses(d, g0, 5, 2, g); // produce 5 offspring
    g.will_track_pedigree = GSC_TRUE;
    g.will_allocate_ids = GSC_FALSE;
    g.offspring_name_prefix = "s";
    GroupNum f2 = self_n_times(d,2,f1,g); // produce 5 offspring that know their parents but will be anonymous parents themselves
    g.will_name_offspring = GSC_FALSE;
    g.will_allocate_ids = GSC_TRUE;
    GroupNum f3 = make_doubled_haploids(d,f2,g); // produce 5 offspring that don't have names or know their parents.

    size_t toprint[] = {0,1,//2,3,4,5 from g0
                     6,7,//8,9,10 from f1
                     11,15,//11,12,13,14,15 from f2
                     16,17//18,19,20 from f3
                    };
    GroupNum printingGroup = make_group_from(d, 8, toprint);

    /*GroupNum gout[5];
    int gs[5];
    get_existing_group_counts(d,5,gout,gs);*/

    // try saving genotypes save_allele_matrix save_transposed_allele_matrix save_group_genotypes save_transposed_group_genotypes
    FILE* fp;
    if ((fp = fopen("test1_save_allele_matrix.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_names_header(fp,d->n_markers,(const char**) d->markers);
    save_allele_matrix(fp,d->m);
    fclose(fp);
    //assert(compareFileToString("test1_save_allele_matrix.txt", TEST1_TRUTH_save_allele_matrix)==0);
    //remove("test1_save_allele_matrix.txt");


    if ((fp = fopen("test1_save_transposed_allele_matrix.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_transposed_allele_matrix(fp,d->m,(const char**) d->markers);
    fclose(fp);
    //assert(compareFileToString("test1_save_transposed_allele_matrix.txt", TEST1_TRUTH_save_transposed_allele_matrix)==0);
    //remove("test1_save_transposed_allele_matrix.txt");


    if ((fp = fopen("test1_save_group_genotypes.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_group_genotypes(fp,d,printingGroup);
    fclose(fp);
    //assert(compareFileToString("test1_save_group_genotypes.txt", TEST1_TRUTH_save_group_genotypes)==0);
    //remove("test1_save_group_genotypes.txt");


    if ((fp = fopen("test1_save_transposed_group_genotypes.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_transposed_group_genotypes(fp,d,printingGroup);
    fclose(fp);
    //assert(compareFileToString("test1_save_transposed_group_genotypes.txt", TEST1_TRUTH_save_transposed_group_genotypes)==0);
    //remove("test1_save_transposed_group_genotypes.txt");


    // try saving counts save_count_matrix save_group_count_matrix
    if ((fp = fopen("test1_save_count_matrix.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_count_matrix(fp,d,'A');
    fclose(fp);
    //assert(compareFileToString("test1_save_count_matrix.txt", TEST1_TRUTH_save_count_matrix)==0);
    //remove("test1_save_count_matrix.txt");


    if ((fp = fopen("test1_save_group_count_matrix.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_group_count_matrix(fp,d,'T',printingGroup);
    fclose(fp);
    //assert(compareFileToString("test1_save_group_count_matrix.txt", TEST1_TRUTH_save_group_count_matrix)==0);
    //remove("test1_save_group_count_matrix.txt");


    printf("...genotype matrix file savers produce the expected output formats\n");


    // try saving bvs save_bvs save_group_bvs
    EffectID effSet1 = {.id=1};

    if ((fp = fopen("test1_save_bvs.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_bvs(fp,d,effSet1);
    fclose(fp);
    //assert(compareFileToString("test1_save_bvs.txt", TEST1_TRUTH_save_bvs)==0);
    //remove("test1_save_bvs.txt");


    if ((fp = fopen("test1_save_group_bvs.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_group_bvs(fp,d,printingGroup,effSet1);
    fclose(fp);
    //assert(compareFileToString("test1_save_group_bvs.txt", TEST1_TRUTH_save_group_bvs)==0);
    //remove("test1_save_group_bvs.txt");


    // try saving local gebvs save_marker_blocks calculate_local_bvs
    MarkerBlocks exampleMB = create_evenlength_blocks_each_chr(d,1);

    if ((fp = fopen("test1_save_marker_blocks.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_marker_blocks(fp,d,exampleMB);
    fclose(fp);
    //assert(compareFileToString("test1_save_marker_blocks.txt", TEST1_TRUTH_save_marker_blocks)==0);
    //remove("test1_save_marker_blocks.txt");


    calculate_local_bvs(d,exampleMB,effSet1,"test1_save_local_bvs.txt");
    //assert(compareFileToString("test1_save_local_bvs.txt", TEST1_TRUTH_save_local_bvs)==0);
    //remove("test1_save_local_bvs.txt");


    calculate_group_local_bvs(d,exampleMB,effSet1,"test1_save_group_local_bvs.txt",printingGroup);
    assert(compareFileToString("test1_save_group_local_bvs.txt", TEST1_TRUTH_save_group_local_bvs)==0);
    remove("test1_save_group_local_bvs.txt");


    printf("...breeding value file savers produce the expected output formats\n");


    // try saving one-step pedigrees save_group_one_step_pedigree save_one_step_pedigree
    g = (GenOptions){.will_name_offspring=GSC_TRUE,
            .offspring_name_prefix="F2",
            .family_size=1,.will_track_pedigree=GSC_TRUE,
            .will_allocate_ids=GSC_TRUE,
            .filename_prefix="",
            .will_save_pedigree_to_file=GSC_FALSE,
            .will_save_bvs_to_file=0,
            .will_save_alleles_to_file=GSC_FALSE,
            .will_save_to_simdata=GSC_TRUE};
    GroupNum f2b = make_random_crosses(d,f1,1,1,g);

    if ((fp = fopen("test1_save_one_step_pedigree.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_one_step_pedigree(fp,d);
    fclose(fp);
    assert(compareFileToString("test1_save_one_step_pedigree.txt", TEST1_TRUTH_save_one_step_pedigrees)==0);
    remove("test1_save_one_step_pedigree.txt");


    if ((fp = fopen("test1_save_group_one_step_pedigree.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_group_one_step_pedigree(fp,d,printingGroup);
    fclose(fp);
    assert(compareFileToString("test1_save_group_one_step_pedigree.txt", TEST1_TRUTH_save_group_one_step_pedigrees)==0);
    remove("test1_save_group_one_step_pedigree.txt");


    // try saving full pedigrees save_group_full_pedigree save_full_pedigree
    if ((fp = fopen("test1_save_full_pedigree.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_full_pedigree(fp,d);
    fclose(fp);
    assert(compareFileToString("test1_save_full_pedigree.txt", TEST1_TRUTH_save_full_pedigrees)==0);
    remove("test1_save_full_pedigree.txt");


    if ((fp = fopen("test1_save_group_full_pedigree.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_group_full_pedigree(fp,d,printingGroup);
    fclose(fp);
    assert(compareFileToString("test1_save_group_full_pedigree.txt", TEST1_TRUTH_save_group_full_pedigrees)==0);
    remove("test1_save_group_full_pedigree.txt");


    // try save-as-you-go savers (ideally you'd run this with a very low CONTIG_WIDTH, and also for all crossing funcs.
    size_t parentix[] = {20};
    GroupNum parent = make_group_from(d,1,parentix);
    g = (GenOptions){.will_name_offspring=GSC_FALSE,
            .offspring_name_prefix="",
            .family_size=(CONTIG_WIDTH+5),.will_track_pedigree=GSC_TRUE,
            .will_allocate_ids=GSC_FALSE,
            .filename_prefix="test1_save_as_you_go",
            .will_save_pedigree_to_file=GSC_TRUE,
            .will_save_bvs_to_file=effSet1.id,
            .will_save_alleles_to_file=GSC_TRUE,
            .will_save_to_simdata=GSC_FALSE};
    make_doubled_haploids(d,parent,g);

    assert(compareRepeatingFileToTable("test1_save_as_you_go-genotype.txt", CONTIG_WIDTH+5,
                                       TEST1_TRUTH_sayg_genotype_header, TEST1_TRUTH_sayg_genotype_bodyrow)==0);
    remove("test1_save_as_you_go-genotype.txt");
    assert(compareRepeatingFileToTable("test1_save_as_you_go-bv.txt", CONTIG_WIDTH+5,
                                       NULL, TEST1_TRUTH_sayg_bv_bodyrow)==0);
    remove("test1_save_as_you_go-bv.txt");
    assert(compareRepeatingFileToTable("test1_save_as_you_go-pedigree.txt", CONTIG_WIDTH+5,
                                       NULL, TEST1_TRUTH_sayg_pedigree_bodyrow)==0);
    remove("test1_save_as_you_go-pedigree.txt");


    return 0;
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

GroupNum just_load(SimData* d) {
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
    if ((fp = fopen("a-test-eff2.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    fwrite(HELPER_EFF2, sizeof(char), strlen(HELPER_EFF2), fp);
    fclose(fp);

    struct GroupAndEffectSet gande = load_all_data(d, "a-test.txt", "a-test-map.txt", "a-test-eff.txt");
    GroupNum g0 = gande.group;

    remove("a-test.txt");
    remove("a-test-map.txt");
    remove("a-test-eff.txt");
    return g0;
}


GroupNum test_loaders(SimData* d) {
    GroupNum g0 = just_load(d);

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

    assert(d->n_eff_sets == 1);
    assert(d->e[0].effects.rows == 2);
    assert(d->e[0].effects.cols == 3);
	// don't mind which order the effects are in so just figure it out, then check based on that.
	int apos = 0;
    if (d->e[0].effect_names[0] == 'A') {
        assert(d->e[0].effect_names[0] == 'A');
        assert(d->e[0].effect_names[1] == 'T');
	} else {
		apos = 1;
        assert(d->e[0].effect_names[0] == 'T');
        assert(d->e[0].effect_names[1] == 'A');
	}
    assert(fabs(d->e[0].effects.matrix[apos][0] - (-0.8)) < TOL);
    assert(fabs(d->e[0].effects.matrix[apos][1] - (-0.1)) < TOL);
    assert(fabs(d->e[0].effects.matrix[apos][2] - (0.1)) < TOL);
	int tpos = 1 - apos;
    assert(fabs(d->e[0].effects.matrix[tpos][0] - (0.9)) < TOL);
    assert(fabs(d->e[0].effects.matrix[tpos][1] - (-0.5)) < TOL);
    assert(fabs(d->e[0].effects.matrix[tpos][2] - (-0.1)) < TOL);
    printf("...marker effects loaded correctly\n");

    assert(load_effects(d, "a-test-eff2.txt").id==2);
    assert(d->n_eff_sets == 2);
    assert(d->e[0].effects.rows == 2);
    assert(d->e[0].effects.cols == 3);
    // don't mind which order the effects are in so just figure it out, then check based on that.
    if (d->e[0].effect_names[0] == 'A') {
        assert(d->e[0].effect_names[0] == 'A');
        assert(d->e[0].effect_names[1] == 'T');
    } else {
        apos = 1;
        assert(d->e[0].effect_names[0] == 'T');
        assert(d->e[0].effect_names[1] == 'A');
    }
    assert(fabs(d->e[0].effects.matrix[apos][0] - (-0.8)) < TOL);
    assert(fabs(d->e[0].effects.matrix[apos][1] - (-0.1)) < TOL);
    assert(fabs(d->e[0].effects.matrix[apos][2] - (0.1)) < TOL);
    assert(fabs(d->e[0].effects.matrix[tpos][0] - (0.9)) < TOL);
    assert(fabs(d->e[0].effects.matrix[tpos][1] - (-0.5)) < TOL);
    assert(fabs(d->e[0].effects.matrix[tpos][2] - (-0.1)) < TOL);

    assert(d->e[1].effects.rows == 1);
    assert(d->e[1].effects.cols == 3);
    assert(d->e[1].effect_names[0] == 'A');
    assert(fabs(d->e[1].effects.matrix[apos][0] - 1) < TOL);
    assert(fabs(d->e[1].effects.matrix[apos][1] - 0) < TOL);
    assert(fabs(d->e[1].effects.matrix[apos][2] - 0) < TOL);

    printf("...second set of marker effects loaded correctly\n");
    remove("a-test-eff2.txt");

    assert(d->current_id.id == 6);
    assert(d->m); // != NULL
	assert(d->m->n_markers == 3);
	assert(d->m->n_genotypes == 6);
    assert(d->m->ids[0].id == 1);
    assert(d->m->ids[1].id == 2);
    assert(d->m->ids[2].id == 3);
    assert(d->m->ids[3].id == 4);
    assert(d->m->ids[4].id == 5);
    assert(d->m->ids[5].id == 6);
    assert(d->m->pedigrees[0][0].id == 0 && d->m->pedigrees[1][0].id == 0);
    assert(d->m->pedigrees[0][1].id == 0 && d->m->pedigrees[1][1].id == 0);
    assert(d->m->pedigrees[0][2].id == 0 && d->m->pedigrees[1][2].id == 0);
    assert(d->m->pedigrees[0][3].id== 0 && d->m->pedigrees[1][3].id == 0);
    assert(d->m->pedigrees[0][4].id == 0 && d->m->pedigrees[1][4].id == 0);
    assert(d->m->pedigrees[0][5].id == 0 && d->m->pedigrees[1][5].id == 0);
    assert(d->m->groups[0].num == g0.num);
    assert(d->m->groups[1].num == g0.num);
    assert(d->m->groups[2].num == g0.num);
    assert(d->m->groups[3].num == g0.num);
    assert(d->m->groups[4].num == g0.num);
    assert(d->m->groups[5].num == g0.num);
    assert(g0.num > 0);

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

GroupNum test_labels(SimData *d, GroupNum g0) {
    // g0: contents of "a-test.txt", contains 6 genotypes

    // Can create a first label.
    assert(d->n_labels == 0);
    const int label1value = 0;
    const LabelID label1 = create_new_label(d, label1value);
    const int label1index = gsc_get_index_of_label(d, label1);
    assert(d->n_labels == 1);
    assert(label1.id == 1); // expected
    assert(label1index == 0); // expected
    assert(d->label_ids != NULL);
    assert(d->label_ids[label1index].id == label1.id);
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
    const LabelID label2 = create_new_label(d, label2value);
    const int label2index = gsc_get_index_of_label(d, label2);
    assert(d->n_labels == 2);
    assert(label2.id == 2); // expected
    assert(label2index == 1); // expected
    assert(d->label_ids != NULL);
    assert(d->label_ids[label1index].id == label1.id);
    assert(d->label_ids[label2index].id == label2.id);
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
    change_label_default(d,label2,newlabel2default);
    assert(d->n_labels == 2);
    assert(d->label_ids != NULL);
    assert(d->label_ids[label1index].id == label1.id);
    assert(d->label_ids[label2index].id == label2.id);
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
    change_label_to(d,g0,label1,newlabel1value);
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
    const int newlabel2values[] = {11, 12, 13, 14, 15, 16, 17};
    change_label_to_values(d, g0, 0, label2, 7, newlabel2values);
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
    const int newerlabel1values[] = {0,1,-1,100};
    change_label_to_values(d, g0, 1, label1, 4, newerlabel1values);
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

    // Test name-setters here too, they use the same procedure as change_label_to_values
    // (Tests when setting whole body not just group)
    const char* newNames[] = {"name0", "aname", "two", "THR33", "7", "name$$"};
    change_names_to_values(d,NO_GROUP,0,6,newNames);
    assert(strncmp(d->m->names[0],newNames[0],strlen(newNames[0])) == 0);
    assert(strncmp(d->m->names[1],newNames[1],strlen(newNames[1])) == 0);
    assert(strncmp(d->m->names[2],newNames[2],strlen(newNames[2])) == 0);
    assert(strncmp(d->m->names[3],newNames[3],strlen(newNames[3])) == 0);
    assert(strncmp(d->m->names[4],newNames[4],strlen(newNames[4])) == 0);
    assert(strncmp(d->m->names[5],newNames[5],strlen(newNames[5])) == 0);
    const char* betterNames[] = {"G01", "G02", "G03", "G04", "G05", "G06"};
    change_names_to_values(d,g0,0,6,betterNames);
    assert(strncmp(d->m->names[0],betterNames[0],strlen(betterNames[0])) == 0);
    assert(strncmp(d->m->names[1],betterNames[1],strlen(betterNames[1])) == 0);
    assert(strncmp(d->m->names[2],betterNames[2],strlen(betterNames[2])) == 0);
    assert(strncmp(d->m->names[3],betterNames[3],strlen(betterNames[3])) == 0);
    assert(strncmp(d->m->names[4],betterNames[4],strlen(betterNames[4])) == 0);
    assert(strncmp(d->m->names[5],betterNames[5],strlen(betterNames[5])) == 0);

    // Can increment a label
    const int increment = -4+1;
    change_label_by_amount(d, g0, label1, 1);
    change_label_by_amount(d, g0, label1, -4);
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
                    .will_track_pedigree = GSC_FALSE,
                    .will_name_offspring = GSC_FALSE,
                    .offspring_name_prefix = NULL,
                    .will_allocate_ids = GSC_FALSE,
                    .filename_prefix = NULL,
                    .will_save_alleles_to_file = GSC_FALSE,
                    .will_save_bvs_to_file = -1,
                    .will_save_pedigree_to_file = GSC_FALSE,
                    .will_save_to_simdata = GSC_TRUE};
    assert(d->n_groups == 1);
    GroupNum f1 = make_random_crosses(d,g0,4,0,g);
    assert(d->n_groups == 2);
    // Check that the labels are set to their defaults
    assert(d->m->labels[label1index][6+0] == label1value);
    assert(d->m->labels[label1index][6+1] == label1value);
    assert(d->m->labels[label1index][6+2] == label1value);
    assert(d->m->labels[label1index][6+3] == label1value);
    assert(d->m->labels[label2index][6+0] == newlabel2default);
    assert(d->m->labels[label2index][6+1] == newlabel2default);
    assert(d->m->labels[label2index][6+2] == newlabel2default);
    assert(d->m->labels[label2index][6+3] == newlabel2default);

    change_label_to(d,f1,label1,newlabel1value + increment);
    change_label_to_values(d,f1,0,label2,4,newerlabel1values);

    // test can split by label (across groups)
    GroupNum groupB = split_by_label_value(d, NO_GROUP, label1, newlabel1value + increment);
    assert(d->n_groups == 3);
    assert(g0.num != groupB.num && f1.num != groupB.num);
    assert(get_group_size(d, groupB) == 2+4);
    size_t Bindexes[2+4];
    get_group_indexes(d, groupB, 2+4, Bindexes);
    assert(Bindexes[0] == 0);
    assert(Bindexes[1] == 5);
    assert(Bindexes[2] == 6 && Bindexes[3] == 7 && Bindexes[4] == 8 && Bindexes[5] == 9);

    size_t outtakes[4] = {6,7,8,9};
    GroupNum f1outtakes = make_group_from(d,4,outtakes);
    assert(d->n_groups == 3); // not 4, because it's corrected by get_new_group_num inside make_group_from
    GroupNum toCombine[2] = {f1, f1outtakes};
    f1 = combine_groups(d, 2, toCombine);
    assert(d->n_groups == 3);
    assert(get_group_size(d, f1) == 4);
    toCombine[0] = g0;
    toCombine[1] = groupB;
    g0 = combine_groups(d, 2, toCombine);
    assert(d->n_groups == 2);
    assert(get_group_size(d, g0) == 6);

    // test can split by label range (across groups)
    groupB = split_by_label_range(d, NO_GROUP, label2, 1, 15);
    assert(d->n_groups == 3);
    assert(g0.num != groupB.num && f1.num != groupB.num);
    assert(get_group_size(d, groupB) == 6);
    get_group_indexes(d, groupB, 6, Bindexes);
    assert(Bindexes[0] == 0);
    assert(Bindexes[1] == 1);
    assert(Bindexes[2] == 2);
    assert(Bindexes[3] == 3);
    assert(Bindexes[4] == 4);
    assert(Bindexes[5] == 7);

    outtakes[0] = 7;
    f1outtakes = make_group_from(d,1,outtakes);
    assert(d->n_groups == 4);
    toCombine[0] = f1;
    toCombine[1] = f1outtakes;
    f1 = combine_groups(d, 2, toCombine);
    assert(d->n_groups == 3);
    assert(get_group_size(d, f1) == 4);
    toCombine[0] = g0;
    toCombine[1] = groupB;
    g0 = combine_groups(d, 2, toCombine);
    assert(d->n_groups == 2);
    assert(get_group_size(d, g0) == 6);

    // and can split from group (within group)
    groupB = split_by_label_value(d, g0, label1, newlabel1value + increment);
    assert(d->n_groups == 3);
    assert(g0.num != groupB.num && f1.num != groupB.num);
    assert(get_group_size(d, groupB) == 2);
    get_group_indexes(d, groupB, 2, Bindexes);
    assert(Bindexes[0] == 0);
    assert(Bindexes[1] == 5);

    toCombine[0] = g0;
    toCombine[1] = groupB;
    g0 = combine_groups(d, 2, toCombine);
    assert(d->n_groups == 2);
    assert(get_group_size(d, g0) == 6);
    assert(get_group_size(d, f1) == 4);


    // and can split from group by label range (within group)
    groupB = split_by_label_range(d, g0, label2, 1, 15);
    assert(d->n_groups == 3);
    assert(g0.num != groupB.num && f1.num != groupB.num);
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
    assert(d->n_groups == 2);
    assert(get_group_size(d, g0) == 6);
    assert(get_group_size(d, f1) == 4);
    delete_group(d,f1);
    assert(d->n_groups == 1);

    // Check label deletion works
    delete_label(d,label1);
    assert(d->n_labels == 1);
    assert(d->m->n_labels == 1);
    assert(d->label_ids[0].id == label2.id);
    assert(d->label_defaults[0] == newlabel2default);
    assert(d->m->labels[0][0] == 11);

    printf("...Label manipulation runs correctly\n");

    return g0;
}

GroupNum test_random_splits(SimData *d, GroupNum g0) {
    // g0: contents of "a-test.txt", contains 6 genotypes

    // Can split into 2; repeat 10x for confidence
    assert(d->n_groups == 1);
    for (int i = 0; i < 10; ++i) {
        GroupNum grpB = split_evenly_into_two(d, g0);
        assert(d->n_groups == 2);
        assert(get_group_size(d,g0) == 3);
        assert(get_group_size(d,grpB) == 3);

        GroupNum toMerge[2] = {g0, grpB};
        g0 = combine_groups(d,2,toMerge);
        assert(d->n_groups == 1);
        assert(get_group_size(d,g0) == 6);
    }

    // Can split into 3; repeat 10x for confidence
    for (int i = 0; i < 10; ++i) {
        GroupNum grpB[3];
        split_evenly_into_n(d,g0,3,grpB);
        assert(d->n_groups == 3);
        assert(get_group_size(d,grpB[0]) == 2);
        assert(get_group_size(d,grpB[1]) == 2);
        assert(get_group_size(d,grpB[2]) == 2);

        g0 = combine_groups(d,3,grpB);
        assert(d->n_groups == 1);
        assert(get_group_size(d,g0) == 6);
    }

    // Can split into bins; repeat 10x for confidence
    for (int i = 0; i < 10; ++i) {
        GroupNum grpB[3];
        int grpBsizes[3] = {1,3,2};
        split_into_buckets(d,g0,3,grpBsizes,grpB);
        assert(d->n_groups == 3);
        assert(get_group_size(d,grpB[0]) == 1);
        assert(get_group_size(d,grpB[1]) == 3);
        assert(get_group_size(d,grpB[2]) == 2);

        g0 = combine_groups(d,3,grpB);
        assert(d->n_groups == 1);
        assert(get_group_size(d,g0) == 6);
    }

    // In future replace these with something that uses the random generator seed
    // Can split by coinflip
    int grpBsizes = 0;
    for (int i = 0; i < 10; ++i) {
        GroupNum grpB = split_randomly_into_two(d, g0);
        assert(d->n_groups == 2);
        int grpBsize = get_group_size(d,grpB);
        assert(get_group_size(d,g0) + grpBsize == 6);
        grpBsizes += grpBsize;

        GroupNum toMerge[2] = {g0, grpB};
        g0 = combine_groups(d,2,toMerge);
        assert(get_group_size(d,g0) == 6);
    }
    assert(grpBsizes > 20 && grpBsizes < 40); // about 1.3% chance of failure by random chance
    assert(get_existing_groups(d,NULL) == 1);
    assert(d->n_groups == 1);

    // Can split by thirds
    grpBsizes = 0;
    int grpBsizes2 = 0;
    for (int i = 0; i < 10; ++i) {
        GroupNum grpB[3];
        split_randomly_into_n(d,g0,3,grpB);
        assert(d->n_groups <= 3); // sometimes, randomly, we may not get all three groups.
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
    assert(get_existing_groups(d,NULL) == 1);
    assert(d->n_groups == 1);

    // Can split by any probability
    grpBsizes = 0;
    for (int i = 0; i < 10; ++i) {
        GroupNum grpB[2];
        double grpBprobs = 3.0/7.0;
        split_by_probabilities(d,g0,2,&grpBprobs,grpB);
        assert(d->n_groups == 2);
        int grpBsize = get_group_size(d,grpB[0]);
        assert(get_group_size(d,grpB[1]) + grpBsize == 6);
        grpBsizes += grpBsize;

        g0 = combine_groups(d,2,grpB);
        assert(get_group_size(d,g0) == 6);
    }
    assert(grpBsizes > 16 && grpBsizes < 39); // about 0.75% chance of failure by random chance
    assert(get_existing_groups(d,NULL) == 1);
    assert(d->n_groups == 1);

    grpBsizes = 0;
    for (int i = 0; i < 10; ++i) {
        GroupNum grpB[2];
        double grpBprobs[2] = {1./12.,11./12.};
        split_by_probabilities(d,g0,2,grpBprobs,grpB);
        assert(d->n_groups == 2);
        int grpBsize = get_group_size(d,grpB[0]);
        assert(get_group_size(d,grpB[1]) + grpBsize == 6);
        grpBsizes += grpBsize;

        g0 = combine_groups(d,2,grpB);
        assert(get_group_size(d,g0) == 6);
    }
    assert(grpBsizes < 11); // about 1% chance of failure by random chance
    assert(get_existing_groups(d,NULL) == 1);
    assert(d->n_groups == 1);

    printf("...Random group splitting runs correctly\n");

    return g0;
}

GroupNum test_specific_splits(SimData *d, GroupNum g0) {
    // g0: contents of "a-test.txt", contains 6 genotypes
    // and there are few enough genotypes they all fit into the first AlleleMatrix

    // Can split into individuals and recombine
    const int g0size = get_group_size(d, g0);
    assert(d->n_groups == 1);
    assert(g0size == 6); // pre-test test validity requirement
    GroupNum g0indivs[10];
    unsigned nindivs = split_into_individuals(d, g0, 10, g0indivs);
    assert(nindivs == 6);
    assert(d->n_groups == 6);
    assert(g0indivs[0].num > 0);
    assert(g0indivs[1].num > 0);
    assert(g0indivs[1].num != g0indivs[0].num);
    assert(g0indivs[2].num > 0);
    assert(g0indivs[2].num != g0indivs[0].num && g0indivs[2].num != g0indivs[1].num);
    assert(g0indivs[3].num > 0);
    assert(g0indivs[3].num != g0indivs[0].num && g0indivs[3].num != g0indivs[1].num && g0indivs[3].num != g0indivs[2].num);
    assert(g0indivs[4].num > 0);
    assert(g0indivs[4].num != g0indivs[0].num && g0indivs[4].num != g0indivs[1].num
            && g0indivs[4].num != g0indivs[2].num && g0indivs[4].num != g0indivs[3].num);
    assert(g0indivs[5].num > 0);
    assert(g0indivs[5].num != g0indivs[0].num && g0indivs[5].num != g0indivs[1].num && g0indivs[5].num != g0indivs[2].num
            && g0indivs[5].num != g0indivs[3].num && g0indivs[5].num != g0indivs[4].num);

    g0 = combine_groups(d,6,g0indivs);
    assert(get_group_size(d,g0) == 6);
    assert(d->n_groups == 1);

    // Can split into families
    GenOptions g = {.family_size = 5,
                    .will_track_pedigree = GSC_TRUE,
                    .will_name_offspring = GSC_FALSE,
                    .offspring_name_prefix = NULL,
                    .will_allocate_ids = GSC_FALSE,
                    .filename_prefix = NULL,
                    .will_save_alleles_to_file = GSC_FALSE,
                    .will_save_bvs_to_file = -1,
                    .will_save_pedigree_to_file = GSC_FALSE,
                    .will_save_to_simdata = GSC_TRUE};
    GroupNum f1 = make_random_crosses(d, g0, 3, 1, g);
    int f1size = 5*3;
    GroupNum families[15];
    assert(d->n_groups == 2);
    int f1nfamilies = split_into_families(d,f1,f1size,families);
    assert(d->n_groups == 4);
    assert(f1nfamilies == 3);
    assert(families[0].num > 0);
    assert(families[1].num > 0);
    assert(families[2].num > 0);
    assert(families[0].num != families[1].num && families[1].num != families[2].num && families[0].num != families[2].num);
    assert(get_group_size(d,families[0]) == 5);
    assert(get_group_size(d,families[1]) == 5);
    assert(get_group_size(d,families[2]) == 5);
    // Check all groups contain families
    size_t f1family1[5];
    assert(get_group_indexes(d,families[0],5,f1family1) == 5);
    assert(d->m->pedigrees[0][f1family1[0]].id == d->m->pedigrees[0][f1family1[1]].id &&
            d->m->pedigrees[0][f1family1[0]].id == d->m->pedigrees[0][f1family1[2]].id &&
            d->m->pedigrees[0][f1family1[0]].id == d->m->pedigrees[0][f1family1[3]].id &&
            d->m->pedigrees[0][f1family1[0]].id == d->m->pedigrees[0][f1family1[4]].id);
    assert(d->m->pedigrees[1][f1family1[0]].id == d->m->pedigrees[1][f1family1[1]].id &&
            d->m->pedigrees[1][f1family1[0]].id == d->m->pedigrees[1][f1family1[2]].id &&
            d->m->pedigrees[1][f1family1[0]].id == d->m->pedigrees[1][f1family1[3]].id &&
            d->m->pedigrees[1][f1family1[0]].id == d->m->pedigrees[1][f1family1[4]].id);

    size_t f1family2[5];
    assert(get_group_indexes(d,families[1],5,f1family2) == 5);
    assert(d->m->pedigrees[0][f1family1[0]].id == d->m->pedigrees[0][f1family1[1]].id &&
            d->m->pedigrees[0][f1family1[0]].id == d->m->pedigrees[0][f1family1[2]].id &&
            d->m->pedigrees[0][f1family1[0]].id == d->m->pedigrees[0][f1family1[3]].id &&
            d->m->pedigrees[0][f1family1[0]].id == d->m->pedigrees[0][f1family1[4]].id);
    assert(d->m->pedigrees[1][f1family1[0]].id == d->m->pedigrees[1][f1family1[1]].id &&
            d->m->pedigrees[1][f1family1[0]].id == d->m->pedigrees[1][f1family1[2]].id &&
            d->m->pedigrees[1][f1family1[0]].id == d->m->pedigrees[1][f1family1[3]].id &&
            d->m->pedigrees[1][f1family1[0]].id == d->m->pedigrees[1][f1family1[4]].id);

    size_t f1family3[5];
    assert(get_group_indexes(d,families[2],5,f1family3) == 5);
    assert(d->m->pedigrees[0][f1family1[0]].id == d->m->pedigrees[0][f1family1[1]].id &&
            d->m->pedigrees[0][f1family1[0]].id == d->m->pedigrees[0][f1family1[2]].id &&
            d->m->pedigrees[0][f1family1[0]].id == d->m->pedigrees[0][f1family1[3]].id &&
            d->m->pedigrees[0][f1family1[0]].id == d->m->pedigrees[0][f1family1[4]].id);
    assert(d->m->pedigrees[1][f1family1[0]].id == d->m->pedigrees[1][f1family1[1]].id &&
            d->m->pedigrees[1][f1family1[0]].id == d->m->pedigrees[1][f1family1[2]].id &&
            d->m->pedigrees[1][f1family1[0]].id == d->m->pedigrees[1][f1family1[3]].id &&
            d->m->pedigrees[1][f1family1[0]].id == d->m->pedigrees[1][f1family1[4]].id);

    // Check separate groups don't represent the same family
    assert(d->m->pedigrees[0][f1family1[0]].id != d->m->pedigrees[0][f1family2[0]].id &&
            d->m->pedigrees[1][f1family1[0]].id != d->m->pedigrees[1][f1family2[0]].id);
    assert(d->m->pedigrees[0][f1family1[0]].id != d->m->pedigrees[0][f1family3[0]].id &&
            d->m->pedigrees[1][f1family1[0]].id != d->m->pedigrees[1][f1family3[0]].id);
    assert(d->m->pedigrees[0][f1family2[0]].id != d->m->pedigrees[0][f1family3[0]].id &&
            d->m->pedigrees[1][f1family2[0]].id != d->m->pedigrees[1][f1family3[0]].id);

    delete_group(d, families[0]);
    delete_group(d, families[1]);
    delete_group(d, families[2]);
    assert(d->n_groups == 1);

    // Can split into halfsib families
    int combinations[2][4];
    combinations[0][0] = 0; combinations[1][0] = 1;
    combinations[0][1] = 0; combinations[1][1] = 2;
    combinations[0][2] = 3; combinations[1][2] = 1;
    combinations[0][3] = 4; combinations[1][3] = 5;
    GroupNum fhs = make_targeted_crosses(d,4,combinations[0], combinations[1],g);
    assert(d->n_groups == 2);

    GroupNum halfsibfamilies[3];
    assert(split_into_halfsib_families(d,fhs,1,3,halfsibfamilies) == 3);
    assert(d->n_groups == 4);
    assert(halfsibfamilies[0].num > 0);
    assert(halfsibfamilies[1].num > 0);
    assert(halfsibfamilies[2].num > 0);

    // Check all halfsib families share the same parent 1
    int fhsfamily1size = get_group_size(d,halfsibfamilies[0]); // assumes the largest hsfamily will be the first one seen
    assert(fhsfamily1size == 10);
    size_t fhsfamily1[10];
    assert(get_group_indexes(d,halfsibfamilies[0],10,fhsfamily1) == fhsfamily1size);
    assert(d->m->pedigrees[0][fhsfamily1[0]].id == d->m->pedigrees[0][fhsfamily1[1]].id &&
            d->m->pedigrees[0][fhsfamily1[0]].id == d->m->pedigrees[0][fhsfamily1[2]].id &&
            d->m->pedigrees[0][fhsfamily1[0]].id == d->m->pedigrees[0][fhsfamily1[3]].id &&
            d->m->pedigrees[0][fhsfamily1[0]].id == d->m->pedigrees[0][fhsfamily1[4]].id &&
            d->m->pedigrees[0][fhsfamily1[0]].id == d->m->pedigrees[0][fhsfamily1[5]].id &&
            d->m->pedigrees[0][fhsfamily1[0]].id == d->m->pedigrees[0][fhsfamily1[6]].id &&
            d->m->pedigrees[0][fhsfamily1[0]].id == d->m->pedigrees[0][fhsfamily1[7]].id &&
            d->m->pedigrees[0][fhsfamily1[0]].id == d->m->pedigrees[0][fhsfamily1[8]].id &&
            d->m->pedigrees[0][fhsfamily1[0]].id == d->m->pedigrees[0][fhsfamily1[9]].id);

    int fhsfamily2size = get_group_size(d,halfsibfamilies[1]);
    assert(fhsfamily2size == 5);
    size_t fhsfamily2[5];
    assert(get_group_indexes(d,halfsibfamilies[1],5,fhsfamily2) == fhsfamily2size);
    assert(d->m->pedigrees[0][fhsfamily2[0]].id == d->m->pedigrees[0][fhsfamily2[1]].id &&
            d->m->pedigrees[0][fhsfamily2[0]].id == d->m->pedigrees[0][fhsfamily2[2]].id &&
            d->m->pedigrees[0][fhsfamily2[0]].id == d->m->pedigrees[0][fhsfamily2[3]].id &&
            d->m->pedigrees[0][fhsfamily2[0]].id == d->m->pedigrees[0][fhsfamily2[4]].id);

    int fhsfamily3size = get_group_size(d,halfsibfamilies[2]);
    assert(fhsfamily3size == 5);
    size_t fhsfamily3[5];
    assert(get_group_indexes(d,halfsibfamilies[2],5,fhsfamily3) == fhsfamily3size);
    assert(d->m->pedigrees[0][fhsfamily3[0]].id == d->m->pedigrees[0][fhsfamily3[1]].id &&
            d->m->pedigrees[0][fhsfamily3[0]].id == d->m->pedigrees[0][fhsfamily3[2]].id &&
            d->m->pedigrees[0][fhsfamily3[0]].id == d->m->pedigrees[0][fhsfamily3[3]].id &&
            d->m->pedigrees[0][fhsfamily3[0]].id == d->m->pedigrees[0][fhsfamily3[4]].id);

    // and that none of the across-groups have the same parent 1
    assert(d->m->pedigrees[0][fhsfamily1[0]].id != d->m->pedigrees[0][fhsfamily2[0]].id);
    assert(d->m->pedigrees[0][fhsfamily1[0]].id != d->m->pedigrees[0][fhsfamily3[0]].id);
    assert(d->m->pedigrees[0][fhsfamily2[0]].id != d->m->pedigrees[0][fhsfamily3[0]].id);

    // Then check halfsibs in the other direction
    fhs = combine_groups(d,3,halfsibfamilies);
    assert(d->n_groups == 2);
    assert(split_into_halfsib_families(d,fhs,2,3,halfsibfamilies) == 3);
    assert(d->n_groups == 4);
    assert(halfsibfamilies[0].num > 0);
    assert(halfsibfamilies[1].num > 0);
    assert(halfsibfamilies[2].num > 0);

    // Check all halfsib families share the same parent 1
    fhsfamily1size = get_group_size(d,halfsibfamilies[0]); // assumes the largest hsfamily will be the first one seen
    assert(fhsfamily1size == 10);
    assert(get_group_indexes(d,halfsibfamilies[0],10,fhsfamily1) == fhsfamily1size);
    assert(d->m->pedigrees[1][fhsfamily1[0]].id == d->m->pedigrees[1][fhsfamily1[1]].id &&
            d->m->pedigrees[1][fhsfamily1[0]].id == d->m->pedigrees[1][fhsfamily1[2]].id &&
            d->m->pedigrees[1][fhsfamily1[0]].id == d->m->pedigrees[1][fhsfamily1[3]].id &&
            d->m->pedigrees[1][fhsfamily1[0]].id == d->m->pedigrees[1][fhsfamily1[4]].id &&
            d->m->pedigrees[1][fhsfamily1[0]].id == d->m->pedigrees[1][fhsfamily1[5]].id &&
            d->m->pedigrees[1][fhsfamily1[0]].id == d->m->pedigrees[1][fhsfamily1[6]].id &&
            d->m->pedigrees[1][fhsfamily1[0]].id == d->m->pedigrees[1][fhsfamily1[7]].id &&
            d->m->pedigrees[1][fhsfamily1[0]].id == d->m->pedigrees[1][fhsfamily1[8]].id &&
            d->m->pedigrees[1][fhsfamily1[0]].id == d->m->pedigrees[1][fhsfamily1[9]].id);

    fhsfamily2size = get_group_size(d,halfsibfamilies[1]);
    assert(fhsfamily2size == 5);
    assert(get_group_indexes(d,halfsibfamilies[1],5,fhsfamily2) == fhsfamily2size);
    assert(d->m->pedigrees[1][fhsfamily2[0]].id == d->m->pedigrees[1][fhsfamily2[1]].id &&
            d->m->pedigrees[1][fhsfamily2[0]].id == d->m->pedigrees[1][fhsfamily2[2]].id &&
            d->m->pedigrees[1][fhsfamily2[0]].id == d->m->pedigrees[1][fhsfamily2[3]].id &&
            d->m->pedigrees[1][fhsfamily2[0]].id == d->m->pedigrees[1][fhsfamily2[4]].id);

    fhsfamily3size = get_group_size(d,halfsibfamilies[2]);
    assert(fhsfamily3size == 5);
    assert(get_group_indexes(d,halfsibfamilies[2],5,fhsfamily3) == fhsfamily3size);
    assert(d->m->pedigrees[1][fhsfamily3[0]].id == d->m->pedigrees[1][fhsfamily3[1]].id &&
            d->m->pedigrees[1][fhsfamily3[0]].id == d->m->pedigrees[1][fhsfamily3[2]].id &&
            d->m->pedigrees[1][fhsfamily3[0]].id== d->m->pedigrees[1][fhsfamily3[3]].id &&
            d->m->pedigrees[1][fhsfamily3[0]].id == d->m->pedigrees[1][fhsfamily3[4]].id);

    // and that none of the across-groups have the same parent 1
    assert(d->m->pedigrees[1][fhsfamily1[0]].id != d->m->pedigrees[1][fhsfamily2[0]].id);
    assert(d->m->pedigrees[1][fhsfamily1[0]].id != d->m->pedigrees[1][fhsfamily3[0]].id);
    assert(d->m->pedigrees[1][fhsfamily2[0]].id != d->m->pedigrees[1][fhsfamily3[0]].id);

    fhs = combine_groups(d,3,halfsibfamilies);
    assert(d->n_groups == 2);

    // Can create a group from indexes
    size_t splitters[4] = {0, 2, 5, 6};
    GroupNum g0b = make_group_from(d, 4, splitters);
    assert(d->n_groups == 3);
    assert(get_group_size(d,g0) == 3 && get_group_size(d,g0b) == 4);
    size_t g0bi[4];
    assert(get_group_indexes(d,g0b,4,g0bi) == 4);
    assert(g0bi[0] == 0 && g0bi[1] == 2 && g0bi[2] == 5 && g0bi[3] == 6);

    GroupNum n6 = make_group_from(d, 1, splitters+3);
    assert(d->n_groups == 4);
    GroupNum combiners[2] = {fhs, n6};
    fhs = combine_groups(d, 2, combiners);
    assert(d->n_groups == 3);
    delete_group(d, fhs);
    assert(d->n_groups == 2);

    size_t g0members[6] = {0,1,2,3,4,5};
    g0 = make_group_from(d,6,g0members);
    assert(d->n_groups >= 1);

    printf("...Specified group splitting runs correctly\n");

    return g0;
}


GroupNum test_grouping(SimData *d, GroupNum g0) {
    g0 = test_labels(d, g0);
    g0 = test_random_splits(d, g0);
    g0 = test_specific_splits(d, g0);

    return g0;
}

int test_effect_calculators(SimData *d, GroupNum g0) {
    DecimalMatrix dec = calculate_group_bvs(d, g0, (EffectID){.id=1});

	assert(dec.rows == 1);
	assert(dec.cols == 6);
    assert(fabs(dec.matrix[0][0] - 1.4) < TOL);
    assert(fabs(dec.matrix[0][1] - 1.4) < TOL);
    assert(fabs(dec.matrix[0][2] - 1.6) < TOL);
    assert(fabs(dec.matrix[0][3] - (-0.1)) < TOL);
    assert(fabs(dec.matrix[0][4] - 0.6) < TOL);
    assert(fabs(dec.matrix[0][5] - (-0.3)) < TOL);

    delete_dmatrix(&dec);
    // and with the second set of effects:
    DecimalMatrix dec2 = calculate_group_bvs(d, g0, (EffectID){.id=2});

    assert(dec2.rows == 1);
    assert(dec2.cols == 6);
    assert(fabs(dec2.matrix[0][0] - 0) < TOL);
    assert(fabs(dec2.matrix[0][1] - 0) < TOL);
    assert(fabs(dec2.matrix[0][2] - 0) < TOL);
    assert(fabs(dec2.matrix[0][3] - 1) < TOL);
    assert(fabs(dec2.matrix[0][4] - 0) < TOL);
    assert(fabs(dec2.matrix[0][5] - 1) < TOL);
    delete_dmatrix(&dec2);

    FILE* fp = fopen("a-test-eff3.txt","w");
    fwrite(HELPER_EFF2, sizeof(char), strlen(HELPER_EFF2), fp);
    fclose(fp);
    EffectID e3 = load_effects(d,"a-test-eff3.txt");
    DecimalMatrix dec3 = calculate_group_bvs(d, g0, e3);
    assert(fabs(dec3.matrix[0][0] - 0) < TOL);
    assert(fabs(dec3.matrix[0][1] - 0) < TOL);
    assert(fabs(dec3.matrix[0][2] - 0) < TOL);
    assert(fabs(dec3.matrix[0][3] - 1) < TOL);
    assert(fabs(dec3.matrix[0][4] - 0) < TOL);
    assert(fabs(dec3.matrix[0][5] - 1) < TOL);
    delete_eff_set(d,e3);
    delete_dmatrix(&dec3);

    remove("a-test-eff3.txt");

    printf("...GEBVs calculated correctly\n");

	return 0;
}

int test_optimal_calculators(SimData *d, GroupNum g0) {
    EffectID eff_set = {.id = 1};
    char* ig = calculate_optimal_haplotype(d, eff_set);
    assert(ig[0] == 'T');
    assert(ig[1] == 'A');
    assert(ig[2] == 'A');
    free(ig);

    double optimal = calculate_optimal_bv(d, eff_set);
    assert(fabs(optimal - 1.8) < TOL);
    assert(fabs(calculate_optimal_bv(d, (EffectID){.id=2}) - 2) < TOL);

    double unoptimal = calculate_minimal_bv(d, eff_set);
    assert(fabs(unoptimal + 2.8) < TOL);
    //assert(fabs(calculate_minimal_bv(d, 1) - 0) < TOL); // this function actually doesn't work when not all alleles have marker effects for somewhere in the genome

    char* founderhaplo = calculate_optimal_possible_haplotype(d, g0, eff_set);
    assert(founderhaplo[0] == 'T');
    assert(founderhaplo[1] == 'A');
    assert(founderhaplo[2] == 'A');
    free(founderhaplo);
    founderhaplo = calculate_optimal_possible_haplotype(d, g0, (EffectID){.id=2});
    assert(founderhaplo[0] == 'A');
    free(founderhaplo);

    double founderoptimal = calculate_optimal_possible_bv(d, g0, eff_set);
    assert(fabs(founderoptimal - 1.8) < TOL);
    assert(fabs(calculate_optimal_possible_bv(d, g0, (EffectID){.id=2}) - 2) < TOL);

    size_t factorout[2] = {4,5};
    GroupNum g0partial = make_group_from(d,2, factorout);
    char* founderhaplo2 = calculate_optimal_possible_haplotype(d, g0partial, eff_set);
    assert(founderhaplo2[0] == 'T');
    assert(founderhaplo2[1] == 'A');
    assert(founderhaplo2[2] == 'T');
    free(founderhaplo2);

    double founderoptimal2 = calculate_optimal_possible_bv(d, g0partial, eff_set);
    assert(fabs(founderoptimal2 - 1.4) < TOL);

    GroupNum recombine[2] = {g0,g0partial};
    GroupNum newg0 = combine_groups(d,2,recombine);
    assert(g0.num == newg0.num); // for validity of further tests

    printf("...Optimal genotype and GEBV calculated correctly\n");

    return 0;
}

int test_crossing(SimData *d, GroupNum g0) {
    GroupNum gall = test_crossing_unidirectional(d, g0);

    test_crossing_randomly(d, g0);

	FILE* fp;
	if ((fp = fopen("a-test-plan.txt", "w")) == NULL) {
		fprintf(stderr, "Failed to create file.\n");
		exit(1);
	}
	fwrite(HELPER_PLAN, sizeof(char), strlen(HELPER_PLAN), fp);
	fclose(fp);
    GroupNum gfile = test_crossing_from_file(d, "a-test-plan.txt");
	remove("a-test-plan.txt");

    GroupNum gselfed = test_crossing_selfing(d, g0);

    assert(gselfed.num != g0.num && gfile.num != gall.num && gfile.num != g0.num && gfile.num != gselfed.num);
	printf("...group number allocations are correct\n");

	return 0;
}

GroupNum test_crossing_unidirectional(SimData *d, GroupNum g0) {
    GenOptions g = {.will_name_offspring=GSC_TRUE,
                    .offspring_name_prefix="F1",
                    .family_size=1,.will_track_pedigree=GSC_TRUE,
                    .will_allocate_ids=GSC_TRUE,
                    .filename_prefix="atestF1",
                    .will_save_pedigree_to_file=GSC_TRUE,
                    .will_save_bvs_to_file=0,
                    .will_save_alleles_to_file=GSC_TRUE,
                    .will_save_to_simdata=GSC_TRUE};
	//AlleleMatrix* a = make_all_unidirectional_crosses(&sd, 0, g);
	//sd.m->next_gen = a;
    GroupNum g1 = make_all_unidirectional_crosses(d, g0, g);

    assert(g1.num != g0.num);
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
    assert(d->m->pedigrees[0][6].id == 1 && d->m->pedigrees[1][6].id == 2); //assumes you're doing the top triangle of crosses
    assert(d->m->pedigrees[0][7].id == 1 && d->m->pedigrees[1][7].id == 3);
    assert(d->m->pedigrees[0][8].id == 1 && d->m->pedigrees[1][8].id == 4);
    assert(d->m->pedigrees[0][9].id == 1 && d->m->pedigrees[1][9].id == 5);
    assert(d->m->pedigrees[0][10].id == 1 && d->m->pedigrees[1][10].id == 6);
    assert(d->m->pedigrees[0][11].id == 2 && d->m->pedigrees[1][11].id == 3);
    assert(d->m->pedigrees[0][12].id == 2 && d->m->pedigrees[1][12].id == 4);
    assert(d->m->pedigrees[0][13].id == 2 && d->m->pedigrees[1][13].id == 5);
    assert(d->m->pedigrees[0][14].id == 2 && d->m->pedigrees[1][14].id == 6);
    assert(d->m->pedigrees[0][15].id == 3 && d->m->pedigrees[1][15].id == 4);
    assert(d->m->pedigrees[0][16].id == 3 && d->m->pedigrees[1][16].id == 5);
    assert(d->m->pedigrees[0][17].id == 3 && d->m->pedigrees[1][17].id == 6);
    assert(d->m->pedigrees[0][18].id == 4 && d->m->pedigrees[1][18].id == 5);
    assert(d->m->pedigrees[0][19].id == 4 && d->m->pedigrees[1][19].id == 6);
    assert(d->m->pedigrees[0][20].id == 5 && d->m->pedigrees[1][20].id == 6);
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

GroupNum test_crossing_from_file(SimData *d, char* fname) {
	// check we can load a plan from a file.
	GenOptions g2 = BASIC_OPT;
	g2.will_track_pedigree = GSC_TRUE;
	g2.family_size = 2;
	//g2.will_save_pedigree_to_file = GSC_TRUE;
	//2.filename_prefix = "atest-dc";
    GroupNum bp = make_double_crosses_from_file(d, fname, g2);
    assert(d->m->pedigrees[0][21].id == d->m->ids[6].id && d->m->pedigrees[1][21].id == 17);
    assert(d->m->pedigrees[0][23].id == d->m->ids[7].id && d->m->pedigrees[1][23].id == 21);
    assert(d->m->pedigrees[0][25].id == d->m->ids[20].id && d->m->pedigrees[1][25].id == 9);
    assert(d->m->pedigrees[0][22].id == 7 && d->m->pedigrees[1][22].id == 17);
    assert(d->m->pedigrees[0][24].id == 8 && d->m->pedigrees[1][24].id == 21);
    assert(d->m->pedigrees[0][26].id == 21 && d->m->pedigrees[1][26].id == 9);
	assert(d->m->n_genotypes == 27);
	assert(d->m->n_markers == 3);
	printf("...crossed combinations from file correctly\n");

	return bp;
}

GroupNum test_crossing_selfing(SimData *d, GroupNum g1) {
    int oldsize = d->m->n_genotypes;
    GenOptions opt = BASIC_OPT;
    opt.will_track_pedigree = GSC_TRUE;
	float h1 = calculate_heterozygosity(d, g1);
    GroupNum g1selfed = self_n_times(d, 5, g1, opt);
	float h2 = calculate_heterozygosity(d,  g1selfed);

    assert(g1selfed.num != g1.num);
	//printf("Heterozygousity reduction from selfing: %f %f\n", h2, h1);
	assert(h1 - h2 > 0);
    assert(d->m->n_genotypes == oldsize + 6);
	assert(d->m->n_markers == 3);
    assert(d->m->groups[oldsize].num == g1selfed.num && d->m->groups[oldsize + 5].num == g1selfed.num
            && d->m->groups[oldsize + 6].num != g1selfed.num);
    assert(d->m->pedigrees[0][oldsize + 0].id == d->m->ids[0].id && d->m->pedigrees[1][oldsize + 0].id == d->m->ids[0].id);
    assert(d->m->pedigrees[0][oldsize + 4].id == d->m->ids[4].id && d->m->pedigrees[1][oldsize + 4].id == d->m->ids[4].id);
	printf("...selfing function correctly reduced heterozygosity by %f%%\n", (h2-h1)*100);

    // test doubled haploids
    GroupNum g1dhap = make_doubled_haploids(d, g1, opt);
    assert(g1dhap.num != g1.num);
    assert(d->m->n_genotypes == oldsize + 2*6);
    assert(d->m->n_markers == 3);
    assert(calculate_heterozygosity(d,  g1dhap) == 0);
    assert(d->m->groups[oldsize + 6].num == g1dhap.num && d->m->groups[oldsize + 11].num == g1dhap.num
            && d->m->groups[oldsize + 12].num != g1dhap.num);
    assert(d->m->pedigrees[0][oldsize + 6].id == d->m->ids[0].id && d->m->pedigrees[1][oldsize + 6].id == d->m->ids[0].id);
    assert(d->m->pedigrees[0][oldsize + 10].id == d->m->ids[4].id && d->m->pedigrees[1][oldsize + 10].id == d->m->ids[4].id);
    //printf("%s -> %s\n", d->m->alleles[0], d->m->alleles[oldsize + 6]);
    assert(strncmp(d->m->alleles[oldsize + 6],d->m->alleles[0],sizeof(char)*6) == 0);
    assert(strncmp(d->m->alleles[oldsize + 11],"AAAATT",sizeof(char)*6) == 0 ||
            strncmp(d->m->alleles[oldsize + 11],"TTAATT",sizeof(char)*6) == 0);
    printf("...doubled haploid function works correctly\n");

    delete_group(d, g1dhap);

    // test cloning
    opt.will_allocate_ids = GSC_FALSE;
    GroupNum g1clones = make_clones(d, g1, GSC_TRUE, opt);
    assert(g1clones.num != g1.num);
    assert(d->m->n_genotypes == oldsize + 2*6);
    assert(d->m->n_markers == 3);
    assert(fabsf(calculate_heterozygosity(d,  g1clones) - h1) < TOL);
    assert(d->m->groups[oldsize + 6].num == g1clones.num && d->m->groups[oldsize + 11].num == g1clones.num
            && d->m->groups[oldsize + 12].num != g1clones.num);
    assert(d->m->pedigrees[0][oldsize + 6].id == d->m->ids[0].id && d->m->pedigrees[1][oldsize + 6].id == d->m->ids[0].id);
    assert(d->m->pedigrees[0][oldsize + 10].id == d->m->ids[4].id && d->m->pedigrees[1][oldsize + 10].id == d->m->ids[4].id);
    assert(strncmp(d->m->alleles[oldsize + 6],d->m->alleles[0],sizeof(char)*6) == 0);
    assert(strncmp(d->m->alleles[oldsize + 11],"ATAATT",sizeof(char)*6) == 0);
    assert(strcmp(d->m->names[1],d->m->names[oldsize + 7]) == 0);
    assert(strcmp(d->m->names[3],d->m->names[oldsize + 9]) == 0);
    printf("...cloning function works correctly\n");

    opt.will_allocate_ids = GSC_TRUE;
    GroupNum gap = make_random_crosses(d , g1, 5, 0, opt);
    GroupNum newparents[2];
    newparents[0] = g1clones;
    newparents[1] = make_random_crosses(d , g1, 5, 0, opt);
    GroupNum cloneparents2 = combine_groups(d, 2, newparents);
    GroupNum g2clones = make_clones(d, cloneparents2, GSC_FALSE, opt);

    delete_group(d, cloneparents2);
    delete_group(d, gap);
    delete_group(d, g2clones);

	return g1selfed;
}

int test_crossing_randomly(SimData *d, GroupNum g1) {
    // we test it correctly does its crossing randomly (requiring a bit of human input)
    // we do not test that the genes were correctly crossed
    // created 7 Apr 2022

    GenOptions gopt = BASIC_OPT;
    gopt.will_track_pedigree = GSC_TRUE;
    // Test random crossing seems about right
    GroupNum g2 = make_random_crosses( d , g1, 4, 0, gopt);
    size_t g2ixs[4];
    assert(get_group_indexes(d, g2, -1, g2ixs) == 4);

    assert(get_group_size(d, g2) == 4);
    int g2minid = d->m->ids[g2ixs[0]].id;
    int g2maxid = d->m->ids[g2ixs[3]].id;
    fprintf(stdout, "Should be random parents: %d, %d, %d, %d\n",
            d->m->pedigrees[0][g2ixs[0]].id, d->m->pedigrees[0][g2ixs[1]].id,
            d->m->pedigrees[0][g2ixs[2]].id, d->m->pedigrees[0][g2ixs[3]].id);
    fprintf(stdout, "Should be random parents: %d, %d, %d, %d\n\n",
            d->m->pedigrees[1][g2ixs[0]].id, d->m->pedigrees[1][g2ixs[1]].id,
            d->m->pedigrees[1][g2ixs[2]].id, d->m->pedigrees[1][g2ixs[3]].id);

    printf("...crossed randomly within a group\n");

    // Test random crossing between two groups seems about right.
    gopt.family_size = 2;
    GroupNum g3 = make_random_crosses_between( d, g1, g2, 3, 0, 0, gopt);

    assert(get_group_size(d, g3) == 6);
    assert(d->m->pedigrees[0][g2ixs[3] + 1].id == d->m->pedigrees[0][g2ixs[3] + 2].id); // family size works
    assert(d->m->pedigrees[1][g2ixs[3] + 1].id == d->m->pedigrees[1][g2ixs[3] + 2].id); // family size works
    assert(d->m->pedigrees[1][g2ixs[3] + 3].id >= g2minid && d->m->pedigrees[1][g2ixs[3] + 3].id <= g2maxid ); //right parent groupings
    assert(d->m->pedigrees[0][g2ixs[3] + 3].id < g2minid); //right parent groupings
    fprintf(stdout, "Should be random parents: %d, %d, %d\n",
            d->m->pedigrees[0][g2ixs[3] + 1].id, d->m->pedigrees[0][g2ixs[3] + 3].id,
            d->m->pedigrees[0][g2ixs[3] + 5].id);
    fprintf(stdout, "Should be random parents: %d, %d, %d\n",
            d->m->pedigrees[1][g2ixs[3] + 1].id, d->m->pedigrees[1][g2ixs[3] + 3].id,
            d->m->pedigrees[1][g2ixs[3] + 5].id);

    delete_group(d, g3);

    gopt.family_size = 1;
    GroupNum tempg = make_group_from( d, 1, g2ixs + 1);
    assert(get_group_size(d, tempg) == 1);
    int g4sizetobe = 6;
    GroupNum g4 = make_random_crosses_between( d, g1, tempg, g4sizetobe, 1, 0, gopt );
    assert(get_group_size(d, g4) == g4sizetobe);
    GroupNum combine[2] = {tempg, g2};
    g2 = combine_groups( d, 2, combine );

    assert(get_group_size(d, g4) == g4sizetobe);
    for (int i = 1; i < g4sizetobe; ++i) {
        for (int j = i + 1; j <= g4sizetobe; ++j) {
            assert(d->m->pedigrees[1][g2ixs[3] + i].id == d->m->pedigrees[1][g2ixs[3] + j].id);// right parent repetition
            assert(d->m->pedigrees[0][g2ixs[3] + i].id != d->m->pedigrees[0][g2ixs[3] + j].id); // should be different first parents
        }
    }
    delete_group( d, g4 );

    //int g5 = make_random_crosses_between( d, g1, g2, 2000, 1, 1, gopt); // expect error
    //delete_group( d, g5 );
    printf("...crossed randomly between two groups\n");
    delete_group( d, g2 );
    return 0;
}

int test_deletors(SimData *d, GroupNum g0) {
    GroupNum groups1[1000];
    assert(d->n_groups < 1000);
    int ngroups1 = get_existing_groups(d, groups1);
    assert(d->n_groups == ngroups1);

    GroupNum groups1b[1000];
    size_t groupcounts1b[1000];
    int ngroups1b = get_existing_group_counts(d, groups1b, groupcounts1b);
    assert(d->n_groups == ngroups1b);
    assert(ngroups1b == ngroups1);
    for (int i = 0; i < ngroups1; ++i) {
        assert(groups1[i].num == groups1b[i].num);
    }

	delete_group(d, g0);
    GroupNum groups2[1000];
    int ngroups2 = get_existing_groups(d, groups2);
    assert(d->n_groups == ngroups2);

	assert(ngroups1 - ngroups2 == 1);
	for (int i = 0; i < ngroups1; ++i) {
        if (i == ngroups1 - 1 || groups1[i].num != groups2[i].num) {
            assert(groups1[i].num == g0.num);
            assert(groups1b[i].num == g0.num);
            assert(groupcounts1b[i] == 6);
			break;
		}
	}
	printf("...group of genotypes cleared correctly\n");

    delete_eff_set(d, (EffectID){.id=1});
    assert(d->n_eff_sets == 1);
    assert(d->e != NULL);
    assert(d->e[0].effects.rows == 1); // check identity
    assert((d->e[0].effects.matrix[0][0] - 1) < TOL); // check identity

    printf("...marker effects set cleared correctly\n");

	delete_simdata(d);
	printf("...SimData cleared correctly\n");

	return 0;
}

int test_block_generator(SimData *d) {
    MarkerBlocks b = create_evenlength_blocks_each_chr(d, 2);

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

    b = create_evenlength_blocks_each_chr(d, 4);

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

int test_iterators(SimData* d, GroupNum gp) {
    GenOptions g = {
        .will_name_offspring = GSC_FALSE,
        .offspring_name_prefix = NULL,
        .family_size = 1,
        .will_track_pedigree = GSC_TRUE,
        .will_allocate_ids = GSC_FALSE,
        .filename_prefix = NULL,
        .will_save_pedigree_to_file = GSC_FALSE,
        .will_save_bvs_to_file = -1,
        .will_save_alleles_to_file = GSC_FALSE,
        .will_save_to_simdata = GSC_TRUE
    };
    int combos[2][3];
    combos[0][0] = 0; combos[1][0] = 0;
    combos[0][1] = 1; combos[1][1] = 2;
    combos[0][2] = 1; combos[1][2] = 5;
    GroupNum f1 = make_targeted_crosses(d, 3, combos[0], combos[1], g);
    const LabelID label1 = create_new_label(d, 0);
    change_label_by_amount(d,f1,label1,1);

    // Bidirectional global iterator; starting forwards
    BidirectionalIterator it1 = create_bidirectional_iter(d, NO_GROUP);
    int i = 0;
    for (; i < 6; ++i) {
        GenoLocation gl = next_forwards(&it1);
        assert(IS_VALID_LOCATION(gl));
        assert(get_id(gl).id == i+1);
        assert(get_label_value(gl,gsc_get_index_of_label(d,label1)) == 0);
        assert(get_group(gl).num == gp.num);
    }
    for (; i < 6+3; ++i) {
        GenoLocation gl = next_forwards(&it1);
        assert(IS_VALID_LOCATION(gl));
        assert(get_id(gl).id == 0);
        assert(get_label_value(gl,gsc_get_index_of_label(d,label1)) == 1);
        assert(get_group(gl).num == f1.num);
    }
    assert(!IS_VALID_LOCATION(next_forwards(&it1)));
    for (int j = 0; j < 3-1; ++j) {
        assert(get_id(next_backwards(&it1)).id == 0);
    }
    assert(get_id(next_backwards(&it1)).id == 6);

    delete_bidirectional_iter(&it1);

    // Bidirectional group iterator; starting backwards
    BidirectionalIterator it2 = create_bidirectional_iter(d, gp);
    assert(get_id(next_backwards(&it2)).id == 6);
    assert(!IS_VALID_LOCATION(next_forwards(&it2)));
    assert(get_id(next_backwards(&it2)).id == 5);
    assert(get_id(next_backwards(&it2)).id == 4);
    assert(get_id(next_backwards(&it2)).id == 3);
    assert(get_id(next_backwards(&it2)).id == 2);
    assert(get_id(next_forwards(&it2)).id == 3);
    assert(get_id(next_backwards(&it2)).id == 2);
    assert(get_id(next_backwards(&it2)).id == 1);
    assert(!IS_VALID_LOCATION(next_backwards(&it2)));

    // RandomAccess global iterator
    RandomAccessIterator it3 = create_randomaccess_iter(d, NO_GROUP);
    assert(get_group(next_get_nth(&it3, 6)).num == f1.num);
    assert(get_group(next_get_nth(&it3, 2)).num == gp.num);
    assert(get_id(next_get_nth(&it3, 2)).id == 3);
    assert(get_id(next_get_nth(&it3, 5)).id == 6);
    assert(!IS_VALID_LOCATION(next_get_nth(&it3,-1)));
    assert(!IS_VALID_LOCATION(next_get_nth(&it3,10)));

    delete_randomaccess_iter(&it3);

    // RandomAccess group iterator
    const char* newnames[4] = {"f1a", "f1b", "f1c", "f1d"};
    change_names_to_values(d,f1,0,4,newnames);
    RandomAccessIterator it4 = create_randomaccess_iter(d, f1);
    assert(get_group(next_get_nth(&it4, 0)).num == f1.num);
    assert(strncmp(get_name(next_get_nth(&it4, 0)),"f1a",5)==0);
    assert(get_first_parent(next_get_nth(&it4, 0)).id == 1);
    assert(get_first_parent(next_get_nth(&it4, 2)).id == 2);
    assert(get_second_parent(next_get_nth(&it4, 1)).id == 3);

    delete_randomaccess_iter(&it4);

    // Empty group iterator
    RandomAccessIterator it5 = create_randomaccess_iter(d,(GroupNum){.num=4});
    assert(!IS_VALID_LOCATION(next_get_nth(&it5,0)));

    delete_group(d, f1);
    delete_label(d, label1);

    printf("...iterators work correctly.\n");

    return 0;
}

int test_getters(SimData* d, GroupNum gp) {
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

    PedigreeID ids[10];
    assert(get_group_ids(d, gp, -1, ids) == 6);
    for (int i = 0; i < 6; ++i) {
        assert(ids[i].id == i+1);
    }

    size_t indexes[6];
    assert(get_group_indexes(d, gp, 6, indexes) == 6);
    for (int i = 0; i < 6; ++i) {
        assert(indexes[i] == i);
    }

    double bvs[6];
    assert(get_group_bvs(d, gp, (EffectID){.id=1}, 6, bvs) == 6);
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
        .will_name_offspring = GSC_FALSE,
        .offspring_name_prefix = NULL,
        .family_size = 1,
        .will_track_pedigree = GSC_TRUE,
        .will_allocate_ids = GSC_FALSE,
        .filename_prefix = NULL,
        .will_save_pedigree_to_file = GSC_FALSE,
        .will_save_bvs_to_file = -1,
        .will_save_alleles_to_file = GSC_FALSE,
        .will_save_to_simdata = GSC_TRUE
    };
    GroupNum f1 = make_targeted_crosses(d,3,combos[0],combos[1],g);

    PedigreeID p1s[3]; PedigreeID p2s[6];
    assert(get_group_parent_ids(d, f1, GSC_UNINIT, 1, p1s) == 3);
    assert(get_group_parent_ids(d, f1, 3, 2, p2s) == 3);
    assert(p1s[0].id == 1);
    assert(p1s[1].id == 2);
    assert(p1s[2].id == 2);
    assert(p2s[0].id == 1);
    assert(p2s[1].id == 3);
    assert(p2s[2].id == 6);

    char* p1ns[3]; char* p2ns[3];
    assert(get_group_parent_names(d, f1, 3, 1, p1ns) == 3);
    assert(get_group_parent_names(d, f1, GSC_UNINIT, 2, p2ns) == 3);
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

int test_data_access(SimData* d, GroupNum gp) {
    test_iterators(d,gp);

    test_getters(d,gp);

    return 0;
}


int subfunctiontest_meiosis(int seed) {
    SimData* d = create_empty_simdata(seed);
    GroupNum g0 = just_load(d);
    // As of yet these do not check anything about crossovers. Nothing about whether it's
    // following the promised probabilities, only that the results, in isolation, look plausible.

    char gamete[10] = { 0 };
    char expected[10];
    char expected2[10];

    RandomAccessIterator it = create_randomaccess_iter(d,g0);
    GenoLocation homozygote = next_get_nth(&it,0);
    GenoLocation heterozygote = next_get_nth(&it,5);
    const int GENOMESIZE = 6;
    assert(IS_VALID_LOCATION(homozygote));
    assert(strcmp(get_name(homozygote), "G01") == 0);
    assert(strncmp(get_alleles(homozygote),"TTAATT",sizeof(char)*GENOMESIZE) == 0);
    assert(IS_VALID_LOCATION(heterozygote));
    assert(strcmp(get_name(heterozygote), "G06") == 0);
    assert(strncmp(get_alleles(heterozygote),"ATAATT",sizeof(char)*GENOMESIZE) == 0);


    // test of generate_gamete
    // with a totally homozygous parent eg G01, gives the expected result eg (TTA)
    strcpy(expected,"T?A?T?");
    generate_gamete(d,get_alleles(homozygote),gamete);
    for (int i = 0; i < GENOMESIZE; ++i) {
        if (expected[i] != '?') {
            assert(expected[i] == gamete[i]);
        } else {
            assert(gamete[i] == 0); // did not overwrite spaces that don't belong to it.
        }
    }
    assert(gamete[GENOMESIZE] == 0);

    // with a heterozygous parent, eg G06, can give both expected results eg (AAT, TAT)
    strcpy(expected,"A?A?T?");
    strcpy(expected2,"T?A?T?");
    memset(gamete,0,sizeof(char)*10);
    int e1count = 0; int e2count = 0;
    while (1) {
        generate_gamete(d,get_alleles(heterozygote),gamete);
        if (gamete[0] == expected[0] && gamete[2] == expected[2] && gamete[4] == expected[4]) {
            ++e1count;
        } else if (gamete[0] == expected2[0] && gamete[2] == expected2[2] && gamete[4] == expected2[4]) {
            ++e2count;
        } else {
            assert(GSC_FALSE); // immediate fail because gamete did not match any of our expectations.
        }

        if (e1count > 0 && e2count > 0) {
            break;
        }

        assert(e1count < 100 && e2count < 100); // let's just fail here if we never get both outcomes.
        // does have a minuscule (0.5^100) chance of failing even on a valid generate_gamete function.
    }
    assert(gamete[GENOMESIZE] == 0);

    // test of generate_doubled_haploid
    // with a totally homozygous parent
    memset(gamete,0,sizeof(char)*10);
    strcpy(expected,"TTAATT");
    generate_doubled_haploid(d,get_alleles(homozygote),gamete);
    assert(strncmp(gamete,expected,sizeof(char)*GENOMESIZE) == 0);
    assert(gamete[GENOMESIZE] == 0);

    // With a heterozygous parent, can give both expected results.
    strcpy(expected,"AAAATT");
    strcpy(expected2,"TTAATT");
    memset(gamete,0,sizeof(char)*10);
    e1count = 0; e2count = 0;
    while (1) {
        generate_doubled_haploid(d,get_alleles(heterozygote),gamete);
        if (strncmp(gamete,expected,sizeof(char)*GENOMESIZE) == 0) {
            ++e1count;
        } else if (strncmp(gamete,expected2,sizeof(char)*GENOMESIZE) == 0)  {
            ++e2count;
        } else {
            assert(GSC_FALSE); // immediate fail because gamete did not match any of our expectations.
        }

        if (e1count > 0 && e2count > 0) {
            break;
        }

        assert(e1count < 100 && e2count < 100); // let's just fail here if we never get both outcomes.
        // does have a minuscule (0.5^100) chance of failing even on a valid generate_doubled_haploid function.
    }
    assert(gamete[GENOMESIZE] == 0);

    // test of generate_clone
    // With a totally homozygous parent
    memset(gamete,0,sizeof(char)*10);
    strcpy(expected,"TTAATT");
    generate_clone(d,get_alleles(homozygote),gamete);
    assert(strncmp(gamete,expected,sizeof(char)*GENOMESIZE) == 0);
    assert(gamete[GENOMESIZE] == 0);

    // With a heterozygous parent.
    strcpy(expected,"ATAATT");
    memset(gamete,0,sizeof(char)*10);
    generate_clone(d,get_alleles(heterozygote),gamete);
    assert(strncmp(gamete,expected,sizeof(char)*GENOMESIZE) == 0);
    assert(gamete[GENOMESIZE] == 0);

    delete_simdata(d);
    printf("...(subfunction) gamete generators work correctly.\n");

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

int compareFileToString(char* filename, const char* target) {
    FILE* fp = fopen(filename, "r");
    int i = 0;

    char c1, c2;

    do {
        c1 = fgetc(fp);
        c2 = target[i];

        if (c1 != c2) break;
        ++i;

    } while (c1 != EOF && c2 != '\0');

    if (c1 == EOF && c2 == '\0') return 0;
    else return -1;
}

int compareRepeatingFileToTable(char* filename, unsigned int expectedNRows, const char* header, const char* body) {
    FILE* fp = fopen(filename, "r");
    int i = 0;
    unsigned int row = 0;
    char processingheader = (header == NULL) ? GSC_FALSE : GSC_TRUE;

    char c1, c2;

    do {
        c1 = fgetc(fp);
        if (processingheader) {
            c2 = header[i];
        } else {
            c2 = body[i];
        }

        if (c2 == '\0') {
            if (c1 == '\n') {
                i = 0;
                if (processingheader) {
                    processingheader = GSC_FALSE;
                } else {
                    ++row;
                }
            } else {
                break;
            }
        } else if (c1 != c2) {
            break;
        } else {
            ++i;
        }

    } while (c1 != EOF && row < expectedNRows);

    if (row == expectedNRows) {
        c1 = fgetc(fp);
    }

    if (c1 == EOF && i == 0 && row == expectedNRows) return 0;
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
    GroupNum g0 = test_loaders(d);
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
    subfunctiontest_meiosis(randomSeed);
	printf("\t\t-> Crossing functions all clear\n");

	//test blocking
	printf("\nNow testing blocking functions:\n");
	test_block_generator(d);
	printf("\t\t-> Blocking functions all clear\n");

	//test file savers
	printf("\nNow testing saver functions:\n");
    test_savers(randomSeed);
    printf("\t\t-> Saver functions all clear\n");

	// test SimData deletors.
	printf("\nNow testing deletor functions:\n");
	test_deletors(d, g0);
	printf("\t\t-> Deletor functions all clear\n");

    printf("\n------- All tests passed. -------\n");


	//testing new grouping functions
    d = create_empty_simdata(randomSeed);
    struct GroupAndEffectSet gande =
            load_all_data(d, "./gt_parents_mr2_50-trimto-5000.txt",
			 "./genetic-map_5112-trimto5000.txt",
             "./qtl_mr2.eff-processed.txt");
    g0 = gande.group;

    /*g0 = load_all_data(d, "./gt_parents_mr2_3000x30000.txt",
                          "./genetic-map_huge30000.txt",
                          "./eff_huge30000.txt");*/

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
	GenOptions gens =  {.will_name_offspring=GSC_FALSE, .offspring_name_prefix=NULL, .family_size=243,
		.will_track_pedigree=GSC_TRUE, .will_allocate_ids=GSC_TRUE,
		.filename_prefix="testcross", .will_save_pedigree_to_file=GSC_FALSE,
        .will_save_bvs_to_file=-1, .will_save_alleles_to_file=GSC_FALSE,
		.will_save_to_simdata=GSC_TRUE};
    int g1 = make_targeted_crosses(d,nc,crosses[0],crosses[1],gens);

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
    split_by_probabilities(d,g0,3,probs,g0bs);*/

    /*// there's too little for my buckets~
    int buckets[3];
    buckets[0] = 46;
    buckets[1] = 13;
    split_into_buckets(d,g0,3,buckets,g0bs);*/

    //int a = get_group_size(d, g0);
    //a = get_group_size(d,g0b);
    /*// all genotypes transferred check
    int a = 0;
    for (int i = 0; i < 5; ++i) {
        a += get_group_size(d, g0bs[i]);
    }*/

    // This should not segfault. Or do anything else weird.
    GenOptions gens = {.will_name_offspring=GSC_TRUE, .offspring_name_prefix="cr", .family_size=243,
        .will_track_pedigree=GSC_TRUE, .will_allocate_ids=GSC_TRUE,
        .filename_prefix="testcross", .will_save_pedigree_to_file=GSC_FALSE,
        .will_save_bvs_to_file=-1, .will_save_alleles_to_file=GSC_FALSE,
        .will_save_to_simdata=GSC_TRUE};
    int ncrosses = 10;
    size_t ixs[3000];
    assert(get_group_size(d,g0) == 50);

        GroupNum g000 = make_random_crosses(d, g0, ncrosses, 0, gens);
    assert(d->n_groups == 2);
    assert(g000.num == 2);
    assert(get_group_size(d,g000)==gens.family_size*ncrosses);
    get_group_indexes(d,g000,-1,ixs);
    assert(ixs[0] == 50);
    assert(ixs[gens.family_size*ncrosses-1] == gens.family_size*ncrosses + 49);

        g000 = make_random_crosses(d, g0, ncrosses, 0, gens);
    assert(d->n_groups == 3);
    assert(g000.num == 3);
    assert(get_group_size(d,g000)==gens.family_size*ncrosses);
    get_group_indexes(d,g000,-1,ixs);
    assert(ixs[0] == gens.family_size*ncrosses + 50);
    assert(ixs[gens.family_size*ncrosses-1] == gens.family_size*ncrosses*2 + 49);

        g000 = make_random_crosses(d, g0, ncrosses, 0, gens);
    assert(d->n_groups == 4);
    assert(g000.num == 4);
    assert(get_group_size(d,g000)==gens.family_size*ncrosses);
    get_group_indexes(d,g000,-1,ixs);
    assert(ixs[0] == gens.family_size*ncrosses*2 + 50);
    assert(ixs[gens.family_size*ncrosses-1] == gens.family_size*ncrosses*3 + 49);

        g000 = make_random_crosses(d, g0, ncrosses, 0, gens);
    assert(d->n_groups == 5);
    assert(g000.num == 5);
    BidirectionalIterator it000 = create_bidirectional_iter(d,NO_GROUP);
    GenoLocation last = next_backwards(&it000);
    assert(get_group(last).num==5);
    int totalngenos = 0;
    AlleleMatrix* cam = d->m;
    do {
        assert(cam->n_genotypes == CONTIG_WIDTH || cam->n_genotypes == (gens.family_size*ncrosses*4 + 50) % CONTIG_WIDTH);
        totalngenos += cam->n_genotypes;
    } while ((cam = cam->next) != NULL);
    assert(totalngenos == gens.family_size*ncrosses*4 + 50);
    assert(get_group_size(d,(GroupNum){.num=4})==gens.family_size*ncrosses);
    assert(get_group_size(d,(GroupNum){.num=3})==gens.family_size*ncrosses);
    assert(get_group_size(d,(GroupNum){.num=2})==gens.family_size*ncrosses);
    assert(get_group_size(d,(GroupNum){.num=1})==50);
    assert(get_group_size(d,(GroupNum){.num=0})==0);
    assert(get_group_size(d,(GroupNum){.num=6})==0);
    get_group_indexes(d,(GroupNum){.num=4},-1,ixs);
    assert(ixs[0] == gens.family_size*ncrosses*2 + 50);
    assert(ixs[gens.family_size*ncrosses-1] == gens.family_size*ncrosses*3 + 49);
    get_group_indexes(d,g000,-1,ixs); // !! indexes start at 7400 instead of the expected 7250
    assert(ixs[0] == gens.family_size*ncrosses*3 + 50);
    assert(ixs[gens.family_size*ncrosses-1] == gens.family_size*ncrosses*4 + 49);
    assert(get_group_size(d,g000)==gens.family_size*ncrosses);

        g000 = make_random_crosses(d, g0, 10, 0, gens);
    assert(d->n_groups == 6);
    assert(g000.num == 6);
    assert(get_group_size(d,g000)==gens.family_size*ncrosses);

    GroupNum m = split_by_bv(d, g0, (EffectID){.id=1}, 55, 1);
    assert(d->n_groups == 7);
    assert(m.num == 7);
    assert(get_group_size(d, m) == 50);

    delete_simdata(d);

    // Testing 5 allele effects bug? Replicability.
    d = create_empty_simdata(randomSeed);
    gande = load_all_data(d, "./gt_parents_mr2_50-trimto-5000.txt",
             "./genetic-map_5112-trimto5000.txt",
             "./qtl_5test.txt");
    /*for (int i = 0; i < d->map.n_chr; ++i) {
        printf("%f, ", d->map.chr_lengths[i]);
    }
    printf("\n");*/

    // Testing break-evenly bug
    GroupNum g2 = make_random_crosses(d,gande.group,1200,0,BASIC_OPT);
    delete_group(d,gande.group);
    GroupNum split[2000];
    for (int i = 0; i < 10; ++i) {
        split_evenly_into_n(d,g2,2000,split);
        get_existing_groups(d,NULL);  // only crashes when this is between the split and recombine, not other way around, and not when it's not here.
        g2 = combine_groups(d,2000,split);
    }


	printf("\nAll done\n");
	delete_simdata(d);

	/*clock_t c;

	printf("\n--------Timing tests--------------\n");
	c = clock();
	SimData* sd = create_empty_simdata();
    int fg0 = load_all_data(sd, "./gt_parents_mr2_50-trimto-5000.txt",
			 "./genetic-map_5112-trimto5000.txt",
			 "./qtl_mr2.eff-processed.txt");
    c = clock() - c;
	printf("Loading took %f seconds to run\n", (double)c / CLOCKS_PER_SEC);

	c = clock();
	GenOptions g = {.will_name_offspring=GSC_FALSE, .offspring_name_prefix=NULL, .family_size=1,
		.will_track_pedigree=GSC_TRUE, .will_allocate_ids=GSC_TRUE,
		.filename_prefix="testcross", .will_save_pedigree_to_file=GSC_FALSE,
        .will_save_bvs_to_file=-1, .will_save_alleles_to_file=GSC_FALSE,
		.will_save_to_simdata=GSC_TRUE};
    int f = make_random_crosses(sd, fg0, 100000, g);
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
    save_group_genotypes(fp, sd, f);
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

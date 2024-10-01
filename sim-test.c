#include "sim-test.h"

float calculate_heterozygosity(SimData* d, GroupNum group_number) {
	int hetcount = 0;
    int gn = get_group_size(d, group_number);
    assert(gn < 1000);
    char* galleles[1000];
    get_group_genes(d, group_number, gn, galleles);

	// uses subjects as the first index
	for (int i = 0; i < gn; i++) {
        for (int j = 0; j < d->genome.n_markers; j += 2) {
			if (galleles[i][j] != galleles[i][j + 1]) {
				hetcount += 1;
			}
		}
	}

    return (float) hetcount / (gn * d->genome.n_markers);
}

void write_to_file(const char* filename, const char* contents) {
    FILE* fp;
    if ((fp = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    fwrite(contents, sizeof(char), strlen(contents), fp);
    fclose(fp);
}

int test_savers(unsigned int rseed) {
    SimData* d = create_empty_simdata(rseed);
    just_load(d);
    assert(d->genome.n_markers == 3);
    assert(d->n_eff_sets == 1);
    assert(d->m->n_genotypes == 6);

    // create some interesting groups
    GenOptions g = {.will_name_offspring=GSC_TRUE,
            .offspring_name_prefix="F1",
            .family_size=1,.will_track_pedigree=GSC_TRUE,
            .will_allocate_ids=GSC_TRUE,
            .filename_prefix="",
            .will_save_pedigree_to_file=GSC_FALSE,
            .will_save_bvs_to_file=NO_EFFECTSET,
            .will_save_alleles_to_file=GSC_FALSE,
            .will_save_to_simdata=GSC_TRUE};

    int firstparents[] = {2, 3, 1, 4, 0};
    int secondparents[] = {4, 1, 3, 2, 5};
    GroupNum f1 = gsc_make_targeted_crosses(d, 5, firstparents, secondparents, NO_MAP, NO_MAP, g); // produce 5 offspring
    g.will_track_pedigree = GSC_TRUE;
    g.will_allocate_ids = GSC_FALSE;
    g.offspring_name_prefix = "s";
    GroupNum f2 = self_n_times(d,2,f1, NO_MAP,g); // produce 5 offspring that know their parents but will be anonymous parents themselves
    g.will_name_offspring = GSC_FALSE;
    g.will_allocate_ids = GSC_TRUE;
    make_doubled_haploids(d,f2, NO_MAP,g); // produce 5 offspring that don't have names or know their parents.

    size_t toprint[] = {0,1,//2,3,4,5 from g0
                     6,7,//8,9,10 from f1
                     11,15,//11,12,13,14,15 from f2
                     16,17//18,19,20 from f3
                    };
    GroupNum printingGroup = make_group_from(d, 8, toprint);

    /*GroupNum gout[5];
    int gs[5];
    get_existing_group_counts(d,5,gout,gs);*/

    /*int REFSIZE = 5000;
    char ref[REFSIZE];
    memset(ref,0,sizeof(ref));*/

    // try saving genotypes save_allele_matrix save_transposed_allele_matrix save_group_genotypes save_transposed_group_genotypes
    FILE* fp;
    if ((fp = fopen("test1_save_allele_matrix.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_names_header(fp,d->genome.n_markers,(const char**) d->genome.marker_names);
    save_allele_matrix(fp,d->m);
    fclose(fp);
    /*strcat(ref,"\tm1\tm2\tam3\n");
    for (int i = 0; i < d->m->n_genotypes; ++i) {
        if (d->m->names[i] == NULL) {
            snprintf()
            strcat(ref,"\tm1\tm2\tam3\n");
        } else {
            strcat(ref,"\tm1\tm2\tam3\n");
        }
    }*/

    assert(compareFileToString("test1_save_allele_matrix.txt", TEST1_TRUTH_save_allele_matrix)==0);
    remove("test1_save_allele_matrix.txt");


    if ((fp = fopen("test1_save_transposed_allele_matrix.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_transposed_allele_matrix(fp,d->m,(const char**) d->genome.marker_names);
    fclose(fp);
    assert(compareFileToString("test1_save_transposed_allele_matrix.txt", TEST1_TRUTH_save_transposed_allele_matrix)==0);
    remove("test1_save_transposed_allele_matrix.txt");


    if ((fp = fopen("test1_save_group_genotypes.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_group_genotypes(fp,d,printingGroup);
    fclose(fp);
    assert(compareFileToString("test1_save_group_genotypes.txt", TEST1_TRUTH_save_group_alleles)==0);
    remove("test1_save_group_genotypes.txt");


    if ((fp = fopen("test1_save_transposed_group_genotypes.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_transposed_group_genotypes(fp,d,printingGroup);
    fclose(fp);
    assert(compareFileToString("test1_save_transposed_group_genotypes.txt", TEST1_TRUTH_save_transposed_group_alleles)==0);
    remove("test1_save_transposed_group_genotypes.txt");


    // try saving counts save_count_matrix save_group_count_matrix
    if ((fp = fopen("test1_save_count_matrix.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_count_matrix(fp,d,'A');
    fclose(fp);
    assert(compareFileToString("test1_save_count_matrix.txt", TEST1_TRUTH_save_count_matrix)==0);
    remove("test1_save_count_matrix.txt");


    if ((fp = fopen("test1_save_group_count_matrix.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_group_count_matrix(fp,d,'T',printingGroup);
    fclose(fp);
    assert(compareFileToString("test1_save_group_count_matrix.txt", TEST1_TRUTH_save_count_matrix_of_group)==0);
    remove("test1_save_group_count_matrix.txt");


    printf("...genotype matrix file savers produce the expected output formats\n");


    // try saving bvs save_bvs save_group_bvs
    EffectID effSet1 = {.id=1};

    if ((fp = fopen("test1_save_bvs.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_bvs(fp,d,effSet1);
    fclose(fp);
    assert(compareFileToString("test1_save_bvs.txt", TEST1_TRUTH_save_bvs)==0);
    remove("test1_save_bvs.txt");


    if ((fp = fopen("test1_save_group_bvs.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_group_bvs(fp,d,printingGroup,effSet1);
    fclose(fp);
    assert(compareFileToString("test1_save_group_bvs.txt", TEST1_TRUTH_save_group_bvs)==0);
    remove("test1_save_group_bvs.txt");


    // try saving local gebvs save_marker_blocks calculate_local_bvs
    MarkerBlocks exampleMB = create_evenlength_blocks_each_chr(d, NO_MAP,1);

    if ((fp = fopen("test1_save_marker_blocks.txt", "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    save_markerblocks(fp,d,exampleMB);
    fclose(fp);
    assert(compareFileToString("test1_save_marker_blocks.txt", TEST1_TRUTH_save_marker_blocks)==0);
    remove("test1_save_marker_blocks.txt");


    calculate_local_bvs(d,exampleMB,effSet1,"test1_save_local_bvs.txt");
    assert(compareFileToString("test1_save_local_bvs.txt", TEST1_TRUTH_save_local_bvs)==0);
    remove("test1_save_local_bvs.txt");


    calculate_group_local_bvs(d,exampleMB,effSet1,"test1_save_group_local_bvs.txt",printingGroup);
    assert(compareFileToString("test1_save_group_local_bvs.txt", TEST1_TRUTH_save_group_local_bvs)==0);
    remove("test1_save_group_local_bvs.txt");

    delete_markerblocks(&exampleMB);


    printf("...breeding value file savers produce the expected output formats\n");


    // try saving one-step pedigrees save_group_one_step_pedigree save_one_step_pedigree
    g = (GenOptions){.will_name_offspring=GSC_TRUE,
            .offspring_name_prefix="F2",
            .family_size=1,.will_track_pedigree=GSC_TRUE,
            .will_allocate_ids=GSC_TRUE,
            .filename_prefix="",
            .will_save_pedigree_to_file=GSC_FALSE,
            .will_save_bvs_to_file=NO_EFFECTSET,
            .will_save_alleles_to_file=GSC_FALSE,
            .will_save_to_simdata=GSC_TRUE};
    make_random_crosses(d,f1,1,1, NO_MAP,g);

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
            .will_save_bvs_to_file=effSet1,
            .will_save_alleles_to_file=GSC_TRUE,
            .will_save_to_simdata=GSC_FALSE};
    make_doubled_haploids(d,parent, NO_MAP,g);

    assert(compareRepeatingFileToTable("test1_save_as_you_go-genotype.txt", CONTIG_WIDTH+5,
                                       TEST1_TRUTH_sayg_genotype_header, TEST1_TRUTH_sayg_genotype_bodyrow)==0);
    remove("test1_save_as_you_go-genotype.txt");
    assert(compareRepeatingFileToTable("test1_save_as_you_go-bv.txt", CONTIG_WIDTH+5,
                                       NULL, TEST1_TRUTH_sayg_bv_bodyrow)==0);
    remove("test1_save_as_you_go-bv.txt");
    assert(compareRepeatingFileToTable("test1_save_as_you_go-pedigree.txt", CONTIG_WIDTH+5,
                                       NULL, TEST1_TRUTH_sayg_pedigree_bodyrow)==0);
    remove("test1_save_as_you_go-pedigree.txt");


    delete_simdata(d);
    return 0;
}

int test_effloaders2(void) {
    int it = 0;
    SimData* d = create_empty_simdata(1);
    const char* f1; MapID tmp;
    char filename[] = "0test-eff.txt";

    f1 = TEST1_EFFLOADERS_MAP;
    write_to_file(filename, f1);
    tmp = load_mapfile(d,filename);
    assert(tmp.id == 1);

    // ------- start testing ----------

    // Regular file
    ++it; filename[it / 26] = (it % 26) + 'A';
    assert(it == 1);
    f1 = TEST1_EFFLOADERS[it-1];
    write_to_file(filename, f1);
    assert(load_effectfile(d,filename).id == it);
    fflush(stdout);
    assert(d->n_eff_sets == it);
    assert(d->eff_set_ids[it-1].id == it);
    assert(d->e[it-1].effect_names[0] == 'A');
    assert(d->e[it-1].effect_names[1] == 'T');
    assert(d->e[it-1].effects.rows == 2);
    assert(d->e[it-1].effects.cols == d->genome.n_markers);
    assert(fabs(d->e[it-1].effects.matrix[0][0] - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][0] - - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[0][1] - 1e2) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][1] - 1e-2) < TOL);
	remove(filename);

    // tabs, no final newline
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_EFFLOADERS[it-1];
    write_to_file(filename, f1);
    assert(load_effectfile(d,filename).id == it);
    fflush(stdout);
    assert(d->n_eff_sets == it);
    assert(d->eff_set_ids[it-1].id == it);
    assert(d->e[it-1].effect_names[0] == 'A');
    assert(d->e[it-1].effect_names[1] == 'T');
    assert(d->e[it-1].effects.rows == 2);
    assert(d->e[it-1].effects.cols == d->genome.n_markers);
    assert(fabs(d->e[it-1].effects.matrix[0][0] - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][0] - - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[0][1] - 1e2) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][1] - 1e-2) < TOL);
	remove(filename);

    // mixed spacing
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_EFFLOADERS[it-1];
    write_to_file(filename, f1);
    assert(load_effectfile(d,filename).id == it);
    fflush(stdout);
    assert(d->n_eff_sets == it);
    assert(d->eff_set_ids[it-1].id == it);
    assert(d->e[it-1].effect_names[0] == 'A');
    assert(d->e[it-1].effect_names[1] == 'T');
    assert(d->e[it-1].effects.rows == 2);
    assert(d->e[it-1].effects.cols == d->genome.n_markers);
    assert(fabs(d->e[it-1].effects.matrix[0][0] - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][0] - - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[0][1] - 1e2) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][1] - 1e-2) < TOL);
	remove(filename);
	
    // rearranged rows
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_EFFLOADERS[it-1];
    write_to_file(filename, f1);
    assert(load_effectfile(d,filename).id == it);
    fflush(stdout);
    assert(d->n_eff_sets == it);
    assert(d->eff_set_ids[it-1].id == it);
    assert(d->e[it-1].effect_names[0] == 'A');
    assert(d->e[it-1].effect_names[1] == 'T');
    assert(d->e[it-1].effects.rows == 2);
    assert(d->e[it-1].effects.cols == d->genome.n_markers);
    assert(fabs(d->e[it-1].effects.matrix[0][0] - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][0] - - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[0][1] - 1e2) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][1] - 1e-2) < TOL);
	remove(filename);
	
    // no header
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_EFFLOADERS[it-1];
    write_to_file(filename, f1);
    assert(load_effectfile(d,filename).id == it);
    fflush(stdout);
    assert(d->n_eff_sets == it);
    assert(d->eff_set_ids[it-1].id == it);
    assert(d->e[it-1].effect_names[0] == 'A');
    assert(d->e[it-1].effect_names[1] == 'T');
    assert(d->e[it-1].effects.rows == 2);
    assert(d->e[it-1].effects.cols == d->genome.n_markers);
    assert(fabs(d->e[it-1].effects.matrix[0][0] - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][0] - - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[0][1] - 1e2) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][1] - 1e-2) < TOL);
	remove(filename);
	
    // rearranged columns
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_EFFLOADERS[it-1];
    write_to_file(filename, f1);
    assert(load_effectfile(d,filename).id == it);
    fflush(stdout);
    assert(d->n_eff_sets == it);
    assert(d->eff_set_ids[it-1].id == it);
    assert(d->e[it-1].effect_names[0] == 'A');
    assert(d->e[it-1].effect_names[1] == 'T');
    assert(d->e[it-1].effects.rows == 2);
    assert(d->e[it-1].effects.cols == d->genome.n_markers);
    assert(fabs(d->e[it-1].effects.matrix[0][0] - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][0] - - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[0][1] - 1e2) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][1] - 1e-2) < TOL);
	remove(filename);
	
    // assorted alleles
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_EFFLOADERS[it-1];
    write_to_file(filename, f1);
    assert(load_effectfile(d,filename).id == it);
    fflush(stdout);
    assert(d->n_eff_sets == it);
    assert(d->eff_set_ids[it-1].id == it);
    assert(d->e[it-1].effect_names[0] == '8');
    assert(d->e[it-1].effect_names[1] == 'A');
    assert(d->e[it-1].effect_names[2] == 'T');
    assert(d->e[it-1].effects.rows == 3);
    assert(d->e[it-1].effects.cols == d->genome.n_markers);
    assert(fabs(d->e[it-1].effects.matrix[1][0] - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][1] - 0) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[2][0] - - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[2][1] - 1e2) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[0][1] - 1e-2) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[0][0] - 0) < TOL);
	remove(filename);
	
    // only one row
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_EFFLOADERS[it-1];
    write_to_file(filename, f1);
    assert(load_effectfile(d,filename).id == it);
    fflush(stdout);
    assert(d->n_eff_sets == it);
    assert(d->eff_set_ids[it-1].id == it);
    assert(d->e[it-1].effect_names[0] == 'A');
    assert(d->e[it-1].effects.rows == 1);
    assert(d->e[it-1].effects.cols == d->genome.n_markers);
    assert(fabs(d->e[it-1].effects.matrix[0][0] - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[0][1] - 0) < TOL);
	remove(filename);
	
    // only one row no header
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_EFFLOADERS[it-1];
    write_to_file(filename, f1);
    assert(load_effectfile(d,filename).id == it);
    fflush(stdout);
    assert(d->n_eff_sets == it);
    assert(d->eff_set_ids[it-1].id == it);
    assert(d->e[it-1].effect_names[0] == 'a');
    assert(d->e[it-1].effects.rows == 1);
    assert(d->e[it-1].effects.cols == d->genome.n_markers);
    assert(fabs(d->e[it-1].effects.matrix[0][0] - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[0][1] - 0) < TOL);
	remove(filename);
	
    // discard markers not tracked by the simulation
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_EFFLOADERS[it-1];
    write_to_file(filename, f1);
    assert(load_effectfile(d,filename).id == it);
    fflush(stdout);
    assert(d->n_eff_sets == it);
    assert(d->eff_set_ids[it-1].id == it);
    assert(d->e[it-1].effect_names[0] == 'A');
    assert(d->e[it-1].effect_names[1] == 'T');
    assert(d->e[it-1].effects.rows == 2);
    assert(d->e[it-1].effects.cols == d->genome.n_markers);
    assert(fabs(d->e[it-1].effects.matrix[0][0] - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][0] - - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[0][1] - 1e2) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][1] - 0) < TOL);
	remove(filename);
	
    // file with too many columns on one row
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_EFFLOADERS[it-1];
    write_to_file(filename, f1);
    assert(load_effectfile(d,filename).id == it);
    fflush(stdout);
    assert(d->n_eff_sets == it);
    assert(d->eff_set_ids[it-1].id == it);
    assert(d->e[it-1].effect_names[0] == 'T');
    assert(d->e[it-1].effects.rows == 1);
    assert(d->e[it-1].effects.cols == d->genome.n_markers);
    assert(fabs(d->e[it-1].effects.matrix[0][0] - - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[0][1] - 1e-2) < TOL);
	remove(filename);
	
    // file with a duplicated marker/allele pair
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_EFFLOADERS[it-1];
    write_to_file(filename, f1);
    assert(load_effectfile(d,filename).id == it);
    fflush(stdout);
    assert(d->n_eff_sets == it);
    assert(d->eff_set_ids[it-1].id == it);
    assert(d->e[it-1].effect_names[0] == 'A');
    assert(d->e[it-1].effect_names[1] == 'T');
    assert(d->e[it-1].effects.rows == 2);
    assert(d->e[it-1].effects.cols == d->genome.n_markers);
    assert(fabs(d->e[it-1].effects.matrix[0][0] - 0) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[0][1] - 1e2) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][0] - - 0.5) < TOL);
    assert(fabs(d->e[it-1].effects.matrix[1][1] - 1e-2) < TOL);
	remove(filename);
	
    // file with no valid lines
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_EFFLOADERS[it-1];
    write_to_file(filename, f1);
    assert(load_effectfile(d,filename).id == NO_EFFECTSET.id);
    fflush(stdout);
    assert(d->n_eff_sets == it-1);
	remove(filename);
	
    delete_simdata(d);

    return 0;
}

static void check_mapfile(SimData* d, char* filename, int it) {
    MapID mp = load_mapfile(d, filename);
    fflush(stdout);

    assert(d->genome.n_markers == 5);
    assert(strcmp(d->genome.marker_names[0],"first") == 0);
    assert(strcmp(d->genome.marker_names[1],"second") == 0);
    assert(strcmp(d->genome.marker_names[2],"3rd") == 0);
    assert(strcmp(d->genome.marker_names[3],"fourth") == 0);
    assert(strcmp(d->genome.marker_names[4],"5") == 0);
    assert(strcmp(d->genome.names_alphabetical[0][0],"3rd") == 0);
    assert(strcmp(d->genome.names_alphabetical[1][0],"5") == 0);
    assert(strcmp(d->genome.names_alphabetical[2][0],"first") == 0);
    assert(strcmp(d->genome.names_alphabetical[3][0],"fourth") == 0);
    assert(strcmp(d->genome.names_alphabetical[4][0],"second") == 0);
    assert(d->genome.n_maps == it);
    assert(mp.id == it);
    for (int i = 0; i < it; ++i) {
        assert(d->genome.map_ids[i].id == i+1);
    }
    assert(d->genome.maps[it-1].n_chr == 2);
    assert(d->genome.maps[it-1].chrs[0].type == GSC_LINKAGEGROUP_SIMPLE);
    assert(d->genome.maps[it-1].chrs[1].type == GSC_LINKAGEGROUP_REORDER);
    assert(d->genome.maps[it-1].chrs[0].map.simple.n_markers == 2);
    assert(fabs(d->genome.maps[it-1].chrs[0].map.simple.expected_n_crossovers - 0.4) < TOL);
    assert(fabs(d->genome.maps[it-1].chrs[0].map.simple.dists[0] - 0) < TOL);
    assert(fabs(d->genome.maps[it-1].chrs[0].map.simple.dists[1] - 1) < TOL);
    assert(d->genome.maps[it-1].chrs[0].map.simple.first_marker_index == 0);
    assert(d->genome.maps[it-1].chrs[1].map.reorder.n_markers == 3);
    assert(fabs(d->genome.maps[it-1].chrs[1].map.reorder.expected_n_crossovers - 0.1) < TOL);
    assert(fabs(d->genome.maps[it-1].chrs[1].map.reorder.dists[0] - 0) < TOL);
    assert(fabs(d->genome.maps[it-1].chrs[1].map.reorder.dists[1] - 0.83) < TOL);
    assert(fabs(d->genome.maps[it-1].chrs[1].map.reorder.dists[2] - 1) < TOL);
    assert(d->genome.maps[it-1].chrs[1].map.reorder.marker_indexes[0] == 4);
    assert(d->genome.maps[it-1].chrs[1].map.reorder.marker_indexes[1] == 2);
    assert(d->genome.maps[it-1].chrs[1].map.reorder.marker_indexes[2] == 3);

    SimData* d2 = create_empty_simdata(3);
    assert(load_mapfile(d2, filename).id == 1);

    assert(d2->genome.n_markers == 5);
    assert(strcmp(d2->genome.marker_names[0],"first") == 0);
    assert(strcmp(d2->genome.marker_names[1],"second") == 0);
    assert(strcmp(d2->genome.marker_names[3],"3rd") == 0);
    assert(strcmp(d2->genome.marker_names[4],"fourth") == 0);
    assert(strcmp(d2->genome.marker_names[2],"5") == 0);
    assert(strcmp(d2->genome.names_alphabetical[0][0],"3rd") == 0);
    assert(strcmp(d2->genome.names_alphabetical[1][0],"5") == 0);
    assert(strcmp(d2->genome.names_alphabetical[2][0],"first") == 0);
    assert(strcmp(d2->genome.names_alphabetical[3][0],"fourth") == 0);
    assert(strcmp(d2->genome.names_alphabetical[4][0],"second") == 0);
    assert(d2->genome.n_maps == 1);
    assert(d2->genome.maps[0].n_chr == 2);
    assert(d2->genome.maps[0].chrs[0].type == GSC_LINKAGEGROUP_SIMPLE);
    assert(d2->genome.maps[0].chrs[1].type == GSC_LINKAGEGROUP_SIMPLE);
    assert(d2->genome.maps[0].chrs[0].map.simple.n_markers == 2);
    assert(fabs(d2->genome.maps[0].chrs[0].map.simple.expected_n_crossovers - 0.4) < TOL);
    assert(fabs(d2->genome.maps[0].chrs[0].map.simple.dists[0] - 0) < TOL);
    assert(fabs(d2->genome.maps[0].chrs[0].map.simple.dists[1] - 1) < TOL);
    assert(d2->genome.maps[0].chrs[0].map.simple.first_marker_index == 0);
    assert(d2->genome.maps[0].chrs[1].map.simple.n_markers == 3);
    assert(fabs(d2->genome.maps[0].chrs[1].map.simple.expected_n_crossovers - 0.1) < TOL);
    assert(fabs(d2->genome.maps[0].chrs[1].map.simple.dists[0] - 0) < TOL);
    assert(fabs(d2->genome.maps[0].chrs[1].map.simple.dists[1] - 0.83) < TOL);
    assert(fabs(d2->genome.maps[0].chrs[1].map.simple.dists[2] - 1) < TOL);
    assert(d2->genome.maps[0].chrs[1].map.simple.first_marker_index == 2);
    delete_simdata(d2);
}

int test_maploaders2(void) {
    FILE* fp;
    int it = 0;
    SimData* d = create_empty_simdata(1);
    const char* f1; MapID tmp;
    char filename[] = "0test-map.txt";

    f1 = TEST1_MAPLOADERS[it];
    if ((fp = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    fwrite(f1, sizeof(char), strlen(f1), fp);
    fclose(fp);
    tmp = load_mapfile(d,filename);
    assert(tmp.id == 1);
    delete_recombination_map(d, tmp);
	remove(filename);

    // First real test: mixed spacing
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_MAPLOADERS[it];
    if ((fp = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    fwrite(f1, sizeof(char), strlen(f1), fp);
    fclose(fp);
    check_mapfile(d, filename, it);
	remove(filename);
	
    // No final line
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_MAPLOADERS[it];
    if ((fp = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    fwrite(f1, sizeof(char), strlen(f1), fp);
    fclose(fp);
    check_mapfile(d, filename, it);
	remove(filename);
	
    // rearranged column order
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_MAPLOADERS[it];
    if ((fp = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    fwrite(f1, sizeof(char), strlen(f1), fp);
    fclose(fp);
    check_mapfile(d, filename, it);
	remove(filename);

    // No header line
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_MAPLOADERS[it];
    if ((fp = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    fwrite(f1, sizeof(char), strlen(f1), fp);
    fclose(fp);
    check_mapfile(d, filename, it);
	remove(filename);
	
    // Mix order
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_MAPLOADERS[it];
    if ((fp = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    fwrite(f1, sizeof(char), strlen(f1), fp);
    fclose(fp);
    check_mapfile(d, filename, it);
    remove(filename);
	
    // Discard extra markers, leave out a missing marker, and change chromosomes and positions
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_MAPLOADERS[it];
    write_to_file(filename, f1);
    assert(load_mapfile(d,filename).id == it);
    assert(d->genome.n_markers == 5);
    assert(strcmp(d->genome.marker_names[0],"first") == 0);
    assert(strcmp(d->genome.marker_names[1],"second") == 0);
    assert(strcmp(d->genome.marker_names[2],"3rd") == 0);
    assert(strcmp(d->genome.marker_names[3],"fourth") == 0);
    assert(strcmp(d->genome.marker_names[4],"5") == 0);
    assert(strcmp(d->genome.names_alphabetical[0][0],"3rd") == 0);
    assert(strcmp(d->genome.names_alphabetical[1][0],"5") == 0);
    assert(strcmp(d->genome.names_alphabetical[2][0],"first") == 0);
    assert(strcmp(d->genome.names_alphabetical[3][0],"fourth") == 0);
    assert(strcmp(d->genome.names_alphabetical[4][0],"second") == 0);
    assert(d->genome.n_maps == it);
    for (int i = 0; i < it; ++i) {
        assert(d->genome.map_ids[i].id == i+1);
    }
    assert(d->genome.maps[it-1].n_chr == 2);
    assert(d->genome.maps[it-1].chrs[0].type == GSC_LINKAGEGROUP_REORDER);
    assert(d->genome.maps[it-1].chrs[1].type == GSC_LINKAGEGROUP_REORDER);
    assert(d->genome.maps[it-1].chrs[1].map.reorder.n_markers == 2);
    assert(fabs(d->genome.maps[it-1].chrs[1].map.reorder.expected_n_crossovers - 0.4) < TOL);
    assert(fabs(d->genome.maps[it-1].chrs[1].map.reorder.dists[0] - 0) < TOL);
    assert(fabs(d->genome.maps[it-1].chrs[1].map.reorder.dists[1] - 1) < TOL);
    assert(d->genome.maps[it-1].chrs[1].map.reorder.marker_indexes[0] == 1);
    assert(d->genome.maps[it-1].chrs[1].map.reorder.marker_indexes[1] == 0);
    assert(d->genome.maps[it-1].chrs[0].map.reorder.n_markers == 2);
    assert(fabs(d->genome.maps[it-1].chrs[0].map.reorder.expected_n_crossovers - 0.083) < TOL);
    assert(fabs(d->genome.maps[it-1].chrs[0].map.reorder.dists[0] - 0) < TOL);
    assert(fabs(d->genome.maps[it-1].chrs[0].map.reorder.dists[1] - 1) < TOL);
    assert(d->genome.maps[it-1].chrs[0].map.reorder.marker_indexes[0] == 3);
    assert(d->genome.maps[it-1].chrs[0].map.reorder.marker_indexes[1] == 2);
    remove(filename);
	
    // Too many columns in one place
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_MAPLOADERS[it];
    write_to_file(filename, f1);
    assert(load_mapfile(d,filename).id == it);
    assert(d->genome.n_markers == 5);
    assert(d->genome.maps[it-1].n_chr == 2);
    assert(d->genome.maps[it-1].chrs[1].type == GSC_LINKAGEGROUP_REORDER);
    assert(d->genome.maps[it-1].chrs[0].type == GSC_LINKAGEGROUP_SIMPLE);
    assert(d->genome.maps[it-1].chrs[0].map.simple.n_markers == 2);
    assert(d->genome.maps[it-1].chrs[0].map.simple.first_marker_index == 0);
    assert(d->genome.maps[it-1].chrs[1].map.reorder.n_markers == 2);
    assert(fabs(d->genome.maps[it-1].chrs[1].map.reorder.expected_n_crossovers - 0.1) < TOL);
    assert(fabs(d->genome.maps[it-1].chrs[1].map.reorder.dists[0] - 0) < TOL);
    assert(fabs(d->genome.maps[it-1].chrs[1].map.reorder.dists[1] - 1) < TOL);
    assert(d->genome.maps[it-1].chrs[1].map.reorder.marker_indexes[0] == 4);
    assert(d->genome.maps[it-1].chrs[1].map.reorder.marker_indexes[1] == 3);
    remove(filename);
	
    // Single line map with header. A single line map doesn't seem very useful but we may as well check against crashes
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_MAPLOADERS[it];
    write_to_file(filename, f1);
    assert(load_mapfile(d,filename).id == it);
    assert(d->genome.n_markers == 5);
    assert(d->genome.maps[it-1].n_chr == 1);
    assert(d->genome.maps[it-1].chrs[0].type == GSC_LINKAGEGROUP_SIMPLE);
    assert(d->genome.maps[it-1].chrs[0].map.simple.n_markers == 1);
    assert(d->genome.maps[it-1].chrs[0].map.simple.first_marker_index == 1);
    remove(filename);
	
    // single line map without header
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_MAPLOADERS[it];
    write_to_file(filename, f1);
    assert(load_mapfile(d,filename).id == it);
    assert(d->genome.n_markers == 5);
    assert(d->genome.maps[it-1].n_chr == 1);
    assert(d->genome.maps[it-1].chrs[0].type == GSC_LINKAGEGROUP_SIMPLE);
    assert(d->genome.maps[it-1].chrs[0].map.simple.n_markers == 1);
    assert(d->genome.maps[it-1].chrs[0].map.simple.first_marker_index == 1);
	remove(filename);
	
    // Single line inappropriate map
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_MAPLOADERS[it];
    write_to_file(filename, f1);
    assert(load_mapfile(d,filename).id == NO_MAP.id);
    assert(d->genome.n_markers == 5);
	remove(filename);

    delete_simdata(d);
	
    return 0;

}

static void check_genotypes(SimData* d, int it, int skip_names) {
    fflush(stdout);
    const int nm = 2;
    const int nl = 3;

    assert(d->n_groups == it);
    assert(d->current_id.id == nl*it);
    assert(d->genome.n_markers == nm);
    assert(d->genome.n_maps == 1);
    assert(d->n_labels == 0);
    assert(d->n_eff_sets == 0);
    assert(d->m->n_labels == 0);
    assert(d->m->next == NULL);
    assert(d->m->n_genotypes == nl*it);
    assert(d->m->n_markers == nm);
    assert(strncmp(d->m->alleles[nl*(it-1) + 0],"AATA",2*nm) == 0);
    assert(strncmp(d->m->alleles[nl*(it-1) + 1],"ATAA",2*nm) == 0);
    assert(strncmp(d->m->alleles[nl*(it-1) + 2],"AATT",2*nm) == 0);
    if (!skip_names) {
        assert(strcmp(d->m->names[nl*(it-1) + 0],"g1") == 0);
        assert(strcmp(d->m->names[nl*(it-1) + 1],"g2") == 0);
        assert(strcmp(d->m->names[nl*(it-1) + 2],"oog3") == 0);
    }
    for (int i = 0; i < nl; ++i) {
        assert(d->m->ids[nl*(it-1) + i].id == nl*(it-1) + i + 1);
        assert(d->m->groups[nl*(it-1) + i].num == it);
        assert(d->m->pedigrees[0][nl*(it-1) + i].id == 0);
        assert(d->m->pedigrees[1][nl*(it-1) + i].id == 0);
    }
}

static void check_matrix_with_different_specification_levels(int seed, char* filename, const char* mapfile, int only_full,
                                                      int has_header, int markers_as_rows, enum gsc_GenotypeFileCellStyle style) {
    FileFormatSpec spec = {
        .filetype=GSC_GENOTYPEFILE_MATRIX,
        .spec={(struct gsc_GenotypeFile_MatrixFormat){
            .has_header=has_header,
            .markers_as_rows=markers_as_rows,
            .cell_style=style
        }}
    };
    int it = 1;

    // Full spec
    SimData* d2 = create_empty_simdata(seed);
    load_mapfile(d2,mapfile);
    load_genotypefile(d2,filename,spec);
    check_genotypes(d2,it,0);

    if (!only_full) {
        // No spec, detect all
        load_genotypefile(d2,filename,DETECT_FILE_FORMAT);
        ++it; check_genotypes(d2,it,0);

        // only header
        spec.spec.matrix.has_header = has_header;
        spec.spec.matrix.markers_as_rows = GSC_UNINIT;
        spec.spec.matrix.cell_style = GSC_GENOTYPECELLSTYLE_UNKNOWN;
        load_genotypefile(d2,filename,spec);
        ++it; check_genotypes(d2,it,0);

        // only orientation
        spec.spec.matrix.has_header = GSC_UNINIT;
        spec.spec.matrix.markers_as_rows = markers_as_rows;
        spec.spec.matrix.cell_style = GSC_GENOTYPECELLSTYLE_UNKNOWN;
        load_genotypefile(d2,filename,spec);
        ++it; check_genotypes(d2,it,0);

        // only style
        spec.spec.matrix.has_header = GSC_UNINIT;
        spec.spec.matrix.markers_as_rows = GSC_UNINIT;
        spec.spec.matrix.cell_style = style;
        load_genotypefile(d2,filename,spec);
        ++it; check_genotypes(d2,it,0);

        // header + orientation
        spec.spec.matrix.has_header = has_header;
        spec.spec.matrix.markers_as_rows = markers_as_rows;
        spec.spec.matrix.cell_style = GSC_GENOTYPECELLSTYLE_UNKNOWN;
        load_genotypefile(d2,filename,spec);
        ++it; check_genotypes(d2,it,0);

    }
    delete_simdata(d2);
}

int test_genoloaders2(void) {
    // Test genotype matrix loading
    /*
name g1 g2 g3
m1 AA AA AA
m2 AT AA TT
    */
    int it = 0;
    int seed = 1;
    SimData* d = create_empty_simdata(seed);
    const char* f1;
    char filename[] = "0test.txt";

    // Prepare map, for check_with_different_specification_levels
    f1 = TEST1_2MARKER_MAP;
    const char* fmap = "2mtest.txt";
    write_to_file(fmap, f1);

    // vertical first, tabs, with corner, with endline. "B"
    ++it; filename[it / 26] = (it % 26) + 'A';
	assert(it == 1);
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,define_matrix_format_details(GSC_UNINIT,GSC_FALSE,GSC_GENOTYPECELLSTYLE_UNKNOWN));
    check_genotypes(d,it,0);
    check_matrix_with_different_specification_levels(seed, filename, fmap, GSC_FALSE,
                                                     GSC_TRUE,GSC_FALSE,GSC_GENOTYPECELLSTYLE_PAIR);
	remove(filename);

    // tabs, wcorner, wendline C
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
    check_matrix_with_different_specification_levels(seed,filename,fmap,1,
                                                     GSC_TRUE,GSC_TRUE,GSC_GENOTYPECELLSTYLE_PAIR);
	remove(filename);
	
    // spaces, wcorner, wendline D
    ++it; filename[it / 26] = it+ 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // spaces, wcorner, wrendline E
    ++it; filename[it / 26] = it+ 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // mixed endlines and spacings, wcorner, wrendline F
    ++it; filename[it / 26] = it+ 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // tabs, noendline, wcorner G
    ++it; filename[it / 26] = it+ 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // spaces, noendline, wcorner H
    ++it; filename[it / 26] = it+ 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // spaces, noendline, nocorner I
    ++it; filename[it / 26] = it+ 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // tabs, noendline, nocorner J
    ++it; filename[it / 26] = it+ 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // K
    ++it; filename[it / 26] = it+ 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // tabs, nocorner, wendline L
    ++it; filename[it / 26] = it+ 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);

    // vertical, tabs, noendline M
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // vertical, tabs, nocorner, noendline N
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
    check_matrix_with_different_specification_levels(seed,filename,fmap,1,
                                                     GSC_TRUE,GSC_FALSE,GSC_GENOTYPECELLSTYLE_PAIR);
	remove(filename);

    // vertical, tabs, nocorner O
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // vertical, commas, nocorner P
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // vertical, commas, corner, ignore extra Q
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // vertical, commas, corner, ignore extra, noendline R
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1]; 
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    fflush(stdout);
    check_genotypes(d,it,0);
	remove(filename);
	
    // vertical, slashpairs S
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1]; 
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // horizontal, slashpairs T
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1]; 
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // horizontal, rearranged markers U
    ++it; filename[it / 26] = it+ 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];  
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // vertical, rearranged markers V
    ++it; filename[it / 26] = it+ 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // vertical, rearranged markers and extra fake marker W
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1]; 
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // horizontal, extra fake marker X
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
    check_matrix_with_different_specification_levels(seed,filename,fmap,1,
                                                     GSC_TRUE,GSC_FALSE,GSC_GENOTYPECELLSTYLE_PAIR);
	remove(filename);

    // horizontal, extra fake marker centrally Y
    ++it; filename[it / 26] = it+ 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1]; ;
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,0);
	remove(filename);
	
    // vertical, blank names Z
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1]; 
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,1);
	remove(filename);
	
    // Markers as rows, no genotype names
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,1);
	remove(filename);
	
    // Markers as rows, no genotype names v2
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,1);
	remove(filename);
	
    // Markers as columns, no genotype names "ZC"
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,1);
	remove(filename);
	
    // Markers as columns, no genotype names v2 "ZD"
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d,filename,DETECT_FILE_FORMAT);
    check_genotypes(d,it,1);
	remove(filename);
	
    // temporarily create new SimData with map.
    f1 = TEST1_2MARKER_MAP; 
    write_to_file(filename, f1);
    SimData* d2 = create_empty_simdata(1);
    load_mapfile(d2, "2marker-map.txt");
    assert(d2->genome.n_markers == 2);
    int cid = 0;
	remove(filename);
	
    // horizontal, single-line w header
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1]; 
    write_to_file(filename, f1);
    load_genotypefile(d2,filename,DETECT_FILE_FORMAT);
    cid += 3;
    assert(d2->n_groups == 1);
    assert(d2->current_id.id == cid);
    assert(d2->genome.n_markers == 2);
    assert(d2->genome.n_maps == 1);
    assert(d2->n_labels == 0);
    assert(d2->n_eff_sets == 0);
    assert(d2->m->n_labels == 0);
    assert(d2->m->next == NULL);
    assert(d2->m->n_genotypes == 3);
    assert(d2->m->n_markers == 2);
    assert(d2->m->alleles[0][0] == 'A'); assert(d2->m->alleles[0][1] == 'A');
    assert(d2->m->alleles[0][2] == '\0'); assert(d2->m->alleles[0][3] == '\0');
    assert(d2->m->alleles[1][0] == 'A'); assert(d2->m->alleles[1][1] == 'T');
    assert(d2->m->alleles[1][2] == '\0'); assert(d2->m->alleles[1][3] == '\0');
    assert(d2->m->alleles[2][0] == 'A'); assert(d2->m->alleles[2][1] == 'A');
    assert(d2->m->alleles[2][2] == '\0'); assert(d2->m->alleles[2][3] == '\0');
    assert(strcmp(d2->m->names[0],"g1") == 0);
    assert(strcmp(d2->m->names[1],"g2") == 0);
    assert(strcmp(d2->m->names[2],"oog3") == 0);
    for (int i = 0; i < 3; ++i) {
        assert(d2->m->ids[i].id == (cid - 2) + i);
        assert(d2->m->groups[i].num == 1);
        assert(d2->m->pedigrees[0][i].id == 0);
        assert(d2->m->pedigrees[1][i].id == 0);
    }
    delete_group(d2,(GroupNum){.num=1});
	remove(filename);
	
    // horizontal, single line no header
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d2,filename,DETECT_FILE_FORMAT);
    cid += 3;
    change_allele_symbol(d2,NULL,'\0','-');
    fflush(stdout);
    assert(d2->n_groups == 1);
    assert(d2->current_id.id == cid);
    assert(d2->genome.n_markers == 2);
    assert(d2->genome.n_maps == 1);
    assert(d2->n_labels == 0);
    assert(d2->n_eff_sets == 0);
    assert(d2->m->n_labels == 0);
    assert(d2->m->next == NULL);
    assert(d2->m->n_genotypes == 3);
    assert(d2->m->n_markers == 2);
    assert(strncmp(d2->m->alleles[0],"--AA",2*2) == 0);
    assert(strncmp(d2->m->alleles[1],"--TA",2*2) == 0);
    assert(strncmp(d2->m->alleles[2],"--AA",2*2) == 0);
    assert(d2->m->names[0] == NULL);
    assert(d2->m->names[1] == NULL);
    assert(d2->m->names[2] == NULL);
    for (int i = 0; i < 3; ++i) {
        assert(d2->m->ids[i].id == (cid - 2) + i);
        assert(d2->m->groups[i].num == 1);
        assert(d2->m->pedigrees[0][i].id == 0);
        assert(d2->m->pedigrees[1][i].id == 0);
    }
    delete_group(d2,(GroupNum){.num=1});
	remove(filename);

    // horizontal, single column
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d2,filename,DETECT_FILE_FORMAT);
    ++cid;
    fflush(stdout);
    assert(d2->n_groups == 1);
    assert(d2->current_id.id == cid);
    assert(d2->genome.n_markers == 2);
    assert(d2->genome.n_maps == 1);
    assert(d2->n_labels == 0);
    assert(d2->n_eff_sets == 0);
    assert(d2->m->n_labels == 0);
    assert(d2->m->next == NULL);
    assert(d2->m->n_genotypes == 1);
    assert(d2->m->n_markers == 2);
    assert(strncmp(d2->m->alleles[0],"AATA",2*2) == 0);
    assert(strcmp(d2->m->names[0],"g1") == 0);
    for (int i = 0; i < 1; ++i) {
        assert(d2->m->ids[i].id == cid - i);
        assert(d2->m->groups[i].num == 1);
        assert(d2->m->pedigrees[0][i].id == 0);
        assert(d2->m->pedigrees[1][i].id == 0);
    }
    delete_group(d2,(GroupNum){.num=1});
	remove(filename);
	
    // horizontal, single column, noheader
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d2,filename,DETECT_FILE_FORMAT);
    ++cid;
    fflush(stdout);
    assert(d2->n_groups == 1);
    assert(d2->current_id.id == cid);
    assert(d2->genome.n_markers == 2);
    assert(d2->genome.n_maps == 1);
    assert(d2->n_labels == 0);
    assert(d2->n_eff_sets == 0);
    assert(d2->m->n_labels == 0);
    assert(d2->m->next == NULL);
    assert(d2->m->n_genotypes == 1);
    assert(d2->m->n_markers == 2);
    assert(strncmp(d2->m->alleles[0],"AATA",2*2) == 0);
    assert(d2->m->names[0] == NULL);
    for (int i = 0; i < 1; ++i) {
        assert(d2->m->ids[i].id == cid + i);
        assert(d2->m->groups[i].num == 1);
        assert(d2->m->pedigrees[0][i].id == 0);
        assert(d2->m->pedigrees[1][i].id == 0);
    }
    delete_group(d2,(GroupNum){.num=1});
	remove(filename);
	
    // vertical, single line
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1]; 
    write_to_file(filename, f1);
    load_genotypefile(d2,filename,DETECT_FILE_FORMAT);
    ++cid;
    fflush(stdout);
    assert(d2->n_groups == 1);
    assert(d2->current_id.id == cid);
    assert(d2->genome.n_markers == 2);
    assert(d2->genome.n_maps == 1);
    assert(d2->n_labels == 0);
    assert(d2->n_eff_sets == 0);
    assert(d2->m->n_labels == 0);
    assert(d2->m->next == NULL);
    assert(d2->m->n_genotypes == 1);
    assert(d2->m->n_markers == 2);
    assert(strncmp(d2->m->alleles[0],"AATA",2*2) == 0);
    assert(strcmp(d2->m->names[0],"g1") == 0);
    for (int i = 0; i < 1; ++i) {
        assert(d2->m->ids[i].id == cid + i);
        assert(d2->m->groups[i].num == 1);
        assert(d2->m->pedigrees[0][i].id == 0);
        assert(d2->m->pedigrees[1][i].id == 0);
    }
    delete_group(d2,(GroupNum){.num=1});
	remove(filename);
	
    // vertical singleline noheader -> fails, can't guarantee marker order.
    //++it; filename[it / 26] = (it % 26) + 'A';
    /*f1 = TEST1_GENOMATRIX_LOADERS[it-1]; "g1,AA,TA\n";
    if ((fp = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "Failed to create file.\n");
        exit(1);
    }
    fwrite(f1, sizeof(char), strlen(f1), fp);
    fclose(fp);
    load_genotypefile(d2,filename);
    ++cid;
    fflush(stdout);
    assert(d2->n_groups == 1);
    assert(d2->current_id.id == cid);
    assert(d2->genome.n_markers == 2);
    assert(d2->genome.n_maps == 1);
    assert(d2->n_labels == 0);
    assert(d2->n_eff_sets == 0);
    assert(d2->m->n_labels == 0);
    assert(d2->m->next == NULL);
    assert(d2->m->n_genotypes == 1);
    assert(d2->m->n_markers == 2);
    assert(strncmp(d2->m->alleles[0],"AATA",2*2) == 0);
    assert(strcmp(d2->m->names[0],"g1") == 0);
    for (int i = 0; i < 1; ++i) {
        assert(d2->m->ids[i].id == cid + i);
        assert(d2->m->groups[i].num == 1);
        assert(d2->m->pedigrees[0][i].id == 0);
        assert(d2->m->pedigrees[1][i].id == 0);
    }
    delete_group(d2,(GroupNum){.num=1});
	remove(filename);*/

    // vertical singlecolumn
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d2,filename,DETECT_FILE_FORMAT);
    ++cid;
    fflush(stdout);
    assert(d2->n_groups == 1);
    assert(d2->current_id.id == cid);
    assert(d2->genome.n_markers == 2);
    assert(d2->genome.n_maps == 1);
    assert(d2->n_labels == 0);
    assert(d2->n_eff_sets == 0);
    assert(d2->m->n_labels == 0);
    assert(d2->m->next == NULL);
    assert(d2->m->n_genotypes == 1);
    assert(d2->m->n_markers == 2);
    assert(d2->m->alleles[0][0] == '\0'); assert(d2->m->alleles[0][1] == '\0');
    assert(d2->m->alleles[0][2] == 'A'); assert(d2->m->alleles[0][3] == 'A');
    assert(strcmp(d2->m->names[0],"g1") == 0);
    for (int i = 0; i < 1; ++i) {
        assert(d2->m->ids[i].id == cid + i);
        assert(d2->m->groups[i].num == 1);
        assert(d2->m->pedigrees[0][i].id == 0);
        assert(d2->m->pedigrees[1][i].id == 0);
    }
    delete_group(d2,(GroupNum){.num=1});
	remove(filename);

    // Markers rows, single row no header
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d2,filename,define_matrix_format_details(GSC_FALSE,GSC_UNINIT,GSC_GENOTYPECELLSTYLE_UNKNOWN));
    cid += 3;
    fflush(stdout);
    assert(d2->n_groups == 1);
    assert(d2->current_id.id == cid);
    assert(d2->genome.n_markers == 2);
    assert(d2->genome.n_maps == 1);
    assert(d2->n_labels == 0);
    assert(d2->n_eff_sets == 0);
    assert(d2->m->n_labels == 0);
    assert(d2->m->next == NULL);
    assert(d2->m->n_genotypes == 3);
    assert(d2->m->n_markers == 2);
    assert(d2->m->alleles[0][0] == 'T'); assert(d2->m->alleles[0][1] == 'T');
    assert(d2->m->alleles[0][2] == '\0'); assert(d2->m->alleles[0][3] == '\0');
    assert(d2->m->names[0] == NULL);
    assert(d2->m->alleles[1][0] == '\0'); assert(d2->m->alleles[1][1] == '\0');
    assert(d2->m->alleles[1][2] == '\0'); assert(d2->m->alleles[1][3] == '\0');
    assert(d2->m->names[1] == NULL);
    assert(d2->m->alleles[2][0] == 'A'); assert(d2->m->alleles[2][1] == 'A');
    assert(d2->m->alleles[2][2] == '\0'); assert(d2->m->alleles[2][3] == '\0');
    assert(d2->m->names[2] == NULL);
    for (int i = 0; i < 3; ++i) {
        assert(d2->m->ids[i].id == cid - (3-1) + i);
        assert(d2->m->groups[i].num == 1);
        assert(d2->m->pedigrees[0][i].id == 0);
        assert(d2->m->pedigrees[1][i].id == 0);
    }
    delete_group(d2,(GroupNum){.num=1});
	remove(filename);
	
    // Markers rows, single row no body
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d2,filename,DETECT_FILE_FORMAT);
    cid += 3;
    fflush(stdout);
    assert(d2->n_groups == 1);
    assert(d2->current_id.id == cid);
    assert(d2->genome.n_markers == 2);
    assert(d2->genome.n_maps == 1);
    assert(d2->n_labels == 0);
    assert(d2->n_eff_sets == 0);
    assert(d2->m->n_labels == 0);
    assert(d2->m->next == NULL);
    assert(d2->m->n_genotypes == 3);
    assert(d2->m->n_markers == 2);
    assert(d2->m->alleles[0][0] == '\0'); assert(d2->m->alleles[0][1] == '\0');
    assert(d2->m->alleles[0][2] == '\0'); assert(d2->m->alleles[0][3] == '\0');
    assert(strcmp(d2->m->names[0],"g1") == 0);
    assert(d2->m->alleles[1][0] == '\0'); assert(d2->m->alleles[1][1] == '\0');
    assert(d2->m->alleles[1][2] == '\0'); assert(d2->m->alleles[1][3] == '\0');
    assert(strcmp(d2->m->names[1],"g2") == 0);
    assert(d2->m->alleles[2][0] == '\0'); assert(d2->m->alleles[2][1] == '\0');
    assert(d2->m->alleles[2][2] == '\0'); assert(d2->m->alleles[2][3] == '\0');
    assert(strcmp(d2->m->names[2],"oog3") == 0);
    for (int i = 0; i < 3; ++i) {
        assert(d2->m->ids[i].id == cid - (3-1) + i);
        assert(d2->m->groups[i].num == 1);
        assert(d2->m->pedigrees[0][i].id == 0);
        assert(d2->m->pedigrees[1][i].id == 0);
    }
    delete_group(d2,(GroupNum){.num=1});
	remove(filename);
	
    // Markers as columns, single row, no body
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d2,filename,DETECT_FILE_FORMAT);
    fflush(stdout);
    assert(d2->n_groups == 0);
    assert(d2->current_id.id == cid);
    assert(d2->genome.n_markers == 2);
    assert(d2->genome.n_maps == 1);
    assert(d2->n_labels == 0);
    assert(d2->n_eff_sets == 0);
    assert(d2->m->n_labels == 0);
    assert(d2->m->next == NULL);
    assert(d2->m->n_genotypes == 0);
    assert(d2->m->n_markers == 2);
	remove(filename);
	
    // Markers as rows, single column, no body
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d2,filename,DETECT_FILE_FORMAT);
    fflush(stdout);
    assert(d2->n_groups == 0);
    assert(d2->current_id.id == cid);
    assert(d2->genome.n_markers == 2);
    assert(d2->genome.n_maps == 1);
    assert(d2->n_labels == 0);
    assert(d2->n_eff_sets == 0);
    assert(d2->m->n_labels == 0);
    assert(d2->m->next == NULL);
    assert(d2->m->n_genotypes == 0);
    assert(d2->m->n_markers == 2);
	remove(filename);
	
    // Markers as rows, no marker name matches
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d2,filename,DETECT_FILE_FORMAT);
    cid += 3;
    fflush(stdout);
    assert(d2->n_groups == 1);
    assert(d2->current_id.id == cid);
    assert(d2->genome.n_markers == 2);
    assert(d2->genome.n_maps == 1);
    assert(d2->n_labels == 0);
    assert(d2->n_eff_sets == 0);
    assert(d2->m->n_labels == 0);
    assert(d2->m->next == NULL);
    assert(d2->m->n_genotypes == 3);
    assert(d2->m->n_markers == 2);
    assert(d2->m->alleles[0][0] == '\0'); assert(d2->m->alleles[0][1] == '\0');
    assert(d2->m->alleles[0][2] == '\0'); assert(d2->m->alleles[0][3] == '\0');
    assert(strcmp(d2->m->names[0],"g1") == 0);
    assert(d2->m->alleles[1][0] == '\0'); assert(d2->m->alleles[1][1] == '\0');
    assert(d2->m->alleles[1][2] == '\0'); assert(d2->m->alleles[1][3] == '\0');
    assert(strcmp(d2->m->names[1],"g2") == 0);
    assert(d2->m->alleles[2][0] == '\0'); assert(d2->m->alleles[2][1] == '\0');
    assert(d2->m->alleles[2][2] == '\0'); assert(d2->m->alleles[2][3] == '\0');
    assert(strcmp(d2->m->names[2],"oog3") == 0);
    for (int i = 0; i < 3; ++i) {
        assert(d2->m->ids[i].id == cid - (3-1) + i);
        assert(d2->m->groups[i].num == 1);
        assert(d2->m->pedigrees[0][i].id == 0);
        assert(d2->m->pedigrees[1][i].id == 0);
    }
    delete_group(d2,(GroupNum){.num=1});
	remove(filename);
	
    // Markers as columns, no marker name matches
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
    load_genotypefile(d2,filename,define_matrix_format_details(GSC_UNINIT,GSC_FALSE,GSC_GENOTYPECELLSTYLE_UNKNOWN));
    cid += 3;
    fflush(stdout);
    assert(d2->n_groups == 1);
    assert(d2->current_id.id == cid);
    assert(d2->genome.n_markers == 2);
    assert(d2->genome.n_maps == 1);
    assert(d2->n_labels == 0);
    assert(d2->n_eff_sets == 0);
    assert(d2->m->n_labels == 0);
    assert(d2->m->next == NULL);
    assert(d2->m->n_genotypes == 3);
    assert(d2->m->n_markers == 2);
    assert(d2->m->alleles[0][0] == '\0'); assert(d2->m->alleles[0][1] == '\0');
    assert(d2->m->alleles[0][2] == '\0'); assert(d2->m->alleles[0][3] == '\0');
    assert(strcmp(d2->m->names[0],"g1") == 0);
    assert(d2->m->alleles[1][0] == '\0'); assert(d2->m->alleles[1][1] == '\0');
    assert(d2->m->alleles[1][2] == '\0'); assert(d2->m->alleles[1][3] == '\0');
    assert(strcmp(d2->m->names[1],"g2") == 0);
    assert(d2->m->alleles[2][0] == '\0'); assert(d2->m->alleles[2][1] == '\0');
    assert(d2->m->alleles[2][2] == '\0'); assert(d2->m->alleles[2][3] == '\0');
    assert(strcmp(d2->m->names[2],"oog3") == 0);
    for (int i = 0; i < 3; ++i) {
        assert(d2->m->ids[i].id == cid - (3-1) + i);
        assert(d2->m->groups[i].num == 1);
        assert(d2->m->pedigrees[0][i].id == 0);
        assert(d2->m->pedigrees[1][i].id == 0);
    }
    delete_group(d2,(GroupNum){.num=1});
	remove(filename);
	
    // Check counts work
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1]; 
    write_to_file(filename, f1);
    load_genotypefile(d2,filename,DETECT_FILE_FORMAT);
    cid += 3;
    fflush(stdout);
    assert(d2->n_groups == 1);
    assert(d2->current_id.id == cid);
    assert(d2->genome.n_markers == 2);
    assert(d2->genome.n_maps == 1);
    assert(d2->n_labels == 0);
    assert(d2->n_eff_sets == 0);
    assert(d2->m->n_labels == 0);
    assert(d2->m->next == NULL);
    assert(d2->m->n_genotypes == 3);
    assert(d2->m->n_markers == 2);
    assert(strncmp(d2->m->alleles[0],"ATAT",2*2) == 0);
    assert(strncmp(d2->m->alleles[1],"TATT",2*2) == 0);
    assert(strncmp(d2->m->alleles[2],"AAAT",2*2) == 0);
    assert(strcmp(d2->m->names[0],"g1") == 0);
    assert(strcmp(d2->m->names[1],"g2") == 0);
    assert(strcmp(d2->m->names[2],"oog3") == 0);
    for (int i = 0; i < 3; ++i) {
        assert(d2->m->ids[i].id == (cid - 2) + i);
        assert(d2->m->groups[i].num == 1);
        assert(d2->m->pedigrees[0][i].id == 0);
        assert(d2->m->pedigrees[1][i].id == 0);
    }
    delete_group(d2,(GroupNum){.num=1});
	remove(filename);
	
    // Check encodings work
    ++it; filename[it / 26] = (it % 26) + 'A';
    f1 = TEST1_GENOMATRIX_LOADERS[it-1];
    write_to_file(filename, f1);
     load_genotypefile(d2,filename,DETECT_FILE_FORMAT);
    cid += 15;
    fflush(stdout);
    assert(d2->n_groups == 1);
    assert(d2->current_id.id == cid);
    assert(d2->genome.n_markers == 2);
    assert(d2->genome.n_maps == 1);
    assert(d2->n_labels == 0);
    assert(d2->n_eff_sets == 0);
    assert(d2->m->n_labels == 0);
    assert(d2->m->next == NULL);
    assert(d2->m->n_genotypes == 15);
    assert(d2->m->n_markers == 2);
    assert(strncmp(d2->m->alleles[0],"GGAA",2*2) == 0);
    assert(strncmp(d2->m->alleles[1],"TTCC",2*2) == 0);
    assert(strncmp(d2->m->alleles[2],"GAAG",2*2) == 0);
    assert(strncmp(d2->m->alleles[3],"AGGA",2*2) == 0);
    assert(strncmp(d2->m->alleles[4],"TCTC",2*2) == 0);
    assert(strncmp(d2->m->alleles[5],"CTCT",2*2) == 0);
    assert(strncmp(d2->m->alleles[6],"CACA",2*2) == 0);
    assert(strncmp(d2->m->alleles[7],"ACCA",2*2) == 0);
    assert(strncmp(d2->m->alleles[8],"TGTG",2*2) == 0);
    assert(strncmp(d2->m->alleles[9],"TGGT",2*2) == 0);
    assert(strncmp(d2->m->alleles[10],"GCGC",2*2) == 0);
    assert(strncmp(d2->m->alleles[11],"GCGC",2*2) == 0);
    assert(strncmp(d2->m->alleles[12],"ATAT",2*2) == 0);
    assert(strncmp(d2->m->alleles[13],"TAAT",2*2) == 0);
    assert(strncmp(d2->m->alleles[14],"CGCG",2*2) == 0);
    assert(strcmp(d2->m->names[0],"g1") == 0);
    assert(strcmp(d2->m->names[1],"g2") == 0);
    assert(strcmp(d2->m->names[2],"oog3") == 0); // we'll just trust rest of the names
    for (int i = 0; i < 15; ++i) {
        assert(d2->m->ids[i].id == (cid - 14) + i);
        assert(d2->m->groups[i].num == 1);
        assert(d2->m->pedigrees[0][i].id == 0);
        assert(d2->m->pedigrees[1][i].id == 0);
    }
    delete_group(d2,(GroupNum){.num=1});
	remove(filename);
	
    delete_simdata(d2);

    delete_simdata(d);
    return 0;
}

struct MultiIDSet just_load(SimData* d) {
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

    struct MultiIDSet gande = load_data_files(d, "a-test.txt", "a-test-map.txt", "a-test-eff.txt",DETECT_FILE_FORMAT);
    //remove("a-test.txt");
    //remove("a-test-map.txt");
    //remove("a-test-eff.txt");
    return gande;
}


GroupNum test_loaders(SimData* d) {
    struct MultiIDSet loaded = just_load(d);
    GroupNum g0 = loaded.group;

    assert(d->genome.n_markers == 3); // all markers loaded
    assert(strcmp(d->genome.marker_names[0], "m1") == 0); // all markers ordered right
    assert(strcmp(d->genome.marker_names[1], "m2") == 0);
    assert(strcmp(d->genome.marker_names[2], "am3") == 0);
    assert(strcmp(*d->genome.names_alphabetical[0], "am3") == 0); // all markers ordered right
    assert(strcmp(*d->genome.names_alphabetical[1], "m1") == 0);
    assert(strcmp(*d->genome.names_alphabetical[2], "m2") == 0);
    printf("...genome map loaded correctly\n");

    assert(d->genome.n_maps == 1);
    assert(d->genome.map_ids[0].id == loaded.map.id);
    assert(d->genome.maps[0].n_chr == 2); // correct number of chromosomes
    assert(d->genome.maps[0].chrs[0].type == GSC_LINKAGEGROUP_SIMPLE);
    assert(d->genome.maps[0].chrs[0].map.simple.n_markers == 2);
    assert(d->genome.maps[0].chrs[0].map.simple.first_marker_index == 0);
    assert(fabs(d->genome.maps[0].chrs[0].map.simple.expected_n_crossovers - (3.1/100)) < TOL);
    assert(fabs(d->genome.maps[0].chrs[0].map.simple.dists[0] - 0) < TOL);
    assert(fabs(d->genome.maps[0].chrs[0].map.simple.dists[1] - 1) < TOL);
    assert(d->genome.maps[0].chrs[1].type == GSC_LINKAGEGROUP_SIMPLE);
    assert(d->genome.maps[0].chrs[1].map.simple.n_markers == 1);
    assert(d->genome.maps[0].chrs[1].map.simple.first_marker_index == 2);
    // assert(d->genome.maps[0].chrs[1].map.simple.expected_n_crossovers // length doesn't matter
    // assert(d->genome.maps[0].chrs[1].map.simple.dists[0]
    printf("...recombination map loaded correctly\n");

    MapID m2 = load_mapfile(d,"2marker-map.txt");
    assert(d->genome.n_maps == 2);
    assert(d->genome.map_ids[0].id == loaded.map.id);
    assert(d->genome.map_ids[1].id == m2.id);
    assert(loaded.map.id != m2.id);
    assert(d->genome.maps[0].n_chr == 2); // correct number of chromosomes
    assert(d->genome.maps[0].chrs[0].type == GSC_LINKAGEGROUP_SIMPLE);
    assert(d->genome.maps[0].chrs[0].map.simple.n_markers == 2);
    assert(d->genome.maps[0].chrs[0].map.simple.first_marker_index == 0);
    assert(fabs(d->genome.maps[0].chrs[0].map.simple.expected_n_crossovers - (3.1/100)) < TOL);
    assert(fabs(d->genome.maps[0].chrs[0].map.simple.dists[0] - 0) < TOL);
    assert(fabs(d->genome.maps[0].chrs[0].map.simple.dists[1] - 1) < TOL);
    assert(d->genome.maps[0].chrs[1].type == GSC_LINKAGEGROUP_SIMPLE);
    assert(d->genome.maps[0].chrs[1].map.simple.n_markers == 1);
    assert(d->genome.maps[0].chrs[1].map.simple.first_marker_index == 2);
    assert(d->genome.maps[1].n_chr == 1); // correct number of chromosomes
    assert(d->genome.maps[1].chrs[0].type == GSC_LINKAGEGROUP_SIMPLE);
    assert(d->genome.maps[1].chrs[0].map.simple.n_markers == 2);
    assert(d->genome.maps[1].chrs[0].map.simple.first_marker_index == 0);
    assert(fabs(d->genome.maps[1].chrs[0].map.simple.expected_n_crossovers - 0.1) < TOL);
    assert(fabs(d->genome.maps[1].chrs[0].map.simple.dists[0] - 0) < TOL);
    assert(fabs(d->genome.maps[1].chrs[0].map.simple.dists[1] - 1) < TOL);
    printf("...second recombination map loaded correctly\n");
    delete_recombination_map(d,m2);

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

    assert(load_effectfile(d, "a-test-eff2.txt").id==2);
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


    GroupNum g1 = load_genotypefile(d, "a-test.txt",DETECT_FILE_FORMAT);

	assert(d->m != NULL);
    assert(d->current_id.id == 12);
	assert(d->m->n_markers == 3);
	assert(d->m->n_genotypes == 12);
    assert(d->m->groups[0].num == g0.num);
    assert(d->m->groups[1].num == g0.num);
    assert(d->m->groups[2].num == g0.num);
    assert(d->m->groups[3].num == g0.num);
    assert(d->m->groups[4].num == g0.num);
    assert(d->m->groups[5].num == g0.num);
    assert(g0.num == 1);
    assert(d->m->groups[6].num == g1.num);
    assert(d->m->groups[7].num == g1.num);
    assert(d->m->groups[8].num == g1.num);
    assert(d->m->groups[9].num == g1.num);
    assert(d->m->groups[10].num == g1.num);
    assert(d->m->groups[11].num == g1.num);
    assert(g1.num == 2);
    assert(strncmp(d->m->alleles[0],"TTAATT", 6) == 0); // G01
    assert(strncmp(d->m->alleles[1],"TTAATT", 6) == 0); // G02
    assert(strncmp(d->m->alleles[2],"TTAATA", 6) == 0); // G03
    assert(strncmp(d->m->alleles[3],"TAAATA", 6) == 0); // G04
    assert(strncmp(d->m->alleles[4],"TTTTTT", 6) == 0); // G05
    assert(strncmp(d->m->alleles[5],"ATAATT", 6) == 0); // G06
    assert(strncmp(d->m->alleles[6],"TTAATT", 6) == 0); // G01
    assert(strncmp(d->m->alleles[7],"TTAATT", 6) == 0); // G02
    assert(strncmp(d->m->alleles[8],"TTAATA", 6) == 0); // G03
    assert(strncmp(d->m->alleles[9],"TAAATA", 6) == 0); // G04
    assert(strncmp(d->m->alleles[10],"TTTTTT", 6) == 0); // G05
    assert(strncmp(d->m->alleles[11],"ATAATT", 6) == 0); // G06
    printf("...second set of genotypes loaded correctly\n");

    delete_group(d, g1);

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
                    .will_save_bvs_to_file = NO_EFFECTSET,
                    .will_save_pedigree_to_file = GSC_FALSE,
                    .will_save_to_simdata = GSC_TRUE};
    assert(d->n_groups == 1);
    GroupNum f1 = make_random_crosses(d,g0,4,0,NO_MAP,g);
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
                    .will_save_bvs_to_file = NO_EFFECTSET,
                    .will_save_pedigree_to_file = GSC_FALSE,
                    .will_save_to_simdata = GSC_TRUE};
    GroupNum f1 = make_random_crosses(d, g0, 3, 1, NO_MAP, g);
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
    GroupNum fhs = make_targeted_crosses(d,4,combinations[0], combinations[1],NO_MAP,NO_MAP,g);
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
    EffectID e3 = load_effectfile(d,"a-test-eff3.txt");
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
                    .will_save_bvs_to_file={.id=1},
                    .will_save_alleles_to_file=GSC_TRUE,
                    .will_save_to_simdata=GSC_TRUE};
	//AlleleMatrix* a = make_all_unidirectional_crosses(&sd, 0, g);
	//sd.m->next_gen = a;
    GroupNum g1 = make_all_unidirectional_crosses(d, g0, NO_MAP, g);

    assert(g1.num != g0.num);
	assert(d->m->n_genotypes == 21);
	assert(d->m->n_markers == 3);
	//assert(strcmp(sd.m->name, "F1g") == 0);
    assert(strcmp(d->m->names[6], "F113") == 0);
    assert(strcmp(d->m->names[7], "F114") == 0);
    assert(strcmp(d->m->names[8], "F115") == 0);
    assert(strcmp(d->m->names[9], "F116") == 0);
    assert(strcmp(d->m->names[10], "F117") == 0);
    assert(strcmp(d->m->names[11], "F118") == 0);
    assert(strcmp(d->m->names[12], "F119") == 0);
    assert(strcmp(d->m->names[13], "F120") == 0);
    assert(strcmp(d->m->names[14], "F121") == 0);
    assert(strcmp(d->m->names[15], "F122") == 0);
    assert(strcmp(d->m->names[16], "F123") == 0);
    assert(strcmp(d->m->names[17], "F124") == 0);
    assert(strcmp(d->m->names[18], "F125") == 0);
    assert(strcmp(d->m->names[19], "F126") == 0);
    assert(strcmp(d->m->names[20], "F127") == 0);
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
    GroupNum bp = make_double_crosses_from_file(d, fname, NO_MAP, NO_MAP, g2);
    assert(d->m->pedigrees[0][21].id == d->m->ids[6].id && d->m->pedigrees[1][21].id == 23);
    assert(d->m->pedigrees[0][23].id == d->m->ids[7].id && d->m->pedigrees[1][23].id == 27);
    assert(d->m->pedigrees[0][25].id == d->m->ids[20].id && d->m->pedigrees[1][25].id == 15);
    assert(d->m->pedigrees[0][22].id == 13 && d->m->pedigrees[1][22].id == 23);
    assert(d->m->pedigrees[0][24].id == 14 && d->m->pedigrees[1][24].id == 27);
    assert(d->m->pedigrees[0][26].id == 27 && d->m->pedigrees[1][26].id == 15);
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
    GroupNum g1selfed = self_n_times(d, 5, g1, NO_MAP, opt);
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
    GroupNum g1dhap = make_doubled_haploids(d, g1, (MapID){.id=1}, opt);
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
    GroupNum gap = make_random_crosses(d , g1, 5, 0, NO_MAP, opt);
    GroupNum newparents[2];
    newparents[0] = g1clones;
    newparents[1] = make_random_crosses(d , g1, 5, 0, NO_MAP, opt);
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
    GroupNum g2 = make_random_crosses( d , g1, 4, 0, NO_MAP, gopt);
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
    GroupNum g3 = make_random_crosses_between( d, g1, g2, 3, 0, 0, NO_MAP, NO_MAP, gopt);

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
    GroupNum g4 = make_random_crosses_between( d, g1, tempg, g4sizetobe, 1, 0, NO_MAP, NO_MAP, gopt );
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

	printf("...SimData cleared correctly\n");

	return 0;
}

int test_block_generator(SimData *d) {
    MarkerBlocks b = create_evenlength_blocks_each_chr(d, NO_MAP, 2);

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

    b = create_evenlength_blocks_each_chr(d, NO_MAP, 4);

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
        .will_save_bvs_to_file = NO_EFFECTSET,
        .will_save_alleles_to_file = GSC_FALSE,
        .will_save_to_simdata = GSC_TRUE
    };
    int combos[2][3];
    combos[0][0] = 0; combos[1][0] = 0;
    combos[0][1] = 1; combos[1][1] = 2;
    combos[0][2] = 1; combos[1][2] = 5;
    GroupNum f1 = make_targeted_crosses(d, 3, combos[0], combos[1], NO_MAP, NO_MAP, g);
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

    delete_bidirectional_iter(&it2);

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
    
    gsc_delete_randomaccess_iter(&it5);
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
        .will_save_bvs_to_file = NO_EFFECTSET,
        .will_save_alleles_to_file = GSC_FALSE,
        .will_save_to_simdata = GSC_TRUE
    };
    GroupNum f1 = make_targeted_crosses(d,3,combos[0],combos[1], NO_MAP, NO_MAP, g);

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
    struct MultiIDSet loaded = just_load(d);
    GroupNum g0 = loaded.group;
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
    gsc_delete_randomaccess_iter(&it);

    // test of generate_gamete
    // with a totally homozygous parent eg G01, gives the expected result eg (TTA)
    strcpy(expected,"T?A?T?");
    generate_gamete(d,get_alleles(homozygote),gamete,0);
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
        generate_gamete(d,get_alleles(heterozygote),gamete,0);
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
    generate_doubled_haploid(d,get_alleles(homozygote),gamete,0);
    assert(strncmp(gamete,expected,sizeof(char)*GENOMESIZE) == 0);
    assert(gamete[GENOMESIZE] == 0);

    // With a heterozygous parent, can give both expected results.
    strcpy(expected,"AAAATT");
    strcpy(expected2,"TTAATT");
    memset(gamete,0,sizeof(char)*10);
    e1count = 0; e2count = 0;
    while (1) {
        generate_doubled_haploid(d,get_alleles(heterozygote),gamete,0);
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

    fclose(fp);

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

    fclose(fp);

    if (c1 == EOF && i == 0 && row == expectedNRows) return 0;
    else return -1;
}

/* main, for testing. Only uses a small dataset. */
int main(int argc, char* argv[]) {
    unsigned int randomSeed = 123123123; //time(NULL);
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
    test_genoloaders2();
    test_maploaders2();
    test_effloaders2();
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

    delete_simdata(d);

    printf("\n------- All tests passed. -------\n");


	//Small test that it can load a larger file alright.
    d = create_empty_simdata(randomSeed);
    // struct MultiIDSet init =
            load_data_files(d, "./gt_parents_mr2_50-trimto-5000.txt",
			 "./genetic-map_5112-trimto5000.txt",
             "./qtl_mr2.eff-processed.txt",DETECT_FILE_FORMAT);
    delete_simdata(d);

	printf("\nAll done\n");

	return 0;
}

Features and Concordance Guide 			{#concordance}
=================

Function names are not all consistent between R and C counterparts of genomicSimulation, particularly for functions introduced very early in the package's life. This is a table of which functions match up with which. 

This can also be interpreted as the list of user-facing functions in genomicSimulationC (though to use genomicSimulationC, you'll also need to consider the delete_ family and the structs in which data is stored).

| genomicSimulationC (C) | genomicSimulation (R) |
| ------------- | ------------- |
| load_all_simdata * | load.data * |
| load_more_transposed_genes_to_simdata | load.more.genotypes |
| load_effects_to_simdata | load.different.effects |
| delete_eff_set | delete.effect.set |
| calculate_minimum_bv | see.minimum.GEBV |
| calculate_optimum_bv | see.optimal.GEBV |
| calculate_optimal_alleles | see.optimal.haplotype |
| create_new_label | make.label |
| ~ ||
| make_all_unidirectional_crosses | cross.all.pairs |
| make_crosses_from_file | cross.combinations.file |
| make_double_crosses_from_file | cross.dc.combinations.file |
| make_doubled_haploids | make.doubled.haploids |
| make_clones | make.clones |
| self_n_times | self.n.times |
| cross_these_combinations | cross.combinations |
| cross_random_individuals | cross.randomly |
| cross_randomly_between | cross.randomly.between |
| ~ ||
| get_existing_group_counts | see.existing.groups  |
| get_group_names | see.group.data(data.type="N") |
| get_group_ids | see.group.data(data.type="D") |
| get_group_indexes | see.group.data(data.type="X") |
| get_group_genes | see.group.data(data.type="G") |
| get_group_bvs | see.group.data(data.type="B") |
| get_group_parent_names(parent=1) | see.group.data(data.type="P1") |
| get_group_parent_names(parent=2) | see.group.data(data.type="P2") |
| get_group_pedigrees | see.group.data(data.type="ped") |
| ~ ||
| set_names_to_values | change.names.to.values |
| set_labels_to_values | change.label.to.values |
| set_labels_to_const | change.label.to.this |
| increment_labels | change.label.by.amount |
| set_label_default | change.label.default |
| delete_label | delete.label | 
| ~ ||
| split_evenly_into_two | break.group.evenly(into.n=2) |
| split_evenly_into_n   | break.group.evenly |
| split_by_specific_counts_into_n  | break.group.into.buckets |
| split_randomly_into_two  | break.group.randomly(into.n=2) |
| split_randomly_into_n  | break.group.randomly |
| split_by_probabilities_into_n | break.group.with.probabilities |
| split_by_bv | select.by.gebv |
| split_into_families | break.group.into.families |
| split_into_halfsib_families | break.group.into.halfsib.families |
| split_into_individuals | break.group.into.individuals |
| split_from_group | make.group |
| split_by_label_value | make.group.from.label |
| split_by_label_range | make.group.from.label.range |
| combine_groups | combine.groups |
| delete_group | delete.group  |
| ~ ||
| save_count_matrix | save.allele.counts(group=NULL) |
| save_count_matrix_of_group | save.allele.counts |
| save_bvs | save.GEBVs(group=NULL) |
| save_group_bvs | save.GEBVs |
| save_simdata | save.genome.model |
| save_allele_matrix | save.genotypes(group=NULL, type="R") |
| save_transposed_allele_matrix | save.genotypes(group=NULL,  type="T") |
| save_group_alleles | save.genotypes(type="R") |
| save_transposed_group_alleles | save.genotypes(type="T") |
| create_n_blocks_by_chr; save_marker_blocks | save.local.GEBVs.by.chr |
| read_block_file; save_marker_blocks | save.local.GEBVs.by.file |
| save_full_pedigree | save.pedigrees(group=NULL, type="R") |
| save_one_step_pedigree | save.pedigrees(group=NULL, type="P") |
| save_group_full_pedigree | save.pedigrees(type="R") |
| save_group_one_step_pedigree | save.pedigrees(type="P") |
| calculate_recombinations_from_file | find.crossovers |
 
*: The R function does not directly call the C function, but they have the same effects.



Features and Concordance Guide 			{#concordance}
=================

Function names are not all consistent between R and C counterparts of genomicSimulation, particularly for functions introduced very early in the package's life. This is a table of which functions match up with which.

This can also be interpreted as the list of user-facing functions in genomicSimulationC (though to use genomicSimulationC, you'll also need to consider the delete_ family and the structs in which data is stored).

| genomicSimulationC (C) | genomicSimulation (R) |
| ------------- | ------------- |
| load_data_files | load.data |
| load_genotypefile | load.genotypes |
| load_mapfile | load.map |
| load_effectfile | load.effects |
| define_matrix_format_details | define.matrix.format.details |
| delete_group | delete.group |
| delete_recombination_map | delete.recombination.map |
| delete_eff_set | delete.effect.set |
| calculate_minimal_bv | see.minimal.GEBV |
| calculate_optimal_bv | see.optimal.GEBV |
| calculate_optimal_haplotype | see.optimal.haplotype |
| create_new_label | create.new.label |
| ~ ||
| make_all_unidirectional_crosses | make.all.unidirectional.crosses |
| make_crosses_from_file | make.crosses.from.file |
| make_double_crosses_from_file | make.double.crosses.from.file |
| make_doubled_haploids | make.doubled.haploids |
| make_clones | make.clones |
| self_n_times | self.n.times |
| make_targeted_crosses | make.targeted.crosses |
| make_random_crosses | make.random.crosses |
| make_random_crosses_between | make.random.crosses.between |
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
| change_names_to_values | change.names.to.values |
| change_allele_symbol | change.allele.symbol |
| change_labels_to_values | change.label.to.values |
| change_labels_to | change.label.to.this |
| change_labels_by_amount | change.label.by.amount |
| change_label_default | change.label.default |
| delete_label | delete.label |
| ~ ||
| split_evenly_into_two | break.group.evenly(into.n=2) |
| split_evenly_into_n   | break.group.evenly |
| split_into_buckets  | break.group.into.buckets |
| split_randomly_into_two  | break.group.randomly(into.n=2) |
| split_randomly_into_n  | break.group.randomly |
| split_by_probabilities | break.group.with.probabilities |
| split_by_bv | break.group.by.gebv |
| split_into_families | break.group.into.families |
| split_into_halfsib_families | break.group.into.halfsib.families |
| split_into_individuals | break.group.into.individuals |
| split_by_label_value | break.group.by.label.value |
| split_by_label_range | break.group.by.label.range |
| make_group_from | make.group |
| combine_groups | combine.groups |
| ~ ||
| save_allele_counts | save.allele.counts |
| save_bvs | save.GEBVs |
| save_genotypes | save.genotypes |
| save_pedigrees | save.pedigrees |
| create_evenlength_blocks_each_chr; calculate_group_local_bvs | save.local.GEBVs.blocks.from.chrsplit |
| read_block_file; calculate_group_local_bvs | save.local.GEBVs.blocks.from.file |
| calculate_recombinations_from_file | find.crossovers |

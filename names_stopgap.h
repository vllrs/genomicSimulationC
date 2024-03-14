/*
NAME UPDATE (between v0.2.4 and v0.2.5 of genomicSimulation)

- In this update, all functions have been prefixed with the module code 'gsc_'.
User-facing functions are also available without the module code
if the GSC_NO_SHORT_NAMES guard is not defined.

- Some functions, additionally, have had their names changed 
for better consistency between C and R versions of genomicSimulation.

Since this is a big interface change, the file you are reading has been provided
to help out/to patch over gaps.

Include this file:
	#include "names_stopgap.h" 
if you wish to still use the old names of functions, but also want to use 
genomicSimulation v>=0.2.5 . Cannot guarantee this will work perfectly.

	#define GSC_DEPRECATED_VERBOSE
Insert the above directive if you would like "Deprecated function call: [name]" 
printed to stderr on every call to a function in this file. To assist with
replacement of deprecated functions.
*/
#include "sim-operations.h"


#define DEPRECATION__PREFIXED(a) fprintf(stderr,"Deprecated function call: %s. Drop-in replacement: gsc_%s\n",a,a)
#define DEPRECATION__STATIC(a) fprintf(stderr,"Deprecated function call: %s. Function is now defined as internal/static. Contact maintainer if this is not suitable.\n",a)
#define DEPRECATION__CHANGED(a,b) fprintf(stderr,"Deprecated function call: %s. Replacement: (gsc_)%s\n",a,b)
#define DEPRECATION__FULL(a) fprintf(stderr,"Deprecated function call: %s. Function is now fully removed.\n",a)

// ---- Functions with changed names (for convergence with R version) ----
void set_label_default(SimData* d, const LabelID whichLabel, const int newDefault) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("set_label_default", "change_label_default");
	#endif
    gsc_change_label_default(d,whichLabel,newDefault);
}

void set_labels_to_const(SimData* d, const GroupNum whichGroup, const LabelID whichLabel, const int setTo) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("set_labels_to_const", "change_label_to");
	#endif
    gsc_change_label_to(d,whichGroup,whichLabel,setTo);
}

void increment_labels(SimData* d, const GroupNum whichGroup, const LabelID whichLabel, const int byValue) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("increment_labels", "change_label_by_amount");
	#endif
    gsc_change_label_by_amount(d,whichGroup,whichLabel,byValue);
}

void set_labels_to_values(SimData* d, const GroupNum whichGroup, const int startIndex, const LabelID whichLabel, const int n_values, const int values[n_values]) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("set_labels_to_values", "change_label_to_values");
	#endif	
    gsc_change_label_to_values(d,whichGroup,startIndex,whichLabel,n_values,values);
}

void set_names_to_values(SimData* d, const GroupNum whichGroup, const int startIndex, const int n_values, const char* values[n_values]) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("set_names_to_values", "change_names_to_values");
	#endif		
    gsc_change_names_to_values(d,whichGroup,startIndex,n_values,values);
}

void get_malloc(const size_t size) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("get_malloc", "gsc_malloc_wrap");
		DEPRECATION__STATIC("gsc_malloc_wrap");
	#endif
}

GroupNum load_transposed_genes_to_simdata(SimData* d, const char* filename) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("load_transposed_genes_to_simdata", "load_genotypes_transposed");
	#endif	
	return gsc_load_genotypes_transposed(d,filename);
}

GroupNum load_more_transposed_genes_to_simdata(SimData* d, const char* filename) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("load_more_transposed_genes_to_simdata", "load_genotypes_transposed");
	#endif	
	return gsc_load_genotypes_transposed(d,filename);
}

GroupNum load_transposed_encoded_genes_to_simdata(SimData* d, const char* filename) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("load_transposed_encoded_genes_to_simdata", "load_genotypes_encoded_and_transposed");
	#endif	
	return gsc_load_genotypes_encoded_and_transposed(d,filename);
}

void load_genmap_to_simdata(SimData* d, const char* filename) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("load_genmap_to_simdata", "load_genmap");
	#endif	
	gsc_load_genmap(d,filename);
}

EffectID load_effects_to_simdata(SimData* d, const char* filename) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("load_effects_to_simdata", "load_effects");
	#endif	
	return gsc_load_effects(d,filename);
}

struct GroupAndEffectSet load_all_simdata(SimData* d, const char* data_file, const char* map_file, const char* effect_file) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("load_all_simdata", "load_all_data");
	#endif	
	return gsc_load_all_data(d,data_file,map_file,effect_file);	
}

void save_group_alleles(FILE* f, SimData* d, const GroupNum group_id) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("save_group_alleles", "save_group_genotypes");
	#endif	
	gsc_save_group_genotypes(f,d,group_id);		
}

void save_transposed_group_alleles(FILE* f, const SimData* d, const GroupNum group_id) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("save_transposed_group_alleles", "save_transposed_group_genotypes");
	#endif	
	gsc_save_transposed_group_genotypes(f,d,group_id);		
}

void save_AM_pedigree(FILE* f, const AlleleMatrix* m, const SimData* parents) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("save_AM_pedigree", "save_allelematrix_full_pedigree");
	#endif	
	gsc_save_allelematrix_full_pedigree(f,m,parents);		
}

void save_count_matrix_of_group(FILE* f, const SimData* d, const char allele, const GroupNum group) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("save_count_matrix_of_group", "save_group_count_matrix");
	#endif	
	gsc_save_group_count_matrix(f,d,allele,group);
}

GroupNum split_from_group( SimData* d, const int n, const size_t indexes_to_split[n]) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("split_from_group", "make_group_from");
	#endif	
	return gsc_make_group_from(d,n,indexes_to_split);
}

unsigned int _split_by_somequality( SimData* d, const GroupNum group_id,
        void* somequality_data,
        GroupNum (*somequality_tester)(GenoLocation, void*, unsigned int, unsigned int, GroupNum*),
        size_t maxentries_results, GroupNum* results) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("_split_by_somequality", "scaffold_split_by_somequality");
	#endif	
	return gsc_scaffold_split_by_somequality(d,group_id,somequality_data,somequality_tester,maxentries_results,results);
}

unsigned int _split_by_someallocation( SimData* d, const GroupNum group_id, void* someallocator_data,
        GroupNum (*someallocator)(GenoLocation, SimData*, void*, unsigned int, unsigned int*, GroupNum*),
        size_t n_outgroups, GroupNum outgroups[n_outgroups]) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("_split_by_someallocation", "scaffold_split_by_someallocation");
	#endif	
	return gsc_scaffold_split_by_someallocation(d,group_id,someallocator_data,someallocator,n_outgroups,outgroups);
}

unsigned int split_by_specific_counts_into_n(SimData* d, const GroupNum group_id, const int n, const int* counts, GroupNum* results) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("split_by_specific_counts_into_n", "split_into_buckets");
	#endif	
	return gsc_split_into_buckets(d,group_id,n,counts,results);
}

unsigned int split_by_probabilities_into_n(SimData* d, const GroupNum group_id, const int n, const double* probs, GroupNum* results) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("split_by_probabilities_into_n", "split_by_probabilities");
	#endif	
	return gsc_split_by_probabilities(d,group_id,n,probs,results);
}

void generate_cross(SimData* d, const char* parent1_genome, const char* parent2_genome, char* output) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("generate_cross", "generate_gamete (called twice)");
	#endif
	gsc_generate_gamete(d,parent1_genome,output);
	gsc_generate_gamete(d,parent2_genome,output+1);
}

GroupNum _make_new_genotypes(SimData* d, const GenOptions g,
        void* parentIterator, void* datastore,
        int (*parentChooser)(void*, void*, unsigned int*, GenoLocation[static 2]),
        void (*offspringGenerator)(SimData*, void*, GenoLocation[static 2], GenoLocation) ) { 
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("_make_new_genotypes", "scaffold_make_new_genotypes");
	#endif
	return gsc_scaffold_make_new_genotypes(d,g,parentIterator,datastore,parentChooser,offspringGenerator);
}

GroupNum cross_random_individuals(SimData* d, const GroupNum from_group, const int n_crosses, const int cap, const GenOptions g) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("cross_random_individuals", "make_random_crosses");
	#endif	
	return gsc_make_random_crosses(d,from_group,n_crosses,cap,g);
}

GroupNum cross_randomly_between(SimData*d, const GroupNum group1, const GroupNum group2, const int n_crosses, const int cap1, const int cap2, const GenOptions g) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("cross_randomly_between", "make_random_crosses_between");
	#endif	
	return gsc_make_random_crosses_between(d,group1,group2,n_crosses,cap1,cap2,g);
}

unsigned int _specialised_random_draw(SimData* d, unsigned int max, unsigned int cap, unsigned int* member_uses, unsigned int noCollision) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("_specialised_random_draw", "randomdraw_replacementrules");
	#endif	
	return gsc_randomdraw_replacementrules(d,max,cap,member_uses,noCollision);
}

GroupNum cross_these_combinations(SimData* d, const int n_combinations, const int* firstParents, const int* secondParents, const GenOptions g) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("cross_these_combinations", "make_targeted_crosses");
	#endif	
	return gsc_make_targeted_crosses(d,n_combinations,firstParents,secondParents,g);
}

int calculate_group_count_matrix_of_allele( const SimData* d, const GroupNum group, const char allele, DecimalMatrix* counts) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("calculate_group_count_matrix_of_allele", "calculate_group_count_matrix");
	#endif
	return gsc_calculate_group_count_matrix(d,group,allele,counts);
}

int calculate_count_matrix_of_allele( const AlleleMatrix* m, const char allele, DecimalMatrix* counts) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("calculate_count_matrix_of_allele", "calculate_count_matrix");
	#endif
	return gsc_calculate_count_matrix(m,allele,counts);
}

DecimalMatrix calculate_full_count_matrix_of_allele( const AlleleMatrix* m, const char allele) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("calculate_full_count_matrix_of_allele", "calculate_full_count_matrix");
	#endif
	return gsc_calculate_full_count_matrix(m,allele);
}

int calculate_group_doublecount_matrix_of_allele( const SimData* d, const GroupNum group, const char allele, DecimalMatrix* counts, const char allele2, DecimalMatrix* counts2) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("calculate_group_doublecount_matrix_of_allele", "calculate_group_count_matrix_pair");
	#endif
	return gsc_calculate_group_count_matrix_pair(d,group,allele,counts,allele2,counts2);
}

int calculate_doublecount_matrix_of_allele( const AlleleMatrix* m , const char allele, DecimalMatrix* counts, const char allele2, DecimalMatrix* counts2) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("calculate_doublecount_matrix_of_allele", "calculate_count_matrix_pair");
	#endif
	return gsc_calculate_count_matrix_pair(m,allele,counts,allele2,counts2);
}

MarkerBlocks create_n_blocks_by_chr(const SimData* d, const int n) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("create_n_blocks_by_chr", "create_evenlength_blocks_each_chr");
	#endif
	return gsc_create_evenlength_blocks_each_chr(d,n);
}

MarkerBlocks read_block_file(const SimData* d, const char* block_file) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("read_block_file", "load_blocks");
	#endif
	return gsc_load_blocks(d,block_file);
}

char* calculate_optimal_alleles(const SimData* d, const EffectID effID) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("calculate_optimal_alleles", "calculate_optimal_haplotype");
	#endif
	return gsc_calculate_optimal_haplotype(d,effID);
}

char* calculate_optimal_available_alleles(const SimData* d, const GroupNum group, const EffectID effID) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("calculate_optimal_available_alleles", "calculate_optimal_possible_haplotype");
	#endif
	return gsc_calculate_optimal_possible_haplotype(d,group,effID);
}

double calculate_optimum_bv(const SimData* d, const EffectID effID) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("calculate_optimum_bv", "calculate_optimal_bv");
	#endif
	return gsc_calculate_optimal_bv(d,effID);
}

double calculate_optimal_available_bv(const SimData* d, const GroupNum group, const EffectID effID) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("calculate_optimal_available_bv", "calculate_optimal_possible_bv");
	#endif
	return gsc_calculate_optimal_possible_bv(d,group,effID);
}

double calculate_minimum_bv(const SimData* d, const EffectID effID) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__CHANGED("calculate_minimum_bv", "calculate_minimal_bv");
	#endif
	return gsc_calculate_minimal_bv(d,effID);
}

// ---- Functions that are now static ----
void get_sorted_markers(SimData* d, int actual_n_markers) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__STATIC("get_sorted_markers");
	#endif		
	//gsc_get_sorted_markers(d,actual_n_markers);
}

void get_chromosome_locations(SimData* d) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__STATIC("get_chromosome_locations");
	#endif		
    //gsc_get_chromosome_locations(d);
}

void set_names(AlleleMatrix* a, const char* prefix, const int suffix, const int from_index) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__STATIC("set_names");
	#endif	
	//gsc_set_names(a,prefix,suffix,from_index);
}

void set_ids(SimData* d, const int from_index, const int to_index) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__STATIC("set_ids");
	#endif	
	//gsc_set_ids(d,from_index,to_index);
}


// ---- Other functions (gsc_ prefix added and no short name option available ----
int randpoi(rnd_pcg_t* rng, double lambda) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("randpoi");
	#endif
	return gsc_randpoi(rng,lambda);
}

DecimalMatrix generate_zero_dmatrix(const int r, const int c) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("generate_zero_dmatrix");
	#endif
	return gsc_generate_zero_dmatrix(r,c);
}

int add_matrixvector_product_to_dmatrix(DecimalMatrix* result, const DecimalMatrix* a, const double* b) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("add_matrixvector_product_to_dmatrix");
	#endif
	return gsc_add_matrixvector_product_to_dmatrix(result,a,b);
}

int	add_doublematrixvector_product_to_dmatrix(DecimalMatrix* result, const DecimalMatrix* amat, const double* avec, const DecimalMatrix* bmat, const double* bvec) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("add_doublematrixvector_product_to_dmatrix");
	#endif
	return add_doublematrixvector_product_to_dmatrix(result,amat,avec,bmat,bvec);
}

struct TableSize get_file_dimensions(const char* filename, const char sep) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_file_dimensions");
	#endif	
	return gsc_get_file_dimensions(filename, sep);
}

int get_from_ordered_uint_list(const unsigned int target, const unsigned int listLen, const unsigned int list[listLen]) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_from_ordered_uint_list");
	#endif	
	return gsc_get_from_ordered_uint_list(target,listLen,list);
}

int get_from_ordered_pedigree_list(const PedigreeID target, const unsigned int listLen, const PedigreeID list[listLen]) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_from_ordered_pedigree_list");
	#endif	
	return gsc_get_from_ordered_pedigree_list(target,listLen,list);
}

int get_from_unordered_str_list(const char* target, const int listLen, const char* list[listLen]) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_from_unordered_str_list");
	#endif	
	return gsc_get_from_unordered_str_list(target,listLen,list);
}

void shuffle_up_to(rnd_pcg_t* rng, size_t* sequence, const size_t total_n, const size_t n_to_shuffle) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("shuffle_up_to");
	#endif		
	gsc_shuffle_up_to(rng, sequence, total_n, n_to_shuffle);
}

int get_integer_digits(const int i) { 
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_integer_digits");
	#endif	
	return gsc_get_integer_digits(i);
}

int get_index_of_label( const SimData* d, const LabelID label ) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_index_of_label");
	#endif		
	return gsc_get_index_of_label(d, label);
}

int get_index_of_eff_set( const SimData* d, const EffectID eff_set_id ) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_index_of_eff_set");
	#endif		
	return gsc_get_index_of_eff_set(d, eff_set_id);
}

LabelID get_new_label_id( const SimData* d ) { 
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_new_label_id");
	#endif	
	return gsc_get_new_label_id(d);
}

EffectID get_new_eff_set_id( const SimData* d ) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_new_eff_set_id");
	#endif	
	return gsc_get_new_eff_set_id(d);
}

GroupNum get_next_free_group_num( const int n_existing_groups, const GroupNum* existing_groups, int* cursor,  GroupNum previous) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_next_free_group_num");
	#endif		
	return gsc_get_next_free_group_num(n_existing_groups, existing_groups, cursor, previous);
}

GroupNum get_new_group_num( SimData* d ) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_new_group_num");
	#endif	
	return gsc_get_new_group_num(d);
}

void get_n_new_group_nums( SimData* d, const int n, GroupNum* result) { 
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_n_new_group_nums");
	#endif	
	gsc_get_n_new_group_nums(d,n,result);
}

void condense_allele_matrix( SimData* d) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("condense_allele_matrix");
	#endif	
	gsc_condense_allele_matrix(d);
}

AlleleMatrix* create_empty_allelematrix(const int n_markers, const int n_labels, const int labelDefaults[n_labels], const int n_genotypes) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("create_empty_allelematrix");
	#endif	
	return gsc_create_empty_allelematrix(n_markers,n_labels,labelDefaults,n_genotypes);
}

AlleleMatrix* get_nth_AlleleMatrix( AlleleMatrix* listStart, const unsigned int n) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_nth_AlleleMatrix");
	#endif	
	return gsc_get_nth_AlleleMatrix(listStart, n);
}

char* get_name_of_id( const AlleleMatrix* start, const PedigreeID id) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_name_of_id");
	#endif	
	return gsc_get_name_of_id(start, id);
}

int get_parents_of_id( const AlleleMatrix* start, const PedigreeID id, PedigreeID output[2]) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_parents_of_id");
	#endif	
	return gsc_get_parents_of_id(start, id,output);
}

void get_ids_of_names( const AlleleMatrix* start, const int n_names, const char* names[n_names], PedigreeID* output) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_ids_of_names");
	#endif	
	gsc_get_ids_of_names(start, n_names,names,output);
}

int get_index_of_child( const AlleleMatrix* start, const PedigreeID parent1id, const PedigreeID parent2id) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_index_of_child");
	#endif	
	return gsc_get_index_of_child(start, parent1id,parent2id);
}

int get_index_of_name( const AlleleMatrix* start, const char* name) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_index_of_name");
	#endif	
	return gsc_get_index_of_name(start, name);
}

PedigreeID get_id_of_index( const AlleleMatrix* start, const int index) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_id_of_index");
	#endif	
	return gsc_get_id_of_index(start,index);
}

char* get_genes_of_index( const AlleleMatrix* start, const int index) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("get_genes_of_index");
	#endif	
	return gsc_get_genes_of_index(start,index);	
}

void delete_allele_matrix(AlleleMatrix* m) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("delete_allele_matrix");
	#endif	
	gsc_delete_allele_matrix(m);	
}

void delete_effect_matrix(EffectMatrix* m) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("delete_effect_matrix");
	#endif	
	gsc_delete_effect_matrix(m);	
}

void save_parents_of(FILE* f, const AlleleMatrix* m, PedigreeID p1, PedigreeID p2) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("save_parents_of");
	#endif		
	gsc_save_parents_of(f,m,p1,p2);
}

void save_manual_bvs(FILE* f, const DecimalMatrix* e, const PedigreeID* ids, const char** names) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("save_manual_bvs");
	#endif
	gsc_save_manual_bvs(f,e,ids,names);
}

int* calculate_min_recombinations_fw1(SimData* d, char* parent1, unsigned int p1num, char* parent2, unsigned int p2num, char* offspring, int certain) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("calculate_min_recombinations_fw1");
	#endif
	return gsc_calculate_min_recombinations_fw1(d,parent1,p1num,parent2,p2num,offspring,certain);
}

int* calculate_min_recombinations_fwn(SimData* d, char* parent1, unsigned int p1num, char* parent2, unsigned int p2num, char* offspring, int window_size, int certain) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("calculate_min_recombinations_fwn");
	#endif
	return gsc_calculate_min_recombinations_fwn(d,parent1,p1num,parent2,p2num,offspring,window_size,certain);
}

static inline int has_same_alleles(const char* p1, const char* p2, const int i) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("has_same_alleles");
	#endif
	return (p1[i<<1] == p2[i<<1] || p1[(i<<1) + 1] == p2[i] || p1[i] == p2[(i<<1) + 1]);
}

static inline int has_same_alleles_window(const char* g1, const char* g2, const int start, const int w) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("has_same_alleles_window");
	#endif	
	int same = GSC_TRUE;
    int i;
    for (int j = 0; j < w; ++j) {
        i = start + j;
        same = same && (g1[i<<1] == g2[i<<1] || g1[(i<<1) + 1] == g2[i] || g1[i] == g2[(i<<1) + 1]);
    }
    return same;
}

int calculate_recombinations_from_file(SimData* d, const char* input_file, const char* output_file, int window_len, int certain) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__PREFIXED("calculate_recombinations_from_file");
	#endif
	return gsc_calculate_recombinations_from_file(d,input_file,output_file,window_len,certain);
}

// ------ Functions that are now completely removed (already deprecated, and not called anywhere) --------

void get_genes_of_id( const AlleleMatrix* start, const PedigreeID id) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__FULL("get_genes_of_id");
	#endif	
}

void get_id_of_child( const AlleleMatrix* start, const PedigreeID parent1id, const PedigreeID parent2id) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__FULL("get_id_of_child");
	#endif	
}

void save_simdata(FILE* f, const SimData* m) {
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__FULL("save_simdata");
	#endif	
}

void validate_bidirectional_cache(gsc_BidirectionalIterator* it) { 
	// could not be implemented without a heap of false positives and 
	// false negatives in the use cases I wanted to use an iterator for.
	// Would give false confidence or false worry, so instead it has been 
	// removed.
	#ifdef GSC_DEPRECATED_VERBOSE
		DEPRECATION__FULL("validate_bidirectional_cache");
	#endif	
}
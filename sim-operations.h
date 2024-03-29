#ifndef SIM_OPERATIONS_H
#define SIM_OPERATIONS_H
/* genomicSimulationC v0.2.4.004 - last edit 8 Mar 2024 */

#ifdef SIM_OPERATIONS
    #define RND_IMPLEMENTATION
#endif
#include "lib/rnd.h"
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#define PI 3.1415926535897932384626433832795028841971693993751
#define TRUE 1
#define FALSE 0
#define UNINITIALISED -1


 /* This section contains settings for the simulation that users can modify if they have the need.
 * To apply the modified settings the simulation tool must be re-compiled.
 * To modify a setting, replace the number (eg 1000) with the new value without modifying the name
 */

/** The largest contiguous block of memory that could be requested in the
 * process of simulation is CONTIG_WIDTH integers. This setting's default value
 * is 1000.
 *
 * This could be decreased to help long simulations or lower-end machines.
 * Increasing this may provide some speed gain.
 */
#define CONTIG_WIDTH 1000

 /** The maximum number of characters allowed in a name field.
 * These include names of SNPs, names of genotypes loaded from files,
 * names of generated genotypes, and save-as-you-go filenames. Default is 30.
 *
 * Increase this if there is a risk some names may be longer than
 * this value.
 */
#define NAME_LENGTH 45

/** @defgroup structs Data Structures
 *
 * How the simulation stores data.
 *
 * genomicSimulation is a state-based package. These data structures
 * are used to store the library's data/state. Many use dynamically
 * allocated memory so call the relevant delete_ function if one exists
 * when you are finished with them.
 *
 * SimData is the central state/data
 * storage struct, and so is a required parameter to most user-facing functions.
 * It contains pointers to a GeneticMap (storing the loaded genome map), an
 * EffectMatrix (storing the loaded allele effects), and an AlleleMatrix
 * (storing metadata and genotypes of founders and simulated offspring).
 *
 * Other structs in this group (TableSize, MarkerBlocks, GenOptions) represented
 * specific types of data and are used as parameters and return value of certain
 * functions.
 *
 * @{
 */

/** A struct representing a single marker location. the attribute
 * `chromosome` represents the chromosome number, and `position` the
 * position on the chromosome in centiMorgans.
 */
typedef struct {
	int chromosome; /**< The chromosome number */
	float position; /**< The distance in centiMorgans along the chromosome */
} MarkerPosition;

/** A simple struct used for returning the dimensions of a matrix or table.*/
struct TableSize {
	int num_columns;
	int num_rows;
};

/** A struct used to store a set of blocks of markers.
 */
typedef struct {
    /** The number of blocks whose details are stored here. */
	int num_blocks;

    /** Pointer to a heap array of length num_blocks
      * containing the number of markers that make up each block */
	int* num_markers_in_block;

    /** Pointer to a heap array of length num_blocks, each
     * entry in which is a pointer to a heap array with length corresponding to
     * the value of the corresponding entry in num_markers_in_block whose values
     * are the indexes in the SimData of the markers that make up that block. */
	int** markers_in_block;
} MarkerBlocks;

/** A row-major heap matrix that contains floating point numbers. `dmatrix` functions
 * are designed to act on this matrix.
 *
 * Rows make the first index of the matrix and columns the second.
 */
typedef struct {
	double** matrix; /**< The actual matrix and contents */
	int rows;        /**< Number of rows in the matrix */
	int cols;        /**< number of columns in the matrix */
} DecimalMatrix;


/** A type representing a program-lifetime-unique identifier for a genotype,
 *  to be used in tracking pedigree.
 */
typedef struct {
    unsigned int id;
} PedigreeID;
#define NO_PEDIGREE (PedigreeID){.id=0}

/** A type representing the identifier of a group of genotypes
 */
typedef struct {
    unsigned int num;
} GroupNum;
#define NO_GROUP (GroupNum){.num=0}

/** A type representing a particular loaded set of marker effects
 */
typedef struct {
    int id;
} EffectID;
#define NOT_AN_EFFECT_SET (EffectID){.id=0}

/** A type representing a particular integer label
 */
typedef struct {
    int id;
} LabelID;
#define NOT_A_LABEL (LabelID){.id=0}

struct GroupAndEffectSet {
    GroupNum group;
    EffectID effectSet;
};

/** A type that contains choices of settings for SimData functions that create a
 * new AlleleMatrix/generation.
 *
 * The family_size parameter will affect how many offspring are produced.
 *
 * The will_name_offspring, will_track_pedigree, and will_allocate_ids parameters
 * affect how much extra detail about the offspring is generated/saved.
 *
 * The will_save_to_simdata toggle allows you the option of generating offspring without
 * saving them in memory. This may be useful in combination with save-as-you-go toggles
 * will_save_pedigree_to_file, will_save_bvs_to_file, and will_save_alleles_to_file,
 * to generate a larger number of offspring than will fit in memory.
*/
typedef struct {
	int will_name_offspring; /**< A boolean: whether generated offspring should be given names. */
    const char* offspring_name_prefix; /**< If `will_name_offspring` is true, generated
                           * offspring are named [offspring_name_prefix][index]. */

	int family_size; /**< The number of offspring to produce from each cross.*/

	int will_track_pedigree; /**< A boolean: whether to track parentage of generated offspring.*/
	int will_allocate_ids; /**< A boolean: whether to allocate generated offspring session-
                            * unique IDs. IDs are used for pedigree tracking. The
                            * offspring of an anonymous individual (one without an ID)
                            * cannot identify that individual as their parent. */

    const char* filename_prefix; /**< A string used in save-as-you-go file names. */
	int will_save_pedigree_to_file; /**< A boolean. If true, the full/recursive
                            * pedigrees of every offspring generated in the cross
                            * are saved to "[filename_prefix}-pedigree.txt", even
                            * if `will_save_to_simdata` is false.
                            * Pedigrees are saved in the format of save_full_pedigree()*/
    EffectID will_save_bvs_to_file; /**< If equal to NOT_AN_EFFECT_SET, no bvs are calculated or saved.
                            * Otherwise, for each offspring in the cross,
                            * the breeding values according
                            * to the marker effect set with this EffectID
                            * are saved to "[filename_prefix}-bv.txt", even
                            * if `will_save_to_simdata` is false.
                            * BVs are saved in the format of save_bvs() */
	int will_save_alleles_to_file; /**< A boolean. If true, the set of alleles
                            * of every offspring generated in the cross
                            * are saved to "[filename_prefix}-genotype.txt", even
                            * if `will_save_to_simdata` is false.
                            * Genotypes are saved in the format of save_group_alleles()*/
	int will_save_to_simdata; /**< A boolean. If true, the generated offspring exist
                            * in the SimData struct after the function executes.
                            * If false, they are discarded after creation. */
} GenOptions;


/** A type that stores the genetic map for a set of markers.
 *
 * To get all markers belonging to a particular chromosome, use the following rule:
 * Chr n includes all markers in `positions` starting at index chr_ends[n-1] up
 * but not including the marker at index chr_ends[n]
 *
 * Chromosomes must be numbered. All chromosomes from 1 up to and including the
 * highest chromosome number found in the loaded map are represented in these
 * arrays.
 *
*/
typedef struct {
	int n_chr; /**< The number of chromosomes represented in the map. This
                * corresponds to the highest numbered chromosome with a tracked
                * marker (some chromosomes in between may be empty) */
	int* chr_ends; /**< An array of ints. The entry at index i is the index in
                * `positions` of the first marker that belongs to Chr(i + 1).
                * The array is n_chr + 1 integers long.*/
	float* chr_lengths; /**< An array of floats. The entry at index i is the length
                * of Chr(i + 1), calculated by
                * `position of last marker - position of first marker`.
                * The array is n_chr entries long. */

	MarkerPosition* positions; /**< An array of MarkerPositions, ordered from lowest to highest. */
} GeneticMap;

/** A linked list entry that stores a matrix of alleles for a set of SNP markers
 * and genotypes.
*/
typedef struct AlleleMatrix AlleleMatrix;
struct AlleleMatrix {
    /** A matrix of SNP markers by lines/genotypes containing pairs of alleles
     * eg TT, TA. Use `alleles[line index][marker index * 2]` to get the
     * first allele and `alleles[lines index][marker index * 2 + 1]` to
     * get the second. If CONTIG_WIDTH lines are saved here, another
     * AlleleMatrix is added to the linked list when there's a need to save more.*/
	char* alleles[CONTIG_WIDTH];

	int n_genotypes; /**< Number of genotypes currently loaded in this matrix.*/
	int n_markers; /**< Number of markers across which genotypes are tracked. This has
                    * redundancy with SimData and other members of its linked list
                    * but it's good to know how big your own `alleles` array is.*/

    char* names[CONTIG_WIDTH]; /**< Array of dynamically allocated strings
                    * containing the names of the lines/genotypes in this matrix.
                    * Guaranteed to be NULL if they do not have names. */
    PedigreeID ids[CONTIG_WIDTH]; /**< Unique ID for each genotype. */
    PedigreeID pedigrees[2][CONTIG_WIDTH]; /**< Two lists of integer IDs of the
                    * parents of this genotype (if tracked), or 0 if we don't know/care.*/
    GroupNum groups[CONTIG_WIDTH]; /**< Group allocation of each genotype. */

    int n_labels; /**< Number of custom labels currently available to this AlleleMatrix. This has
                    * redundancy with SimData and other members of its linked list
                    * but it's good to know how big your own `labels` array is.*/
    int** labels; /**< Pointer to list of labels. Size of first dimension is n_labels,
                           * of second dimension is arrays of labels of length CONTIG_WIDTH*/

	AlleleMatrix* next; /**< Pointer to the next AlleleMatrix in the linked list,
                         * or NULL if this entry is the last. */
};

/** A type that stores a matrix of effect values and their names.
 */
typedef struct {
	DecimalMatrix effects; /**< Effect on breeding value of alleles at markers.
        * Rows correspond to `effect_names`/alleles, columns to markers. */
	char* effect_names; /**< Character array containing allele characters ordered
        * to match rows of `effects`. */
} EffectMatrix;

/** Composite type that is used to run crossing simulations.
 *
 * The core of this type is a list of markers. These are used to index the rows
 * of the allele matrix and the position map, and the columns of the effect matrix.
 */
typedef struct {
	int n_markers;  /**< The number of markers/length of `markers`. */
	char** markers; /**< Array of strings containing the names of markers. */

    int n_labels; /**< The number of custom labels in the simulation.*/
    LabelID* label_ids; /**< The identifier number of each label in the simulation, in order
                     * of their lookup index. */
    int* label_defaults; /**< Array containing the default (birth) value of each
                          * custom label. */

	GeneticMap map; /**< A GeneticMap. If this is set, then `markers`
                     * will be ordered and all markers have a known position.*/
	AlleleMatrix* m; /**< Pointer to an AlleleMatrix, which stores data and
                      * metadata of founders and simulated offspring. The
                      * AlleleMatrix is start of a linked list if there are
                      * many genotypes. */

    int n_eff_sets; /**< The number of sets of allele effects in the simulation **/
    EffectID* eff_set_ids; /**< The identifier number of each set of allele effects in the simulation,
                     * ordered by their lookup index. */
    EffectMatrix* e; /**< Array of n_eff_sets EffectMatrix, optional for the use of the simulation.
                     * Used for calculating breeding values from which alleles
                     * a genotype has at each marker.*/

    rnd_pcg_t rng; /**< Random number generator working memory. */
    PedigreeID current_id; /**< Highest SimData-unique ID that has been generated
                              * so far. Used to track which IDs have already been
                              * given out.*/
    size_t n_groups; /**< Number of groups currently existing in simulation. It is
                        * guaranteed to never be less than the number of groups in simulation
                        * even if not perfectly accurate. */
} SimData;

extern const GenOptions BASIC_OPT;
/** @} */

/** @defgroup maths Mathematical functions
 *
 * For mathematical and statistical operations as required by the package.
 *
 * Includes matrix operations defined on a DecimalMatrix struct, and
 * draws from certain random distributions.
 *
 * @{
 */

int randpoi(rnd_pcg_t* rng, double lambda);

DecimalMatrix generate_zero_dmatrix(const int r, const int c);
int add_matrixvector_product_to_dmatrix(DecimalMatrix* result, const DecimalMatrix* a, const double* b);
int add_doublematrixvector_product_to_dmatrix(DecimalMatrix* result, const DecimalMatrix* amat, const double* avec,
                                              const DecimalMatrix* bmat, const double* bvec);

/** @} */

/** @defgroup supporters Utils/Supporting Functions
 *
 * @{
 */
struct TableSize get_file_dimensions(const char* filename, const char sep);
int get_from_ordered_uint_list(const unsigned int target, const unsigned int listLen, const unsigned int list[listLen]);
int get_from_ordered_pedigree_list(const PedigreeID target, const unsigned int listLen, const PedigreeID list[listLen]);
int get_from_unordered_str_list(const char* target, const int listLen, const char* list[listLen]);
void shuffle_up_to(rnd_pcg_t* rng, size_t* sequence, const size_t total_n, const size_t n_to_shuffle);

LabelID create_new_label(SimData* d, const int setTo);
void set_label_default(SimData* d, const LabelID whichLabel, const int newDefault);
void set_labels_to_const(SimData* d, const GroupNum whichGroup, const LabelID whichLabel, const int setTo);
void increment_labels(SimData* d, const GroupNum whichGroup, const LabelID whichLabel, const int byValue);
void set_labels_to_values(SimData* d, const GroupNum whichGroup, const int startIndex, const LabelID whichLabel,
                          const int n_values, const int values[n_values]);
void set_names_to_values(SimData* d, const GroupNum whichGroup, const int startIndex, const int n_values, const char* values[n_values]);

void get_sorted_markers(SimData* d, int actual_n_markers);
void get_chromosome_locations(SimData *d);

void set_names(AlleleMatrix* a, const char* prefix, const int suffix, const int from_index);
void set_ids(SimData* d, const int from_index, const int to_index);
int get_integer_digits(const int i);
int get_index_of_label( const SimData* d, const LabelID label );
int get_index_of_eff_set( const SimData* d, const EffectID eff_set_id );
LabelID get_new_label_id( const SimData* d );
EffectID get_new_eff_set_id( const SimData* d );
GroupNum get_next_free_group_num( const int n_existing_groups, const GroupNum* existing_groups, int* cursor,  GroupNum previous);
GroupNum get_new_group_num( SimData* d );
void get_n_new_group_nums( SimData* d, const int n, GroupNum* result);
void condense_allele_matrix( SimData* d);
void* get_malloc(const size_t size);

int _simdata_pos_compare(const void *pp0, const void *pp1);
int _descending_double_comparer(const void* pp0, const void* pp1);
int _ascending_double_comparer(const void* pp0, const void* pp1);
int _ascending_float_comparer(const void* p0, const void* p1);
int _ascending_index_comparer(const void* p0, const void* p1);
int _ascending_GroupNum_comparer(const void* p0, const void* p1);
/**@}*/


/** @defgroup loaders Setup Functions
 *
 * For setup of the simulation (loading founders, genetic maps, and optionally allele effects).
 *
 * @{
 */
AlleleMatrix* create_empty_allelematrix(const int n_markers, const int n_labels, const int labelDefaults[n_labels], const int n_genotypes);
SimData* create_empty_simdata(RND_U32 RNGseed);
void clear_simdata(SimData* d);

GroupNum load_transposed_genes_to_simdata(SimData* d, const char* filename);
GroupNum load_more_transposed_genes_to_simdata(SimData* d, const char* filename);
GroupNum load_genes_to_simdata(SimData* d, const char* filename); //@ add
GroupNum load_more_genes_to_simdata(SimData* d, const char* filename); //@ add
GroupNum load_transposed_encoded_genes_to_simdata(SimData* d, const char* filename);
void load_genmap_to_simdata(SimData* d, const char* filename);
EffectID load_effects_to_simdata(SimData* d, const char* filename);
struct GroupAndEffectSet load_all_simdata(SimData* d, const char* data_file, const char* map_file, const char* effect_file);
/** @} */


/** @defgroup getters Data Access and Search Functions
 *
 * For non-persistent access to simulation results and contents.
 *
 * @{
 */

    /** @defgroup iterators Genotype Iterators
     *
     * For iterating through the genotypes in the simulation.
     * It is possible to iterate through one group or through
     * every genotype in the simulation.
     *
     * @{
     */

/** An AlleleMatrix/AlleleMatrix index coordinate of a particular
 *  genotype in the simulation. To be used to look up details of
 *  that genotype using the `get_` family of functions.
*/
typedef struct {
    AlleleMatrix* localAM; /**< Pointer to the AlleleMatrix in which
                            * the genotype can be found. */
    int localPos; /**< Index in the localAM where the genotype can be
                   * found (min value: 0. Max value: CONTIG_WIDTH-1). */
} GenoLocation;

/** Constant representing a nonexistent location in the simulation.
 */
//extern const GenoLocation INVALID_GENO_LOCATION;
#define INVALID_GENO_LOCATION (GenoLocation){.localAM=0,.localPos=-1}
#define IS_VALID_LOCATION(g) (g.localAM != 0 && g.localPos != -1)

/** Identify whether a GenoLocation is INVALID_GENO_LOCATION
 *
 * @param g location to check.
 * @return FALSE if g has either of the attributes of
 * INVALID_GENO_LOCATION, TRUE otherwise
 */
static inline int isValidLocation(const GenoLocation g) {
    // Either entry of INVALID_GENO_LOCATION is inappropriate in a valid GenoLocation
    return (g.localAM != INVALID_GENO_LOCATION.localAM &&
            g.localPos != INVALID_GENO_LOCATION.localPos);
}

/** A structure to iterate forwards and backwards through all
 *  genotypes in a SimData or through only the members of a group.
 *  @see create_bidirectional_iter
 */
typedef struct {
    SimData* d; /**< Simulation data through which to iterate */
    const GroupNum group; /**< Group through which to iterate. If it is 0,
                          * then iterate through all genotypes in the simulation.
                          * Otherwise, iterate through members of the group with
                          * this as their group number. */
    size_t localPos; /**< Local index (index within the cachedAM) of the genotype in the linked list
                       * of AlleleMatrix beginning at `d->m` where the
                       * iterator's 'cursor' currently sits. */

    AlleleMatrix* cachedAM; /**< Pointer to the AlleleMatrix from the linked list
                              * of AlleleMatrix beginning at `d->m` where the
                              * iterator's 'cursor' currently sits. Contains
                              * the genotype at `localPos`. */
    unsigned int cachedAMIndex; /**< Index of `cachedAM` in the linked list of
                                  * AlleleMatrix beginning at `d->m`. `d->m`
                                  * is considered to be index 0. */

    char atEnd; /**< Boolean that is TRUE if the iterator's 'cursor' is on
                  * the last genotype (genotype with the highest index in the
                  * SimData) that fulfils the `group` critera of this iterator. */
    char atStart; /**< Boolean that is TRUE if the iterator's 'cursor' is on
                    * the first genotype (genotype with the lowest index in the
                    * SimData) that fulfils the `group` critera of this iterator. */

} BidirectionalIterator;

/** A structure to iterate forwards through all
 *  positions in the AlleleMatrix linked list in SimData. Used
 *  in @see condense_allele_matrix. Internal, not recommended for end users.
 */
struct GappyIterator {
    GenoLocation cursor;
    unsigned int cursorAMIndex;
};

/** A structure to search and cache indexes of all
 *  genotypes in a SimData or of all the members of a group.
 *  @see create_randomaccess_iter
 */
typedef struct {
    SimData* d; /**< Simulation data through which to iterate */
    const GroupNum group; /**< Group through which to iterate. If it is 0,
                          * then iterate through all genotypes in the simulation.
                          * Otherwise, iterate through members of the group with
                          * this as their group number. */

    unsigned int cacheSize; /**< Length in GenoLocations of `cache` */
    GenoLocation* cache; /**< Array iteratively updated with the known
                           * genotypes in the simulation that fulfil the
                           * `group` criteria of the iterator as they
                           * are discovered during calls to next_ functions */

    int largestCached; /**< Local/group index (that is, index in `cache`) of the
                         * highest cell in `cache` that has been filled. */
    int groupSize; /**< If the number of genotypes in the simulation that fulfil
                     * the iterator's `group` criteria is known, it is saved here.
                     * This value is left uninitialised until then. */
} RandomAccessIterator;

BidirectionalIterator create_bidirectional_iter( SimData* d, const GroupNum group);
RandomAccessIterator create_randomaccess_iter( SimData* d, const GroupNum group);

AlleleMatrix* get_nth_AlleleMatrix( AlleleMatrix* listStart, const unsigned int n);

GenoLocation set_bidirectional_iter_to_start(BidirectionalIterator* it);
GenoLocation set_bidirectional_iter_to_end(BidirectionalIterator* it);
GenoLocation next_forwards(BidirectionalIterator* it);
GenoLocation next_backwards(BidirectionalIterator* it);
GenoLocation next_get_nth(RandomAccessIterator* it, const size_t n);

//GenoLocation nextgappy_get_gap(struct GappyIterator* it);
//GenoLocation nextgappy_get_nongap(struct GappyIterator* it);
//GenoLocation nextgappy_valid_pos(struct GappyIterator* it);
//void _move_genotype(GenoLocation from, GenoLocation to, int* label_defaults);
    /**@}*/

    /** @defgroup liteget Getting data from an Iterator
     *
     * These functions take as a parameter a GenoLocation (the output
     * of iterators), and access the data for the single genotype that
     * the iterator has found.
     *
     * They are very lightweight: they have no error-checking and
     * implemented inline.
     *
     * @{
     */
/** Get the name of a genotype
 *
 * @param loc location of the relevant genotype
 * @return shallow copy of the name of the
 * genotype at location `loc`
 */
static inline char* get_name(const GenoLocation loc) {
    return loc.localAM->names[loc.localPos];
}

/** Get the name of a genotype
 *
 * @param loc location of the relevant genotype
 * @return shallow copy of the name of the
 * genotype at location `loc`
 */
/** Set the name of a genotype
 *
 * @param loc location of the relevant genotype
 * @param name name of the genotype
 * at location `loc` after this call
 */
static inline void set_name(const GenoLocation loc, char* name) {
    char* oldname = loc.localAM->names[loc.localPos];
    if (oldname != NULL) free(oldname);
    loc.localAM->names[loc.localPos] = name;
}

/** Get the alleles of a genotype
 *
 * @param loc location of the relevant genotype
 * @return shallow copy of the allele string of the
 * genotype at location `loc`. The loci are ordered
 * according to the genetic map (by chromosome, then
 * by location). The entire string is 2*n characters
 * long, where n is the number of loci. The two
 * alleles of each locus are presented side-by-side,
 * at positions 2i and 2i+1.
 */
static inline char* get_alleles(const GenoLocation loc) {
    return loc.localAM->alleles[loc.localPos];
}

/** Get the first/left parent of a genotype
 *
 * @param loc location of the relevant genotype
 * @return id of the left parent of the genotype
 * at location `loc`
 */
static inline PedigreeID get_first_parent(const GenoLocation loc) {
    return loc.localAM->pedigrees[0][loc.localPos];
}

/** Get the second/right parent of a genotype
 *
 * @param loc location of the relevant genotype
 * @return id of the right parent of the genotype
 * at location `loc`
 */
static inline PedigreeID get_second_parent(const GenoLocation loc) {
    return loc.localAM->pedigrees[1][loc.localPos];
}

/** Get the persistent id of a genotype
 *
 *  The persistent id is the number used in pedigree
 *  tracing, and is unique for the lifetime of the simulation.
 *
 * @param loc location of the relevant genotype
 * @return id of the genotype
 * at location `loc`
 */
static inline PedigreeID get_id(const GenoLocation loc) {
    return loc.localAM->ids[loc.localPos];
}

/** Get the current group membership of a genotype
 *
 * @param loc location of the relevant genotype
 * @return group number of the group affiliation
 * of the genotype at location `loc`
 */
static inline GroupNum get_group(const GenoLocation loc) {
    return loc.localAM->groups[loc.localPos];
}

/** Set the current group membership of a genotype
 *
 * @param loc location of the relevant genotype
 * @param group GroupNum of the group affiliation
 * of the genotype at location `loc` after this call
 */
static inline void set_group(const GenoLocation loc, const GroupNum group) {
    loc.localAM->groups[loc.localPos] = group;
}

//static inline int get_bv(GenotypeLocation loc) {}

/** Get the value of a specific label of a genotype
 *
 * @param loc location of the relevant genotype
 * @param labelIndex index of the relevant label.
 * @return value of the `labelIndex`th label
 * of the genotype at location `loc`
 */
static inline int get_label_value(const GenoLocation loc, const int labelIndex) {
    return loc.localAM->labels[labelIndex][loc.localPos];
}
    /**@}*/

    /** @defgroup search Data Searching Functions
     *
     * These functions search the simulation data for the genotype
     * that matches a particular known piece of information, eg a
     * name, id, global index, or set of parents. Depending on
     * the size of the simulation, they may not be fast.
     *
     * @{
     */
char* get_name_of_id( const AlleleMatrix* start, const PedigreeID id);
char* get_genes_of_id ( const AlleleMatrix* start, const PedigreeID id);
int get_parents_of_id( const AlleleMatrix* start, const PedigreeID id, PedigreeID output[2]);
void get_ids_of_names( const AlleleMatrix* start, const int n_names, const char* names[n_names], PedigreeID* output);
PedigreeID get_id_of_child( const AlleleMatrix* start, const PedigreeID parent1id, const PedigreeID parent2id);
int get_index_of_child( const AlleleMatrix* start, const PedigreeID parent1id, const PedigreeID parent2id);
int get_index_of_name( const AlleleMatrix* start, const char* name);
PedigreeID get_id_of_index( const AlleleMatrix* start, const int index);
char* get_genes_of_index( const AlleleMatrix* start, const int index);
    /**@}*/

    /** @defgroup collgetters Collective Data Access Functions
     *
     * These functions return vector data, rather than data on
     * a single genotype or single group.
     *
     * The `get_group_` family are left here for legacy purposes: using an
     * iterator is the new and less memory intensive way to do the tasks
     * these were used for.
     *
     * @{
     */
int get_group_size( const SimData* d, const GroupNum group_id);
int get_group_genes( const SimData* d, const GroupNum group_id, int group_size, char** output);
int get_group_names( const SimData* d, const GroupNum group_id, int group_size, char** output);
int get_group_ids( const SimData* d, const GroupNum group_id, int group_size, PedigreeID* output);
int get_group_indexes( const SimData* d, const GroupNum group_id, int group_size, size_t* output);
int get_group_bvs( const SimData* d, const GroupNum group_id, const EffectID effID, int group_size, double* output);
int get_group_parent_ids( const SimData* d, const GroupNum group_id, int group_size, const int whichParent, PedigreeID* output);
int get_group_parent_names( const SimData* d, const GroupNum group_id, int group_size, const int whichParent, char** output);
int get_group_pedigrees( const SimData* d, const GroupNum group_id, int group_size, char** output);

int get_existing_groups( SimData* d, GroupNum* output);
int get_existing_group_counts( SimData* d, GroupNum* out_groups, size_t* out_sizes);
    /**@}*/
/**@}*/


/** @defgroup groupmod Seletion/Group Modification Functions
 *
 * For simulation of selection or structure in breeding programs.
 *
 * @{
 */
GroupNum combine_groups( SimData* d, const int list_len, const GroupNum group_ids[list_len]);
GroupNum split_from_group( SimData* d, const int n, const size_t indexes_to_split[n]);
GroupNum split_by_label_value( SimData* d, const GroupNum group, const LabelID whichLabel, const int valueToSplit);
GroupNum split_by_label_range( SimData* d, const GroupNum group, const LabelID whichLabel, const int valueLowBound, const int valueHighBound);

// GENERIC
unsigned int _split_by_somequality( SimData* d, const GroupNum group_id,
        void* somequality_data,
        GroupNum (*somequality_tester)(GenoLocation, void*, unsigned int, unsigned int, GroupNum*),
        size_t maxentries_results, GroupNum* results);
// APPLICATIONS
unsigned int split_into_individuals( SimData* d, const GroupNum group_id, size_t maxentries_results, GroupNum results[maxentries_results]);
    //GroupNum __split_by_quality__individuate(GenoLocation loc, void** datastore, unsigned int maxgroups, unsigned int groupsfound, GroupNum** results);
unsigned int split_into_families(SimData* d, const GroupNum group_id, size_t maxentries_results, GroupNum results[maxentries_results]);
    //GroupNum __split_by_quality__family(GenoLocation loc, void** datastore, unsigned int maxgroups, unsigned int groupsfound, GroupNum** results);
unsigned int split_into_halfsib_families( SimData* d, const GroupNum group_id, const int parent, size_t maxentries_results, GroupNum results[maxentries_results]);
    //GroupNum __split_by_quality__halfsib1(GenoLocation loc, void** datastore, unsigned int maxgroups, unsigned int groupsfound, GroupNum* results);
    //GroupNum __split_by_quality__halfsib2(GenoLocation loc, void** datastore, unsigned int maxgroups, unsigned int groupsfound, GroupNum* results);
    //    GroupNum __split_by_quality__halfsibtemplate(GenoLocation loc, void** datastore, unsigned int maxgroups, unsigned int groupsfound, GroupNum* results, PedigreeID (*getparent)(GenoLocation));

// GENERIC
unsigned int _split_by_someallocation( SimData* d, const GroupNum group_id, void* someallocator_data,
        GroupNum (*someallocator)(GenoLocation, SimData*, void*, unsigned int, unsigned int*, GroupNum*),
        size_t n_outgroups, GroupNum outgroups[n_outgroups]);
// APPLICATIONS
GroupNum split_evenly_into_two(SimData* d, const GroupNum group_id);
    // GroupNum __split_by_allocator__knowncounts(GenoLocation loc, SimData* d, void* datastore, unsigned int n_outgroups, unsigned int* subgroupsfound, GroupNum* outgroups)
unsigned int split_evenly_into_n(SimData* d, const GroupNum group_id, const int n, GroupNum* results);
unsigned int split_by_specific_counts_into_n(SimData* d, const GroupNum group_id, const int n, const int* counts, GroupNum* results);
GroupNum split_randomly_into_two(SimData* d, const GroupNum group_id);
unsigned int split_randomly_into_n(SimData* d, const GroupNum group_id, const int n, GroupNum* results);
    // GroupNum __split_by_allocator__equalprob(GenoLocation loc, SimData* d, void* datastore, unsigned int n_outgroups, unsigned int* subgroupsfound, GroupNum* outgroups);
unsigned int split_by_probabilities_into_n(SimData* d, const GroupNum group_id, const int n, const double* probs, GroupNum* results);
    // GroupNum __split_by_allocator__unequalprobs(GenoLocation loc, SimData* d, void* datastore, unsigned int n_outgroups, unsigned int* subgroupsfound, GroupNum* outgroups);

/**@}*/


/** @defgroup crossers Crossing and Progression Simulation Functions
 *
 * For simulation of progression steps in breeding programs.
 *
 * @{
 */

    /** @defgroup meiosis Meiosis Simulation Functions
     *
     * For simulation of meiosis.
     *
     * @{
     */
void generate_gamete(SimData* d, const char* parent_genome, char* output);
void generate_cross(SimData* d, const char* parent1_genome, const char* parent2_genome, char* output);
void generate_doubled_haploid(SimData* d, const char* parent_genome, char* output);
void generate_clone(SimData* d, const char* parent_genome, char* output);
    /**@}*/

// SUPPORTER FUNCS FOR GENOPTIONS OPTIONS
FILE* _genoptions_setup_save_pedigrees(const GenOptions g);
FILE* _genoptions_setup_save_bvs(const SimData* d, const GenOptions g, int* effIndexp);
FILE* _genoptions_setup_save_genotypes(const SimData* d, const GenOptions g);
void _genoptions_save_pedigrees(FILE* fp, SimData* d, AlleleMatrix* tosave);
void _genoptions_save_bvs(FILE* fe, EffectMatrix* effMatrices, int effIndex, AlleleMatrix* tosave);
void _genoptions_save_genotypes(FILE* fg, AlleleMatrix* tosave);
void _genoptions_give_names_and_ids(AlleleMatrix* am, SimData* d, const GenOptions g);

// Double underscores in front of name represent parameter functions for a generic function.
// PARAMETER FUNCTIONS FOR THE FOLLOWING GENERIC
void __make_offspring_cross(SimData* d, void* datastore, GenoLocation parents[static 2], GenoLocation putHere);
void __make_offspring_self_n_times(SimData* d, void* datastore, GenoLocation parents[static 2], GenoLocation putHere);
void __make_offspring_doubled_haploids(SimData* d, void* datastore, GenoLocation parents[static 2], GenoLocation putHere);
void __make_offspring_clones(SimData* d, void* datastore, GenoLocation parents[static 2], GenoLocation putHere);

int __parentChooser_cross_randomly(void* parentIterator, void* datastore, unsigned int* counter, GenoLocation parents[static 2]);
int __parentChooser_cross_randomly_between(void* parentIterator, void* datastore, unsigned int* counter, GenoLocation parents[static 2]);
int __parentChooser_cross_combos(void* parentIterator, void* datastore, unsigned int* counter, GenoLocation parents[static 2]);
int __parentChooser_selfing(void* parentIterator, void* datastore, unsigned int* counter, GenoLocation parents[static 2]);
int __parentChooser_cloning(void* parentIterator, void* datastore, unsigned int* counter, GenoLocation parents[static 2]);

// GENERIC
GroupNum _make_new_genotypes(SimData* d, const GenOptions g,
        void* parentIterator, void* datastore,
        int (*parentChooser)(void*, void*, unsigned int*, GenoLocation[static 2]),
        void (*offspringGenerator)(SimData*, void*, GenoLocation[static 2], GenoLocation) );
// APPLICATIONS
GroupNum cross_random_individuals(SimData* d, const GroupNum from_group, const int n_crosses, const int cap, const GenOptions g);
	int _random_cross_checks(SimData* d, const GroupNum from_group, const int n_crosses, const int cap);
	unsigned int _specialised_random_draw(SimData* d, unsigned int max, unsigned int cap, unsigned int* member_uses, unsigned int noCollision);
GroupNum cross_randomly_between(SimData*d, const GroupNum group1, const GroupNum group2, const int n_crosses, const int cap1, const int cap2, const GenOptions g);
GroupNum cross_these_combinations(SimData* d, const int n_combinations, const int* firstParents, const int* secondParents, const GenOptions g);
GroupNum self_n_times(SimData* d, const unsigned int n, const GroupNum group, const GenOptions g);
GroupNum make_doubled_haploids(SimData* d, const GroupNum group, const GenOptions g);
GroupNum make_clones(SimData* d, const GroupNum group, const int inherit_names, const GenOptions g);

GroupNum make_all_unidirectional_crosses(SimData* d, const GroupNum from_group, const GenOptions g);
GroupNum make_n_crosses_from_top_m_percent(SimData* d, const int n, const int m, const GroupNum group, const EffectID effID, const GenOptions g);
GroupNum make_crosses_from_file(SimData* d, const char* input_file, const GenOptions g);
GroupNum make_double_crosses_from_file(SimData* d, const char* input_file, const GenOptions g);
/**@}*/


/** @defgroup calculators Breeding Value and Allele Count Calculators
 *
 * For calculations related to the loaded allele effects and internal additive breeding value model.
 *
 * @{
 */
GroupNum split_by_bv(SimData* d, const GroupNum group, const EffectID effID, const int top_n, const int lowIsBest);
DecimalMatrix calculate_group_bvs(const SimData* d, const GroupNum group, const EffectID effID);
DecimalMatrix calculate_bvs( const AlleleMatrix* m, const EffectMatrix* e);
int calculate_group_count_matrix_of_allele( const SimData* d, const GroupNum group, const char allele, DecimalMatrix* counts);
int calculate_group_doublecount_matrix_of_allele( const SimData* d, const GroupNum group, const char allele, DecimalMatrix* counts, const char allele2, DecimalMatrix* counts2);
int calculate_count_matrix_of_allele( const AlleleMatrix* m, const char allele, DecimalMatrix* counts);
int calculate_doublecount_matrix_of_allele( const AlleleMatrix* m , const char allele, DecimalMatrix* counts, const char allele2, DecimalMatrix* counts2);
DecimalMatrix calculate_full_count_matrix_of_allele( const AlleleMatrix* m, const char allele);

MarkerBlocks create_n_blocks_by_chr(const SimData* d, const int n);
MarkerBlocks read_block_file(const SimData* d, const char* block_file);
void calculate_group_local_bvs(const SimData* d, const MarkerBlocks b, const EffectID effID, const char* output_file, const GroupNum group);
void calculate_local_bvs(const SimData* d, const MarkerBlocks b, const EffectID effID, const char* output_file);

char* calculate_optimal_alleles(const SimData* d, const EffectID effID);
char* calculate_optimal_available_alleles(const SimData* d, const GroupNum group, const EffectID effID);
double calculate_optimum_bv(const SimData* d, const EffectID effID);
double calculate_optimal_available_bv(const SimData* d, const GroupNum group, const EffectID effID);
double calculate_minimum_bv(const SimData* d, const EffectID effID);
/**@}*/


/** @defgroup deletors Deletor Functions
 *
 * For deleting and free associated memory of data structures.
 *
 * @ingroup structs
 * @{
 */
void delete_group(SimData* d, const GroupNum group_id);
void delete_label(SimData* d, const LabelID whichLabel);
void delete_genmap(GeneticMap* m);
void delete_allele_matrix(AlleleMatrix* m);
void delete_effect_matrix(EffectMatrix* m);
void delete_eff_set(SimData* d, EffectID whichID);
void delete_simdata(SimData* m);
void delete_markerblocks(MarkerBlocks* b);
void delete_dmatrix(DecimalMatrix* m);
void delete_bidirectional_iter(BidirectionalIterator* it);
void delete_randomaccess_iter(RandomAccessIterator* it);
/**@}*/


/** @defgroup savers Saving Functions
 *
 * For saving persistent simulation results.
 *
 * @{
 */
void save_simdata(FILE* f, const SimData* m);

void save_marker_blocks(FILE* f, const SimData* d, const MarkerBlocks b);

void save_allele_matrix(FILE* f, const AlleleMatrix* m);
void save_names_header(FILE* f, size_t n, const char* names[n]);
void save_transposed_allele_matrix(FILE* f, const AlleleMatrix* m, const char** markers);
void save_group_alleles(FILE* f, SimData* d, const GroupNum group_id);
void save_transposed_group_alleles(FILE* f, const SimData* d, const GroupNum group_id);

void save_group_one_step_pedigree(FILE* f, const SimData* d, const GroupNum group);
void save_one_step_pedigree(FILE* f, const SimData* d);
void save_group_full_pedigree(FILE* f, const SimData* d, const GroupNum group);
void save_full_pedigree(FILE* f, const SimData* d);
void save_AM_pedigree(FILE* f, const AlleleMatrix* m, const SimData* parents);
void save_parents_of(FILE* f, const AlleleMatrix* m, PedigreeID p1, PedigreeID p2);

void save_group_bvs(FILE* f, const SimData* d, const GroupNum group, const EffectID effID);
void save_manual_bvs(FILE* f, const DecimalMatrix* e, const PedigreeID* ids, const char** names);
void save_bvs(FILE* f, const SimData* d, const EffectID effID);

void save_count_matrix(FILE* f, const SimData* d, const char allele);
void save_count_matrix_of_group(FILE* f, const SimData* d, const char allele, const GroupNum group);

/**@}*/


/** @defgroup recomb Recombination Calculators
 *
 * Experimental functions for retroactively calculating number of recombinations.
 *
 * This functionality is for interest only. It is not clear, or tidy,
 * or checked against real data.
 *
 * @{
 */
int* calculate_min_recombinations_fw1(SimData* d, char* parent1, unsigned int p1num, char* parent2,
        unsigned int p2num, char* offspring, int certain); // forward filling, window size 1
int* calculate_min_recombinations_fwn(SimData* d, char* parent1, unsigned int p1num, char* parent2,
        unsigned int p2num, char* offspring, int window_size, int certain); // forward filling, window size n


// int has_same_alleles(char* p1, char* p2, int i);
// int has_same_alleles_window(char* g1, char* g2, int start, int w);
/** Simple operator to determine if at marker i, two genotypes share at least
 * one allele. Checks only 3 of four possible permutations because assumes
 * there cannot be more than two alleles at a given marker.
 *
 * @param p1 pointer to a character array genotype of the type stored in an AlleleMatrix
 * (2*n_markers long, representing the two alleles at a marker consecutively) for the first
 * of the genotypes to compare.
 * @param p2 pointer to a character array genotype for the second of the genotypes to compare.
 * @param i index of the marker at which to perform the check
 * @returns boolean result of the check
 */
static inline int has_same_alleles(const char* p1, const char* p2, const int i) {
    return (p1[i<<1] == p2[i<<1] || p1[(i<<1) + 1] == p2[i] || p1[i] == p2[(i<<1) + 1]);
}
// w is window length, i is start value
/** Simple operator to determine if at markers with indexes i to i+w inclusive, two genotypes
 * share at least one allele. Checks only 3 of four possible permutations at each marker
 * because assumes there cannot be more than two alleles at a given marker. For the return value
 * to be true, there must be at least one match at every one of the markers in the window.
 *
 * @param g1 pointer to a character array genotype of the type stored in an AlleleMatrix
 * (2*n_markers long, representing the two alleles at a marker consecutively) for the first
 * of the genotypes to compare.
 * @param g2 pointer to a character array genotype for the second of the genotypes to compare.
 * @param start index of the first marker in the window over which to perform the check
 * @param w length of the window over which to perform the check
 * @returns boolean result of the check
 */
static inline int has_same_alleles_window(const char* g1, const char* g2, const int start, const int w) {
    int same = TRUE;
    int i;
    for (int j = 0; j < w; ++j) {
        i = start + j;
        same = same && (g1[i<<1] == g2[i<<1] || g1[(i<<1) + 1] == g2[i] || g1[i] == g2[(i<<1) + 1]);
    }
    return same;
}

int calculate_recombinations_from_file(SimData* d, const char* input_file, const char* output_file,
        int window_len, int certain);
/**@}*/

#endif

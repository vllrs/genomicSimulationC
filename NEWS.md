Latest News       {#news}
===========

# genomicSimulation (development version)

## New Features

- It is now possible to add "centre" values to a set of marker effects. The centre values are subtracted from the breeding value contribution of each marker when calculating breeding values. Centres can be added using new functions `change_eff_set_centres_to_values`, `change_eff_set_centres_of_markers`, and `change_eff_set_centres_of_allele_count`, or by adding an extra column titled "centre" to the input effect file. These centre values can be useful, for example, for marker effects derived from GBLUP, where proper use of these marker effects requires subtracting a factor based on allele frequency from each marker's score.
- CHANGE: when initialising the simulation with a genotype file but no map file, the default behaviour is now to generate a map where each genetic marker from the genotype file belongs to a separate chromosome/linkage group, and therefore each marker's inheritance is unlinked/independent of the inheritance of alleles at other markers. The old default (generated map where all markers are on the same chromosome, 1cM apart) is still available via function `gsc_create_uniformspaced_recombmap`. The new behaviour is available via function `gsc_create_unlinked_recombmap`.
- A warning is now displayed if a map is loaded that has a very large chromosome length to number of markers ratio, which might correspond to loading a map with (invalid) base pair positions instead of (correct) centimorgan positions.
- A warning is now displayed if the primary map contains duplicate marker names, because data from further input files will only be loaded into one of the duplicates.
- The simulation internally now does store chromosome names, not just chromosome orderings.

## Improvements

- Generalisation of header-parsing logic of `load_mapfile` and `load_effectfile`, to prepare for accepting further columns in future.
- BREAKING CHANGE: `calculate_local_bvs` now returns local BVs in a DecimalMatrix and can be called with a group number or NO_GROUP, like other `calculate_` functions. Removed `calculate_group_local_bvs` as its function was subsumed. To calculate and save local BVs to a file, use the function `save_local_bvs`.
- `EffectMatrix` has been replaced with new struct `MarkerEffects`. To match, `gsc_delete_effect_matrix` has been renamed to `gsc_delete_effect_set`.
- Makefile has a new target "sharedlib" (run `make sharedlib`) to compile genomicSimulation as a shared library, rather than a standalone executable that runs the test suite.
- The conversion script that creates the library for the genomicSimulation R package now converts `fprintf(stderr, "msg")` to `Rprintf("NOTE! msg")`, instead of `warning("msg")`, so that non-stopping errors like these are still printed if a segfault occurs. (By default, R hides `warning` messages until the end of execution scope, and so loses them if a segfault occurs). 
- The conversion script that creates the library for the genomicSimulation R package now adds a tag to converted files below the version number, noting that they have been through the conversion process, and what version of the conversion script converted them. The first version-numbered conversion script is version "v2".

## Bug Fixes

- Fix a segfault that occured when calling `load_data_files` with `genotype_file = NULL`.
- Fix a segfault that occured when trying to create iterators for a simulation object with no loaded genotypes.
- Fix a bug where marker names had trailing spaces (and therefore, when loading more input files later, could not be matched) when automatically generating a genetic map from a genotype file. 
- Fix a bug in `create_evenlength_blocks_each_chr` where too many markers would be allocated to the last block because genetic distances were being cumulatively summed even though they were already cumulative values. This bug would date back to v0.2.5, when ability to load multiple recombination maps was added.
- Memory allocations of size 0 are no longer requested when attempting to use a MarkerBlocks struct that has been deleted.
- `SimData.n_groups` is no longer decremented if a call is made to `delete_group` for a group number with zero members (a nonexistent group). 
- Fix a bug during genotype file loading where alleles would not be saved to the SimData if the program was creating a genetic map from the markers in the file. 

# genomicSimulation v0.2.6

## New Features

- BREAKING CHANGE: `load_genotypefile` and `load_data_files` now take an extra final parameter of type `FileFormatSpec`. This parameters allows users to specify some or all of the layout and formatting details of an input genotype file. During the loading process, the sections tasked with detecting the already-specified features of the layout are then skipped.
- BREAKING CHANGE: Introduced new type aliases: `GSC_ID_T`, `GSC_GLOBALX_T`, `GSC_LOCALX_T`, `GSC_GENOLEN_T`. Some function type signatures may have changed, as they are now using these type aliases, and the type aliases have consistent definitions across all functions in the library.
- `make_targeted_crosses` will now skip invalid pairings, where one or more parents does not exist, instead of halting execution at the first invalid pairing. If any pairings were skipped, it will print a debug message listing the number of invalid pairings detected.
- Improved Makefile. It now supports both debug and release builds.

## Bug Fixes

- Log messages printed during genotype file loading now accurately reflect the number of markers per genotype that were successfully loaded. Previously, it incorrectly printed out the total number of markers in the stored genetic map, while claiming it was the number of markers in the genotype file that had been accurately matched to the map.
- Input genotype files now have `fclose` called on them after being parsed.
- Fixed memory leaks in `load_genotypefile` and `delete_simdata` introduced in v0.2.5.
- Fixed segmentation fault when calling `make_clones` with `inherit_names = GSC_TRUE` on genotypes with no names.
- More edge cases covered in `condense_allele_matrix`, to patch a rare memory leak.
- Fixed memory access error in `make_random_crosses_between` when the function was called with breeding caps on both groups.
- `get_index_of_name` no longer tries to compare its target name to the names of nameless genotypes. Also, it now prints a warning and returns an invalid value when it cannot find the target name, instead of calling exit().
- Group RandomAccessIterators no longer request memory for their cache structures before checking if there are actually that many members in the group. This fixes a bug where passing a large invalid `n` to `next_get_nth` would cause the iterator to hang, trying to request incredibly large blocks of memory.

## Improvements 

- Polish and improve the Templates page of the documentation.
- Improvements to the genotype matrix layout detection processes, to cover more edge cases and use a more comprehensible decision tree.
- Improvements to input file loading tests and allele count matrix calculation tests.
- Helper function `shuffle_up_to` is now size-agnostic and can shuffle arrays of any type, instead of only integer arrays.
- BREAKING CHANGE: The file output interface (functions with prefix `save_`) has undergone a simplification and rewrite. There is now one function per data type, which offers whole-simulation vs single-group output and output file format choices as parameters, much like the R version of genomicSimulation (as opposed to separate function whether the output was for a single group or all genotypes, and also different functions per ouput formats). All user-facing functions now take a file name instead of a file pointer. See below for details.
	- BREAKING CHANGE: `save_markerblocks` now offers two file output formats, but neither exactly matches the old file format. To create a marker blocks file with the old format, call `save_markerblocks` without a map, then run `sed 's/^/0\t0\tb0\tb\t/' f.txt | sed '1i\Chrom\tPos\tName\tClass\tMarkers'`.
	- BREAKING CHANGE: `save_names_header`, `save_allele_matrix`, `save_transposed_allele_matrix`, `save_group_genotypes`, and `save_transposed_group_genotypes` have been replaced with the single function `save_genotypes`. 
	- BREAKING CHANGE: `save_group_count_matrix` and `save_count_matrix` have been replaced with the single function `save_allele_counts`. This new combined function also allows for printing a "transposed" count matrix (where genetic markers are rows, rather than columns), as was possible for a genotype matrix prior to this update. 
	- BREAKING CHANGE: `save_one_step_pedigree`, `save_group_one_step_pedigree`, `save_full_pedigree`, and `save_group_full_pedigree` have been replaced with the single function `save_pedigrees`.
	- BREAKING CHANGE: `save_bvs`, `save_group_bvs`, and `save_manual_bvs` have been combined into the single function (with updated function signature) `save_bvs`.
- BREAKING CHANGE: Several features among the calculation functions interface (functions with prefix `calculate_`) have undergone a restructure and rewrite.
	- BREAKING CHANGE: `calculate_group_count_matrix`, `calculate_count_matrix`, and `calculate_full_count_matrix` have been combined into the single function (with updated function signature) `calculate_allele_counts`.
	- BREAKING CHANGE: `calculate_bvs` and `calculate_group_bvs` have been combined into the single function (with updated function signature) `calculate_bvs`.  
- BREAKING CHANGE: Now using _Bool type for true/false values. The GSC_LOGICVAL type is used only for variables that need three values: true/false/invalid. GSC_UNINIT has been renamed to GSC_NA.
- The dimensions of a DecimalMatrix are now stored in larger size_t variables.
- BREAKING CHANGE: Removed function `make_n_crosses_from_top_m_percent` for disobeying standard genomicSimulation rules of division of functionality. A section has been added to the Templates section of the documentation with a drop-in replacement of this function.
- More informative error message when calling `split_by_bv` on a nonexistent group.
- Removed unused struct MarkerPosition.


# genomicSimulation v0.2.5

## New Features 

- Add ability to load multiple recombination maps into simulation at a time. Loading a genetic map file will return an MapID that represents that particular map. All crossing functions now have additional parameters that allow you to select which recombination map to use when generating gametes from each group of parents. The crossing functions default to using the first loaded recombination map if not otherwise specified. Maps can be removed from memory using delete_recombination_map. The list of genetic markers changed by the simulation, however, is still immutable, so genetic markers not present in the first map loaded will not be tracked or simulated.
- New function `change_allele_symbol` changes the internal representation of a particular allele. It does not modify corresponding marker effects. 

## Bug Fixes

- Fix a bug in large simulations where gaps in the AlleleMatrix linked list failed to be patched in `condense_allele_matrix`, and some genotypes lost their group numbers (group numbers reverted to 0). 

## Improvements

- Smarter and more robust genotype matrix, genetic map, and effect file loaders. See R version vignette or Templates page of C version documentation for more information. 
- All C code now has C-style "namespace guards" (the prefix `gsc_` in front of every name). Data structures and user-friendly functions are still exported with non-prefixed "short names" unless `GSC_NO_SHORT_NAMES` is defined before import. The excellent stb libraries (https://github.com/nothings/stb) served as guides for how to do this.
- Several function names have been changed, for better consistency within functional categories and better consistency with the R version of genomicSimulation. To keep using the old names (at least while their function signatures remain the same), import the `names_stopgap.h` header file. Optionally, `#define GSC_DEPRECATED_VERBOSE` before importing `names_stopgap.h`, and it will tell you exactly what names have changed to what, for a nice and easy find-and-replace way to update to the new names. 
- CONTIG_WIDTH and NAME_LENGTH (the "global settings") can now be altered without needing to modify the genomicSimulation source files, using macro definitions. See instructions at the top of sim-operations.h
- genomicSimulation's malloc and free can be replaced with alternate implementations, using macro definitions. See instructions at the top of sim-operations.h. 
- Significant under-the-hood changes to crossing functions. They now call `_make_genotypes`, a generic function, to reduce code duplication. 
- The same script and same random seed re-run post-0.2.4.003 will now produce slightly different genotypes. This is because the two gametes that make a cross are now generated successively (first one, then the other) rather than simultaneously.
- Under-the-hood changes to `condense_allele_matrix`, for easier debugging. It now takes advantage of a `GappyIterator` to help it shuffle the positions of genotypes in the AlleleMatrix linked list.
- Remove all Variable Length Arrays (VLAs). All buffers used internally in genomicSimulation now allocate themselves sizeof(int)\*CONTIG_WIDTH bytes of space on the stack, or use heap memory if that is insufficient.  

# genomicSimulation 0.2.4.003

## Bug Fixes 

- Group modification functions (`split_` prefix) no longer fail if there are more than 10 groups in existence. 

## Improvements 

- Significant under-the-hood changes to `split_`-prefixed functions. They now use BidirectionalIterator and RandomAccessIterator to find members of the group to split, and some call `_split`-prefixed generic functions to reduce code duplication.
- `split_`-prefixed functions now return the number of groups they created, rather than being `void` functions. Thanks to the GroupNum type (introduced v0.2.4) we can avoid confusion over whether the return value of these functions represents a group identifier or a count of groups. 
	- For users, this will mean it is no longer necessary to zero the output array passed to these functions so that you can identify how many entries are the newly-created groups after the function finishes execution. Instead, the return value lets you know how many entries have been saved to the output array.
- BREAKING CHANGE: `split_into_families` and `split_into_halfsib_families` now take a parameter for length of the output array, and will truncate output if they find more families than that.
- Added a cache `n_groups` inside SimData to track the number of groups present in simulation at the current point in time. This allows users to know how much space to allocate before calling functions like `get_existing_groups` and `get_existing_group_counts`, and the program to know how much space to allocate inside `split_into_families` and `split_into_halfsib_families`. This will help avoid issues around truncated output or memory overruns.
	- Specifically, this cache is an upper bound on the number of groups currently present in simulation, because `split_from_group` is not aware if it has destroyed a group by moving the entire group into its new one. 
	- `get_existing_groups` and `get_existing_group_counts` will update `n_groups` to the true value of the number of groups in simulation as a side-effect of their calls. 
- BREAKING CHANGE: Function signatures of `get_existing_groups` and `get_existing_group_counts` are slightly different. They can now be passed NULL to output parameters that you don't need, and they assume that their not-NULL output parameters are of length `SimData.n_groups` at least, rather than asking for the length of their output parameters separately.
- Gradual transformation of genotype indexes from `unsigned int` type to a `size_t` type. This conversion is now implemented in the iterators and `get_group_indexes`.
- New iterator helper function `set_group`, to go alongside `get_group`, `get_id`, etc. It is currently the only iterator helper function that modifies data.
- BidirectionalIterators now use local position rather than global position. Just in case of hypothetical really huge simulations.
- Removed some qsort comparator functions that are no longer called by anything.


# genomicSimulation 0.2.4

## New Features

- Add ability to load multiple effect files. Multiple sets of marker effects can now be available at one time. Loading an effect file will return an EffectID that represents that particular set of marker effects. Sets of marker effects can be removed from memory using delete_eff_set.

## Improvements

- BREAKING CHANGE: The four different types of identifiers have been tidied into four struct types (GroupNum, PedigreeID, LabelID and EffectID), and all functions have been updated to use these struct types. Each of these types also has a named invalid value that may be returned if a function fails: NO_GROUP, NO_PEDIGREE, NOT_A_LABEL, and NOT_AN_EFFECT_SET. 
- BREAKING CHANGE: Because sets of marker effects now have identifiers, load_all_simdata now returns a struct GroupAndEffectSet, which is a wrapper around a GroupNum and EffectID pair, instead of just a group number.
- Added test cases for functions that cover saving data to output files.
- The save-as-you-go setting for saving breeding values in GenOptions now takes the value of the id from an EffectID, if you wish to save the breeding values from that set of marker effects, or 0/NOT_AN_EFFECT_SET, to not use the save-as-you-go functionality.
 
## Bug Fixes
- Standardised file-output functions so that genotype names are consistently substituted with their PedigreeIDs, if the genotype does not have a name. 
- Stopped save-as-you-go genotype saving repeating the header row (of genetic marker names) every 1000 rows. The header row now only appears in the first row.
- Semicolon separators appear correctly in the output of save_marker_blocks
- Pedigrees are printed consistently between save_group_[]_pedigree and save_[]_pedigree variants of the same function. Parents are no longer skipped in favour of starting with grandparents in group-specific functions.


# genomicSimulation 0.2.3.002

- genomicSimulation had a chance of a segmentation fault in load_all_simdata for certain marker effect files. The chance was higher for effect files listing few markers or listing many alleles. This release is a quick-fix for this bug.
- Also fixed a flaw in Rconversion.sh, the script that converts genomicSimulationC's code to genomicSimulation (R version)

# genomicSimulation 0.2.3

## New Features

- Add functions calculate_optimal_available_alleles and calculate_optimal_available_bv, to calculate the best combination of alleles and best possible breeding value score given the pool of alleles available in a particular group.
- Added custom integer labels, eg. for tracking age:
	- Create a new label: create_new_label (every genotype has a value for every label; multiple labels can exist at a time).
	- Change the values of a label for some or all genotypes: set_labels_to_const; increment_labels; set_labels_to_values
	- Group manipulation using the values of a label: split_by_label_value; split_by_label_range
- Added two iterators: BidirectionalIterator and RandomAccessIterator. They can iterate through every genotype in the simulation, or through just the members of a specific group. Iterators need to be reset or destroyed and re-created when genotypes are created or destroyed.
	- Creating and initialising iterators: create_bidirectional_iter + set_bidirectional_iter_to_start or set_bidirectional_iter_to_end;  create_randomaccess_iter
	- Iterating: next_forwards; next_backwards; next_get_nth
	- Deleting: delete_bidirectional_iter; delete_randomaccess_iter
- Iterators return GenoLocation values. A new family of get_ functions has been added to access the data of the genotype at a particular GenoLocation.
	- get_name; get_alleles; get_first_parent; get_second_parent; get_id; get_group; get_label_value
	
## Improvements 

- BREAKING CHANGE: Removed automatic re-seeding in get_chromosome_locations. Users should now set the random number generator seed in their own programs (see the example sim.c). This removes the unexpected behaviour where genomicSimulation overwrites your own random seeds.
- BREAKING CHANGE: Collective data access functions (get_group_\* family, get_existing_groups and get_existing_group_counts) now use parameter modification (save their results to an array passed as a parameter) rather than by returning a heap array.
- BREAKING CHANGE: Changed function signature of cross_these_combinations to allow for more flexible usage: instead of one `combinations[2][n_crosses]` parameter, there are a pair of `firstParents[n_crosses]` and `secondParents[n_crosses]` parameters.
- Added 'const' to all readonly parameters. Compilation should no longer produce warnings.
- Swapped to using Mattias Gustavsson's "rnd" implementation of a Permuted Congruential Generator (https://github.com/mattiasgustavsson/libs) for random number generation, rather than the builtin and subpar "rand()".

## Bug Fixes 

- Up NAME_LENGTH (the buffer size for name loading) to 45, as a temporary fix.
- Allow uint ordered search to ignore 0s mid-list, to solve a bug when not all genotypes are allocated IDs

# genomicSimulationC 0.2.2

- Swapped out a custom randlim() in shuffle_up_to for the basic call
- Removed a "negative ID" warning for get_group_parent_names when the group contains members that have no known parents.
- Fixed a bug in get_n_new_group_nums which incorrectly identified some existing groups as empty.
- Fixed some memory leaks in the tests
- Fixed an inappropriate memory access bug in crossing functions when will_allocate_ids = FALSE
- Fixed an inappropriate memory access bug in make_clones when trying to inherit the name of a parent with no name.

# genomicSimulationC 0.2.1

- Added function make_clones (and its genotype-generating subfunction generate_clone) to make identical duplicates of members of a group.
- Added function to randomly cross between two groups, with the option of a cap on the number of uses of each group member as a parent of a cross:
    - cross_randomly_between
- Functions split_into_individuals and split_into_families now return the groups they create.
- Added more functions for splitting groups:
    - split_into_halfsib_families
    - split_evenly_into_two
    - split_evenly_into_n
    - split_by_specific_counts_into_n
    - split_randomly_into_two
    - split_randomly_into_n
    - split_by_probabilities_into_n

Minor changes:

- Risk of losing reference to SimData object in load_all_simdata fixed.
- Fixed read from freed memory error in self_n_times, fixed some dereference null pointer risks in other crossing functions.
- Improve parameter checking in crossing and saving functions.
- Fixed get_group_[data] functions requesting a 0-length block of heap space when they are provided with a nonexistent group. Instead they now stop and return NULL
- Fix segfault when trying to select more individuals by GEBV than exist in the group (e.g. asking for the best 5 members of a group of 2). Now, it just moves all group members to the new selected group, and doesn't worry about the missing requested remainder.
- Added two new helper functions: get_n_new_group_nums and shuffle_up_to


# genomicSimulationC 0.2

- State corresponding to [genomicSimulation R package v0.2](https://github.com/vllrs/genomicSimulation/releases/tag/v0.2) release.
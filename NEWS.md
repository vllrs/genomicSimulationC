Latest News       {#news}
===========

## Improvements

- Significant under-the-hood changes to crossing functions. They now call `_make_genotypes`, a generic function, to reduce code duplication. 
- The same script and same random seed re-run post-0.2.4.003 will now produce slightly different genotypes. This is because the two gametes that make a cross are now generated successively (first one, then the other) rather than simultaneously.

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
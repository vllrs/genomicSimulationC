Latest News       {#news}
===========

# Since last release

- Swapped out a custom randlim() in shuffle_up_to for the basic call

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
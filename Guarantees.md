# Guarantees

A list of understandings about the simulation tool that you can trust will hold. C function names are used, but the assumptions still hold in the R version:

- No function can change genetic marker indexes or order, or the genetic map. The genetic map is set and static from the point it is loaded. A new SimData must be created to change it
- No function except delete_group can change genotype indexes or invalidate iterators. delete_group deletes all genotypes in the named groups, and shuffles all genotypes after the now-deleted ones backwards (retaining order) to fill in gaps where now-deleted genotypes sat. 
- No function will change the relative order of genotypes (follows from above).
- Progeny from the same cross will be given sequential indexes (and IDs, if IDs are allocated). That is, when GenOptions.family_size > 1 when calling a crossing function, that fullsib family are stored side-by-side to each other.
- Genotypes loaded from an input file will be stored (and given IDs, and indexes) in the same order as they appeared in the input file. 
- The order of genetic markers may not match the order of genetic markers in the input file: they will be re-ordered to match the genetic map (eg sorted by lower chromosome number and earlier position).


## Valid values

- Genotype indexes will be integer values greater than or equal to 0. These correspond to the position the genotype is currently stored at inside the simulation's data structures.
- Genotype IDs will be integer values greater than or equal to 1. 0 denotes no/unknown ID. This is used to track pedigree, and is unique to each genotype for the lifetime of the simulation.
- Group IDs/group numbers will be integer values greater than or equal to 1. All genotypes will be a member of a valid group; 0 is an error code. This is used to group together genotypes and run commands on groups of genotypes.
- Label IDs will be integer values greater than or equal to 1. Labels can be used to save information on each genotype (eg age, sex).


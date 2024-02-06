Features and Templates Guide        {#templates}
============================

While [the documentation pages of genomicSimulationC](@ref modules) and its counterpart [genomicSimulation](https://github.com/vllrs/genomicSimulation) do describe all the features the simulation tool offers, they don't do so in a way that is necessarily easy to browse or that's helpful for new users to understand the simulation tool.

This page provides a set of templates for implementing common breeding project actions in genomicSimulationC and its R counterpart.

[TOC]

# Setting Up

## Set the random number generator
In C, you should manually set the seed for the random number generator when you create the SimData object, before loading data. The R version borrows R's own random number generator so has no need for setting a seed, or for manually creating the simulation data structure before loading data.

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>
<td>
```{C}
SimData* sd = create_empty_simdata(time(NULL));
```
<td>
</table>

## Load a genetic map and a set of founder genotypes

<table>
<tr><th>Task <th>Input Files <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>
<td>Genotype file/Allele file:
- File location: genotype-file.txt
- File contents:
```
name	G01	G02	G03	G04	G05	G06
m1	TT	TT	TT	TA	TT	AT
m3	TT	TT	TA	TA	TT	TT
m2	AA	AA	AA	AA	TT	AA
```

Map file:
- File location: map-file.txt
- File contents:
```
marker chr pos
m3 3 15
m2 1 8.3
m1 1 5.2
```
<td>
```{C}
SimData* d = create_empty_simdata(time(NULL));
struct GroupAndEffectSet init = load_all_simdata(d, "genotype-file.txt", "map-file.txt", NULL);
GroupNum founders = init.group;
```
<td>
```{R}
founders <- load.data(allele.file="genotype-file.txt", map.file="map-file.txt")
```
</table>

Consider using the get_group_ family of functions (genomicSimulationC) and see.group.data function (genomicSimulation) to confirm data is correctly loaded.

## Load a genetic map and several sets of founder genotypes

<table>
<tr><th>Task <th>Input Files <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>
<td>Genotype file/Allele file:
- File location: genotype-file.txt
- File contents:
```
name	G01	G02	G03	G04	G05
m1	TT	TT	TT	TA	TT
m3	TT	TT	TA	TA	TT
m2	AA	AA	AA	AA	TT
```

Genotype file 2:
- File location: genotype-file2.txt
- File contents:
```
name	G06
m1	AT
m2	AA
m3	TT
```

Map file:
- File location: map-file.txt
- File contents:
```
marker chr pos
m3 3 15
m2 1 8.3
m1 1 5.2
```
<td>
```{C}
SimData* d = create_empty_simdata(time(NULL));
struct GroupAndEffectSet init = load_all_simdata(d, "genotype-file.txt", "map-file.txt", NULL);
GroupNum founders_a = init.group;
GroupNum founders_b = load_more_transposed_genes_to_simdata(d, "genotype-file2.txt");
```
<td>
```{R}
founders_a <- load.data(allele.file="genotype-file.txt", map.file="map-file.txt")
founders_b <- load.more.genotypes("genotype-file2.txt")
```
</table>

## Load a genetic map, a set of founder genotypes, and a set of additive trait effects

<table>
<tr><th>Task <th>Input Files <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>
<td>Genotype file/Allele file:
- File location: genotype-file.txt
- File contents:
```
name	G01	G02	G03	G04	G05	G06
m1	TT	TT	TT	TA	TT	AT
m3	TT	TT	TA	TA	TT	TT
m2	AA	AA	AA	AA	TT	AA
```

Map file:
- File location: map-file.txt
- File contents:
```
marker chr pos
m3 3 15
m2 1 8.3
m1 1 5.2
```

Effect file:
- File location: eff-file.txt
- File contents: (line order does not matter)
```
m1 A -0.8
m2 A -0.1
m3 A 0.1
m1 T 0.9
m3 T -0.1
```
<td>
```{C}
SimData* d = create_empty_simdata(time(NULL));
struct GroupAndEffectSet init = load_all_simdata(d, "genotype-file.txt", "map-file.txt", "eff-file.txt");
GroupNum founders = init.group;
EffectID eff1 = init.effectSet;
```
<td>
```{R}
# Note that the output is slightly different when an effect file is loaded.
init <- load.data(allele.file="genotype-file.txt", map.file="map-file.txt", effect.file="eff-file.txt")
founders <- init$groupNum
eff1 <- init$effectID
```
</table>

## Load another set of additive trait effects

This assumes the simulation is set up, that is, one of **Load a genetic map and a set of founder genotypes**, **Load a genetic map and several sets of founder genotypes**, **Load a genetic map, a set of founder genotypes, and a set of additive trait effects** or equivalent has been carried out. If the simulation is not yet set up and you wish to load trait effect values, see the section above, titled **Load a genetic map, a set of founder genotypes, and a set of additive trait effects**.

<table>
<tr><th>Task <th>Input Files <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>
<td>Effect file:
- File location: eff-file2.txt
<td>
```{C}
EffectID eff2 = load_effects_to_simdata(d, "eff-file2.txt");
```
<td>
```{R}
eff2 <- load.different.effects("eff-file2.txt")
```
</table>

## Input file formats:
The file format expected for the genotype file is something like:
```
name	G01	G02	G03	G04	G05	G06
m1	TT	TT	TT	TA	TT	AT
m3	TT	TT	TA	TA	TT	TT
m2	AA	AA	AA	AA	TT	AA
```
where G01, G02, ..., are names of the starter set of genotypes; m1, m2, ..., are the markers; and entries in the matrix are two-character pairs listing the two alleles that genotype has at that marker. Alleles can be any non-space character and the order of the two alleles is saved.

Cells may be space-separated or tab-separated. The order of rows and columns does not matter. The value in the first cell ("name" in this example) is ignored.

The C version of the package has additional functions for loading differently-formatted genotype files. This functionality is on track to be improved and shared between versions in some upcoming update.

The map file should be formatted as follows:
```
marker chr pos
m3 3 15
m2 1 8.3
m1 1 5.2
```

The first (header) line's values are not checked. After that, all rows must have three space-or-tab-separated values. The first should be the marker name, the second an integer representing the chromosome number, and the third a decimal representing the position of the SNP along the chromosome in centiMorgans (cM). The order of the rows in the file does not matter.

Only markers that appear both in the genotype file and the map file are used in simulation. The simulation tool will print as output the number of markers that were loaded and the number that were discarded.

Loading effect file(s) is optional for running the simulation. Effect files should should be formatted as follows:
```
m1 A -0.8
m2 A -0.1
m3 A 0.1
m1 T 0.9
m3 T -0.1
```

The first column should be a marker name, corresponding to a name used in the map and genotype files. The second should be the allele (non-space character, eg A) this effect value corresponds to. The third should be a decimal representing the additive effect value of that allele for that marker.

A particular marker/allele combination not being included in the file is equivalent to that combination having an effect value of 0. If a particular marker/allele combination is included multiple times in a file, only the last occurence is saved. If a marker name in the file does not match any marker tracked by the simulation, that row will be ignored. The simulation will print out the number of marker/allele pairs for which effects are loaded: for the sample file above, that would be 5.

### Loading from other file formats
Updates to expand the range of allowed input formats are coming soon.


# Crossing & Other Ways to Generate New Genotypes: Plant-themed
Templates in this section assume you have loaded a set of founders whose group number is saved in the variable `founders`.

## Single seed descent

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>For six generations grow a single seed from each plant to maturity.
<td>
```{C}
GroupNum f6 = self_n_times(d, 6, founders, BASIC_OPT);
```
<td>
```{R}
f6 <- self.n.times(founders, n=6)
```
</table>

## Creating halfsib or fullsib families

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Make 20 random crosses and collect families of 6 full siblings from the result of each cross.
<td>
```{C}
GenOptions opt = {.will_name_offspring=FALSE, .offspring_name_prefix=NULL, .family_size=6,
		.will_track_pedigree=TRUE, .will_allocate_ids=TRUE,
		.filename_prefix=NULL, .will_save_pedigree_to_file=FALSE,
		.will_save_bvs_to_file=NOT_AN_EFFECT_ID, .will_save_alleles_to_file=FALSE,
		.will_save_to_simdata=TRUE};
GroupNum crosses = cross_random_individuals(d, founders, 20, 0, opt);
GroupNum families[20];
split_into_families(d, crosses, families);
```
<td>
```{R}
crosses <- cross.randomly(founders, n.crosses=20, offspring=6)
families <- break.group.into.families(crosses)
```
</table>

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Make 10 random crosses with the second founder genotype, and collect halfsib families of 6 half-siblings from the result of each cross.
<td>
```{C}
int targetparent = 1;
GroupNum targetparent_group = split_from_group(d, 1, &targetparent);

GenOptions opt = {.will_name_offspring=FALSE, .offspring_name_prefix=NULL, .family_size=6,
		.will_track_pedigree=TRUE, .will_allocate_ids=TRUE,
		.filename_prefix=NULL, .will_save_pedigree_to_file=FALSE,
		.will_save_bvs_to_file=NOT_AN_EFFECT_ID, .will_save_alleles_to_file=FALSE,
		.will_save_to_simdata=TRUE};
GroupNum crosses = cross_randomly_between(d, targetparent_group, founders, 10, 0, 0, opt);

GroupNum families[10];
split_into_halfsib_families(d, crosses, 1, families);
```
<td>
```{R}
targetparent_group <- make.group(c(1L))

crosses <- cross.randomly.between(targetparent_group, founders, n.crosses=10, offspring=6)

families <- break.group.into.halfsib.families(crosses)
```
</table>

## Updating marker effects

Updating the marker effect estimates is a task that must be done outside of genomicSimulation. Export required data from the simulation, run the marker effect re-estimation, then re-import marker effects as in **Swap out the set of additive trait effects for another**.

## Trials and locations

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose you want to simulate the success of the founder population in different environments.
<td>
```{C}
GroupNum location1 = make_clones(d, founders, TRUE, BASIC_OPT);
GroupNum location2 = make_clones(d, founders, TRUE, BASIC_OPT);
```
<td>
```{R}
location1 <- make.clones(founders)
location2 <- make.clones(founders)
```
</table>

Once separate copies exist for trials at different locations, simulate selection with different intensities or heritabilities for each different location group. Eg **Select on phenotype (simulated with a given heritability)**.

## Backcrossing

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Pull out a chosen founder and repeatedly (for 20 generations) cross back to it.
<td>
```{C}
int targetparent = 1;
GroupNum targetparent_group = split_from_group(d, 1, &targetparent);

GroupNum backcross_generations[20];
backcross_generation[0] = founders;
for (int i = 1; i < 20; ++i) {
	backcross_generations[i] = cross_randomly_between(d, targetparent_group, backcross_generations[i-1], 10, 0, 0, BASIC_OPT);
}
```
<td>
```{R}
targetparent_group <- make.group(c(1L))

backcross_generations <- rep(0L, times=20);
backcross_generations[1] <- founders
for (ii in 1:20) {
	backcross_generations[ii] <- cross.randomly.between(targetparent_group, backcross_generations[ii-1], n.crosses=10)
}
```
</table>

For marker-assisted backcrossing, select each loop on presence of the marker and/or the estimated breeding value. See **Select on a qualitative trait** and **Select on a qualitative trait, then a quantitative trait.**


# Crossing & Other Ways to Generate New Genotypes: Animal-themed
Templates in this section assume you have loaded two sets of founders whose group id are saved in the variables `cows` and `bulls`. See **Load a genetic map and several sets of founder genotypes**.

## Split offspring into male and female

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose you randomly cross your founders, then want to identify the male and female calves among the offspring.
<td>
```{C}
GroupNum offspring = cross_randomly_between(d, cows, bulls, 10, 0, 0, BASIC_OPT);

GroupNum offspring_f = split_randomly_into_two(d, offspring);
GroupNum offspring_m = offspring;
```
<td>
```{R}
offspring <- cross.randomly.between(cows, bulls, n.crosses=10)

temporary <- break.group.randomly(offspring, into.n = 2)
offspring_f <- temporary[1]
offspring_m <- temporary[2]
rm(temporary)
```
</table>

split_randomly_into_two or break.group.randomly can be used to split a group into two sub-groups by flipping a coin on each group member. This may result in two groups that are not quite the same size. To split a group of eg. 10 offspring into two groups of exactly 5 members, the cousin functions, split_evenly_into_two or break.group.evenly, can be used instead.

## Introducing fresh blood to male and female kernels
Suppose the `offspring_f` and `offspring_m` groups exist, as created in the previous section **Split offspring into male and female**.

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose you want to add your new calves to the breeding kernels from which their parents were chosen.
<td>
```{C}
GroupNum cowGroups[2];
cowGroups[0] = cows;
cowGroups[1] = offspring_f;

cows = combine_groups(d, 2, cowGroups);

GroupNum bullGroups[2];
bullGroups[0] = bulls;
bullGroups[1] = offspring_m;

bulls = combine_groups(d, 2, bullGroups);
```
<td>
```{R}
cows <- combine.groups(c(cows,offspring_f))
bulls <- combine.groups(c(bulls,offspring_m))
```
</table>

## Random mating with caps on number of offspring per animal

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose each cow should only have one calf this generation.
<td>
```{C}
GroupNum offspring = cross_randomly_between(d, cows, bulls, 10, 1, 0, BASIC_OPT);
```
<td>
```{R}
offspring <- cross.randomly.between(cows, bulls, cap1=1, n.crosses=10)
```
</table>

## Mating all females to a good male

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose the best bull in the population is the only one that will father calves under the current breeding strategy.
<td>
```{C}
GroupNum bestbull = split_by_bv(d, bulls, eff1, 1, FALSE); # where eff1 is an EffectID representing the marker effect set to use to calculate bvs

GroupNum offspring = cross_randomly_between(d, cows, bestbull, 10, 1, 0, BASIC_OPT);
```
<td>
```{R}
bestbull <- select.by.gebv(bulls, number=1, eff.set=eff1) # by default, this function uses the first effect set, so `eff.set=eff1` is optional
offspring <- cross.randomly.between(cows, bestbull, cap1=1, n.crosses=10)
```
</table>

## Making specific chosen matings

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose you want to cross Cow1 to Bull1, Cow2 to Bull2, and Daisy to Bull1
<td>
```{C}
cow1_index = get_index_of_name(d->m, "Cow1");
cow2_index = get_index_of_name(d->m, "Cow2");
cow3_index = get_index_of_name(d->m, "Daisy");
bull1_index = get_index_of_name(d->m, "Bull1");
bull2_index = get_index_of_name(d->m, "Bull2");

int crossingPlan[2][3];
crossingPlan[0][0] = cow1_index;
crossingPlan[0][1] = cow2_index;
crossingPlan[0][2] = cow3_index;
crossingPlan[1][0] = bull1_index;
crossingPlan[1][1] = bull2_index;
crossingPlan[1][2] = bull1_index;

GroupNum offspring = cross_these_combinations(d, 3, crossingPlan[0], crossingPlan[1], BASIC_OPT);
```
<td>
```{R}
offspring <- cross.combinations(c("Cow1", "Cow2", "Daisy"), c("Bull1", "Bull2", "Bull1"))
```
</table>

The other option is to create an input file of the following format:
```
Cow1 Bull1
Cow2 Bull2
Daisy Bull1
```
and call `make_crosses_from_file` or `cross.combinations.file`.

## Three-Way Crosses

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose you want to make a single Breed1/Breed2//Breed3 cross. That is, cross the F1 of a mating between Breed1 and Breed2 to Breed3. Assumes the genotypes for Breed1, Breed2, and Breed3 have those names.
<td>
```{C}
breed1_index = get_index_of_name(d->m, "Breed1");
breed2_index = get_index_of_name(d->m, "Breed2");
breed3_index = get_index_of_name(d->m, "Breed3");

int crossingPlan[2][1];
crossingPlan[0][0] = breed1_index;
crossingPlan[1][0] = breed2_index;

GroupNum f1 = cross_these_combinations(d, 1, crossingPlan[0], crossingPlan[1] BASIC_OPT);
int* f1_index = malloc(sizeof(int) * 1);
get_group_indexes(d, f1, 1, f1_index); //we know this group has only one member

// Re-use crossingPlan
crossingPlan[0][0] = breed3_index;
crossingPlan[1][0] = f1_index[0];
free(f1_index);

GroupNum f3way = cross_these_combinations(d, 1, crossingPlan[0], crossingPlan[1], BASIC_OPT);
```
<td>
```{R}
f1 <- cross.combinations(c("Breed1"), c("Breed2"))
f1_index <- see.group.data(f1, "X")
f3way <- cross.combinations(c("Breed3"), f1_index)
```
</table>

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose you want 25 offspring of a Breed1/Breed2//Breed3 cross. Assumes the genotypes for Breed1, Breed2, and Breed3 have those names.
<td>
```{C}
breed1_index = get_index_of_name(d->m, "Breed1");
breed2_index = get_index_of_name(d->m, "Breed2");
breed3_index = get_index_of_name(d->m, "Breed3");

int crossingPlan[2][1];
crossingPlan[0][0] = breed1_index;
crossingPlan[1][0] = breed2_index;

GenOptions opt = {.family_size=5,
		.will_name_offspring=FALSE, .offspring_name_prefix=NULL,
		.will_track_pedigree=TRUE, .will_allocate_ids=TRUE,
		.filename_prefix=NULL, .will_save_pedigree_to_file=FALSE,
		.will_save_bvs_to_file=FALSE, .will_save_alleles_to_file=FALSE,
		.will_save_to_simdata=TRUE};

GroupNum f1 = cross_these_combinations(d, 1, crossingPlan[0], crossingPlan[1], opt);
int* f1_indexes = malloc(sizeof(int) * 5);
get_group_indexes(d, f1, 5, f1_indexes);

int crossingPlanb[2][5];
crossingPlanb[0][0] = breed3_index; crossingPlanb[0][1] = breed3_index;
crossingPlanb[0][2] = breed3_index; crossingPlanb[0][3] = breed3_index;
crossingPlanb[0][4] = breed3_index;
crossingPlanb[1][0] = f1_indexes[0]; crossingPlanb[1][1] = f1_indexes[1];
crossingPlanb[1][2] = f1_indexes[2]; crossingPlanb[1][3] = f1_indexes[3];
crossingPlanb[1][4] = f1_indexes[4];
free(f1_indexes);

GroupNum f3way = cross_these_combinations(d, 5, crossingPlanb[0], crossingPlanb[1], opt);
```
<td>
```{R}
f1 <- cross.combinations(c("Breed1"), c("Breed2"), offspring=5)
# This time, f1 is a group of multiple individuals
f1_indexes <- see.group.data(f1, "X")
f3way <- cross.combinations(rep("Breed3", times=length(f1_indexes)), f1_indexes, offspring=5)

```
</table>


## Tracking age and culling individuals that are too old
genomicSimulation's custom labels can be used to track age (or some other known trait that can be represented with an integer for each individual).

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose you only want to breed 3-year-old animals. 
<td>
```{C}
SimData* d = create_empty_simdata(1234567);
struct GroupAndEffectSet init = load_all_simdata(d, "genotype-file.txt", "map-file.txt", NULL);
GroupNum animals = init.group;

// Create a new label to represent age, with default/at-birth value of 0.
LabelID ageLabel = create_new_label(d, 0);

// Founders are 3 years old at the beginning.
set_labels_to_const(d, animals, ageLabel, 3);

for (int year = 0; year < 10; ++year) {
	GroupNum breedingGroup = split_by_label_value(d, animals, ageLabel, 3);
	
	// Do some breeding/selection steps as appropriate for the breeding program, eg:
	GruopNum offspring = cross_random_individuals(d, breedingGroup, 50, 0, BASIC_OPT);
	// Offspring will have the default value for the label i.e. ageLabel = 0
	
	GroupNum toCombine[3] = {animals, breedingGroup, offspring};
	animals = combine_groups(d, 3, toCombine);
	
	// Increase age of all by 1
	increment_labels(d, animals, ageLabel, 1);
}
```
<td>
```{R}
animals <- load.data("genotype-file.txt", "map-file.txt")

# Create a new label to represent age, with default/at-birth value of 0.
ageLabel <- make.label(0L)

# Founders are 3 years old at the beginning.
change.label.to.this(ageLabel, 3L, animals);

for (year in 0:10) {
	breedingGroup <- make.group.from.label(ageLabel, 3L, animals);
	
	# Do some breeding/selection steps as appropriate for the breeding program, eg:
	offspring <- cross.randomly(breedingGroup, 50, 0)
	# Offspring will have the default value for the label i.e. ageLabel = 0
	
	animals <- combine.groups(c(animals, breedingGroup, offspring))
	
	# Increase age of all by 1
	change.label.by.amount(ageLabel, 1, animals);
}
```
</table>


<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose you want to cull animals that are 12 years old or older.
<td>
```{C}
GroupNum toCull = split_by_label_range(d, animals, ageLabel, 12, 1000000); 
// some very large number for upper bound. User is responsible for the values
// of the custom label, so you can work out a number your script should never 
// set the label above.

if (toCull > 0) { // split_by_label_range did find some animals to cull
	delete_group(d, toCull);
}
```
or
```{C}
GroupNum toKeep = split_by_label_range(d, animals, ageLabel, 0, 11); 

delete_group(d, animals);
```
<td>
```{R}
toCull <- make.group.from.label.range(ageLabel, 12, 1000000, animals)

if (toCull > 0) { // split_by_label_range did find some animals to cull
	delete.group(toCull) // delete from memory
	rm(toCull) // remove `toCull` label from R environment.
}
```
or
```{R}
toKeep <- make.group.from.label.range(ageLabel, 0, 11, animals)

delete.group(animals)
rm(animals)
```
</table>


# Selection
Templates in this section assume you have a group `target` from which you want to select a subset of 'good' genotypes. Scripts are provided for selecting in different ways using the 'custom selection interface'. This interface is basically the existence the function `make.group`, allowing users to use any kind of scripting/data manipulation to choose the indexes of the genotypes they want to select.

Templates will be formatted as functions, so that selecting the 'good' genotypes out of `target` is simply a matter of calling that section's selection function.

Note the concept of the 'custom selection interface' is mainly appropriate to the R version of genomicSimulation. Due to C's lack of vector manipulation/sorting functionality compared to R, implementing these selection methods in C is less about scripting and more about knowing how to sort/loop through datasets efficiently in C. The writer of this page isn't going to bother to write up templates for all of these selection tasks in C. If the task is popular enough it will probably be written up as an easier-to-use library function instead of a template, anyway.


## Manually select on true breeding value

The template selection template. This does the exact same thing as `split_by_bv` / `select.by.gebv`, but using the 'custom selection method interface'. Hopefully this aids in understanding the concept.

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose you want to select the top 10 genotypes in `target` according to their true breeding value as calculated by the internal breeding value calculator.
<td>
```{C}
GroupNum select_10_with_best_GEBV(GroupNum group) {
  unsigned int group_size = get_group_size( d, group );
  int group_indexes[group_size];
  memset(group_indexes, 0, sizeof(int)*group_size);
  get_group_indexes( d, group, group_size, group_indexes);
  double group_bvs[group_size];
  memset(group_bvs, 0, sizeof(double)*group_size);
  get_group_bvs( d, group, group_size, group_bvs);

  # A pretty inefficient way of finding the top 10 scores and their corresponding
  # indexes, but it will do for this example.
  double min_score = group_bvs[0];
  for (int i = 1; i < group_size; ++i) {
    if (group_bvs[i] < min_score) {
	  min_score = group_bvs[i];
	}
  }

  int top_individuals[10];
  int currentTopIndex = 0;
  for (int j = 0; j < 10; ++j) {
	currentTopIndex = 0;
	for (int i = 1; i < group_size; ++i) {
	  if (group_bvs[i] > group_bvs[currentTopIndex]) {
	    currentTopIndex = i;
	  }
	}
	top_individuals[j] = group_indexes[currentTopIndex];
	group_bvs[currentTopIndex] = min_score;
  }

  return split_from_group(d, 10, top_individuals);
}
```
<td>
```{R}
select.10.best.GEBV <- function(group) {
  info <- data.frame(Index=see.group.data(group,"X"),
                     GEBV=see.group.data(group,"BV"))

  selected_group <- make.group(info[order(info$GEBV, decreasing=TRUE),]$Index[1:10])
  return( selected_group )
}
```
</table>


## Select on phenotype (simulated with a given heritability)

<table>
<tr><th>Task <th>genomicSimulation (R)
<tr><td>Suppose you want to select the top 10 genotypes in `target` according to a phenotype simulated with a certain $h^2$ heritability value.
<td>
```{R}
select.top.10.phenotypes <- function(group, heritability) {
  info <- data.frame(Index=see.group.data(group,"X"),
                     GEBV=see.group.data(group,"BV"))

  # simulate phenotype = genotype + environmental variation
  # using normally distributed Ve and heritability H^2 = (Ve + Vg)/Vg
  Vg <- var(info$GEBV)
  Ve <- Vg/heritability - Vg
  info$Pheno <- info$GEBV + rnorm(length(info$GEBV), mean=0, sd = sqrt(Ve))

  # Select those with the top phenotype
  selected <- make.group(info[order(info$Pheno, decreasing=TRUE),]$Index[1:10])
  return( selected )

}
```
</table>

## Select on a qualitative trait

If you have a trait with only a small number of possible outcomes - perhaps a marker for which each individual can be homozygous, heterozygous, or carry no copies.

<table>
<tr><th>Task <th>genomicSimulation (R)
<tr><td>Suppose you want to select on a trait that has possible values
Your effect file looks like:
```
POLLED T 1
```
<td>
```{R}

```
</table>

## Select on a qualitative trait, then a quantitative trait.

Suppose you have an effect file for a quantitative trait, called `eff-large.txt`, as well as a one-line effect file for a qualitative trait, as described in the above section **Select on a qualitative trait**, named `eff-small.txt`.


## Select on an index weighting two traits.



(this list of examples will continue to be expanded...)

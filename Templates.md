Features and Templates Guide        {#templates}
============================

While [the documentation pages of genomicSimulationC](@ref modules) and its counterpart [genomicSimulation](https://github.com/vllrs/genomicSimulation) do describe all the features the simulation tool offers, they don't do so in a way that is necessarily easy to browse or that's helpful for new users to understand the simulation tool. 

This page provides a set of templates for implementing common breeding project actions in genomicSimulationC and its R counterpart. 

[TOC]

# Setting Up

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
SimData* d = create_empty_simdata();
int founders = load_all_simdata(d, "genotype-file.txt", "map-file.txt", NULL);
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
SimData* d = create_empty_simdata();
int founders_a = load_all_simdata(d, "genotype-file.txt", "map-file.txt", NULL);
int founders_b = load_more_transposed_genes_to_simdata(d, "genotype-file2.txt");
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
SimData* d = create_empty_simdata();
int founders = load_all_simdata(d, "genotype-file.txt", "map-file.txt", "eff-file.txt");
```
<td>
```{R}
founders <- load.data(allele.file="genotype-file.txt", map.file="map-file.txt", effect.file="eff-file.txt")
```
</table>

## Swap out the set of additive trait effects for another

This assumes the simulation is set up, that is, one of **Load a genetic map and a set of founder genotypes**, **Load a genetic map and several sets of founder genotypes**, **Load a genetic map, a set of founder genotypes, and a set of additive trait effects** or equivalent has been carried out. If the simulation is not yet set up and you with to load trait effect values, see the section above, titled **Load a genetic map, a set of founder genotypes, and a set of additive trait effects**.

<table>
<tr><th>Task <th>Input Files <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td> 
<td>Effect file:
- File location: eff-file2.txt
<td> 
```{C}
load_effects_to_simdata(d, "eff-file2.txt");
```
<td>
```{R}
load.different.effects("eff-file2.txt")
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

Loading an effect file is optional for running the simulation. Its format should be:
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
Templates in this section assume you have loaded a set of founders whose group id is saved in the variable `founders`. 

## Single seed descent

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>For six generations grow a single seed from each plant to maturity.
<td> 
```{C}
int f6 = self_n_times(d, 6, founders, BASIC_OPT);
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
		.will_save_bvs_to_file=FALSE, .will_save_alleles_to_file=FALSE,
		.will_save_to_simdata=TRUE};
int crosses = cross_random_individuals(d, founders, 20, 0, opt);
int families[20];
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
int targetparent_group = split_from_group(d, 1, &targetparent);

GenOptions opt = {.will_name_offspring=FALSE, .offspring_name_prefix=NULL, .family_size=6,
		.will_track_pedigree=TRUE, .will_allocate_ids=TRUE,
		.filename_prefix=NULL, .will_save_pedigree_to_file=FALSE,
		.will_save_bvs_to_file=FALSE, .will_save_alleles_to_file=FALSE,
		.will_save_to_simdata=TRUE};
int crosses = cross_randomly_between(d, targetparent_group, founders, 10, 0, 0, opt);		

int families[10];
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
int location1 = make_clones(d, founders, TRUE, BASIC_OPT);
int location2 = make_clones(d, founders, TRUE, BASIC_OPT);
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
int targetparent_group = split_from_group(d, 1, &targetparent);

int backcross_generations[20];
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
int offspring = cross_randomly_between(d, cows, bulls, 10, 0, 0, BASIC_OPT);

int offspring_f = split_randomly_into_two(d, offspring);
int offspring_m = offspring;
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

split_randomly_into_two or break.group.randomly can be used to split a group into two sub-groups by flipping a coin on each group member. This may result in two groups that are not quite the same size. To split a group of eg. 10 offspring into two groups of exactly 5 members, the cousin functions of split_evenly_into_two or break.group.evenly can be used instead.

## Introducing fresh blood to male and female kernels
Suppose the `offspring_f` and `offspring_m` groups exist, as created in the previous section **Split offspring into male and female**. 

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose you want to add your new calves to the breeding kernels from which their parents were chosen. 
<td> 
```{C}
int cowGroups[2];
cowGroups[0] = cows;
cowGroups[1] = offspring_f;

cows = combine_groups(d, 2, cowGroups);

int bullGroups[2];
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
int offspring = cross_randomly_between(d, cows, bulls, 10, 1, 0, BASIC_OPT);
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
int bestbull = split_by_bv(d, bulls, 1, FALSE);

int offspring = cross_randomly_between(d, cows, bestbull, 10, 1, 0, BASIC_OPT);
```
<td>
```{R}
bestbull <- select.by.gebv(bulls, number=1)
offspring <- cross.randomly.between(cows, bestbull, cap1=1, n.crosses=10)
```
</table>

## Making specific chosen matings

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose you want to cross Cow1 to Bull1, Cow2 to Bull2, and Daisy to Bull1
<td> 
```{C}
int namesNeeded[5];
int IDsNeeded[5];
namesNeeded[0] = "Cow1";
namesNeeded[1] = "Cow2";
namesNeeded[2] = "Daisy";
namesNeeded[3] = "Bull1";
namesNeeded[4] = "Bull2";
get_ids_of_names(d->m, 5, namesNeeded, IDsNeeded);

int crossingPlan[2][3];
crossingPlan[0][0] = IDsNeeded[0];
crossingPlan[0][1] = IDsNeeded[1];
crossingPlan[0][2] = IDsNeeded[2];
crossingPlan[1][0] = IDsNeeded[3];
crossingPlan[1][1] = IDsNeeded[4];
crossingPlan[1][2] = IDsNeeded[3];

int offspring = cross_these_combinations(d, 3, crossingPlan, BASIC_OPT);
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

## Tracking age and discarding individuals that are 'too old'
Updates to add an age/custom flag to each individual is coming soon.


# Selection 

## Manually select on true breeding value

The template selection template. This does the exact same thing as split_by_bv / select.by.gebv, but using the 'custom selection method interface'. Hopefully this aids in understanding the concept.


## Select on phenotype (simulated with a given heritability)

## Select on a qualitative trait

## Select on a qualitative trait, then a quantitative trait.

## Select on an index weighting two traits.



 











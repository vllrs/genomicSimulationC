Features and Templates Guide        {#templates}
============================

While [the documentation pages of genomicSimulationC](@ref modules) and its counterpart [genomicSimulation](https://github.com/vllrs/genomicSimulation) do describe all the features the simulation tool offers, they don't do so in a way that facilitates quick prototyping or help new users get insight into setting up simulations.

This page provides a set of templates for implementing common breeding project actions in genomicSimulationC and its R counterpart.

[TOC]

# Setting Up

## Input file formats:
genomicSimulation has three types of input files: a genotype matrix, a genetic map, and a list of marker effects. If multiple files are provided at once in a call to `load_data_files` (C)/`load.data` (R), the map file will be loaded first, followed by the genotype matrix and marker effect file.

### Genetic map files

The simplest genetic map file is formatted as follows:
```
marker chr pos
m3 3 15
m2 1 8.3
m1 1 5.2
```
Other valid genetic map files might include:
```
chr marker pos
1A 1243509 173.2
1A 2350898 462.2
1B 4360992 32.009
2A 1243556 243.5
```
or
```
gene 10 3.24
othergene 10 8.3e-1
etc 15 1.203e2
```

The header line is optional. If there is a header line, the three columns "marker" "chr" and "pos" may be rearranged into any ordering. If the header is not provided, the order is assumed to be marker name in the first column, followed by chromosome in the second column, followed by position in the third.

Regarding marker names: For all genomicSimulation features to work, genetic markers must have names. There is no issue with purely numeric names.

Regarding chromosomes/linkage group: any alphanumeric combination may be used to denote a chromosome/linkage group, for example '9' or '1A'. There is no limit on the number of unique chromosomes/linkage groups in the map.

Regarding marker position: this represents the position in centimorgans of each marker along its linkage group. (Distance of 1 cM = expected probability of 0.01 chromosomal crossovers in that range.)

The cells in the map file may be separated by spaces, tabs, commas, or any combination thereof. Cell spacers do not need to be consistent across the file, but therefore marker names/linkage group names cannot contain spaces, tabs, or commas. 

The order in which genetic markers are presented in the file does not matter. If a marker is duplicated in the file, it will be recorded twice as two separate markers in the same position and with the same name.  However, because loading genotypes depends on searching for markers by name, if two markers have the same name, genotypes may not be loaded correctly. Therefore it is suggested that all markers have unique names.

### Genotype matrix files

The simplest genotype matrix file is formatted as follows:
```
name	G01	G02	G03	G04	G05	G06
m1	TT	TT	TT	TA	TT	AT
m3	TT	TT	TA	TA	TT	TT
m2	AA	AA	AA	AA	CC	AA
```
where G01, G02, ..., are names of the founder genotypes; m1, m2, ..., are the genetic markers; and entries in the matrix are the alleles that the founder genotypes have at those markers. 

Other valid genotype matrix files might include:
```
 m100, m101, m102
cand1,0,0,1
cand2,1,2,2
cand4,2,1,2
```
or
```
marker1	T/T	T/T	T/T	T/A	T/T	A/T
marker3	T/T	T/T	T/A	T/A	T/T	T/T
marker2	A/A	A/A	A/A	A/A	T/T	A/A
```

The genotype matrix can be row-major or column-major (that is, the genetic markers may be rows, or columns). The program will determine the orientation by attempting to match row and column headers with names of markers tracked by the simulation (extracted from a genetic map file). If the simulation does not yet have a list of tracked markers (that is, no genetic map has been simultaneously or previously provided), then it defaults to assuming that rows represent genetic markers.  

The order in which genetic markers are presented in the file does not matter. Candidate genotypes, however, will be saved internally in the simulation in the order that they appear in the file. Candidate genotype names do not need to be unique. Candidate genotype names are optional, if candidates are stored as columns. genomicSimulation cannot read a gentoype matrix with no row headers, so if candidates are stored as rows, then they are required to have names.

Genetic markers' names are required. 

When there are both row and column headers, the first/corner cell (the one that contained the value "name" in the first example above) can be deleted or can contain any text. Its value will be ignored. 

The table cells in the genotype matrix file may be separated by spaces, tabs, commas, or any combination thereof. Cell spacers do not need to be consistent across the file. 

There are several options for how the alleles at each marker can be presented. All allele pair cells in the same genotype matrix must be in the same format. The format of the allele pairs in the genotype matrix is automatically detected. There are four acceptable formats for allele pairs in the genotype matrix:

1. Any pair of characters (eg. "AA", "TA", "nW"). Each character is an allele. This is a format that specifies allele phase (ie. "AT" and "TA" are different genotypes).
2. Any pair of characters, separated by a forwards slash "/" character (eg. "A/A", "T/A", "n/W"). The two characters either side of the slash are the alleles. This is a format that specifies allele phase (ie. "AT" and "TA" are different genotypes).
3. Alternate allele counts ("0", "1", "2"). This is a format that does not specify allele phase, so phase of heterozygotes will be randomised when loaded. The counts represent the number of copies of the alternate allele (stored inside genomicSimulation as "A", while the reference allele is stored as "T"). "0" then corresponds to "TT", "2" to "AA", and "1" will be randomised as either "TA" or "AT". Corresponding marker effect files for calculating breeding values must use allele "A" consistently to represent the alternate allele, and "T" to represent the reference. (alternatively, `change_allele_symbol` (C)/`change.allele.symbol` (R) can be used after loading to change T and/or A to other symbols.)
4. (subset of) IUPAC encodings of DNA bases ("G", "A", "T", "C", "R", "Y", "M", "K", "S", "W", "N"). This is a format that does not specify allele phase, so phase of heterozygotes will be randomised when loaded. See below for meanings of the IUPAC encoding symbols that genomicSimulation can parse.

<table>
<tr><th>IUPAC symbol <th>genomicSimulation genotype
<tr><td>G<td>GG
<tr><td>A<td>AA
<tr><td>T<td>TT
<tr><td>C<td>CC
<tr><td>R<td>GA or AG
<tr><td>Y<td>TC or CT
<tr><td>M<td>AC or CA
<tr><td>K<td>GT or TG
<tr><td>S<td>GC or CG
<tr><td>W<td>AT or TA
<tr><td>N<td>\\0\\0 (nulls represent unknown genotype)
</table>

**Note you might have a genotype matrix of only "AA", "AT", and "TT". This uses "alternate allele counts"-style encoding (like format 3) but presents it in a format that looks like pairs of alleles (format 1).** genomicSimulation expects allele pair encodings to include haplotype phase, (that is, to have four possible values for genotypes of two alleles, not three: eg. "AA", "AT", "TA", and "TT" instead of just "AA", "AT", "TT"). 

Two options for loading a dataset with non-phased "AA"/"AT"/"TT" allele pairs are:

- Use “haplotyping”/”haplotype phasing”/”haplotype inference” software to infer whether heterozygotes are "AT" or "TA", before loading into genomicSimulation.
- Find-and-replace "TT" with "0", "AT" with "1", and "AA" with "2" before loading into genomicSimulation. genomicSimulation will then randomise the phase of each haplotype.

### Marker effect files

Loading marker effect file(s) is optional for running the simulation. The simplest marker effect file is formatted as follows:
```
marker allele eff
m1 A -0.8
m2 A -0.1
m3 A 0.1
m1 T 0.9
m3 T -0.1
```
Other valid marker effect files might include:
```
marker eff allele
m1243509 0.1 A
m1243509 -0.1 T
m2350898 0.15 T
m2350898 -0.1 A
```
or
```
specialgene G 1.0
```

A single marker effect file represents the allele effects for one trait. Multiple files like this can be loaded in order to calculate breeding values for multiple traits separately.

The header line is optional. If there is a header line, the three columns "marker" "allele" and "eff" may be rearranged into any ordering. If the header is not provided, the order is assumed to be marker name in the first column, followed by allele, followed by additive effect value.

The order of the rows in the file does not matter. If multiple rows exist in the file for the same marker name/allele combination, only the last row's additive effect value will be saved. 

Regarding marker names: Only rows in this file whose marker name matches a marker name from a previously- or simultaneously-loaded map file or genotype matrix file will be loaded. 

Regarding allele: This column should be the allele (non-space character, eg "A") that this effect value corresponds to. 

Regarding allele effect: This column should be a decimal representing the additive effect value of the same-row allele for the same-row marker. It can be positive or negative or zero, and can be represented by an integer, a decimal, or scientific notation. 

Cells may be separated by spaces, tabs, commas, or any combination thereof. Cell spacers do not need to be consistent across the file. 

### Loading from other file formats
Updates to expand the range of allowed input formats are coming soon.



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
struct MultiIDSet init = load_data_files(d, "genotype-file.txt", "map-file.txt", NULL, DETECT_FILE_FORMAT);
GroupNum founders = init.group;
```
<td>
```{R}
init <- load.data(genotype.file="genotype-file.txt", map.file="map-file.txt")
founders = init$groupNum
```
</table>

The simulation will track genotypes across all markers present in the map file. Which markers are tracked by the simulation cannot currently be changed. Line order does not matter in any file.

Consider using the get_group_ family of functions (genomicSimulationC) and see.group.data function (genomicSimulation) to confirm data is correctly loaded.

## Load a genetic map and genotype file, specifying file format

By default, genomicSimulation attempts to detect the format of input files as it loads them. When a file loading command is run, genomicSimulation prints out its assumptions about the format of the file.

If the printed logs are incorrect about the layout of your map file or effect file, the most likely cause is a mis-spelling in the header row of the file. If the printed logs are incorrect about the size of the table in your map file or effect file, the issue may be that the numbers cannot be detected as numbers, or that there are issues with the column spacers/line breaks that cause there to be more than three entries to a row.

If the printed logs are incorrect about the layout or format of your genotype matrix file, you can manually specify the file format for genomicSimulation, as follows. 

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
FileFormatSpec mformat = define_matrix_format_details(GSC_TRUE,GSC_TRUE,GSC_GENOTYPECELLSTYLE_PAIR);
struct MultiIDSet init = load_data_files(d, "genotype-file.txt", "map-file.txt", NULL, mformat);
```
<td>
```{R}
mformat <- define.matrix.format.details(has.header=TRUE, markers.as.rows=TRUE, cell.style="Pair")
init <- load.data(genotype.file="genotype-file.txt", map.file="map-file.txt",format=mformat)
```
</table>

It is not required that you manually specify all details of the file format. Whatever details of the file format are not manually specified will be automatically detected.

See the documentation for the function `define_matrix_format_details` (C)/`define.matrix.format.details` (R) for more information.

## Load a genetic map and several sets of founder genotypes

<table>
<tr><th>Task <th>Input Files <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>
<td>Genotype file/Allele file:
- File location: genotype-file.txt
- File contents:
```
	G01	G02	G03	G04	G05
m1	TT	TT	TT	TA	TT
m3	TT	TT	TA	TA	TT
m2	AA	AA	AA	AA	TT
```

Genotype file 2:
- File location: genotype-file2.txt
- File contents:
```
 G06
m1,AT
m2,AA
m3,TT
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
struct MultiIDSet init = load_data_files(d, "genotype-file.txt", "map-file.txt", NULL, DETECT_FILE_FORMAT);
GroupNum founders_a = init.group;
GroupNum founders_b = load_genotypefile(d, "genotype-file2.txt", DETECT_FILE_FORMAT);
```
<td>
```{R}
init <- load.data(genotype.file="genotype-file.txt", map.file="map-file.txt")
founders_a <- init$groupNum
founders_b <- load.genotypes("genotype-file2.txt")
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
- File contents:
```
marker allele eff
m1 A -0.8
m2 A -0.1
m3 A 0.1
m1 T 0.9
m3 T -0.1
```
<td>
```{C}
SimData* d = create_empty_simdata(time(NULL));
struct MultiIDSet init = load_data_files(d, "genotype-file.txt", "map-file.txt", "eff-file.txt", DETECT_FILE_FORMAT);
GroupNum founders = init.group;
EffectID eff1 = init.effectSet;
MapID map1 = init.map;
```
<td>
```{R}
init <- load.data(genotype.file="genotype-file.txt", map.file="map-file.txt", effect.file="eff-file.txt")
founders <- init$groupNum
eff1 <- init$effectID
map1 <- init$mapID
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
EffectID eff2 = load_effectfile(d, "eff-file2.txt");
```
<td>
```{R}
eff2 <- load.effects("eff-file2.txt")
```
</table>

## Load another linkage map

Markers in this map that are not present in the primary linkage map will be ignored. (That is, the first genetic map loaded must contain all markers you wish to be tracked by the simulation. Later genetic maps can only provide different recombination frequencies or chromosome allocations).

<table>
<tr><th>Task <th>Input Files <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>
<td>Map file:
- File location: mapfile2.txt
<td>
```{C}
MapID map2 = load_mapfile(d, "mapfile2.txt");
```
<td>
```{R}
map2 <- load.map("mapfile2.txt")
```
</table>

# Crossing & Other Ways to Generate New Genotypes: Plant-themed
Templates in this section assume you have loaded a set of founders whose group number is saved in the variable `founders`.

## Single seed descent

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>For six generations grow a single seed from each plant to maturity.
<td>
```{C}
GroupNum f6 = self_n_times(d, 6, founders, NO_MAP, BASIC_OPT);
```
<td>
```{R}
f6 <- self.n.times(founders, n=6)
```
</table>

Using NO_MAP for the recombination map parameter in any crossing function in the C version causes the tool to default to the first/primary recombination map. This is the same default as in the R version.

## Creating halfsib or fullsib families

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Make 20 random crosses and collect families of 6 full siblings from the result of each cross.
<td>
```{C}
GenOptions opt = {.will_name_offspring=FALSE, .offspring_name_prefix=NULL, .family_size=6,
		.will_track_pedigree=TRUE, .will_allocate_ids=TRUE,
		.filename_prefix=NULL, .will_save_pedigree_to_file=FALSE,
		.will_save_bvs_to_file=NO_EFFECTSET, .will_save_alleles_to_file=FALSE,
		.will_save_to_simdata=TRUE};
GroupNum crosses = make_random_crosses(d, founders, 20, 0, NO_MAP, opt);
GroupNum families[20];
split_into_families(d, crosses, families);
```
<td>
```{R}
crosses <- make.random.crosses(founders, n.crosses=20, offspring=6)
families <- break.group.into.families(crosses)
```
</table>

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Make 10 random crosses with the second founder genotype, and collect halfsib families of 6 half-siblings from the result of each cross.
<td>
```{C}
int targetparent = 1;
GroupNum targetparent_group = make_group_from(d, 1, &targetparent);

GenOptions opt = {.will_name_offspring=FALSE, .offspring_name_prefix=NULL, .family_size=6,
		.will_track_pedigree=TRUE, .will_allocate_ids=TRUE,
		.filename_prefix=NULL, .will_save_pedigree_to_file=FALSE,
		.will_save_bvs_to_file=NO_EFFECTSET, .will_save_alleles_to_file=FALSE,
		.will_save_to_simdata=TRUE};
GroupNum crosses = make_random_crosses_between(d, targetparent_group, founders, 10, 0, 0, NO_MAP, NO_MAP, opt);

GroupNum families[10];
split_into_halfsib_families(d, crosses, 1, families);
```
<td>
```{R}
targetparent_group <- make.group(c(1L))

crosses <- make.random.crosses.between(targetparent_group, founders, n.crosses=10, offspring=6)

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
GroupNum targetparent_group = make_group_from(d, 1, &targetparent);

GroupNum backcross_generations[20];
backcross_generation[0] = founders;
for (int i = 1; i < 20; ++i) {
	backcross_generations[i] = make_random_crosses_between(d, targetparent_group, backcross_generations[i-1], 10, 0, 0, NO_MAP, NO_MAP, BASIC_OPT);
}
```
<td>
```{R}
targetparent_group <- make.group(c(1L))

backcross_generations <- rep(0L, times=20);
backcross_generations[1] <- founders
for (ii in 1:20) {
	backcross_generations[ii] <- make.random.crosses.between(targetparent_group, backcross_generations[ii-1], n.crosses=10)
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
GroupNum offspring = make_random_crosses_between(d, cows, bulls, 10, 0, 0, NO_MAP, NO_MAP, BASIC_OPT);

GroupNum offspring_f = split_randomly_into_two(d, offspring);
GroupNum offspring_m = offspring;
```
<td>
```{R}
offspring <- make.random.crosses.between(cows, bulls, n.crosses=10)

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
GroupNum offspring = make_random_crosses_between(d, cows, bulls, 10, 1, 0, NO_MAP, NO_MAP, BASIC_OPT);
```
<td>
```{R}
offspring <- make.random.crosses.between(cows, bulls, cap1=1, n.crosses=10)
```
</table>

## Random mating between the best performing individuals in a group

genomicSimulation used to have a function `make_n_crosses_from_top_m_percent`, which would identify a highest-performing slice of a group (of a chosen size), make random crosses between them, but leave the high performers as members of their original group. This function was removed because . Here is an implementation that can be used to replicated its functionality:

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Function to make n crosses from the m percent of candidates in a group that have the highest breeding values.
<td>
```{C}
GroupNum make_n_crosses_from_top_m_percent(SimData* d, GroupNum* group, const int n, const float m) {
	# Split off the best candidates 
	int mcount = (int)(get_group_size(d,*group) * m / 100); # close enough, you could round differently if you preferred
	GroupNum best_candidates = split_by_bv(d, *group, NO_EFFECTSET, mcount, 0);
	
	# Make random crosses to get offspring 
	GroupNum offspring = make_random_crosses(d, best_candidates, n, 0, NO_MAP, BASIC_OPT);
	
	# Put the best candidates back into their original group
	# Note: you cannot assume that the group number representing "group"
	# will not change during a combine_groups operation
	GroupNum tmp[2] = {*group, best_candidates};
	*group = combine_groups(d, 2, tmp);
	
	return offspring;
}
```
<td>
```{R}
make.n.crosses.from.top.m.percent <- function(group, n, m) {
	# Split off the best candidates 
	best.candidates <- break.group.by.gebv(group, percentage=m)
	
	# Make random crosses to get offspring 
	offspring <- make.random.crosses(best.candidates, n.crosses=n)
	
	# Put the best candidates back into their original group
	# Note: you cannot assume that the group number representing "group"
	# will not change as a result of this line
	group <<- combine.groups(c(group, best.candidates))
	
	return(offspring)
}
```
</table>

## Mating all females to a good male

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose the best bull in the population is the only one that will father calves under the current breeding strategy.
<td>
```{C}
GroupNum bestbull = split_by_bv(d, bulls, eff1, 1, FALSE); # where eff1 is an EffectID representing the marker effect set to use to calculate bvs

GroupNum offspring = make_random_crosses_between(d, cows, bestbull, 10, 1, 0, NO_MAP, NO_MAP, BASIC_OPT);
```
<td>
```{R}
bestbull <- break.group.by.gebv(bulls, number=1, eff.set=eff1) # by default, this function uses the first effect set, so `eff.set=eff1` is optional
offspring <- make.random.crosses.between(cows, bestbull, cap1=1, n.crosses=10)
```
</table>

## Making specific chosen matings

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose you want to cross Cow1 to Bull1, Cow2 to Bull2, and Daisy to Bull1
<td>
```{C}
cow1_index = gsc_get_index_of_name(d->m, "Cow1");
cow2_index = gsc_get_index_of_name(d->m, "Cow2");
cow3_index = gsc_get_index_of_name(d->m, "Daisy");
bull1_index = gsc_get_index_of_name(d->m, "Bull1");
bull2_index = gsc_get_index_of_name(d->m, "Bull2");

int crossingPlan[2][3];
crossingPlan[0][0] = cow1_index;
crossingPlan[0][1] = cow2_index;
crossingPlan[0][2] = cow3_index;
crossingPlan[1][0] = bull1_index;
crossingPlan[1][1] = bull2_index;
crossingPlan[1][2] = bull1_index;

GroupNum offspring = make_targeted_crosses(d, 3, crossingPlan[0], crossingPlan[1], NO_MAP, NO_MAP, BASIC_OPT);
```
<td>
```{R}
offspring <- make.targeted.crosses(c("Cow1", "Cow2", "Daisy"), c("Bull1", "Bull2", "Bull1"))
```
</table>

The other option is to create an input file of the following format:
```
Cow1 Bull1
Cow2 Bull2
Daisy Bull1
```
and call `make_crosses_from_file` or `make.crosses.from.file`.

## Three-Way Crosses

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose you want to make a single Breed1/Breed2//Breed3 cross. That is, cross the F1 of a mating between Breed1 and Breed2 to Breed3. Assumes the genotypes for Breed1, Breed2, and Breed3 have those names.
<td>
```{C}
breed1_index = gsc_get_index_of_name(d->m, "Breed1");
breed2_index = gsc_get_index_of_name(d->m, "Breed2");
breed3_index = gsc_get_index_of_name(d->m, "Breed3");

int crossingPlan[2][1];
crossingPlan[0][0] = breed1_index;
crossingPlan[1][0] = breed2_index;

GroupNum f1 = make_targeted_crosses(d, 1, crossingPlan[0], crossingPlan[1], NO_MAP, NO_MAP, BASIC_OPT);
int* f1_index = malloc(sizeof(int) * 1);
get_group_indexes(d, f1, 1, f1_index); //we know this group has only one member

// Re-use crossingPlan
crossingPlan[0][0] = breed3_index;
crossingPlan[1][0] = f1_index[0];
free(f1_index);

GroupNum f3way = make_targeted_crosses(d, 1, crossingPlan[0], crossingPlan[1], NO_MAP, NO_MAP, BASIC_OPT);
```
<td>
```{R}
f1 <- make.targeted.crosses(c("Breed1"), c("Breed2"))
f1_index <- see.group.data(f1, "X")
f3way <- make.targeted.crosses(c("Breed3"), f1_index)
```
</table>

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose you want 25 offspring of a Breed1/Breed2//Breed3 cross. Assumes the genotypes for Breed1, Breed2, and Breed3 have those names.
<td>
```{C}
breed1_index = gsc_get_index_of_name(d->m, "Breed1");
breed2_index = gsc_get_index_of_name(d->m, "Breed2");
breed3_index = gsc_get_index_of_name(d->m, "Breed3");

int crossingPlan[2][1];
crossingPlan[0][0] = breed1_index;
crossingPlan[1][0] = breed2_index;

GenOptions opt = {.family_size=5,
		.will_name_offspring=FALSE, .offspring_name_prefix=NULL,
		.will_track_pedigree=TRUE, .will_allocate_ids=TRUE,
		.filename_prefix=NULL, .will_save_pedigree_to_file=FALSE,
		.will_save_bvs_to_file=NO_EFFECTSET, .will_save_alleles_to_file=FALSE,
		.will_save_to_simdata=TRUE};

GroupNum f1 = make_random_crosses(d, 1, crossingPlan[0], crossingPlan[1], NO_MAP, NO_MAP, opt);
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

GroupNum f3way = make_targeted_crosses(d, 5, crossingPlanb[0], crossingPlanb[1], NO_MAP, NO_MAP, opt);
```
<td>
```{R}
f1 <- make.targeted.crosses(c("Breed1"), c("Breed2"), offspring=5)
# This time, f1 is a group of multiple individuals
f1_indexes <- see.group.data(f1, "X")
f3way <- make.targeted.crosses(rep("Breed3", times=length(f1_indexes)), f1_indexes, offspring=5)

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
struct MultiIDSet init = load_data_files(d, "genotype-file.txt", "map-file.txt", NULL, DETECT_FILE_FORMAT);
GroupNum animals = init.group;

// Create a new label to represent age, with default/at-birth value of 0.
LabelID ageLabel = create_new_label(d, 0);

// Founders are 3 years old at the beginning.
change_label_to(d, animals, ageLabel, 3);

for (int year = 0; year < 10; ++year) {
	GroupNum breedingGroup = split_by_label_value(d, animals, ageLabel, 3);
	
	// Do some breeding/selection steps as appropriate for the breeding program, eg:
	GruopNum offspring = make_random_crosses(d, breedingGroup, 50, 0, NO_MAP, BASIC_OPT);
	// Offspring will have the default value for the label i.e. ageLabel = 0
	
	GroupNum toCombine[3] = {animals, breedingGroup, offspring};
	animals = combine_groups(d, 3, toCombine);
	
	// Increase age of all by 1
	change_label_by_amount(d, animals, ageLabel, 1);
}
```
<td>
```{R}
init <- load.data("genotype-file.txt", "map-file.txt")
animals <- init$groupNum

# Create a new label to represent age, with default/at-birth value of 0.
ageLabel <- create.new.label(0L)

# Founders are 3 years old at the beginning.
change.label.to.this(ageLabel, 3L, animals);

for (year in 0:10) {
	breedingGroup <- break.group.by.label.value(ageLabel, 3L, animals);
	
	# Do some breeding/selection steps as appropriate for the breeding program, eg:
	offspring <- make.random.crosses(breedingGroup, 50, 0)
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
toCull <- break.group.by.label.range(ageLabel, 12, 1000000, animals)

if (toCull > 0) { // split_by_label_range did find some animals to cull
	delete.group(toCull) // delete from memory
	rm(toCull) // remove `toCull` label from R environment.
}
```
or
```{R}
toKeep <- break.group.by.label.range(ageLabel, 0, 11, animals)

delete.group(animals)
rm(animals)
```
</table>

## Species with differing male and female recombination rates

<table>
<tr><th>Task <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>Suppose you only want to simulate a cross in a species with differing male and female recombination maps.
<td>
```{C}
SimData* d = create_empty_simdata(7654321);
MapID female_map = load_mapfile("fmap.txt");
MapID male_map = load_mapfile("mmap.txt");
GroupNum female_pop = load_genotypefile("fgenos.txt", DETECT_FILE_FORMAT);
GroupNum male_pop = load_genotypefile("mgenos.txt", DETECT_FILE_FORMAT);

GroupNum offspring = make_random_crosses_between(d, female_pop, male_pop, 1, 0, 0, female_map, male_map, BASIC_OPT);
```
<td>
```{R}
female_map <- load.map("fmap.txt")
male_map <- load.map("mmap.txt")
female_pop <- load.genotypes("fgenos.txt")
male_pop <- load.genotypes("mgenos.txt") 

offspring <- make.random.crosses.between(female_pop, male_pop, 1, map1=female_map, map2=male_map)
```
</table>


# Selection
Templates in this section assume you have a group `target` from which you want to select a subset of 'good' genotypes. Scripts are provided for selecting in different ways using the 'custom selection interface'. This interface is basically the existence the function `make.group`, allowing users to use any kind of scripting/data manipulation to choose the indexes of the genotypes they want to select.

Templates will be formatted as functions, so that selecting the 'good' genotypes out of `target` is simply a matter of calling that section's selection function.

Note the concept of the 'custom selection interface' is mainly appropriate to the R version of genomicSimulation. Due to C's lack of vector manipulation/sorting functionality compared to R, implementing these selection methods in C is less about scripting and more about knowing how to sort/loop through datasets efficiently in C. The writer of this page isn't going to bother to write up templates for all of these selection tasks in C. If the task is popular enough it will probably be written up as an easier-to-use library function instead of a template, anyway.


## Manually select on true breeding value

The template selection template. This does the exact same thing as `split_by_bv` / `break.group.by.gebv`, but using the 'custom selection method interface' of `make_group`/`make.group`. The following sections will show how to simulate different selection strategies using this procedure.

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

  return make_group_from(d, 10, top_individuals);
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
  # using normally distributed Ve and heritability H^2 = Vg/(Vg + Ve)  
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
TARGETGENE T 1
```
<td>

</table>

## Select on a qualitative trait, then a quantitative trait.

Suppose you have an effect file for a quantitative trait, called `eff-large.txt`, as well as a one-line effect file for a qualitative trait, as described in the above section **Select on a qualitative trait**, named `eff-small.txt`.


## Select on an index weighting two traits.



(this list of examples will continue to be expanded...)

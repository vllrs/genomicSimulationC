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

## Creating halfsib or fullsib families

## Updating marker effects

## Trials and locations

## Backcrossing


# Crossing & Other Ways to Generate New Genotypes: Animal-themed

## Split offspring into male and female

## Introducing fresh blood to separate male and female kernels

## Random mating with caps on number of offspring per animal

## Mating all females to a good male

## Making specific chosen matings

## Tracking age and discarding individuals that are 'too old'
Updates to add an age/custom flag to each individual is coming soon.


# Selection 

## Manually select on true breeding value

The template selection template. This does the exact same thing as split_by_bv / select.by.gebv, but using the 'custom selection method interface'. Hopefully this aids in understanding the concept.


## Select on phenotype (simulated with a given heritability)

## Select on a qualitative trait

## Select on a qualitative trait, then a quantitative trait.

## Select on an index weighting two traits.



 











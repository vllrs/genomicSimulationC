Features and Templates Guide        {#templates}
============================

This page provides a set of templates for implementing common breeding project actions in genomicSimulationC and its R counterpart.

[TOC]

# Setting Up

## Input file formats:
genomicSimulation has three types of input files: a genotype matrix, a genetic map, and a list of marker effects. At minimum, a genotype matrix is required to run simulations, because genomicSimulation does not generate you founder genotypes.

The list of genetic markers tracked by the simulation is fixed after the simulation is set up (until you reset the simulation (R)/create another SimData object (C)). No other markers can be added and none can be removed by input files loaded after the simulation is set up. If a map file is the first file loaded, or the simulation is set up using the `load_data_files` (C)/`load.data` (R) function, then the markers listed in that genetic map file are the ones tracked by the simulation. If `load_genotypefile` (C)/`load.genotypes` (R) is the first function run, before loading a map, then the simulation tracks only the list of genetic markers from the genotype matrix, and the first genetic map in the simulation is an invented map with equally-spaced markers along a single chromosome.

The input file formats are described below.

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



## Setting up the simulation and the random number generator
In C, you should manually set the seed for the random number generator when you create the SimData object, before loading data. The R version borrows R's own random number generator so has no need for setting a seed.

<table>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
```{C}
SimData* d = create_empty_simdata(time(NULL));
```
<td>
N/A
</table>

In the remainder of this page, the variable `d` in genomicSimulationC scripts will always represent a pointer to the simulation object holding the current simulation's data.

In the R version, where only a single simulation data object exists at a time, the simulation is automatically set up at the first call to `load.data` or `load.genotypes` or `load.map`, and is reset whenever the functions `clear.simdata` or `load.data` are called. (Unlike the other loading functions, the function `load.data` can only be used as the first call to set up a simulation, as it deletes and clears any previous simulation data when it is run).

## Load a genetic map and a set of founder genotypes

<table>
<tr><th>Example Input Files <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
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
init <- load.data("genotype-file.txt", "map-file.txt")
founders = init$groupNum
```
</table>

As these functions run, they will print out debug information, including whether they detect header rows in the file, how many rows they read from the file, and how many rows had to be skipped for malformed/malformatted data. You can read this information to check that the input files are being read as expected.

The first list of genetic markers loaded into a simulation object, whether that list comes from a genetic map file or from the headers of a genotype matrix file, is the total list of all genetic markers that will be tracked by this simulation object. During the lifetime of the simulation, the list of markers that are being tracked cannot be changed or added to. 

If, as in the table above, a genetic map and a genotype matrix are loaded in the same command (eg. `load_data_files` (C)/`load.data` (R)), then the genetic map takes precedence for setting the genetic markers that will be tracked. All markers that appear in the genetic map file will be tracked, including those that are missing from the genotype matrix (the founders from the genotype matrix will have blank/unknown alleles for those missing markers). Any markers that appear in the genotype matrix but not the genetic map will be ignored.

The line order of the genotype matrix and the genetic map does not need to match. genomicSimulation will match genetic markers using their names.

You can use the `get_group_` family of functions (C)/the `see.group.data` function (R) to confirm data is correctly loaded:

<table>
<caption>Accessing simulation data to check information is correctly loaded</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
For users confident enough to use an iterator, the following method is recommended:
```{C}
// Say we want to save some information about each founder to a logging file
FILE* logfile = fopen("log.txt", "w");

// Set up an iterator that will step through all members of the "founders" group
BidirectionalIterator it = create_bidirectional_iter(d, founders); // <-- Iterator setup
GenoLocation loc = set_bidirectional_iter_to_start(&it);           // <--/
while (IS_VALID_LOCATION(loc)) {                                   // <-- Iterator end condition
	// access some data on this founder. A range of `get_` functions are available
	char* this_founders_name = get_name(loc);
	char* this_founders_alleles = get_alleles(loc);
	// save some of that data, perhaps
	fprintf(logfile, "%i %s\n", get_id(loc).id, this_founders_name);
	
	loc = next_forwards(&it);                                      // <-- Iterator advances to next founder
}
fclose(logfile);
```
Otherwise, the alternate and slightly slower method uses the `get_group_` family of functions, which are more similar to the R version's `see.group.data`:
```{C}
unsigned int num_founders = get_group_size(d, founders);

// Get the names of every founder
char* founder_names[num_founders]; // VLA, can replace with a known maximum length or a heap allocation if your C version does not support VLAs
get_group_names(d, founders, num_founders, founder_names);
// names will now be stored in the founder_names array

// Get the allele sequence (a character vector of length 2*{num markers}) of every founder
char* founder_alleles[num_founders];
get_group_genes(d, founders, num_founders, founder_alleles);
// genotypes will now be stored in the founder_alleles array. Pairs of alleles for the same genetic marker are consecutive. 
```
In the C version, you can also peek at internal data structures, like the list of stored genetic markers. This is not possible in the R version.
```{C}
// (The order of the markers in the genotypes from `get_alleles` & `get_group_genes` is:
char** marker_names_in_order = d->genome.marker_names;
// and there are
unsined int num_markers = d->genome.n_markers;
// of the markers.)
```
<td>
```{R}
# See the stored genotypes, as a (markers x founders) matrix
# The markers may be reordered compared to the input genotype matrix file.
# The markers will be ordered according to their order in the internal genetic map
genomatrix <- see.group.gene.data(founders)

# Other potential checks:
founder.names <- see.group.data(founders, "N")
founder.genotypes.as.strings <- see.group.data(founders, "G") # <-- like concatenating the columns of "genomatrix"
```
</table>

## Load a genetic map and genotype file, specifying file format

By default, genomicSimulation attempts to detect the format of input files as it loads them. When a file loading command is run, genomicSimulation prints out its assumptions about the format of the file.

If the printed logs are incorrect about the layout of your map file or effect file, the most likely cause is a mis-spelling in the header row of the file. If the printed logs are incorrect about the size of the table in your map file or effect file, the issue may be that the numbers cannot be detected as numbers, or that there are issues with the column spacers/line breaks that cause there to be more than three entries to a row.

If the printed logs are incorrect about the layout or format of your genotype matrix file, you can manually specify the file format for genomicSimulation, as follows. 

<table>
<tr><th>Example Input Files <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
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
init <- load.data(allele.file="genotype-file.txt", map.file="map-file.txt",format=mformat)
```
</table>

It is not required that you manually specify all details of the file format. Whatever details of the file format are not manually specified will be automatically detected.

See the documentation for the function `define_matrix_format_details` (C)/`define.matrix.format.details` (R) for more information.

## Load a genetic map and several sets of founder genotypes

<table>
<tr><th>Example Input Files <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
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
init <- load.data("genotype-file.txt", "map-file.txt")
founders_a <- init$groupNum
founders_b <- load.genotypes("genotype-file2.txt")
```
</table>

## Load a genetic map, a set of founder genotypes, and a set of additive trait effects

<table>
<tr><th>Example Input Files <th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
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
init <- load.data(allele.file="genotype-file.txt", map.file="map-file.txt", effect.file="eff-file.txt")
founders <- init$groupNum
eff1 <- init$effectID
map1 <- init$mapID
```
</table>

This provides an introduction to three types of identifiers you will use in genomicSimulation. They are easier to distinguish in the C version, since they are defined as different types (`GroupNum` and `EffectID` and `MapID`). In the R environment, all genomicSimulation identifiers are integer values. An R user has to remember (or use variable names to help them remember) which one is which. 
- `GroupNum`: represents a group of candidates, eg your founder population or your F1 crosses.
- `EffectID`: represents a set of additive trait effects. You can have multiple of these sets of marker effects loaded to represent different traits. Each time you calculate breeding values, you choose which set of marker effects to use. By default it will use the first loaded set of marker effects.
- `MapID`: represents a particular arrangement of your simulation's tracked genetic markers in terms of chromosome positions and recombination frequencies. You can have multiple maps loaded at a time. Each time you generate offspring, you can choose which map to use to generate gametes from the parents. By default it will use the first loaded map.

## Load another set of additive trait effects

This assumes the simulation is set up, that is, one of **Load a genetic map and a set of founder genotypes**, **Load a genetic map and several sets of founder genotypes**, **Load a genetic map, a set of founder genotypes, and a set of additive trait effects** or equivalent has been carried out. If the simulation is not yet set up and you wish to load trait effect values, see the section above, titled **Load a genetic map, a set of founder genotypes, and a set of additive trait effects**.


<table>
<caption>To load an alternative set of marker effects stored in the file "eff-file2.txt":</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
```{C}
EffectID eff2 = load_effectfile(d, "eff-file2.txt");
```
<td>
```{R}
eff2 <- load.effects("eff-file2.txt")
```
</table>

If multiple sets of marker effects are loaded, you can choose which one to use each time you tell genomicSimulation to calculate breeding values, and so can calculate breeding values for multiple traits.

## Load another linkage map

Markers in this map that are not present in the primary linkage map will be ignored. (That is, the first genetic map loaded must contain all markers you wish to be tracked by the simulation. Later genetic maps can only provide different recombination frequencies or chromosome allocations).

<table>
<caption>To load an alternative recombination map stored in the file "mapfile2.txt":</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
```{C}
MapID map2 = load_mapfile(d, "mapfile2.txt");
```
<td>
```{R}
map2 <- load.map("mapfile2.txt")
```
</table>

## Simulating a breeding program

<table>
<caption>To simulate recurrent truncation selection on a trait over 10 generations</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
For the readability of the C simulation, we assume there exists a pair of helpful functions `double mean(int length, double* values)` and `double var(int length, double* values)` that can be called to calculate the mean and variance, respectively, of an array of doubles. 
```{C}
// Set up simulation object
SimData* d = create_empty_simdata(time(NULL));
struct MultiIDSet init = load_data_files(d, "genotype-file.txt", "map-file.txt", "eff-file.txt", DETECT_FILE_FORMAT);
GroupNum parents = init.group;
EffectID effset1 = init.effectSet;
MapID map1 = init.map;

// Set up a file where we will output results
FILE* results = fopen("simresults.txt", "w");
fprintf(results, "Generation,MeanBV,VarBV\n");

double bvs[10]; // assumes even in the founder population, there are only 10 candidates
get_group_bvs(d, parents, effset1, 10, bvs);
fprintf(results, "%i,%d,%d\n", 0, mean(10, bvs), var(10, bvs));

for (int gen = 1; gen <= 10; ++gen) {
	// Generate 25 new genotypes by randomly crossing among last generation's parents
	// (BASIC_OPT is an alias for a default set of progeny generation options. 
	// To use different progeny generation options, create a replacement GenOptions struct.
	// There are examples of this in later sections of this page.)
	GroupNum f1 = make_random_crosses(d, parents, 0, 25, map1, BASIC_OPT); 
	
	// Now, perform truncation selection by taking the 10 candidates with 
    // the best breeding values out of group "f1", and put them in group "f1selection".
	GroupNum f1selection = split_by_bv(d, f1, effset1, 10, 0);
	
	// (Delete last generation's parents and the F1s that were not selected, to free 
    // some memory, since we will not use them again.)
    delete_group(d,f1);
	delete_group(d,parents);

	// Log some information about the current generation. Here, we log some statistics about
    // the breeding values measured in this generation.
	get_group_bvs(d, f1selection, effset1, 10, bvs);
	fprintf(results, "%i,%d,%d\n", gen, mean(10, bvs), var(10, bvs));
	
	// Define the current generation as parents of the next generation for the loop
	parents = f1selection; // Redefining the variable "parents"
}

fclose(results);
```
<td>
```{R}
# Set up simulation
init <- load.data("genotype-file.txt", "map-file.txt", "eff-file.txt")
parents <- init$groupNum
effset1 <- init$effectID
map1 <- init$mapID

# Set up dataframe to save results
founders.bv <- see.group.data(parents, "BV")
results <- data.frame(Generation=0, MeanBV=mean(founders.bv), VarBV=var(founders.bv))

for (gen in 1:10) {
  # Generate 25 new genotypes by randomly crossing among last generation's parents.
  # ("map=map1" is optional, because we have only loaded one map, and it will default to the first map that was loaded)
  f1 <- make.random.crosses(parents, n.crosses=25, map=map1)
  
  # Now, perform truncation selection by taking the 10 candidates with 
  # the best breeding values out of group "f1", and put them in group "f1.selection".
  # ("eff.set=effset1" is optional, because we have only loaded one set of
  # marker effects, and it will default to using the first loaded set of marker effects)
  f1.selection <- break.group.by.GEBV(f1, number=10, eff.set=effset1)
  
  # (Delete last generation's parents and the F1s that were not selected, to free 
  # some memory, since we will not use them again.)
  delete.group(c(parents, f1))

  # Log some information about the current generation. Here, we log some statistics about
  # the breeding values measured in this generation.
  f1.selection.info <- see.group.data(f1.selection, "B")
  results <- rbind(results, c(gen, mean(f1.selection.info), var(f1.selection.info))
  
  # Define the current generation as parents of the next generation for the loop
  parents <- f1.selection
  
  rm(f1, f1.selection, f1.selection.info) # (optional tidying of R environment)
}
```
</table>


The remainder of this document will provide sample scripts for simulating features that may appear in real breeding programs. The scripts for these breeding program components can be substituted or added into the structure above to develop a simulation of a more complex breeding program.

# Crossing & Generating New Genotypes: Plant-themed
Templates in this section assume you have loaded a set of founders whose group number is saved in the variable `founders`.

## Single seed descent and selfing

<table>
<caption>For six generations, to grow a single seed from each plant to maturity:</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
```{C}
GroupNum f6 = self_n_times(d, 6, founders, NO_MAP, BASIC_OPT);
```
<td>
```{R}
f6 <- self.n.times(founders, n=6)
```
</table>

Using NO_MAP for the recombination map parameter in any crossing function in the C version causes the tool to default to the first/primary recombination map. This default behaviour is the same in the R version.

## Creating halfsib or fullsib families

<table>
<caption>Make 20 random crosses and create families of 6 full siblings from each cross:</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr><td>
```{C}
GenOptions opt = {.will_name_offspring=FALSE, 
                  .offspring_name_prefix=NULL, 
				  .family_size=6,     // <---- Non-default value
		          .will_track_pedigree=TRUE, 
				  .will_allocate_ids=TRUE,
		          .filename_prefix=NULL, 
				  .will_save_pedigree_to_file=FALSE,
		          .will_save_bvs_to_file=NO_EFFECTSET, 
				  .will_save_alleles_to_file=FALSE,
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
<caption>Make 10 random crosses with the candidate stored at index 1, and create halfsib families of 6 half siblings from each cross:</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
```{C}
int targetparent = 1;
GroupNum targetparent_group = make_group_from(d, 1, &targetparent);

GenOptions opt = {.will_name_offspring=FALSE, 
                  .offspring_name_prefix=NULL, 
				  .family_size=6,     // <---- Non-default value
		          .will_track_pedigree=TRUE, 
				  .will_allocate_ids=TRUE,
		          .filename_prefix=NULL, 
				  .will_save_pedigree_to_file=FALSE,
		          .will_save_bvs_to_file=NO_EFFECTSET, 
				  .will_save_alleles_to_file=FALSE,
		          .will_save_to_simdata=TRUE};
GroupNum crosses = make_random_crosses_between(d, targetparent_group, founders, 10, 0, 0, NO_MAP, NO_MAP, opt);

GroupNum families[10];
split_into_halfsib_families(d, crosses, 2, families); // <-- 2 for "group by second parent"
```
<td>
```{R}
targetparent_group <- make.group(c(1))

crosses <- make.random.crosses.between(targetparent_group, founders, n.crosses=10, offspring=6)

families <- break.group.into.halfsib.families(crosses, 2) # <-- 2 for "group by second parent"
```
</table>

## Updating marker effects

Updating the marker effect estimates is a task that must be done outside of genomicSimulation. Export required data from the simulation, run the marker effect re-estimation, then re-import marker effects as in preceding section **Swap out the set of additive trait effects for another**.

## Trials and locations

Suppose you want to simulate the success of the founder population in different environments.

<table>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
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
<caption>Pull out a chosen founder and repeatedly (for 20 generations) cross back to it</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
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
targetparent_group <- make.group(c(1))

backcross_generations <- rep(0, times=20) # initialise vector of correct size
backcross_generations[1] <- founders
for (ii in 1:20) {
	backcross_generations[ii] <- make.random.crosses.between(targetparent_group, backcross_generations[ii-1], n.crosses=10)
}
```
</table>

For marker-assisted backcrossing, select each loop on presence of the marker and/or the estimated breeding value. See the sections **Select on a qualitative trait** and **Select on a qualitative trait, then a quantitative trait** further down this page.


# Crossing & Generating New Genotypes: Animal-themed
Templates in this section assume you have loaded two sets of founders whose group numbers are saved in the variables `cows` and `bulls`. See preceding section **Load a genetic map and several sets of founder genotypes**.

## Split offspring into male and female

<table>
<caption>Out of 10 offspring of random crosses, randomly define some as male and some as female</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
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

split_randomly_into_two or break.group.randomly can be used to split a group into two sub-groups by flipping a coin on each group member. This may result in two groups that are not quite the same size. To split a group of eg. 10 offspring into two groups of exactly 5 members, the cousin functions, `split_evenly_into_two` (C)/`break.group.evenly` (R), can be used instead:

<table>
<caption>Out of 10 offspring of random crosses, get exactly 5 male calves and 5 females</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
```{C}
GroupNum offspring = make_random_crosses_between(d, cows, bulls, 10, 0, 0, NO_MAP, NO_MAP, BASIC_OPT);

GroupNum offspring_f = split_randomly_into_two(d, offspring);
GroupNum offspring_m = offspring;
```
<td>
```{R}
offspring <- make.random.crosses.between(cows, bulls, n.crosses=10)

temporary <- break.group.evenly(offspring, into.n = 2)
offspring_f <- temporary[1]
offspring_m <- temporary[2]
rm(temporary)
```
</table>

Or, for uneven proportions, use the functions 
- `split_by_probabilities` (C)/`break.group.with.probabilities` (R), corresponding to `split_randomly_into_n` (C)/`break.group.randomly` (R), to randomly allocate sex by performing a random draw/coin toss for each candidate to determine which group it is placed into, or 
- `split_into_buckets` (C)/`break.group.into.buckets` (R), corresponding to `split_evenly_into_n` (C)/`break.group.evenly` (R), to randomly sort the list of candidates and allocate specific numbers of them to each group.

<table>
<caption>Out of 10 offspring of random crosses, get exactly 1 male calf and 9 females</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
```{C}
GroupNum offspring = make_random_crosses_between(d, cows, bulls, 10, 0, 0, NO_MAP, NO_MAP, BASIC_OPT);

// split_into_buckets needs to be passed n-1 "bucket sizes" to make n groups
// so we tell it we want 9 females and the remainder (n=10-9=1) will be allocated male 
int bucketsize = 9; 
GroupNum temporary[2];
split_into_buckets(d, offspring, 2, &bucketsize, &temporary);
GroupNum offspring_f = temporary[0];
GroupNum offspring_m = temporary[1];
```
<td>
```{R}
offspring <- make.random.crosses.between(cows, bulls, n.crosses=10)

# break.group.into.buckets takes one less "bucket size" than the number of groups it will make
# so we tell it we want 9 females and the remainder (n=10-9=1) will be allocated male.
temporary <- break.group.into.buckets(offspring, buckets=c(9))
offspring_f <- temporary[1]
offspring_m <- temporary[2]
rm(temporary)
```
</table>

As an alternative, you could create a custom label to track the animal's sex (eg. the label having value 1 for a female, 2 for a male). Sex could be initially set by passing a vector of 1s and 2s (random or otherwise) to `change_label_to_values` (C)/`change.label.to.values` (R). Afterwards, `split_by_label_value` (C)/`break.group.by.label.value` (R) can be used to separate the males and females as needed, and as needed `combine_groups` (C)/`combine.groups` (R) can merge the populations back together. 

However, for most simulations of breeding programs (it is likely every generation would require separating out the males and females before crossing them), it would be more computationally efficient to maintain two separate "male" and "female" groups, rather than using a custom label and repeatedly merging and unmerging the groups. 

## Introducing fresh blood to male and female kernels
Suppose the `offspring_f` and `offspring_m` groups exist, as created in the previous section **Split offspring into male and female**.

<table>
<caption>Introduce new calves to the breeding kernels</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
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

It may be useful to use age tracking in these kernels. See the following section, **Tracking age and culling individuals that are too old**.

## Tracking age and culling individuals that are too old
genomicSimulation's custom labels can be used to track age (or some other known trait that can be represented with an integer for each individual).

Using custom labels involves three main operations:
- Creating the label and giving it a default value with `create_new_label` (C)/`create.new.label` (R)
- Modifying the values of certain or all candidate genotypes, eg. to increase the age of the animals. Relevant functions include `change_label_by_amount` (C)/`change.label.by.amount` (R), `change.label.to` (C)/`change.label.to.this` (R).
- Splitting a group based on the values of a custom label, so that you can carry out some crossing or selection operations on only a subset of the population (eg. those of a certain age). Relevant functions include `split_by_label_range` (C)/`break.group.by.label.range` (R), `split_by_label_value` (C)/`break.group.by.label.value` (R).

<table>
<caption>A population where only animals 3 years or older can breed</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
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
	// Ages 3-20 can breed.
	GroupNum breedingGroup = split_by_label_range(d, animals, ageLabel, 3, 20); 
	
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
	// Ages 3-20 can breed.
	breedingGroup <- break.group.by.label.range(ageLabel, 3, 20, animals);
	
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
<caption>Cull animals that are 12 years or older</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
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

## Random mortality

In simulating a breeding program, there may often be stages where some of the population is lost in a way that can be best simulated as a random event. 

<table>
<caption>Some proportion of the animals will be lost this generation</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
```{C}
double mortality = 0.06; // 6% chance of mortality

GroupNum temporary[2];
split_by_probabilities(d, animals, 1, &mortality, &temporary); // <-- 6% chance to be in first group, (100-6)% chance to survive.
GroupNum lost_animals = temporary[0];
GroupNum animals_survived = temporary[1];

delete_group(d, lost_animals);
```
<td>
```{R}
mortality <- 0.06  # 6% chance of mortality

temporary <- break.group.with.probabilities(animals, mortality)  # <-- 6% chance to be in first group, (100-6)% chance to survive.
lost_animals <- temporary[1]
animals_survived <- temporary[2] 
rm(temporary)

delete.group(lost_animals)   # <-- to free the memory being used to store the candidates from "lost_animals". The group "lost_animals" can no longer be used after this line
rm(lost_animals)             # <-- optional: to tidy up the R environment by removing the stale variable "lost_animals" so it does not confuse you later
```
</table>

Because all offspring are generated independently in genomicSimulation, if you have a random culling/offspring loss step immediately after the command in which those genotypes were generated, you could get speedups in simulation by simply only generating the surviving offspring:

<table>
<caption>Some proportion of offspring are immediately lost</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
```{C}
// Produce offspring
GroupNum offspring = make_random_crosses(d, breedingGroup, 50000, 0, NO_MAP, BASIC_OPT);

// Immediate mortality
double mortality = 0.06; // 6% chance of mortality
GroupNum temporary[2];
split_by_probabilities(d, offspring, 1, &mortality, &temporary);
GroupNum lost_offspring = temporary[0];
GroupNum offspring_survived = temporary[1];
delete_group(d, lost_offspring);
```
or, to combine those steps and potentially speed up simulation:
```{C}
// 50000 offspring at 6% mortality rate. Calculate number of survivors:
int n_survivors = 0;
for (int i = 0; i < 50000; ++i) {
	// this line, run 50000 times, adds one survivor to the count if a random draw from a uniform distribution 0-1 is greater than 0.06
	n_survivors += ((double)rand() / (double)RAND_MAX) > 0.06;
	// If you have a proper C random number generation library included, use that instead.
	// (if you have a proper random number generation library, you can probably use its random binomial distribution 
	// function, as show in R version opposite, instead of this loop)
}

GroupNum offspring = make_random_crosses(d, breedingGroup, n_survivors, 0, NO_MAP, BASIC_OPT);
```
<td>
```{R}
// Produce offspring
offspring <- make.random.crosses(breedingGroup, 50000)

// Immediate mortality
mortality <- 0.06  # 6% chance of mortality
temporary <- break.group.with.probabilities(offspring, mortality)
lost_offspring <- temporary[1]
offspring_survived <- temporary[2] 
rm(temporary)

delete.group(lost_offspring)
rm(lost_offspring)
```
or, to combine those steps and potentially speed up simulation:
```
// 50000 offspring at 6% mortality rate. Calculate number of survivors:
n_survivors <- rbinom(1,50000,1-0.06)
offspring <- make.random.crosses(breedingGroup, n_survivors)
```
</table>

## Culling or fitness-based animal loss

Alternatively, rather than a random mortality system, you could simulate loss of candidates according to their estimated breeding values. Perhaps those with the lowest breeding values for some fitness trait are the ones that are lost, or perhaps the ones with the highest for an economic value index are sold. 

<table>
<caption>The worst animals will be lost this generation</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
```{C}
double mortality = 0.06; // 6% mortality proportion
int num_will_lose = round(mortality * get_group_size(d, animals));

GroupNum lost_animals = split_by_bv(d, animals, fitness_effid, num_will_lose, 1);
delete_group(d, lost_animals);
```
If the situation is more complex than just deleting the animals with the lowest breeding values, the following equivalent template will be more useful:
```{C}
double mortality = 0.06; // 6% mortality proportion
unsigned int num_animals = get_group_size(d, animals);
unsigned int num_will_lose = round(mortality * num_animals);
unsigned int index_of_lost_animals[num_will_lose];

// You would then access some group data and choose the indexes of the animals that
// will be lost however you would like. Here I pretend I have created a function 
// "void find_worst_animals(uint num_animals, uint* animal_indexes, uint* animal_bvs, uint find_how_many, uint* save_worst_animal_indexes_here);"
// and I need to pass animal indexes and breeding values to that function.
unsigned int index_animals[num_animals];
get_group_indexes(d, animals, num_animals, index_animals);
double bv_animals[num_animals];
get_group_bvs(d, animals, fitness_effid, num_animals, bv_animals);
find_worst_animals(num_animals, index_animals, bv_animals, num_will_lose, index_of_lost_animals);

GroupNum lost_animals = make_group(d, num_will_lose, index_of_lost_animals);
delete_group(d, lost_animals);
```
<td>
```{R}
mortality <- 0.06

lost_animals <- break.group.by.gebv(animals, low.score.best=TRUE, percentage=mortality*100, eff.set=fitness.effset)
delete.group(lost_animals)
rm(lost_animals)
```
If the situation is more complex than just deleting the animals with the lowest breeding values, the following equivalent template will be more useful:
```
mortality <- 0.06

info <- data.frame(Index=see.group.data(animals,"X"),
                   Fitness=see.group.data(animals,"BV",eff.set=fitness.effset))
	 
# Use the "info" table to choose which will be lost. We show a simple example here
# of just taking the animals with the lowest breeding values, but more complex
# choices could be made.
num_will_lose <- round(mortality * length(info$Index))
index_of_lost_animals <- info[order(info$GEBV, decreasing=FALSE),]$Index[1:num_will_lose]

lost_animals <- make.group(index_of_lost_animals)
delete.group(lost_animals)
rm(lost_animals)
```
</table>

The important thing to remember when using a system like this is that any `delete_group` (C)/`delete.group` (R) call can potentially invalidate any position indexes calculated before it. There is a potential pitfall in a situation when you calculate scores for multiple groups, or multiple subgroups of a single group, and can become tempted to delete the unwanted groups as you go (see table below. The examples shown are for the R version of genomicSimulation, but the same pitfall exists in the C version).

<table>
<caption>Pitfalls when deleting groups</caption>
<tr><th>INCORRECT SEQUENCE<th>CORRECT SEQUENCE
<tr>
<td>
In this scenario, we are performing the same culling strategy as in the previous table. Some proportion of animals with the lowest breeding values will be lost/culled. However, in this case, we have different loss rates by sex, and a custom label that tracks the sex of each animal.
```{R}
female.mortality <- 0.06
male.mortality <- 0.4

info <- data.frame(Index=see.group.data(animals,"X"), # <-- position indexes accessed here
                   Fitness=see.group.data(animals,"BV",eff.set=fitness.effset),
				   Sex=see.group.data(animals,"Label",label=label_sex))        
info.females <- info[Sex==1,]
info.males <- info[Sex==2,]

num_females_will_lose <- round(female.mortality * length(info.females$Index))
index_of_lost_females <- info.females[order(info.females$GEBV, decreasing=FALSE),]$Index[1:num_females_will_lose]
lost_females <- make.group(index_of_lost_females)     # <-- make.group uses position indexes
delete.group(lost_females)                            # <-- all previously calculated position indexes are now potentially invalid
rm(lost_females)

num_males_will_lose <- round(male.mortality * length(info.males$Index))
index_of_lost_males <- info.males[order(info.males$GEBV, decreasing=FALSE),]$Index[1:num_males_will_lose]
													 # <-- info tables contain invalid position indexes now
lost_males <- make.group(index_of_lost_males)        # <-- !! make.group uses invalid position indexes, so the animals selected to be lost_males may not correspond to the ones we intended to select. 
delete.group(lost_males)
rm(lost_males)
```
<td>
To fix this, you could simply rearrange the script so that you only run `delete.group` once all interactions with the `info` tables are done:
```{R}
female.mortality <- 0.06
male.mortality <- 0.4

info <- data.frame(Index=see.group.data(animals,"X"), # <-- position indexes accessed here
                   Fitness=see.group.data(animals,"BV",eff.set=fitness.effset),
				   Sex=see.group.data(animals,"Label",label=label_sex))        
info.females <- info[Sex==1,]
info.males <- info[Sex==2,]

num_females_will_lose <- round(female.mortality * length(info.females$Index))
index_of_lost_females <- info.females[order(info.females$GEBV, decreasing=FALSE),]$Index[1:num_females_will_lose]
lost_females <- make.group(index_of_lost_females)     # <-- make.group uses position indexes

num_males_will_lose <- round(male.mortality * length(info.males$Index))
index_of_lost_males <- info.males[order(info.males$GEBV, decreasing=FALSE),]$Index[1:num_males_will_lose]
lost_males <- make.group(index_of_lost_males)         # <-- make.group uses position indexes that are still valid

# As a helpful check, you may want to remove the tables containing position indexes
# from the R environment so that there is no temptation to access them after groups are
# deleted and the position indexes stored in them are no longer valid
rm(info, info.females, info.males)

delete.group(c(lost_females, lost_males))             # <-- position indexes from info tables now invalid (but the tables are already deleted)
rm(lost_females,lost_males)
```
Be thoughtful about where in your script it is safe to place `delete.group` commands.
</table>


## Random mating with caps on number of offspring per animal

<table>
<caption>Any cow should have at most one calf this generation</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
```{C}
GroupNum offspring = make_random_crosses_between(d, cows, bulls, 10, 1, 0, NO_MAP, NO_MAP, BASIC_OPT);
```
<td>
```{R}
offspring <- make.random.crosses.between(cows, bulls, cap1=1, n.crosses=10)
```
</table>

## Mating all females to a good male

<table>
<caption>A breeding strategy where only the best bull in the population will father calves</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
```{C}
// Assuming there exists a variable "eff1" which is an EffectID for the set of marker effects to use to calculate breeding values here
GroupNum bestbull = split_by_bv(d, bulls, eff1, 1, FALSE);

GroupNum offspring = make_random_crosses_between(d, cows, bestbull, 10, 1, 0, NO_MAP, NO_MAP, BASIC_OPT);

// optional: put him back with the other bulls
GroupNum toCombine[2] = {bulls, bestbull};
bulls = combine_groups(d, 2, toCombine);
```
<td>
```{R}
# by default, this function uses the first effect set loaded, so if "eff1" is 
# the first effect set loaded, `eff.set=eff1` is optional:
bestbull <- break.group.by.gebv(bulls, number=1, eff.set=eff1) 
offspring <- make.random.crosses.between(cows, bestbull, cap1=1, n.crosses=10)

# optional: put him back with the other bulls
bulls <- combine.groups(c(bulls, bestbull))
```
</table>

## Random mating between the best performing individuals in a group

genomicSimulation used to have a function `make_n_crosses_from_top_m_percent`, which would identify a highest-performing slice of a group and make random crosses between them, but leave the high performers as members of their original group. Here are code snippets that replicate its functionality:

<table>
<caption>A function to make n random crosses between the m percent of candidates in a group that have the highest breeding values</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
```{C}
GroupNum make_n_crosses_from_top_m_percent(SimData* d, GroupNum* group, const int n, const float m, EffectID effset) {
	# Split off the best candidates 
	int mcount = (int)(get_group_size(d,*group) * m / 100); # close enough, you could round differently if you preferred
	GroupNum best_candidates = split_by_bv(d, *group, effset, mcount, 0);
	
	# Make random crosses to get offspring, using default map
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
make.n.crosses.from.top.m.percent <- function(group, n, m, effset=1) {
	# Split off the best candidates 
	best.candidates <- break.group.by.gebv(group, percentage=m, eff.set=effset)
	
	# Make random crosses to get offspring 
	offspring <- make.random.crosses(best.candidates, n.crosses=n)
	
	# Put the best candidates back into their original group
	# Note: you cannot assume that the group number representing "group"
	# will not change as a result of this line, thus the double arrow
	# to make sure any changes to the group number representing "group"
	# are also saved in this function's calling environment. 
	# (alternatively you could return "group" alongside the new "offspring" group)
	group <<- combine.groups(c(group, best.candidates))
	
	return(offspring)
}
```
</table>

## Making specific chosen matings

<table>
<caption>Suppose you want to cross Cow1 to Bull1, Cow2 to Bull2, and Daisy to Bull1</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
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
# The R version of the function (make.targeted.crosses) can be passed the indexes of these parent animals, like in the C version, 
# or it can use the names of the parent animals.
offspring <- make.targeted.crosses(c("Cow1", "Cow2", "Daisy"), c("Bull1", "Bull2", "Bull1"))
```
</table>

The other option is to create an input file of the following format:
```
Cow1 Bull1
Cow2 Bull2
Daisy Bull1
```
and call `make_crosses_from_file` (C)/ `make.crosses.from.file` (R).

## Three-Way Crosses

Suppose you want to make a single Breed1/Breed2//Breed3 cross. That is, cross the F1 of a mating between Breed1 and Breed2 to Breed3. Assumes the parent you've chosen to represent Breed1 is named "Breed1Parent", and likewise for Breed2 and Breed3:

<table>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
```{C}
breed1_index = gsc_get_index_of_name(d->m, "Breed1Parent");
breed2_index = gsc_get_index_of_name(d->m, "Breed2Parent");
breed3_index = gsc_get_index_of_name(d->m, "Breed3Parent");

int crossingPlan[2][1];
crossingPlan[0][0] = breed1_index;
crossingPlan[1][0] = breed2_index;

GroupNum f1 = make_targeted_crosses(d, 1, crossingPlan[0], crossingPlan[1], NO_MAP, NO_MAP, BASIC_OPT);
int f1_index[1];
get_group_indexes(d, f1, 1, f1_index); //we know this group has only one member

// Re-use crossingPlan
crossingPlan[0][0] = breed3_index;
crossingPlan[1][0] = f1_index[0];

GroupNum f3way = make_targeted_crosses(d, 1, crossingPlan[0], crossingPlan[1], NO_MAP, NO_MAP, BASIC_OPT);
```
<td>
```{R}
f1 <- make.targeted.crosses(c("Breed1Parent"), c("Breed2Parent"))
f1_index <- see.group.data(f1, "X")
f3way <- make.targeted.crosses(c("Breed3Parent"), f1_index)
```
</table>

And if you wanted more than one offspring from this three way cross:

<table>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
```{C}
breed1_index = gsc_get_index_of_name(d->m, "Breed1Parent");
breed2_index = gsc_get_index_of_name(d->m, "Breed2Parent");
breed3_index = gsc_get_index_of_name(d->m, "Breed3Parent");

int crossingPlan[2][1];
crossingPlan[0][0] = breed1_index;
crossingPlan[1][0] = breed2_index;

GenOptions opt = {.family_size=5,  // <-- changed from default
		.will_name_offspring=FALSE, 
		.offspring_name_prefix=NULL,
		.will_track_pedigree=TRUE, 
		.will_allocate_ids=TRUE,
		.filename_prefix=NULL, 
		.will_save_pedigree_to_file=FALSE,
		.will_save_bvs_to_file=NO_EFFECTSET, 
		.will_save_alleles_to_file=FALSE,
		.will_save_to_simdata=TRUE};

GroupNum f1 = make_random_crosses(d, 1, crossingPlan[0], crossingPlan[1], NO_MAP, NO_MAP, opt);
int f1_indexes[5];
get_group_indexes(d, f1, 5, f1_indexes);

int crossingPlanb[2][5];
crossingPlanb[0][0] = breed3_index; crossingPlanb[0][1] = breed3_index;
crossingPlanb[0][2] = breed3_index; crossingPlanb[0][3] = breed3_index;
crossingPlanb[0][4] = breed3_index;
crossingPlanb[1][0] = f1_indexes[0]; crossingPlanb[1][1] = f1_indexes[1];
crossingPlanb[1][2] = f1_indexes[2]; crossingPlanb[1][3] = f1_indexes[3];
crossingPlanb[1][4] = f1_indexes[4];

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



## Species with differing male and female recombination rates

<table>
<caption>Simulate crossing in a species with differing male and female recombination maps</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
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
Templates in this section assume you have a group `target` from which you want to select a subset of 'good' genotypes. Scripts are provided for selecting in different ways, showing some ways to extract information about the genotypes in the group, and showing how to select using this information using the function `make_group` (C)/`make.group` (R).

The C implementations in this section will be very naive/inefficient. The focus is on showing what kind of custom selection criteria can be implemented in genomicSimulation, not how best to sort/search through data efficiently in C. Sometimes, only R examples will be shown, because R's built-in table-manipulation functionalities make these tasks much simpler.

## Manually select on true breeding value

This is a selection function template. You can use the builtin functions `split_by_bv` (C)/ `break.group.by.gebv` (R) to carry out this exact task. It provides however a simple example of genomicSimulations "custom selection method interface", which involves accessing group data, choosing certain candidates based on that data, then calling `make_group` (C)/`make.group` (R) to move the selected candidates into a new group. This is the overall structure that the following sections will also use.

<table>
<caption>Select the top 10 candidates in `target` according to their true breeding value as calculated by the internal breeding value calculator</caption>
<tr><th>genomicSimulationC (C) <th>genomicSimulation (R)
<tr>
<td>
```{C}
GroupNum select_10_with_best_GEBV(GroupNum group, EffectID eff1) {
  unsigned int group_size = get_group_size( d, group );
  
  int group_indexes[group_size]; // VLA usage: swap for a heap array or a maximum group size if you don't like that
  memset(group_indexes, 0, sizeof(int)*group_size);
  get_group_indexes( d, group, group_size, group_indexes);
  
  double group_bvs[group_size];
  memset(group_bvs, 0, sizeof(double)*group_size);
  get_group_bvs( d, group, eff1, group_size, group_bvs);

  # A pretty inefficient way of finding the top 10 scores and their corresponding
  # indexes, but it will do for this example.
  double min_score = group_bvs[0];
  for (unsigned int i = 1; i < group_size; ++i) {
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
                     GEBV=see.group.data(group,"BV",eff.set=eff1))

  selected_group <- make.group(info[order(info$GEBV, decreasing=TRUE),]$Index[1:10])
  return( selected_group )
}
```
</table>


## Select on phenotype (simulated with a given heritability)

<table>
<caption>Select the top 10 candidates in `target` according to a phenotype simulated with a certain $H^2$ broad-sense heritability value.</caption>
<tr><th>genomicSimulation (R)
<tr>
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

## Select on an index weighting two traits.

<table>
<caption>Select the top 10 candidates in `target` according to an index weighting two traits</caption>
<tr><th>genomicSimulation (R)
<tr>
<td>
```{R}
select.top.10.by.selectionindex <- function(group, trait1.effset, trait2.effset, weighting=0.7) {
  info <- data.frame(PositionIndex=see.group.data(group,"X"),
                     BV1=see.group.data(group,"BV", eff.set=trait1.effset),
					 BV2=see.group.data(group,"BV", eff.set=trait2.effset))

  # you could alter this formula as suits to calculate the selection index				 
  info$SelectionIndex <- info$BV1 * weighting + info$BV2 * (1-weighting)
  
  selected <- make.group(info[order(info$SelectionIndex, decreasing=TRUE),]$PositionIndex[1:10])
  return( selected )
}
```
</table>

## Select on a qualitative trait

If you have a trait with only a small number of possible outcomes - perhaps a marker for which each individual can be homozygous, heterozygous, or carry no copies - and you want to select for individuals that carry that trait.

One option is to write a very small marker effect file, that looks something like:
```
TARGETGENE T 1
```

You can load this effect file as in the earlier section **Load another set of additive trait effects**, and then use it with the function `split_by_bv` (C)/ `break.group.by.gebv` (R). Homozygous carriers will have a breeding value calculated as 2, heterozygous carriers will score 1, and non-carriers will score 0. 

<table>
<caption>Use a custom marker effect file to select on a qualitative trait</caption>
<tr><th>genomicSimulation (R)
<tr>
<td>
```{R}
qualtrait <- load.effects("one_line_effectfile.txt")

select.top.10.qual.genotypes <- function(group) {
	return(break.group.by.GEBV(group, number=10, eff.set=qualtrait))
}
# or perhaps
select.homozygous.carriers <- function(group) {
	info <- data.frame(Index=see.group.data(group,"X"),
                       GEBV=see.group.data(group,"BV",eff.set=qualtrait))
	
	return(make.group(info$Index[info$GEBV == 2]))
}
```
</table>

Alternatively, you could read the genotype matrix manually: 

<table>
<caption>Access the genotype matrix to select on a qualitative trait</caption>
<tr><th>genomicSimulation (R)
<tr>
<td>
```{R}
select.homozygous.carriers <- function(group, allele='T') {
	counts <- see.group.gene.data(group, count.allele=allele)
	selected_names <- colnames(counts)[counts["TARGETGENE",]==2]
	rm(counts)
	
	info <- data.frame(Index=see.group.data(group,"X"),
                       Names=see.group.data(group,"N"))
	return(make.group(info$Index[info$Names %in% selected_names]))
}

# With this method, you can also access phase information.
# Say you wanted to select carriers of that "T" allele
# who inherited the allele from their mother/first parent
select.carriers.inheriting.from.parent1 <- function(group) {
	genotypes <- see.group.gene.data(group)
	selected_names <- colnames(genotypes)[genotypes["TARGETGENE",] %in% c("TA", "TT")]
	rm(genotypes)
	
	info <- data.frame(Index=see.group.data(group,"X"),
                       Names=see.group.data(group,"N"))
	return(make.group(info$Index[info$Names %in% selected_names]))
}
```
</table>

## Select on a qualitative trait, then a quantitative trait.

Suppose you have an effect file for a quantitative trait, called `eff-large.txt`, as well as a one-line effect file for a qualitative trait, as described in the above section **Select on a qualitative trait**, named `eff-small.txt`.

You would like to select primarily on the allele listed in `eff-small.txt`, and secondarily select on the quantitative trait. This means you would like your selection to be wholly homozygous for the allele from `eff-small.txt`, if possible, and if you have enough homozygous candidates in the population, the selection should be the homozygous candidates with the highest scores for the quantitative trait. If your population does not have enough homozygous individuals to fill up your selection, then you will take all the homozygous individuals and the best of the heterozygous individuals. And only if you do not have enough homozygotes and heterozygotes will you add the best of the non-carriers to your selection.

<table>
<caption>Use two marker effect files, select first for carriers of a specific allele, then on a quantitative trait</caption>
<tr><th>genomicSimulation (R)
<tr>
<td>
```{R}
quanttrait <- load.effects("eff-large.txt")
qualtrait <- load.effects("eff-small.txt") # has one line and will return scores of 0, 1, or 2
# Instead of the "eff-small.txt" file, you could also use a see.group.gene.data alternative method to calculate info$Copies, as discussed in the previous section.

select.n.best.carriers <- function(population, n) {
	info <- data.frame(Index=see.group.data(group,"X"),
                       BV=see.group.data(group,"BV",eff.set=quanttrait),
					   Copies=see.group.data(group,"BV",eff.set=qualtrait))

	take.all.or.best.n <- function(infosubset, n) {
		if (n <= 0) { return(c()) }
		if (n >= dim(infosubset)[1]) {
			return(infosubset$Index)
		} else { # rank and return best n
			return(info[order(-info$BV),]$Index[1:n]
		}
	}
	
	will_be_selected <- take.all.or.best.n(info[info$Copies == 2,], n)
	will_be_selected <- take.all.or.best.n(info[info$Copies == 1,], n-length(will_be_selected))
	will_be_selected <- take.all.or.best.n(info[info$Copies == 0,], n-length(will_be_selected))

	return(make.group(will_be_selected))
}
```
</table>



(this list of examples will continue to be expanded...)

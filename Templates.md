Features and Templates Guide        {#templates}
============================

While [the documentation pages of genomicSimulationC](@ref modules) and its counterpart [genomicSimulation](https://github.com/vllrs/genomicSimulation) do describe all the features the simulation tool offers, they don't do so in a way that is necessarily easy to browse or that's helpful for new users to understand the simulation tool. 

This page provides a set of templates for implementing common breeding project actions in genomicSimulationC and its R counterpart. 

[TOC]

# Setting Up

Consider using the get_group_ family of functions (genomicSimulationC) and see.group.data function (genomicSimulation) to confirm data is correctly loaded.

## Load a genetic map and a set of founder genotypes

## Load a genetic map and several sets of founder genotypes

## Load a genetic map, a set of founder genotypes, and a set of additive trait effects

## Swap out the set of additive trait effects for another

## Input file formats: 

### Loading from other file formats
Updates to expand the range of allowed input formats are coming soon.


# Crossing & Other Ways to Generate New Genotypes: Plant-themed
Templates in this section assume you have loaded a set of founders whose group id is saved in the variable `founders`. 

## Single seed descent

In genomicSimulationC (C):
```

```

In genomicSimulation (R):
```

```

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



 











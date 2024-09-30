# genomicSimulation, C version

This project is a set of C functions that run stochastic simulations of genetics in breeding schemes. [An overview of the tool (and its R package counterpart) is available at this link.](https://doi.org/10.1093/g3journal/jkac216).

All functions are located in `sim-operations.c` and its header file `sim-operations.h`. A Makefile for the `gcc` compiler [is provided](https://github.com/vllrs/genomicSimulationC/blob/main/Makefile).

[Latest news updates can be found on the page linked here.](https://vllrs.github.io/genomicSimulationC/html/news.html)

### Guide
Documentation (generated with doxygen) is available at [the GitHub pages site](https://vllrs.github.io/genomicSimulationC/html/index.html) of [the repository](https://github.com/vllrs/genomicSimulationC). Also available there are [guides to setting up and simulating particular crossing designs](https://vllrs.github.io/genomicSimulationC/html/templates.html), and [lists of equivalent functions in the R and C versions of the package](https://vllrs.github.io/genomicSimulationC/html/concordance.html). The file `sim.c` in this repository provides a simplistic simulation template.

### R package version
[An R package allowing these functions to be called in R is also provided under the name genomicSimulation.](https://github.com/vllrs/genomicSimulation)

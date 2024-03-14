# genomicSimulation, C version

This project is a set of C functions that run stochastic simulations of genetics in breeding schemes. [An overview of the tool (and its R package counterpart) is available at this link.](https://doi.org/10.1093/g3journal/jkac216).

All functions are located in `sim-operations.c` and its header file `sim-operations.h`. A Makefile for the `gcc` compiler [is provided](https://github.com/vllrs/genomicSimulationC/blob/main/Makefile).

[Latest news updates can be found on the page linked here.](https://vllrs.github.io/genomicSimulationC/html/news.html)

### Guide
Documentation (generated with doxygen) is available at [the GitHub pages site](https://vllrs.github.io/genomicSimulationC/html/index.html) of [the repository](https://github.com/vllrs/genomicSimulationC). Also available there is an explanation of simulation methods and assumptions, a template for setting up a simulation (in the file `sim.c`), and guides to simulating particular crossing designs using the package's functions.

### R package version
[An R package allowing these functions to be called in R is also provided under the name genomicSimulation.](https://github.com/vllrs/genomicSimulation)

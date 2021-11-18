CFLAGS = -Wall -std=c99

.PHONY: clean all
.DEFAULT_GOAL := all

all: sim sim-operations.o
debug: gsim gsim-operations.o sim-test
test: sim-test

osim: sim.c sim-operations.c sim-operations.h
	gcc -std=c99 -lm sim.c sim-operations.c -O3 -o sim

sim: sim.o sim-operations.o
	gcc -lm sim.o sim-operations.o -o sim
	
gsim: gsim.o gsim-operations.o
	gcc -lm sim.o sim-operations.o -g -o sim
	
sim.o: sim.c
	gcc $(CFLAGS) -c sim.c
	
gsim.o: sim.c
	gcc $(CFLAGS) -g -c sim.c

sim-operations.o: sim-operations.c sim-operations.h
	gcc $(CFLAGS) -c sim-operations.c

gsim-operations.o: sim-operations.c sim-operations.h
	gcc $(CFLAGS) -g -c sim-operations.c
	
sim-test: sim-test.o gsim-operations.o
	gcc -lm sim-test.o sim-operations.o -o sim-test
	
sim-test.o: sim-test.c sim-test.h
	gcc $(CFLAGS) -g -c sim-test.c

clean:
	rm sim-operations.o sim.o sim sim-test sim-test.o
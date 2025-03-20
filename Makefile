CFLAGS = -std=c99 -Wvla -Wall
RFLAGS = -O2 -Wall -mfpmath=sse -msse2 -mstackrealign
SHARELIBFLAGS = -shared -fPIC

LIBS = -lm
DEPS = sim-operations.h
TESTDEPS = $(DEPS) sim-test.h
CSOURCE = $(DEPS:.h=.c)
TESTCSOURCE = $(TESTDEPS:.h=.c)

PATH_REG99_OBJ = $(TESTDEPS:.h=-99.o)
PATH_DEBUG99_OBJ = $(TESTDEPS:.h=-g.o)
PATH_RSTYLE_OBJ = $(TESTDEPS:.h=-rstyle.o)
PATH_SHAREDLIB_OBJ = $(DEPS:.h=-so.o)

PATH_SHAREDLIB = $(DEPS:.h=.so)

all: sim-test-g
	
# Different compilations of the test suite
$(PATH_REG99_OBJ): $(TESTCSOURCE) $(TESTDEPS)
	gcc -c -o $@ $(@:-99.o=.c) $(CFLAGS)

$(PATH_DEBUG99_OBJ): $(TESTCSOURCE) $(TESTDEPS)
	gcc -g -c -o $@ $(@:-g.o=.c) $(CFLAGS)

# RSTYLE compiles and runs with the flags that R uses when compiling sources for an R library
$(PATH_RSTYLE_OBJ): $(TESTCSOURCE) $(TESTDEPS) 
	gcc -c -o $@ $(@:-rstyle.o=.c) $(RFLAGS)
	
sim-test-99: $(PATH_REG99_OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

sim-test-g: $(PATH_DEBUG99_OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

sim-test-rstyle: $(PATH_RSTYLE_OBJ)
	gcc -o $@ $^ $(RFLAGS) $(LIBS)

# Compilation to a shared library
sharedlib: lib
lib: $(PATH_SHAREDLIB_OBJ)
	gcc -o $(PATH_SHAREDLIB) $^ $(LIBS) $(SHARELIBFLAGS)

$(PATH_SHAREDLIB_OBJ): $(CSOURCE) $(DEPS)
	gcc -c -o $@ $(@:-so.o=.c) $(CFLAGS) $(SHARELIBFLAGS)

# clean 
clean:
	rm $(PATH_REG99_OBJ) $(PATH_DEBUG99_OBJ) $(PATH_RSTYLE_OBJ) $(PATH_SHAREDLIB_OBJ)
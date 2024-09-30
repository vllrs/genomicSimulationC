CFLAGS = -std=c99 -Wvla -Wall
RFLAGS = -O2 -Wall -mfpmath=sse -msse2 -mstackrealign
LIBS = -lm
DEPS = sim-operations.h sim-test.h

CSOURCE = $(DEPS:.h=.c)
PATH_REG99_OBJ = $(DEPS:.h=-99.o)
PATH_DEBUG99_OBJ = $(DEPS:.h=-g.o)
PATH_RSTYLE_OBJ = $(DEPS:.h=-rstyle.o)

all: sim-test-g
	
$(PATH_REG99_OBJ): $(CSOURCE) $(DEPS)
	gcc -c -o $@ $(@:-99.o=.c) $(CFLAGS)

$(PATH_DEBUG99_OBJ): $(CSOURCE) $(DEPS)
	gcc -g -c -o $@ $(@:-g.o=.c) $(CFLAGS)
	
$(PATH_RSTYLE_OBJ): $(CSOURCE) $(DEPS)
	gcc -c -o $@ $(@:-rstyle.o=.c) $(RFLAGS)
	
sim-test-99: $(PATH_REG99_OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)
	
sim-test-g: $(PATH_DEBUG99_OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)
	
sim-test-rstyle: $(PATH_RSTYLE_OBJ)
	gcc -o $@ $^ $(RFLAGS) $(LIBS)
	
clean:
	rm $(PATH_REG99_OBJ) $(PATH_DEBUG99_OBJ) $(PATH_RSTYLE_OBJ)
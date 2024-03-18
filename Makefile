CFLAGS = -Wall -std=c99 -Wvla
LIBS = -lm
DEPS = sim-operations.h sim-test.h
OBJ = sim-operations.o sim-test.o

all: sim-test

%.o: %.c $(DEPS)
	gcc -g -c -o $@ $< $(CFLAGS)
	
sim-test: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)
	
clean:
	rm $(OBJ) sim-test.exe
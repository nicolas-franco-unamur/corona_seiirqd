#makefile for corona_seiirqd program
CC=gcc
CFLAGS=-Wall -Wno-unused-variable
LDFLAGS=-lm -lgsl -lgslcblas

#uncomment this if need of additional memory
LDFLAGS+= -Wl,-stack_size,0xc0000000,-stack_addr,0xc0000000 -O3 -mfpmath=sse -msse2


EXEC=corona

all: $(EXEC)

corona: corona_seiirad.o struct_fct.o mcmc.o simulation.o scenario.o output.o
	$(CC) -o $@ $^ $(LDFLAGS)

corona_seiirad.o: corona_seiirqd.c corona.h command.h
	$(CC) -o $@ -c $< $(CFLAGS)

mcmc.o: mcmc.c corona.h command.h
	$(CC) -o $@ -c $< $(CFLAGS)

simulation.o: simulation.c corona.h command.h
	$(CC) -o $@ -c $< $(CFLAGS)

struct_fct.o: struct_fct.c corona.h command.h
	$(CC) -o $@ -c $< $(CFLAGS)

scenario.o: scenario.c corona.h command.h
	$(CC) -o $@ -c $< $(CFLAGS)

output.o: output.c corona.h command.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
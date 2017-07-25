# Makefile for LibBinary examples

# C Compiler (HDF5 MPI compiler must be used)
CC = h5pcc

# Directory for Compiled objects
ODIR = ./obj

# LibBinary details

# Libbinary (requires static linking)
LIBBINARY = -lbinary
# Path to libbinary
LDIR = ../
# Path to libbinary header
IDIR = ../include

# Compiler flags
CFLAGS = -I$(IDIR) -L$(LDIR) -Wall -O3 -std=c11

# Dynamic libraries (require dynamics linking)
DYLIBS = -lfftw3_mpi -lfftw3 -lm -lgsl -lgslcblas

all: benchmark restart input

benchmark: $(ODIR)/benchmark.o
	$(CC) $< -Bstatic $(LIBBINARY) $(CFLAGS) -Bdynamic $(DYLIBS) -o $@ 

restart: $(ODIR)/restart.o
	$(CC) $< -Bstatic $(LIBBINARY) $(CFLAGS) -Bdynamic $(DYLIBS) -o $@

input: $(ODIR)/input.o
	$(CC) $< -Bstatic $(LIBBINARY) $(CFLAGS) -Bdynamic $(DYLIBS) -o $@

$(ODIR)/%.o: ./src/%.c
	$(CC) -c $^ $(CFLAGS) -o $@

$(ODIR):
	mkdir -p obj

.PHONY: clean

clean:
	rm -rf obj
	rm -f benchmark
	rm -f restart
	rm -f newsim
	rm -f input

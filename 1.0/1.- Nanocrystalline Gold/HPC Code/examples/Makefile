# C Compiler
CC = h5pcc
IDIR = ../include
ODIR = ./obj
LDIR = ../lib

CFLAGS = -I$(IDIR) -L$(LDIR) -Wall -O3 -std=c11
DYLIBS = -lfftw3_mpi -lfftw3 -lm -lgsl -lgslcblas
STLIBS = -lbinary

all: benchmark restart input

benchmark: $(ODIR)/benchmark.o
	$(CC) $< -Bstatic $(STLIBS) $(CFLAGS) -Bdynamic $(DYLIBS) -o $@ 

restart: $(ODIR)/restart.o
	$(CC) $< -Bstatic $(STLIBS) $(CFLAGS) -Bdynamic $(DYLIBS) -o $@

input: $(ODIR)/input.o
	$(CC) $< -Bstatic $(STLIBS) $(CFLAGS) -Bdynamic $(DYLIBS) -o $@

$(ODIR)/%.o: ./src/%.c
	$(CC) -c $^ $(CFLAGS) -o $@

.PHONY: clean

clean:
	rm -f obj/*
	rm -f benchmark
	rm -f restart
	rm -f newsim
	rm -f input

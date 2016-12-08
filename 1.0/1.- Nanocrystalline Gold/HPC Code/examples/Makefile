# C Compiler
CC = h5pcc
IDIR = ../include
ODIR = ./obj
LDIR = ../lib

CFLAGS = -I$(IDIR) -L$(LDIR) -Wall -O3 -ffast-math -std=c11
DYLIBS = -lfftw3_mpi -lfftw3 -lm -lgsl -lgslcblas
STLIBS = -lbinary

benchmark: $(ODIR)/benchmark.o
	$(CC) $< $(CFLAGS) $(DYLIBS) -Bstatic $(STLIBS) -o $@

newsim: $(ODIR)/newsim.o
	$(CC) $< $(CFLAGS) $(DYLIBS) -Bstatic $(STLIBS) -o $@

restart: $(ODIR)/restart.o
	$(CC) $< $(CFLAGS) $(DYLIBS) -Bstatic $(STLIBS) -o $@

$(ODIR)/%.o: ./src/%.c
	$(CC) -c $^ $(CFLAGS) -o $@

.PHONY: clean

clean:
	rm -f obj/*
	rm -f benchmark
	rm -f restart
	rm -f newsim

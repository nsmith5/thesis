# Make file for LibBinary

# Must use the HDF5 C compiler
CC = h5pcc

# Targets
IDIR = ./include
ODIR = ./obj
SDIR = ./src

CFLAGS = -I$(IDIR) -Wall -O3 -std=c11
LIBS = -lfftw3_mpi -lfftw3 -lm -lgsl -lgslcblas

# Headers
_DEPS = binary.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

# Objects we'd like to compile
_OBJ = io.o state.o dynamics.o error.o setup.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# Use the archiver tool to make a static library
# (bundle all the objects together)
libbinary.a: $(OBJ)
	ar rc $@ $^

# Compiling C objects and dumping them in the /obj
# directory
$(ODIR)/%.o: $(SDIR)/%.c $(DEPS) $(ODIR)
	$(CC) -c $< $(CFLAGS) -o $@

# Make the object directory *before* you compile 
# things and put them there
$(ODIR):
	mkdir obj

# Fake target (clean everything up)
.PHONY: clean

clean:
	rm -rf obj
	rm -f libbinary.a

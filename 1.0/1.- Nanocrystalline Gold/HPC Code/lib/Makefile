# C Compiler
CC = h5pcc
IDIR = ../include
ODIR = ./obj
SDIR = ../src

CFLAGS = -I$(IDIR) -Wall -O3 -ffast-math -std=c11
LIBS = -lfftw3_mpi -lfftw3 -lm -lgsl -lgslcblas

# Headers
_DEPS = state.h io.h setup.h error.h dynamics.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

# Objects we'd like to compile
_OBJ = io.o state.o dynamics.o error.o setup.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

libbinary.so: $(OBJ)
	$(CC) -shared $^ $(CFLAGS) $(LIBS) -o $@

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -fPIC -c $< $(CFLAGS) -o $@


.PHONY: clean

clean:
	rm -f obj/*
	rm -f libbinary.so

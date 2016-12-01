CC = h5pcc
CFLAGS = -I./include -Wall -g -std=c11
LIBS = -lfftw3_mpi -lfftw3 -lm -lgsl -lgslcblas

main: obj/state.o obj/main.o obj/io.o obj/dynamics.o obj/error.o obj/setup.o
	$(CC) obj/* $(CFLAGS) $(LIBS) -o main

obj/main.o: src/main.c
	$(CC) -c src/main.c $(CFLAGS) -o obj/main.o

obj/state.o: src/state.c include/state.h
	$(CC) -c src/state.c $(CFLAGS) -o obj/state.o

obj/io.o: src/io.c include/io.h
	$(CC) -c src/io.c $(CFLAGS) -o obj/io.o

obj/dynamics.o: src/dynamics.c include/dynamics.h
	$(CC) -c src/dynamics.c $(CFLAGS) -o obj/dynamics.o

obj/error.o: src/error.c include/error.h
	$(CC) -c src/error.c $(CFLAGS) -o obj/error.o

obj/setup.o: src/setup.c include/setup.h
	$(CC) -c src/setup.c $(CFLAGS) -o obj/setup.o

.PHONY: clean

clean:
	rm -f data/Data.h5
	rm -f obj/*
	rm -f main

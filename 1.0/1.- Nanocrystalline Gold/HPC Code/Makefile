CC = h5pcc
CFLAGS = -I./include -Wall -O3 -std=c99
LIBS = -lfftw3_mpi -lfftw3 -lm

main: obj/state.o obj/main.o obj/io.o obj/dynamics.o obj/error.o obj/random.o
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

obj/random.o: src/random.c include/random.h
	$(CC) -c src/random.c $(CFLAGS) -o obj/random.o

.PHONY: clean

clean:
	rm -f data/Data.h5
	rm -f obj/*
	rm -f main

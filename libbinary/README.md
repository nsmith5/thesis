# LibBinary

LibBinary is a high performance C library for simulating the binary XPFC model.
It uses MPI for large scale distributed parallelism. 

# Building and Dependencies

LibBinary is an MPI library and as such it requires an MPI Compiler. Both
OpenMPI and MPICH have been tested to work, but once you chose one you must be
sure that you compile all relevant dependencies with the same thing. Speaking of
dependencies, LibBinary requires:

 - `gsl` the GNU Scientific Library
 - `HDF5` a high performance parallel I/O library for saving simulation results
   (compiled with your MPI compiler)
 - `FFTW3` the 'fastest fourier transform in the west` for distributed parallel
   transforms (compiled with your MPI compiler)

Additionally, building the library required GNUMake to run the relevant build
scripts. In the simpliest case the library can be built with a simple

```bash
$ make
mkdir obj
h5pcc -c src/io.c -I./include -Wall -O3 -std=c11 -o obj/io.o
h5pcc -c src/state.c -I./include -Wall -O3 -std=c11 -o obj/state.o
h5pcc -c src/dynamics.c -I./include -Wall -O3 -std=c11 -o obj/dynamics.o
h5pcc -c src/error.c -I./include -Wall -O3 -std=c11 -o obj/error.o
h5pcc -c src/setup.c -I./include -Wall -O3 -std=c11 -o obj/setup.o
ar rc libbinary.a obj/io.o obj/state.o obj/dynamics.o obj/error.o obj/setup.o
```

To build examples, head to the `/examples` directory and use `make` again.

```bash
$ cd examples
$ make
mkdir -p obj
h5pcc -c src/benchmark.c -I../include -L../ -Wall -O3 -std=c11 -o obj/benchmark.o
h5pcc obj/benchmark.o -Bstatic -lbinary -I../include -L../ -Wall -O3 -std=c11 -Bdynamic -lfftw3_mpi -lfftw3 -lm -lgsl -lgslcblas -o benchmark
h5pcc -c src/restart.c -I../include -L../ -Wall -O3 -std=c11 -o obj/restart.o
h5pcc obj/restart.o -Bstatic -lbinary -I../include -L../ -Wall -O3 -std=c11 -Bdynamic -lfftw3_mpi -lfftw3 -lm -lgsl -lgslcblas -o restart
h5pcc -c src/input.c -I../include -L../ -Wall -O3 -std=c11 -o obj/input.o
h5pcc obj/input.o -Bstatic -lbinary -I../include -L../ -Wall -O3 -std=c11 -Bdynamic -lfftw3_mpi -lfftw3 -lm -lgsl -lgslcblas -o input
rm obj/restart.o obj/input.o obj/benchmark.o
```

# Usage

# Questions
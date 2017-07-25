# LibBinary

LibBinary is a high performance C library for simulating the binary XPFC model. It uses
MPI for large scale distributed parallelism. 

# Building and Dependencies

LibBinary is an MPI library and as such it requires an MPI Compiler. Both OpenMPI and
MPICH have been tested to work, but once you chose one you must be sure that you compile
all relevant dependencies with the same thing. Speaking of dependencies, LibBinary requires:

 - `gsl` the GNU Scientific Library
 - `HDF5` a high performance parallel I/O library for saving simulation results (compiled with your MPI compiler)
 - `FFTW3` the 'fastest fourier transform in the west` for distributed parallel transforms (compiled with your MPI compiler)

Additionally, building the library required GNUMake to run the relevant build scripts. In the
simpliest case the library can be built with a simple

```bash
$ make
```

# Usage

# Questions
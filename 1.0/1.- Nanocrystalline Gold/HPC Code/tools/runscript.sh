#!/bin/bash
#PBS -l nodes=2:ppn=12
#PBS -l walltime=4:00:00
#PBS -A fyw-423-ad
#PBS -o output
#PBS -e error
#PBS -V
#PBS -N Gold_Run

module load HDF5
module load FFTW
module load GSL

cd "$PBS_O_WORKDIR" 

cd lib
make clean
make
cd ../examples
make clean
make 
mv restart ..
cd ..

timeout -s SIGUSR1 3.98h mpiexec -np 24 ./restart output.h5 00057000

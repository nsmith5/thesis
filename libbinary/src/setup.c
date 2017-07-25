/*  Setup state and teardown state */

#include <mpi.h>
#include <hdf5.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>     // Catch signals for killing the simulation gracefully
#include <stdbool.h>    // Booleans...just because they help reading

#include "binary.h"

#define WISDOM_FILE "fft_wisdom"    // Save fftw plans to this file

double start_time = 0.0;

void init (int    argc,
           char **argv)
{
    int rank, err;

    /* Initialize MPI Runtime */
    MPI_Init(&argc, &argv);

    /* Initialize the clock */
    start_time = MPI_Wtime();

    /* Initialize fftw */
    fftw_mpi_init ();

    /* Import wisdom from file and broadcast */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* If the wisdom file exists, try to load it */
    if (rank == 0 && access (WISDOM_FILE, F_OK) != -1)
    {
        err = fftw_import_wisdom_from_filename (WISDOM_FILE);
        if (err == 0) my_error("Importing FFTW wisdom failed!");
    }

    /* Broadcast wisdom if any exists */
    fftw_mpi_broadcast_wisdom (MPI_COMM_WORLD);

    /* Wait for everyone */
    MPI_Barrier (MPI_COMM_WORLD);
    return;
}

void finalize (void)
{
    int rank, err;

    /* Gather wisdom from procs and save to file */
    fftw_mpi_gather_wisdom (MPI_COMM_WORLD);

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    if (rank == 0) // Only rank=0 proc will save the wisdom
    {
        err = fftw_export_wisdom_to_filename (WISDOM_FILE);
        if (err == 0)
        {
            remove (WISDOM_FILE);
            my_error ("Failed to correctly export FFTW wisdom");
        }
    }

    /* Wait for everyone */
    MPI_Barrier (MPI_COMM_WORLD);

    /* Clean up threads, fftw and finalize MPI runtime */
    fftw_mpi_cleanup ();
    MPI_Finalize ();

    return;
}

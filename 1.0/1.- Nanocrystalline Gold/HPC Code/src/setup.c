#include <mpi.h>
#include <hdf5.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#include <stdbool.h>

#include "setup.h"
#include "state.h"
#include "error.h"
#include "io.h"

bool time_to_leave = false;

static void sig_handler(int signo)
{
	/* If we receive SIGUSR1 its time to get the hell out */
    if (signo == SIGUSR1)
    {
		MPI_Barrier (MPI_COMM_WORLD);
    	mpi_print("Caught SIGUSR1, time to leave\n");
		time_to_leave = true;
	}
}

void init (int    argc,
           char **argv)
{
    int rank, err;

    /* Initialize MPI Runtime */
    MPI_Init(&argc, &argv);

    /* Register signal handler for SIGUSR1 */
    signal(SIGUSR1, sig_handler);

    /* Initialize fftw */
    fftw_mpi_init ();

    /* Import wisdom from file and broadcast */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0 && access ("data/plans.wisdom", F_OK) != -1)
    {
        err = fftw_import_wisdom_from_filename ("data/plans.wisdom");
        if (err == 0) my_error("Importing FFTW wisdom failed!");
    }

    fftw_mpi_broadcast_wisdom (MPI_COMM_WORLD);

    MPI_Barrier (MPI_COMM_WORLD);
    return;
}

void finalize (void)
{
    int rank, err;

    /* Gather wisdom from procs and save to file */
    fftw_mpi_gather_wisdom (MPI_COMM_WORLD);

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        err = fftw_export_wisdom_to_filename ("data/plans.wisdom");
        if (err == 0)
        {
            remove ("data/plans.wisdom");
            my_error ("Failed to correctly export FFTW wisdom");
        }
    }

    MPI_Barrier (MPI_COMM_WORLD);

    // Clean up threads, fftw and finalize MPI runtime
    fftw_mpi_cleanup ();
    MPI_Finalize ();

    return;
}

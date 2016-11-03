#include <mpi.h>
#include <hdf5.h>
#include <fftw3-mpi.h>
#include <stdio.h>
#include <unistd.h>		// access () for checking file existance

#include "error.h"
#include "state.h"
#include "dynamics.h"
#include "io.h"

#define FILENAME "data/Data.h5"

void init (int argc, char **argv);
void finalize (void);

int main (int argc, char **argv)
{
    int N = 4096;
    double dx = 0.1;
    double dt = 0.1;
    double D = 1.0;
    state *s;
    hid_t file_id;

    /*
     * Initialize MPI Runtime and FFTW
     */
    init (argc, argv);
    s = create_state (N, dx, dt, D);
    file_id = io_init (FILENAME);
    /*
     * Make a square initial condition
     */
    make_square (s, 1.0);

    /*
     * - Time Step the state
     * - Save each time step
     * - Measure the time for the whole loop
     */

	MPI_Barrier (MPI_COMM_WORLD);
    double t1 = MPI_Wtime();
    for (int i = 0; i < 10; i++)
    {
      step (s);
      //save_state (s, file_id);
    }

	MPI_Barrier (MPI_COMM_WORLD);
    double t2 = MPI_Wtime();
    /*
     * Print out the results of the time trial
     */
	int rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	if (rank == 0)
    	printf("Elapsed time for loop is %f\n", t2-t1);

    /*
     * Do clean up of the system before exiting
     */
    destroy_state (s);
    io_finalize (file_id);
    MPI_Barrier (MPI_COMM_WORLD);
    finalize ();
    return 0;
}

void init (int    argc,
           char **argv)
{
    MPI_Init(&argc, &argv);

    // Initialize fftw
    fftw_mpi_init ();

    // Import wisdom from file and broadcast
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0 && access ("data/plans.wisdom", F_OK) != -1)
    {
        int err = fftw_import_wisdom_from_filename ("data/plans.wisdom");
        if (err == 0) my_error("Importing FFTW wisdom failed!");
    }
    fftw_mpi_broadcast_wisdom (MPI_COMM_WORLD);

    MPI_Barrier (MPI_COMM_WORLD);
  return;
}

void finalize (void)
{
    // Gather wisdom from procs and save to file
    fftw_mpi_gather_wisdom (MPI_COMM_WORLD);
    int rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        int err = fftw_export_wisdom_to_filename ("data/plans.wisdom");
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

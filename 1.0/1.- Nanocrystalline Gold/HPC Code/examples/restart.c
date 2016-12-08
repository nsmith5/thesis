/*
 *      Restart Simulation
 *
 * Usage:
 *  $ timeout -s SIGUSR1 <walltime> mpiexec -np <procs> newsim <datafile> <groupname>
 *
 * Restart a simulation from hdf5 datafile <datafile> from state in group <groupname>
 * with <procs> processes and run for <walltime>
 */

/* Library headers */
#include <mpi.h>
#include <hdf5.h>
#include <math.h>
#include <assert.h>

/* Local Headers */
#include "error.h"
#include "state.h"
#include "dynamics.h"
#include "io.h"
#include "setup.h"

#define PI 2.0*acos(0.0)

int main (int    argc,
          char **argv)
{
    // TODO: fix these magic numbers, this parameters should be loaded from file
    int N = 256;
    double dx = 0.125;
    double dt = 0.00125;
    int rank, size;
	hid_t file_id;
	state* s;

    /* Initialize the system */
    init (argc, argv);
	s = create_state (N, dx, dt);
    assert (s != NULL);

	mpi_print("Initializing I/O..\n");
	file_id = io_init_from_file (argv[1]);

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

	MPI_Barrier (MPI_COMM_WORLD);
    for (int i = 0; i < size; i++)
    {
        if (rank == i)
        {
            printf("Hi!, Im process %d and I've got %d row of data\n", rank, (int)s->local_n0);
        }
    }
	MPI_Barrier (MPI_COMM_WORLD);

    load_state (s, file_id, argv[2]);

    mpi_print("\nTime Stepping...\n");
    /* Do stuff */
    while (!time_to_leave)
	{
        step (s);
		if (s->step % 1000 == 0)
		{
			mpi_print("\tSaving state");
			save_state (s, file_id);
		}
    }
	MPI_Barrier (MPI_COMM_WORLD);

	/* Shut 'er down */
	mpi_print("Shut er down!");
    io_finalize (file_id);
    destroy_state (s);
    finalize ();

    return 0;
}

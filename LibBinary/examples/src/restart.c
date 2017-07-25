/*
 *      Restart Simulation
 *
 * Usage:
 *  $ mpiexec -np <procs> restart <datafile> <groupname> <walltime>
 *
 * Restart a simulation from hdf5 datafile <datafile> from state in group <groupname>
 * with <procs> processes and run for <walltime> minutes
 */

/* Library headers */
#include <mpi.h>
#include <hdf5.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

/* Local Headers */
#include "binary.h"

#define PI 2.0*acos(0.0)

int main (int    argc,
          char **argv)
{
    int runtime;
    double now;
	hid_t file_id;
	state* s;

    /* Initialize the system */
    init (argc, argv);
    file_id = io_init_from_file (argv[1]);
	s = load_state (file_id, argv[2]);
    runtime = atoi (argv[3]);
    assert (s != NULL);
    now = MPI_Wtime();

    mpi_print("\nTime Stepping...\n");
    /* Do stuff */
    while ((now - start_time)/60.0 < (double)runtime)
	{
        step (s);
		if (s->step % 1000 == 0)
		{
			mpi_print("\tSaving state");
			save_state (s, file_id);
            now = MPI_Wtime();
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

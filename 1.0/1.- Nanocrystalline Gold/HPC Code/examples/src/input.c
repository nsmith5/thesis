/*
 *      Simulation from Input file
 *
 * Usage:
 *  >> mpiexec -np <procs> ./input <inputfile> <outputfile> <walltime>
 *
 * Will start a new simulation that will run for <walltime> minutes using <procs>
 * processes. Parameters will be loaded from an input file <inputfile> and output
 * will be saved to <outputfile>
 */

/* Library headers */
#include <mpi.h>
#include <hdf5.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

/* Local Header */
#include "binary.h"

#define PI 2.0*acos(0.0)

int main (int    argc,
          char **argv)
{
    int runtime;
	hid_t file_id;
	state* s;
    double now;

    /* Initialize the system */
    init (argc, argv);
	s = new_state_from_file (argv[1]);
    file_id = io_init_new_file (argv[2]);
    runtime = atoi (argv[3]);
    assert (s != NULL);

    mpi_print("\nTime Stepping...\n");
    
    now = MPI_Wtime();
    
    /* Do stuff */
    while ((now-start_time)/60.0 < (double)runtime)
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

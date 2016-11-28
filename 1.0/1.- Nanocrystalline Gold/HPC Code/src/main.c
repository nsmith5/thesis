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
#define FILENAME "data/Data.h5"

int main (int    argc,
          char **argv)
{
    int N = 1024;
    double dx = 0.125;
    double dt = 0.00125;
    int rank, size;
	hid_t file_id;
	state* s;

    /* Initialize the system */
    init (argc, argv);
	s = create_state (N, dx, dt);
    assert (s != NULL);
	mpi_print("Loading file\n");
	sleep(1);
	file_id = io_init_from_file (FILENAME);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);
   	mpi_print("Done loading file");
	sleep(1); 

	mpi_print("Loading state from file");
    load_state (s, file_id, "00191000");

	mpi_print("Done loading state from file\n");

    mpi_print("Time Stepping");
    /* Do stuff */
   	while (!time_to_leave)
    {
        step (s);
		if (s->step % 1000 == 0)
		{
			mpi_print("Saving state");
			save_state (s, file_id);
		}
    }

	/* Shut 'er down */
	mpi_print("Shut er down!");
    io_finalize (file_id);
    destroy_state (s);
    finalize ();

    return 0;
}

/*
 *      New Simulation
 *
 * $ timeout -s SIGUSR1 <walltime> mpiexec -np <procs> newsim <outputfile>
 *
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

	mpi_print("Initializing I/O..\n");
	file_id = io_init_new_file (FILENAME);

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

	MPI_Barrier (MPI_COMM_WORLD);
    for (int i = 0; i < size; i++)
    {
        if (rank == i)
        {
            printf("Hi!, Im proc rank %d and my local_n0 is %d\n", rank, (int)s->local_n0);
        }
    }
	MPI_Barrier (MPI_COMM_WORLD);

    //load_state (s, file_id, "00422000");

    s->eta = 2.0;
    s->chi = 1.0;
    s->epsilon0 = 30.0;
    s->sigma0 = 0.15;
    s->sigma = 0.07;
    s->omega = 0.30;
    s->Wc = 1.0;

    s->Mn = 1.0;
    s->Mc = 1.0;

    s->k0 = 2*PI;
    s->alpha = 0.8;
    s->beta = 6.0;
    s->rho = 0.8660254037844386;
    s->alphac = 0.5;

    set_C(s);

    for (int i = 0; i < s->local_n0; i++)
    {
        for (int j = 0; j < N; j++)
        {
            int ij = i*2*((N>>1) + 1) + j;
            s->c[ij] = 0.3;
            s->n[ij] = 0.05;
        }
    }

    MPI_Barrier (MPI_COMM_WORLD);
    double t1 = MPI_Wtime ();

    for (int i = 0; i < 100; i++)
	{
        step (s);
    }

	MPI_Barrier (MPI_COMM_WORLD);
    double t2 = MPI_Wtime ();

    printf("Time for 100 time steps = %gs\n", t2-t1);

	/* Shut 'er down */
	mpi_print("Shut er down!");
    io_finalize (file_id);
    destroy_state (s);
    finalize ();

    return 0;
}

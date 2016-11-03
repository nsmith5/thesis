#include <mpi.h>
#include <hdf5.h>
#include <fftw3-mpi.h>
#include <stdio.h>
#include <unistd.h>		// access () file existance
#include <math.h>
#include <assert.h>

#include "random.h"
#include "error.h"
#include "state.h"
#include "dynamics.h"
#include "io.h"

#define PI 2.0*acos(0.0)
#define FILENAME "data/Data.h5"

void init (int argc, char **argv);
void finalize (void);

int main (int argc, char **argv)
{
    int N = 512;
    double dx = 0.125;
    double dt = 0.00125;
    state *s;
    hid_t file_id;

    /* Initialize the system */
    init (argc, argv);
    s = create_state (N, dx, dt);
    file_id = io_init (FILENAME);
    assert (s != NULL);

    /* Set Free Energies Parameters */
    s->eta      = 2.0;
    s->chi      = 1.0;
    s->epsilon0 = 30.0;
    s->sigma0   = 0.15;
    s->sigma    = 0.07;
    s->omega    = 0.30;
    s->Wc       = 1.0;
    s->kbT      = 0.0028;

    /* Dynamic Parameters */
    s->Mn = 1.0;
    s->Mc = 1.0;

    /* Set Correlation Function */
    s->k0 = 2*PI;
    s->alpha = 0.8;
    s->beta = 6.0;
    s->rho = sqrt(3)/2.0;
    s->alphac = 0.5;
    set_C (s);

    for (int i = 0; i < s->local_n0; i++)
        for (int j = 0; j < s->N; j++)
        {
            s->n[i*2*(s->N/2 + 1) + j] = 0.0;
            s->c[i*2*(s->N/2 + 1) + j] = 0.5;
        }

    /*
     * - Time Step the state
     * - Save each time step
     * - Measure the time for the whole loop
     */
    save_state (s, file_id);

	MPI_Barrier (MPI_COMM_WORLD);
     double t1 = MPI_Wtime();
    //for (int i = 0; i < 10; i++)
    //{
    //  step (s);
    //  save_state (s, file_id);
    //}
    step(s);
    step(s);


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

    io_finalize (file_id);
    destroy_state (s);
    printf ("Got here!");
    MPI_Barrier (MPI_COMM_WORLD);
    finalize ();
    return 0;
}

void init (int    argc,
           char **argv)
{
    MPI_Init(&argc, &argv);
    seed (20);

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

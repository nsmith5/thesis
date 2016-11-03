im#include <stdlib.h>
#include <math.h>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include "state.h"

#define PI 2*acos(0)

double calc_k(int i, int j, int N, double dx)
{
    /*  Compute the wavenumber at [i,j] for N x N system    */

    double L = N*dx;

    double kx2 = i < (N>>1) + 1 ? pow (2*PI*i/L, 2) : pow (2*PI*(N-i)/L, 2);
    double ky2 = j < (N>>1) + 1 ? pow (2*PI*j/L, 2) : pow (2*PI*(N-j)/L, 2);

    return sqrt(kx2 + ky2);
}

state* create_state (int N, double dx, double dt)
{
    /* Make a state pointer */
    state *s = malloc (sizeof (state));
    if (s == NULL) return NULL;

    /* Get numerical parameters sorted */
    s->dx = dx;
    s->dt = dt;
    s->t = 0;
    s->step = 0;

    ptrdiff_t local_alloc = fftw_mpi_local_size_2d (N, N/2 + 1,
                                                    MPI_COMM_WORLD,
                                                    $s->local_n0,
                                                    $s->local_0_start);

    /* Allocate soo much memory */
    s->c = fftw_alloc_real (2 * local_alloc);
    s->n = fftw_alloc_real (2 * local_alloc);
    s->fnnl = fftw_complex (local_alloc);
    s->fcnl = fftw_complex (local_alloc);
    s->fc = fttw_complex (local_alloc);
    s->fn = fftw_complex (local_alloc);
    s->fCn = fftw_complex (local_alloc);
    s->fxin = fftw_complex (local_alloc);
    s->fxic = fftw_complex (local_alloc);
    s->Cn = fftw_alloc_real (2 * local_alloc);
    s->nnl = fftw_alloc_real (2 * local_alloc);
    s->cnl = fftw_alloc_real (2 * local_alloc);
    s->k2 = fftw_alloc_real (local_alloc);
    s->C = fftw_alloc_real (local_alloc);

    /* Check that our mountain of memory was obtained */
    if ( s->c       ||
         s->n       ||
         s->fnnl    ||
         s->fncl    ||
         s->fc      ||
         s->fn      ||
         s->fCn     ||
         s->fxin    ||
         s->fxic    ||
         s->Cn      ||
         s->nnl     ||
         s->cnl     ||
         s->k2      ||
         s->C       == NULL)
    {
        free (s->c);
        free (s->n);
        free (s->fnnl);
        free (s->fncl);
        free (s->fc);
        free (s->fn);
        free (s->fCn);
        free (s->fxin);
        free (s->fxic);
        free (s->Cn);
        free (s->nnl);
        free (s->cnl);
        free (s->k2);
        free (s->C);
        return NULL;
    }

    /* Make FFT Plans */
    s->fft_plan = fftw_mpi_plan_dft_r2c_2d (N, N, s->T, s->fT, MPI_COMM_WORLD,
                                            FFTW_MEASURE);
    s->ifft_plan = fftw_mpi_plan_dft_c2r_2d (N, N, s->fT, s->T, MPI_COMM_WORLD,
                                             FFTW_MEASURE);

    /* Make k2 operator */
    set_k2 (s);

    MPI_Barrier(MPI_COMM_WORLD);
  	return s;
}

void destroy_state (state* s)
{
    // Free Memory of state
    if (s != NULL)
    {
        fftw_destroy_plan (s->fft_plan);
        fftw_destroy_plan (s->ifft_plan);
        free (s->c);
        free (s->n);
        free (s->fnnl);
        free (s->fncl);
        free (s->fc);
        free (s->fn);
        free (s->fCn);
        free (s->fxin);
        free (s->fxic);
        free (s->Cn);
        free (s->nnl);
        free (s->cnl);
        free (s->k2);
        free (s->C);
  	    free (s);
    }
}

void set_C (state* S)
{
    /*
     * Call this function if the temperature has changed and the correlation
     * function must be reset
     */
     double k;
     for (int i = 0; i < s->local_n0; i++)
        for (int j = 0; j < N/2 + 1; j++)
        {
                k = calc_k(i, j, s->N, s->dx);
                s->C[i*(N/2+1) + j] = exp(-s->sigma * s->sigma * s->k0 * s->k0 /
                    (2.0 * s->rho * s->beta)) * exp(-pow(k - s->k0, 2.0) /
                    (2.0 * pow(s->alphac, 2)));
        }
}

void set_k2 (state* S)
{
    double k;

    for (int i = 0; i < s->local_n0; i++)
    {
        for (int j = 0; j < N/2 + 1; j++)
        {
            k = calc_k (i, j, s->N, s->dx);
            s->k2[i(N/2+1) + j] = k*k;
        }
    }

    return;
}

#include <stdlib.h>
#include <math.h>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h>

#include "state.h"

#define RNG_SEED 123
#define PI 2*acos(0)

static double ipow (double x, int p)
{
    /* Unsafe x^p computation */
    if (p == 1) return x;
    else return x*ipow(x, p-1);
}

double calc_k (int i, int j, int N, double dx)
{
    /*  Compute the wavenumber at [i,j] for N x N system    */

    double L = N*dx;

    double kx2 = i < (N>>1) + 1 ? ipow (2*PI*i/L, 2) : ipow (2*PI*(N-i)/L, 2);
    double ky2 = j < (N>>1) + 1 ? ipow (2*PI*j/L, 2) : ipow (2*PI*(N-j)/L, 2);

    return sqrt(kx2 + ky2);
}

void set_k2 (state* s)
{
    double k;
    int ij;

    for (int i = 0; i < s->local_n1; i++)
    {
        for (int j = 0; j < s->N; j++)
        {
            k = calc_k (i + s->local_1_start, j, s->N, s->dx);
            ij = i*s->N + j;
            s->k2[ij] = k*k;
        }
    }

    return;
}

state* create_state (int N, double dx, double dt)
{
    /* Make a state pointer */
    int rank;
    ptrdiff_t local_alloc;

    state *s = malloc (sizeof (state));
    if (s == NULL) return NULL;

    /* Get numerical parameters sorted */
    s->N = N;
    s->dx = dx;
    s->dt = dt;
    s->t = 0.0;
    s->step = 0;

    /* Allocate RNG and seed */
    s->rng = gsl_rng_alloc (gsl_rng_default);
    if (s->rng == NULL)
    {
        gsl_rng_free (s->rng);
        return NULL;
    }
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    gsl_rng_set (s->rng, RNG_SEED * (rank + 1));

    local_alloc =
        fftw_mpi_local_size_2d_transposed (N,
                                           N/2 + 1,
                                           MPI_COMM_WORLD,
                                           &s->local_n0,
                                           &s->local_0_start,
                                           &s->local_n1,
                                           &s->local_1_start);

    /* Allocate soo much memory */
    s->c    = fftw_alloc_real (2 * local_alloc);
    s->n    = fftw_alloc_real (2 * local_alloc);
    s->fnnl = fftw_alloc_complex (local_alloc);
    s->fcnl = fftw_alloc_complex (local_alloc);
    s->fc   = fftw_alloc_complex (local_alloc);
    s->fn   = fftw_alloc_complex (local_alloc);
    s->fCn  = fftw_alloc_complex (local_alloc);
    s->fxin = fftw_alloc_complex (local_alloc);
    s->fxic = fftw_alloc_complex (local_alloc);
    s->Cn   = fftw_alloc_real (2 * local_alloc);
    s->nnl  = fftw_alloc_real (2 * local_alloc);
    s->cnl  = fftw_alloc_real (2 * local_alloc);
    s->k2   = fftw_alloc_real (local_alloc);
    s->C    = fftw_alloc_real (local_alloc);

    /* Check that our mountain of memory was obtained */
    if ( s->c     == NULL  ||
         s->n     == NULL  ||
         s->fnnl  == NULL  ||
         s->fcnl  == NULL  ||
         s->fc    == NULL  ||
         s->fn    == NULL  ||
         s->fCn   == NULL  ||
         s->fxin  == NULL  ||
         s->fxic  == NULL  ||
         s->Cn    == NULL  ||
         s->nnl   == NULL  ||
         s->cnl   == NULL  ||
         s->k2    == NULL  ||
         s->C     == NULL)
    {
        fftw_free (s->c);
        fftw_free (s->n);
        fftw_free (s->fnnl);
        fftw_free (s->fcnl);
        fftw_free (s->fc);
        fftw_free (s->fn);
        fftw_free (s->fCn);
        fftw_free (s->fxin);
        fftw_free (s->fxic);
        fftw_free (s->Cn);
        fftw_free (s->nnl);
        fftw_free (s->cnl);
        fftw_free (s->k2);
        fftw_free (s->C);
        return NULL;
    }

    /* Make FFT Plans */
    s->fft_plan =
        fftw_mpi_plan_dft_r2c_2d (N,
                                  N,
                                  s->c,
                                  s->fc,
                                  MPI_COMM_WORLD,
                                  FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);

    s->ifft_plan =
        fftw_mpi_plan_dft_c2r_2d (N,
                                  N,
                                  s->fc,
                                  s->c,
                                  MPI_COMM_WORLD,
                                  FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN);

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
        fftw_free (s->c);
        fftw_free (s->n);
        fftw_free (s->fnnl);
        fftw_free (s->fcnl);
        fftw_free (s->fc);
        fftw_free (s->fn);
        fftw_free (s->fCn);
        fftw_free (s->fxin);
        fftw_free (s->fxic);
        fftw_free (s->Cn);
        fftw_free (s->nnl);
        fftw_free (s->cnl);
        fftw_free (s->k2);
        fftw_free (s->C);
        gsl_rng_free (s->rng);
  	    free (s);
    }
}

void set_C (state* s)
{
    /*
     * Call this function if the temperature has changed and the correlation
     * function must be reset
     */
    double k;
    int ij;

    for (int i = 0; i < s->local_n1; i++)
    {
        for (int j = 0; j < s->N; j++)
        {
            k = calc_k(i + s->local_1_start, j, s->N, s->dx);
            ij = i * s->N + j;
            s->C[ij]  = exp(-ipow (s->sigma * s->k0, 2.0) / (2.0 * s->rho * s->beta));
            s->C[ij] *= exp(-ipow (k - s->k0, 2.0) / (2.0 * ipow(s->alphac, 2)));
        }
    }

    return;
}

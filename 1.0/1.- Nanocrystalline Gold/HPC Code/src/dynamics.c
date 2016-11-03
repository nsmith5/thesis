#include <math.h>
#include <complex.h>
#include <fftw3-mpi.h>

#include "state.h"
#include "dynamics.h"
#include "random.h"

static void normalize (state *s, fftw_complex *field)
{
    /* Normalize a fourier space field after fft */

    int ij;
    double norm_scale = 1.0 /(s->N * s->N);

    for (int i = 0; i < s->local_n0; i++)
    {
        for (int j = 0; j < (s->N>>1) + 1; j++)
        {
            ij = i*(s->N/2 + 1) + j;

            field[ij] *= norm_scale;
        }
    }

    MPI_Barrier (MPI_COMM_WORLD);

    return;
}

static double F_mix (double c, state *s)
{
    double epsilon = -4.0 + s->epsilon0*(s->sigma - s->sigma0);

    return s->omega*(
                c*log(2.0*c) + (1-c)*log(2.0-2.0*c) +
                0.5 * epsilon * pow(c-0.5, 2)
            );
}

static double dF_mixdc (double c, state *s)
{
    return s->omega*log(c/(1.0-c));
}

static void noise(state *s)
{
    fftw_complex c_scale, n_scale;
    c_scale = I*sqrt(s->kbT * s->Mc / s->dx / s->dx / s->dt);
    n_scale = I*sqrt(s->kbT * s->Mn / s->dx / s->dx / s->dt);

    int ij;

    for (int i = 0; i < s->local_n0; i++)
    {
        for (int j = 0; j < (s->N>>1) + 1; j++)
        {
            ij = i*((s->N>>1) + 1) + j;

            s->fxin[ij] = calc_k (i + s->local_0_start, j, s->N, s->dx) *
                          rand_normal_complex () * n_scale;

            s->fxic[ij] = calc_k (i + s->local_0_start, j, s->N, s->dx) *
                          rand_normal_complex () * c_scale;
        }
    }

    MPI_Barrier (MPI_COMM_WORLD);

    return;
}

static void set_nonlinear (state *s)
{
    double third = 1.0/3.0;
    int ij;

    for (int i = 0; i < s->local_n0; i++)
    {
        for (int j = 0; j < s->N; j++)
        {
            ij = i*2*((s->N>>1) + 1) + j;

            s->nnl[ij]  = pow(s->n[ij], 2)*( -0.5 * s->eta + third * s->n[ij]);
            s->nnl[ij] += F_mix (s->c[ij], s);
            s->nnl[ij] += -exp(-pow(s->c[ij] - 1.0, 2)/(2.0*pow(s->alphac, 2)));

            s->cnl[ij]  = (s->n[ij] + 1.0) *  dF_mixdc (s->c[ij], s);
            s->cnl[ij] += 0.5 * s->n[ij] * (s->c[ij] - 1.0) / pow(s->alphac, 2) *
                          exp(-pow(s->c[ij] - 1.0, 2)/(2.0*pow(s->alphac, 2))) *
                          s->Cn[ij];
        }
    }

    MPI_Barrier (MPI_COMM_WORLD);

    return;
}

static void calccorr (state *s)
{
    int ij;

    for (int i = 0; i < s->local_n0; i++)
    {
        for (int j = 0; j < (s->N>>1) + 1; j++)
        {
            ij = i*((s->N>>1) + 1) + j;

            s->fCn[ij] = s->C[ij]*s->fn[ij];
        }
    }

    fftw_mpi_execute_dft_c2r (s->ifft_plan, s->fCn, s->Cn);

    return;
}

static void propagate (state *s)
{
    double epsilon = -4.0 + s->epsilon0 * (s->sigma - s->sigma0);
    int ij;

    for (int i = 0; i < s->local_n0; i++)
    {
        for (int j = 0; j < (s->N>>1) + 1; j++)
        {
            ij = i*((s->N>>1) + 1) + j;

            s->fn[ij]  = 1.0 / (1.0 + s->dt * s->k2[ij]);
            s->fn[ij] *= s->fn[ij] - s->dt * s->k2[ij] * s->fnnl[ij]; /*+ s->dt * s->fxin[ij] */

            s->fc[ij]  = 1.0 / (1.0 + s->dt * s->k2[ij] * (s->omega * epsilon + s->Wc * s->k2[ij]));
            s->fc[ij] *= s->fc[ij] - s->dt * s->k2[ij] * s->fcnl[ij];  /* + s->dt * s->fxic[ij] */
        }
    }

    MPI_Barrier (MPI_COMM_WORLD);

    return;
}

void step(state *s)
{
    /* Fourier transform fields */
    fftw_mpi_execute_dft_r2c (s->fft_plan, s->n, s->fn);
    fftw_mpi_execute_dft_r2c (s->fft_plan, s->c, s->fc);
    normalize (s, s->fn);
    normalize (s, s->fc);

    /* Calculate (C * n) and the nonlinear term */
    calccorr (s);
    set_nonlinear (s);

    /* Fourier transform the nonlinear terms */
    fftw_mpi_execute_dft_r2c (s->fft_plan, s->nnl, s->fnnl);
    fftw_mpi_execute_dft_r2c (s->fft_plan, s->cnl, s->fcnl);
    normalize (s, s->fnnl);
    normalize (s, s->fcnl);

    /* Make some noise and propagate the fields in fourier space */
    noise (s);
    propagate (s);

    /* Inverse fourier transform the propagated fields */
    fftw_mpi_execute_dft_c2r (s->ifft_plan, s->fn, s->n);
    fftw_mpi_execute_dft_c2r (s->ifft_plan, s->fc, s->c);

    /* Update the time and step #   */
    s->t += s->dt;
    s->step += 1;

    MPI_Barrier (MPI_COMM_WORLD);

    return;
}
